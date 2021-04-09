# Use 'groundhog' package to use versions of packages at time of publication
  groundhog::groundhog.library(c("lme4","MuMIn"),"2020-12-20")

# Read data file of calculated taper values
  TaperSample <- read.csv("DataFiles/DataFile_TaperParameterSample.csv")
  
#### 1. Compare models explaining variation in taper parameter ####
    
    # Models WITHOUT family random effect
    
    # Null model
    modelnull <- lmer(b1.iso~ (1|Site), data = TaperSample)
    # With one variable
    model1a <- lmer(b1.iso~log(DBH) + (1|Site), data = TaperSample)
    model1b <- lmer(b1.iso~log(HOM) + (1|Site), data = TaperSample)
    model1c <- lmer(b1.iso~log(WSG) + (1|Site), data = TaperSample)
    # With two variables
    model2a <- lmer(b1.iso~log(DBH) + log(HOM) + (1|Site), data = TaperSample)
    model2b <- lmer(b1.iso~log(DBH) + log(WSG) + (1|Site), data = TaperSample)
    model2c <- lmer(b1.iso~log(HOM) + log(WSG) + (1|Site), data = TaperSample)
    # With three variables
    model3a <- lmer(b1.iso~log(DBH) + log(HOM) + log(WSG) + (1|Site), data = TaperSample)
    

    # Make a table comparing models (for main text)    
    ModelCompare <- anova(modelnull, model1a, model1b, model1c, model2a, model2b, model2c, model3a)
    ModelCompare$R2marginal <- c(r.squaredGLMM(modelnull)[1], 
                                 r.squaredGLMM(model1a)[1], 
                                 r.squaredGLMM(model1b)[1], 
                                 r.squaredGLMM(model1c)[1], 
                                 r.squaredGLMM(model2a)[1], 
                                 r.squaredGLMM(model2b)[1], 
                                 r.squaredGLMM(model2c)[1],
                                 r.squaredGLMM(model3a)[1])
    ModelCompare$R2conditional <- c(r.squaredGLMM(modelnull)[2], 
                                    r.squaredGLMM(model1a)[2], 
                                    r.squaredGLMM(model1b)[2], 
                                    r.squaredGLMM(model1c)[2], 
                                    r.squaredGLMM(model2a)[2], 
                                    r.squaredGLMM(model2b)[2], 
                                    r.squaredGLMM(model2c)[2],
                                    r.squaredGLMM(model3a)[2])
    ModelCompare$Description <- c(paste("b = ",round(summary(modelnull)$coefficients[1],3)," + Site", sep=""),
                                  paste("b = ",round(summary(model1a)$coefficients[1,1],3)," + ",round(summary(model1a)$coefficients[2,1],3),"*log(DAB) + Site", sep=""),
                                  paste("b = ",round(summary(model1b)$coefficients[1,1],3)," + ",round(summary(model1b)$coefficients[2,1],3),"*log(HOM) + Site", sep=""),
                                  paste("b = ",round(summary(model1c)$coefficients[1,1],3)," + ",round(summary(model1c)$coefficients[2,1],3),"*log(WSG) + Site", sep=""),
                                  paste("b = ",round(summary(model2a)$coefficients[1,1],3)," + ",round(summary(model2a)$coefficients[2,1],3),"*log(DAB) ",round(summary(model2a)$coefficients[3,1],3),"*log(HOM) ", "+ Site", sep=""),
                                  paste("b = ",round(summary(model2b)$coefficients[1,1],3)," + ",round(summary(model2b)$coefficients[2,1],3),"*log(DAB) ",round(summary(model2b)$coefficients[3,1],3),"*log(WSG) ", "+ Site", sep=""),
                                  paste("b = ",round(summary(model2c)$coefficients[1,1],3)," + ",round(summary(model2c)$coefficients[2,1],3),"*log(HOM) ",round(summary(model2c)$coefficients[3,1],3),"*log(WSG) ", "+ Site", sep=""),
                                  paste("b = ",round(summary(model3a)$coefficients[1,1],3)," + ",round(summary(model3a)$coefficients[2,1],3),"*log(DAB) ",round(summary(model3a)$coefficients[3,1],3),"*log(HOM) ",round(summary(model3a)$coefficients[4,1],3),"*log(WSG) ", "+ Site", sep=""))
    ModelCompare$Description <- as.character(ModelCompare$Description)
    ModelCompare$RMSE <- c((mean((fitted(modelnull)-TaperSample$b1.iso)^2))^0.5,
                           (mean((fitted(model1a)-TaperSample$b1.iso)^2))^0.5,
                           (mean((fitted(model1b)-TaperSample$b1.iso)^2))^0.5,
                           (mean((fitted(model1c)-TaperSample$b1.iso)^2))^0.5,
                           (mean((fitted(model2a)-TaperSample$b1.iso)^2))^0.5,
                           (mean((fitted(model2b)-TaperSample$b1.iso)^2))^0.5,
                           (mean((fitted(model2c)-TaperSample$b1.iso)^2))^0.5,
                           (mean((fitted(model3a)-TaperSample$b1.iso)^2))^0.5)
    
    ModelResults <- ModelCompare[,c("Description","Df","AIC","R2marginal","R2conditional","RMSE")]
    ModelResults <- ModelResults[order(ModelResults$AIC),]
    ModelResults$dAIC <- ModelResults$AIC-min(ModelResults$AIC)
    
  # Verify that residuals of best model don't vary by site using a 1-way ANOVA
    ModelResiduals <- data.frame(residual = residuals(model3a),
                                 site = TaperSample$Site)
    residualANOVA <- lm(residual~site, data=ModelResiduals)
    summary(residualANOVA)
    
#### 2. Save results in Table 2 ####    
    write.csv(ModelResults, file="ResultsFiles/Table2_TaperMixedModelComparison.csv",
              row.names = F)

#### 3. Taper variation among sites and families ####

     # Krusal-Wallis tests for variation among sites
      b.Test <- kruskal.test(list(TaperSample[TaperSample$Site=='AMA','b1.iso'],
                                  TaperSample[TaperSample$Site=='BCI','b1.iso'],
                                  TaperSample[TaperSample$Site=='BKT','b1.iso'],
                                  TaperSample[TaperSample$Site=='HKK','b1.iso'],
                                  TaperSample[TaperSample$Site=='KCH','b1.iso']))
      
      b.Test30 <- kruskal.test(list(TaperSample[TaperSample$Site=='AMA' & TaperSample$DBH>=30,'b1.iso'],
                                    TaperSample[TaperSample$Site=='BCI' & TaperSample$DBH>=30,'b1.iso'],
                                    TaperSample[TaperSample$Site=='BKT' & TaperSample$DBH>=30,'b1.iso'],
                                    TaperSample[TaperSample$Site=='HKK' & TaperSample$DBH>=30,'b1.iso'],
                                    TaperSample[TaperSample$Site=='KCH' & TaperSample$DBH>=30,'b1.iso']))
      
      # Are there significant differences among families?
        kruskal.test(x = TaperSample$b1.iso, g = TaperSample$Family)
        
      # additional variation explained by family
        model3a_fam <- lmer(b1.iso~log(DBH) + log(HOM) + log(WSG) + (1|Site) + (1|Family), data = TaperSample)
        r.squaredGLMM(model3a_fam)

#### 4. Make Table 3 with taper by site ####
  
  # Calculate EDBH, basal area, and AGB for each tree
    agb.allometry <- function(E,wsg,dbh) {exp(-1.803-0.976*E+0.976*log(wsg)
                                              + 2.673*log(dbh) - 0.0299*(log(dbh)^2))}
    
    taper.eqn <- function(d,h,b1) {d/(exp(-b1*(h-1.3)))}
    
    TaperSample$EDBH <- taper.eqn(d = TaperSample$DBH,
                                  h = TaperSample$HOM,
                                  b1 = TaperSample$b1.iso)
    
    TaperSample$BA <- pi*(TaperSample$EDBH/2)^2
    
    # Define environmental factor "E" used in Chave et al. 2014 allometry for each site
    # Previously run code (takes a while to run) to find values
    # Retrieve values of E for each plot using database from Chave et al. 2014
        # Call in source info:
        #source("http://chave.ups-tlse.fr/pantropical_allometry/readlayers.r")
        # Define coordinates of each plot:
        #AMA.coord <- cbind(-70.27,-3.81); E.AMA <- retrieve_raster("E",AMA.coord)
        #BCI.coord <- cbind(-79.85, 9.15); E.BCI <- retrieve_raster("E",BCI.coord)
        #BKT.coord <- cbind(103.75,1.25); E.BKT <- retrieve_raster("E",BKT.coord)
        #HKK.coord <- cbind(99.22, 15.63); E.HKK <- retrieve_raster("E",HKK.coord)
        #KCH.coord <- cbind(99.80, 7.54); E.KCH <- retrieve_raster("E",KCH.coord)
    TaperSample$E <- NA
      TaperSample[TaperSample$Site=="AMA","E"] <- -0.07928769
      TaperSample[TaperSample$Site=="BCI","E"] <- 0.04944549
      TaperSample[TaperSample$Site=="BKT","E"] <- -0.05956875
      TaperSample[TaperSample$Site=="HKK","E"] <- 0.3194663
      TaperSample[TaperSample$Site=="KCH","E"] <- 0.04786947
      
    TaperSample$AGB <- agb.allometry(E = TaperSample$E,
                                     wsg = TaperSample$WSG,
                                     dbh = TaperSample$DBH)
    
    TaperBySite <- data.frame(Site = c("AMA","BCI","BKT","HKK","KCH"),
                              n = NA,
                              Mean = NA,
                              SD = NA,
                              Mean.BA = NA,
                              Mean.AGB = NA)
    
    for(i in 1:length(TaperBySite$Site)){
      TaperBySite[i,"n"] <- length(TaperSample[TaperSample$Site==TaperBySite$Site[i],"b1.iso"])
      TaperBySite[i,"Mean"] <- mean(TaperSample[TaperSample$Site==TaperBySite$Site[i],"b1.iso"])
      TaperBySite[i,"SD"] <- sd(TaperSample[TaperSample$Site==TaperBySite$Site[i],"b1.iso"])
      TaperBySite[i,"Mean.BA"] <- weighted.mean(x = TaperSample[TaperSample$Site==TaperBySite$Site[i],"b1.iso"],
                                                w = TaperSample[TaperSample$Site==TaperBySite$Site[i],"BA"])
      TaperBySite[i,"Mean.AGB"] <- weighted.mean(x = TaperSample[TaperSample$Site==TaperBySite$Site[i],"b1.iso"],
                                                 w = TaperSample[TaperSample$Site==TaperBySite$Site[i],"AGB"])
    }
    
    write.csv(TaperBySite, file="ResultsFiles/Table3_TaperVariationBySite.csv")
      
#### 5. Circularity variation among sites and families ####
    CircSample <- read.csv("DataFiles/DataFile_CircSample.csv")
    
    # Krusal-Wallis tests for difference among sites
      iso.Test <- kruskal.test(list(CircSample[CircSample$Site=='AMA','iso'],
                                  CircSample[CircSample$Site=='BCI','iso'],
                                  CircSample[CircSample$Site=='BKT','iso'],
                                  CircSample[CircSample$Site=='HKK','iso'],
                                  CircSample[CircSample$Site=='KCH','iso']))
      
      iso.Test30 <- kruskal.test(list(CircSample[CircSample$Site=='AMA' & CircSample$DBH>=30,'iso'],
                                    CircSample[CircSample$Site=='BCI' & CircSample$DBH>=30,'iso'],
                                    CircSample[CircSample$Site=='BKT' & CircSample$DBH>=30,'iso'],
                                    CircSample[CircSample$Site=='HKK' & CircSample$DBH>=30,'iso'],
                                    CircSample[CircSample$Site=='KCH' & CircSample$DBH>=30,'iso']))
      
  # Are there significant differences among families?
    kruskal.test(x = CircSample$iso, g = CircSample$Family)
        
#### 5. Make Table S6 with circularity by site ####
    
  # Calculate EDBH, basal area, and AGB for each tree
 
    CircSample$EDBH <- taper.eqn(d = CircSample$DBH,
                                  h = CircSample$HOM,
                                  b1 = CircSample$b1.iso)
    
    CircSample$BA <- pi*(CircSample$EDBH/2)^2
    
    CircSample$E <- NA
      CircSample[CircSample$Site=="AMA","E"] <- -0.07928769
      CircSample[CircSample$Site=="BCI","E"] <- 0.04944549
      CircSample[CircSample$Site=="BKT","E"] <- -0.05956875
      CircSample[CircSample$Site=="HKK","E"] <- 0.3194663
      CircSample[CircSample$Site=="KCH","E"] <- 0.04786947
      
    CircSample$AGB <- agb.allometry(E = CircSample$E,
                                     wsg = CircSample$WSG,
                                     dbh = CircSample$DBH)
    
    CircBySite <- data.frame(Site = c("AMA","BCI","BKT","HKK","KCH"),
                             n = NA,
                              Mean = NA,
                              SD = NA,
                              Mean.BA = NA,
                              Mean.AGB = NA)
    
    for(i in 1:length(CircBySite$Site)){
      CircBySite[i,"n"] <- length(CircSample[CircSample$Site==CircBySite$Site[i],"iso"])
      CircBySite[i,"Mean"] <- mean(CircSample[CircSample$Site==CircBySite$Site[i],"iso"])
      CircBySite[i,"SD"] <- sd(CircSample[CircSample$Site==CircBySite$Site[i],"iso"])
      CircBySite[i,"Mean.BA"] <- weighted.mean(x = CircSample[CircSample$Site==CircBySite$Site[i],"iso"],
                                                w = CircSample[CircSample$Site==CircBySite$Site[i],"BA"])
      CircBySite[i,"Mean.AGB"] <- weighted.mean(x = CircSample[CircSample$Site==CircBySite$Site[i],"iso"],
                                                 w = CircSample[CircSample$Site==CircBySite$Site[i],"AGB"])
    }
    
    write.csv(CircBySite, file="ResultsFiles/TableS6_CircVariationTable.csv")
