# Read data file of calculated taper values
  TaperSample <- read.csv("DataFiles/DataFile_TaperParameterSample.csv")
  sites <- c("AMA","BCI","BKT","HKK","KCH")

#### Estimate AGB with alternate taper models ####
        
  # Define necessary functions
    agb.allometry <- function(E,wsg,dbh) {exp(-1.803-0.976*E+0.976*log(wsg)
                                              + 2.673*log(dbh) - 0.0299*(log(dbh)^2))}
    
    taper.eqn <- function(d,h,b1) {d/(exp(-b1*(h-1.3)))}
    
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
      
  # Calculate AGB with various taper scenarios
    # Uncorrected AGB
      TaperSample$AGB_Uncorected <- agb.allometry(E = TaperSample$E,
                                                 wsg = TaperSample$WSG,
                                                 dbh = TaperSample$DBH)
      
    # AGB without site and wood density differences to sort by diameter size
      TaperSample$AGB_Sort <- agb.allometry(E = mean(TaperSample$E),
                                            wsg = mean(TaperSample$WSG),
                                            dbh = TaperSample$DBH)  
      
    # AGB with measured taper
      TaperSample$EDBH_MeasTaper <- taper.eqn(d = TaperSample$DBH,
                                              h = TaperSample$HOM,
                                              b1 = TaperSample$b1.iso)
      TaperSample$AGB_MeasTaper <- agb.allometry(E = TaperSample$E,
                                                 wsg = TaperSample$WSG,
                                                 dbh = TaperSample$EDBH_MeasTaper)    
        
  # For each site, use bootstrapping to estimate AGB with CI's under various scenarios
    nboot <- 1000
    
    # Amacayacu
    
        AGB_AMA <- data.frame(sample = 1:nboot,
                          AGB_FullModel = NA,
                          AGB_Cross3Pr = NA,
                          AGB_Cross2Pr = NA,
                          AGB_CrossMed = NA,
                          b_med = NA)
        
        set.seed(1)
        for(i in 1:nboot){
          data.site <- TaperSample[TaperSample$Site=="AMA",]
          
          # Full taper model
          sample.full <- sample(x = 1:dim(TaperSample)[1],
                                size = dim(TaperSample)[1],
                                replace = T)
          data.full <- TaperSample[sample.full,]
          model.full <- lmer(b1.iso~log(DBH) + log(HOM) + log(WSG) + (1|Site), data = data.full)
          data.site$b <- predict(model.full, newdata = data.site)
          data.site$EDBH <- taper.eqn(d = data.site$DBH,
                                      h = data.site$HOM,
                                      b1 = data.site$b)
          data.site$AGB <- agb.allometry(E = data.site$E,
                                         wsg = data.site$WSG,
                                         dbh = data.site$EDBH)
          AGB_AMA$AGB_FullModel[i] <- sum(data.site$AGB)
          
          # 3-parameter cross-validation
          data.LOO <- data.full[!(data.full$Site=="AMA"),]
          model.3pr <- lmer(b1.iso~log(DBH) + log(HOM) + log(WSG) + (1|Site), data = data.LOO)
            model.3pr.int <- summary(model.3pr)$coefficients[1,1]
            model.3pr.dbh <- summary(model.3pr)$coefficients[2,1]
            model.3pr.hom <- summary(model.3pr)$coefficients[3,1]
            model.3pr.wsg <- summary(model.3pr)$coefficients[4,1]
          data.site$b <- model.3pr.int + model.3pr.dbh*log(data.site$DBH) + model.3pr.hom*log(data.site$HOM) + model.3pr.wsg*log(data.site$WSG)

          data.site$EDBH <- taper.eqn(d = data.site$DBH,
                                      h = data.site$HOM,
                                      b1 = data.site$b)
          data.site$AGB <- agb.allometry(E = data.site$E,
                                         wsg = data.site$WSG,
                                         dbh = data.site$EDBH)
          AGB_AMA$AGB_Cross3Pr[i] <- sum(data.site$AGB)
          
          # 2-parameter cross-validation
          data.LOO <- data.full[!(data.full$Site=="AMA"),]
          model.2pr <- lmer(b1.iso~log(DBH) + log(HOM) + (1|Site), data = data.LOO)
            model.2pr.int <- summary(model.2pr)$coefficients[1,1]
            model.2pr.dbh <- summary(model.2pr)$coefficients[2,1]
            model.2pr.hom <- summary(model.2pr)$coefficients[3,1]
          data.site$b <- model.2pr.int + model.2pr.dbh*log(data.site$DBH) + model.2pr.hom*log(data.site$HOM)

          data.site$EDBH <- taper.eqn(d = data.site$DBH,
                                      h = data.site$HOM,
                                      b1 = data.site$b)
          data.site$AGB <- agb.allometry(E = data.site$E,
                                         wsg = data.site$WSG,
                                         dbh = data.site$EDBH)
          AGB_AMA$AGB_Cross2Pr[i] <- sum(data.site$AGB)
          
          # "Median" taper value--tree where half of AGB is from smaller trees
          data.LOO <- data.LOO[order(data.LOO$AGB_Sort),]
          data.LOO$CmAGB <- NA
          for(j in 1:dim(data.LOO)[1]){
            data.LOO$CmAGB[j] <- sum(data.LOO$AGB_Sort[1:j])
          }
          
          medDBH <- data.LOO[data.LOO$CmAGB > sum(data.LOO$AGB_MeasTaper)/2,"DBH"][1]
          model.DAB <- lm(b1.iso~log(DBH), data = data.LOO)
          b_med <- predict(model.DAB, newdata = data.frame(DBH = medDBH))
          data.site$EDBH <- taper.eqn(d = data.site$DBH,
                                      h = data.site$HOM,
                                      b1 = b_med)
          data.site$AGB <- agb.allometry(E = data.site$E,
                                         wsg = data.site$WSG,
                                         dbh = data.site$EDBH)
          AGB_AMA$AGB_CrossMed[i] <- sum(data.site$AGB)
          AGB_AMA$b_med[i] <- b_med
        }

    # Barro Colorado
    
        AGB_BCI <- data.frame(sample = 1:nboot,
                          AGB_FullModel = NA,
                          AGB_Cross3Pr = NA,
                          AGB_Cross2Pr = NA,
                          AGB_CrossMed = NA,
                          b_med = NA)
        
        set.seed(1)
        for(i in 1:nboot){
          data.site <- TaperSample[TaperSample$Site=="BCI",]
          
          # Full taper model
          sample.full <- sample(x = 1:dim(TaperSample)[1],
                                size = dim(TaperSample)[1],
                                replace = T)
          data.full <- TaperSample[sample.full,]
          model.full <- lmer(b1.iso~log(DBH) + log(HOM) + log(WSG) + (1|Site), data = data.full)
          data.site$b <- predict(model.full, newdata = data.site)
          data.site$EDBH <- taper.eqn(d = data.site$DBH,
                                      h = data.site$HOM,
                                      b1 = data.site$b)
          data.site$AGB <- agb.allometry(E = data.site$E,
                                         wsg = data.site$WSG,
                                         dbh = data.site$EDBH)
          AGB_BCI$AGB_FullModel[i] <- sum(data.site$AGB)
          
          # 3-parameter cross-validation
          data.LOO <- data.full[!(data.full$Site=="BCI"),]
          model.3pr <- lmer(b1.iso~log(DBH) + log(HOM) + log(WSG) + (1|Site), data = data.LOO)
            model.3pr.int <- summary(model.3pr)$coefficients[1,1]
            model.3pr.dbh <- summary(model.3pr)$coefficients[2,1]
            model.3pr.hom <- summary(model.3pr)$coefficients[3,1]
            model.3pr.wsg <- summary(model.3pr)$coefficients[4,1]
          data.site$b <- model.3pr.int + model.3pr.dbh*log(data.site$DBH) + model.3pr.hom*log(data.site$HOM) + model.3pr.wsg*log(data.site$WSG)

          data.site$EDBH <- taper.eqn(d = data.site$DBH,
                                      h = data.site$HOM,
                                      b1 = data.site$b)
          data.site$AGB <- agb.allometry(E = data.site$E,
                                         wsg = data.site$WSG,
                                         dbh = data.site$EDBH)
          AGB_BCI$AGB_Cross3Pr[i] <- sum(data.site$AGB)
          
          # 2-parameter cross-validation
          data.LOO <- data.full[!(data.full$Site=="BCI"),]
          model.2pr <- lmer(b1.iso~log(DBH) + log(HOM) + (1|Site), data = data.LOO)
            model.2pr.int <- summary(model.2pr)$coefficients[1,1]
            model.2pr.dbh <- summary(model.2pr)$coefficients[2,1]
            model.2pr.hom <- summary(model.2pr)$coefficients[3,1]
          data.site$b <- model.2pr.int + model.2pr.dbh*log(data.site$DBH) + model.2pr.hom*log(data.site$HOM)

          data.site$EDBH <- taper.eqn(d = data.site$DBH,
                                      h = data.site$HOM,
                                      b1 = data.site$b)
          data.site$AGB <- agb.allometry(E = data.site$E,
                                         wsg = data.site$WSG,
                                         dbh = data.site$EDBH)
          AGB_BCI$AGB_Cross2Pr[i] <- sum(data.site$AGB)
          
          # "Median" taper value--tree where half of AGB is from smaller trees
          data.LOO <- data.LOO[order(data.LOO$AGB_Sort),]
          data.LOO$CmAGB <- NA
          for(j in 1:dim(data.LOO)[1]){
            data.LOO$CmAGB[j] <- sum(data.LOO$AGB_Sort[1:j])
          }
          
          medDBH <- data.LOO[data.LOO$CmAGB > sum(data.LOO$AGB_MeasTaper)/2,"DBH"][1]
          model.DAB <- lm(b1.iso~log(DBH), data = data.LOO)
          b_med <- predict(model.DAB, newdata = data.frame(DBH = medDBH))
          data.site$EDBH <- taper.eqn(d = data.site$DBH,
                                      h = data.site$HOM,
                                      b1 = b_med)
          data.site$AGB <- agb.allometry(E = data.site$E,
                                         wsg = data.site$WSG,
                                         dbh = data.site$EDBH)
          AGB_BCI$AGB_CrossMed[i] <- sum(data.site$AGB)
          AGB_BCI$b_med[i] <- b_med
        }

    # Bukit Timah
    
        AGB_BKT <- data.frame(sample = 1:nboot,
                          AGB_FullModel = NA,
                          AGB_Cross3Pr = NA,
                          AGB_Cross2Pr = NA,
                          AGB_CrossMed = NA,
                          b_med = NA)
        
        set.seed(1)
        for(i in 1:nboot){
          data.site <- TaperSample[TaperSample$Site=="BKT",]
          
          # Full taper model
          sample.full <- sample(x = 1:dim(TaperSample)[1],
                                size = dim(TaperSample)[1],
                                replace = T)
          data.full <- TaperSample[sample.full,]
          model.full <- lmer(b1.iso~log(DBH) + log(HOM) + log(WSG) + (1|Site), data = data.full)
          data.site$b <- predict(model.full, newdata = data.site)
          data.site$EDBH <- taper.eqn(d = data.site$DBH,
                                      h = data.site$HOM,
                                      b1 = data.site$b)
          data.site$AGB <- agb.allometry(E = data.site$E,
                                         wsg = data.site$WSG,
                                         dbh = data.site$EDBH)
          AGB_BKT$AGB_FullModel[i] <- sum(data.site$AGB)
          
          # 3-parameter cross-validation
          data.LOO <- data.full[!(data.full$Site=="BKT"),]
          model.3pr <- lmer(b1.iso~log(DBH) + log(HOM) + log(WSG) + (1|Site), data = data.LOO)
            model.3pr.int <- summary(model.3pr)$coefficients[1,1]
            model.3pr.dbh <- summary(model.3pr)$coefficients[2,1]
            model.3pr.hom <- summary(model.3pr)$coefficients[3,1]
            model.3pr.wsg <- summary(model.3pr)$coefficients[4,1]
          data.site$b <- model.3pr.int + model.3pr.dbh*log(data.site$DBH) + model.3pr.hom*log(data.site$HOM) + model.3pr.wsg*log(data.site$WSG)

          data.site$EDBH <- taper.eqn(d = data.site$DBH,
                                      h = data.site$HOM,
                                      b1 = data.site$b)
          data.site$AGB <- agb.allometry(E = data.site$E,
                                         wsg = data.site$WSG,
                                         dbh = data.site$EDBH)
          AGB_BKT$AGB_Cross3Pr[i] <- sum(data.site$AGB)
          
          # 2-parameter cross-validation
          data.LOO <- data.full[!(data.full$Site=="BKT"),]
          model.2pr <- lmer(b1.iso~log(DBH) + log(HOM) + (1|Site), data = data.LOO)
            model.2pr.int <- summary(model.2pr)$coefficients[1,1]
            model.2pr.dbh <- summary(model.2pr)$coefficients[2,1]
            model.2pr.hom <- summary(model.2pr)$coefficients[3,1]
          data.site$b <- model.2pr.int + model.2pr.dbh*log(data.site$DBH) + model.2pr.hom*log(data.site$HOM)

          data.site$EDBH <- taper.eqn(d = data.site$DBH,
                                      h = data.site$HOM,
                                      b1 = data.site$b)
          data.site$AGB <- agb.allometry(E = data.site$E,
                                         wsg = data.site$WSG,
                                         dbh = data.site$EDBH)
          AGB_BKT$AGB_Cross2Pr[i] <- sum(data.site$AGB)
          
          # "Median" taper value--tree where half of AGB is from smaller trees
          data.LOO <- data.LOO[order(data.LOO$AGB_Sort),]
          data.LOO$CmAGB <- NA
          for(j in 1:dim(data.LOO)[1]){
            data.LOO$CmAGB[j] <- sum(data.LOO$AGB_Sort[1:j])
          }
          
          medDBH <- data.LOO[data.LOO$CmAGB > sum(data.LOO$AGB_MeasTaper)/2,"DBH"][1]
          model.DAB <- lm(b1.iso~log(DBH), data = data.LOO)
          b_med <- predict(model.DAB, newdata = data.frame(DBH = medDBH))
          data.site$EDBH <- taper.eqn(d = data.site$DBH,
                                      h = data.site$HOM,
                                      b1 = b_med)
          data.site$AGB <- agb.allometry(E = data.site$E,
                                         wsg = data.site$WSG,
                                         dbh = data.site$EDBH)
          AGB_BKT$AGB_CrossMed[i] <- sum(data.site$AGB)
          AGB_BKT$b_med[i] <- b_med
        }  
        
    # Huai Kha Khaeng
    
        AGB_HKK <- data.frame(sample = 1:nboot,
                          AGB_FullModel = NA,
                          AGB_Cross3Pr = NA,
                          AGB_Cross2Pr = NA,
                          AGB_CrossMed = NA,
                          b_med = NA)
        
        set.seed(1)
        for(i in 1:nboot){
          data.site <- TaperSample[TaperSample$Site=="HKK",]
          
          # Full taper model
          sample.full <- sample(x = 1:dim(TaperSample)[1],
                                size = dim(TaperSample)[1],
                                replace = T)
          data.full <- TaperSample[sample.full,]
          model.full <- lmer(b1.iso~log(DBH) + log(HOM) + log(WSG) + (1|Site), data = data.full)
          data.site$b <- predict(model.full, newdata = data.site)
          data.site$EDBH <- taper.eqn(d = data.site$DBH,
                                      h = data.site$HOM,
                                      b1 = data.site$b)
          data.site$AGB <- agb.allometry(E = data.site$E,
                                         wsg = data.site$WSG,
                                         dbh = data.site$EDBH)
          AGB_HKK$AGB_FullModel[i] <- sum(data.site$AGB)
          
          # 3-parameter cross-validation
          data.LOO <- data.full[!(data.full$Site=="HKK"),]
          model.3pr <- lmer(b1.iso~log(DBH) + log(HOM) + log(WSG) + (1|Site), data = data.LOO)
            model.3pr.int <- summary(model.3pr)$coefficients[1,1]
            model.3pr.dbh <- summary(model.3pr)$coefficients[2,1]
            model.3pr.hom <- summary(model.3pr)$coefficients[3,1]
            model.3pr.wsg <- summary(model.3pr)$coefficients[4,1]
          data.site$b <- model.3pr.int + model.3pr.dbh*log(data.site$DBH) + model.3pr.hom*log(data.site$HOM) + model.3pr.wsg*log(data.site$WSG)

          data.site$EDBH <- taper.eqn(d = data.site$DBH,
                                      h = data.site$HOM,
                                      b1 = data.site$b)
          data.site$AGB <- agb.allometry(E = data.site$E,
                                         wsg = data.site$WSG,
                                         dbh = data.site$EDBH)
          AGB_HKK$AGB_Cross3Pr[i] <- sum(data.site$AGB)
          
          # 2-parameter cross-validation
          data.LOO <- data.full[!(data.full$Site=="HKK"),]
          model.2pr <- lmer(b1.iso~log(DBH) + log(HOM) + (1|Site), data = data.LOO)
            model.2pr.int <- summary(model.2pr)$coefficients[1,1]
            model.2pr.dbh <- summary(model.2pr)$coefficients[2,1]
            model.2pr.hom <- summary(model.2pr)$coefficients[3,1]
          data.site$b <- model.2pr.int + model.2pr.dbh*log(data.site$DBH) + model.2pr.hom*log(data.site$HOM)

          data.site$EDBH <- taper.eqn(d = data.site$DBH,
                                      h = data.site$HOM,
                                      b1 = data.site$b)
          data.site$AGB <- agb.allometry(E = data.site$E,
                                         wsg = data.site$WSG,
                                         dbh = data.site$EDBH)
          AGB_HKK$AGB_Cross2Pr[i] <- sum(data.site$AGB)
          
          # "Median" taper value--tree where half of AGB is from smaller trees
          data.LOO <- data.LOO[order(data.LOO$AGB_Sort),]
          data.LOO$CmAGB <- NA
          for(j in 1:dim(data.LOO)[1]){
            data.LOO$CmAGB[j] <- sum(data.LOO$AGB_Sort[1:j])
          }
          
          medDBH <- data.LOO[data.LOO$CmAGB > sum(data.LOO$AGB_MeasTaper)/2,"DBH"][1]
          model.DAB <- lm(b1.iso~log(DBH), data = data.LOO)
          b_med <- predict(model.DAB, newdata = data.frame(DBH = medDBH))
          data.site$EDBH <- taper.eqn(d = data.site$DBH,
                                      h = data.site$HOM,
                                      b1 = b_med)
          data.site$AGB <- agb.allometry(E = data.site$E,
                                         wsg = data.site$WSG,
                                         dbh = data.site$EDBH)
          AGB_HKK$AGB_CrossMed[i] <- sum(data.site$AGB)
          AGB_HKK$b_med[i] <- b_med
        }  
        
    # Khao Chong
    
        AGB_KCH <- data.frame(sample = 1:nboot,
                          AGB_FullModel = NA,
                          AGB_Cross3Pr = NA,
                          AGB_Cross2Pr = NA,
                          AGB_CrossMed = NA,
                          b_med = NA)
        
        set.seed(1)
        for(i in 1:nboot){
          data.site <- TaperSample[TaperSample$Site=="KCH",]
          
          # Full taper model
          sample.full <- sample(x = 1:dim(TaperSample)[1],
                                size = dim(TaperSample)[1],
                                replace = T)
          data.full <- TaperSample[sample.full,]
          model.full <- lmer(b1.iso~log(DBH) + log(HOM) + log(WSG) + (1|Site), data = data.full)
          data.site$b <- predict(model.full, newdata = data.site)
          data.site$EDBH <- taper.eqn(d = data.site$DBH,
                                      h = data.site$HOM,
                                      b1 = data.site$b)
          data.site$AGB <- agb.allometry(E = data.site$E,
                                         wsg = data.site$WSG,
                                         dbh = data.site$EDBH)
          AGB_KCH$AGB_FullModel[i] <- sum(data.site$AGB)
          
          # 3-parameter cross-validation
          data.LOO <- data.full[!(data.full$Site=="KCH"),]
          model.3pr <- lmer(b1.iso~log(DBH) + log(HOM) + log(WSG) + (1|Site), data = data.LOO)
            model.3pr.int <- summary(model.3pr)$coefficients[1,1]
            model.3pr.dbh <- summary(model.3pr)$coefficients[2,1]
            model.3pr.hom <- summary(model.3pr)$coefficients[3,1]
            model.3pr.wsg <- summary(model.3pr)$coefficients[4,1]
          data.site$b <- model.3pr.int + model.3pr.dbh*log(data.site$DBH) + model.3pr.hom*log(data.site$HOM) + model.3pr.wsg*log(data.site$WSG)

          data.site$EDBH <- taper.eqn(d = data.site$DBH,
                                      h = data.site$HOM,
                                      b1 = data.site$b)
          data.site$AGB <- agb.allometry(E = data.site$E,
                                         wsg = data.site$WSG,
                                         dbh = data.site$EDBH)
          AGB_KCH$AGB_Cross3Pr[i] <- sum(data.site$AGB)
          
          # 2-parameter cross-validation
          data.LOO <- data.full[!(data.full$Site=="KCH"),]
          model.2pr <- lmer(b1.iso~log(DBH) + log(HOM) + (1|Site), data = data.LOO)
            model.2pr.int <- summary(model.2pr)$coefficients[1,1]
            model.2pr.dbh <- summary(model.2pr)$coefficients[2,1]
            model.2pr.hom <- summary(model.2pr)$coefficients[3,1]
          data.site$b <- model.2pr.int + model.2pr.dbh*log(data.site$DBH) + model.2pr.hom*log(data.site$HOM)

          data.site$EDBH <- taper.eqn(d = data.site$DBH,
                                      h = data.site$HOM,
                                      b1 = data.site$b)
          data.site$AGB <- agb.allometry(E = data.site$E,
                                         wsg = data.site$WSG,
                                         dbh = data.site$EDBH)
          AGB_KCH$AGB_Cross2Pr[i] <- sum(data.site$AGB)
          
          # "Median" taper value--tree where half of AGB is from smaller trees
          data.LOO <- data.LOO[order(data.LOO$AGB_Sort),]
          data.LOO$CmAGB <- NA
          for(j in 1:dim(data.LOO)[1]){
            data.LOO$CmAGB[j] <- sum(data.LOO$AGB_Sort[1:j])
          }
          
          medDBH <- data.LOO[data.LOO$CmAGB > sum(data.LOO$AGB_MeasTaper)/2,"DBH"][1]
          model.DAB <- lm(b1.iso~log(DBH), data = data.LOO)
          b_med <- predict(model.DAB, newdata = data.frame(DBH = medDBH))
          data.site$EDBH <- taper.eqn(d = data.site$DBH,
                                      h = data.site$HOM,
                                      b1 = b_med)
          data.site$AGB <- agb.allometry(E = data.site$E,
                                         wsg = data.site$WSG,
                                         dbh = data.site$EDBH)
          AGB_KCH$AGB_CrossMed[i] <- sum(data.site$AGB)
          AGB_KCH$b_med[i] <- b_med
        }  
        
#### Summarize results ####
    AGB_list <- list(AGB_AMA, AGB_BCI, AGB_BKT, AGB_HKK, AGB_KCH)    
    
    AGB_Results <- data.frame(site=sites,
                              AGB_Uncorrected=NA,
                              AGB_MeasTaper=NA,
                              AGB_FullModel_Mean=NA,
                              AGB_FullModel_Min95=NA,
                              AGB_FullModel_Max95=NA,
                              AGB_Cross3Pr_Mean=NA,
                              AGB_Cross3Pr_Min95=NA,
                              AGB_Cross3Pr_Max95=NA,
                              AGB_Cross2Pr_Mean=NA,
                              AGB_Cross2Pr_Min95=NA,
                              AGB_Cross2Pr_Max95=NA,
                              AGB_CrossMed_Mean=NA,
                              AGB_CrossMed_Min95=NA,
                              AGB_CrossMed_Max95=NA,
                              AGB_MeasTaper_Std=NA,
                              AGB_FullModel_Mean_Std=NA,
                              AGB_FullModel_Min95_Std=NA,
                              AGB_FullModel_Max95_Std=NA,
                              AGB_Cross3Pr_Mean_Std=NA,
                              AGB_Cross3Pr_Min95_Std=NA,
                              AGB_Cross3Pr_Max95_Std=NA,
                              AGB_Cross2Pr_Mean_Std=NA,
                              AGB_Cross2Pr_Min95_Std=NA,
                              AGB_Cross2Pr_Max95_Std=NA,
                              AGB_CrossMed_Mean_Std=NA,
                              AGB_CrossMed_Min95_Std=NA,
                              AGB_CrossMed_Max95_Std=NA)
    
    for(i in 1:length(sites)){
      AGB_Results$AGB_Uncorrected[i] <- sum(TaperSample[TaperSample$Site==sites[i],"AGB_Uncorected"])
      AGB_Results$AGB_MeasTaper[i] <- sum(TaperSample[TaperSample$Site==sites[i],"AGB_MeasTaper"])
      AGB_Results$AGB_FullModel_Mean[i] <- mean(AGB_list[[i]]$AGB_FullModel)
      AGB_Results$AGB_FullModel_Min95[i] <- quantile((AGB_list[[i]]$AGB_FullModel), 0.025)
      AGB_Results$AGB_FullModel_Max95[i] <- quantile((AGB_list[[i]]$AGB_FullModel), 0.975)
      AGB_Results$AGB_Cross3Pr_Mean[i] <- mean(AGB_list[[i]]$AGB_Cross3Pr)
      AGB_Results$AGB_Cross3Pr_Min95[i] <- quantile((AGB_list[[i]]$AGB_Cross3Pr), 0.025)
      AGB_Results$AGB_Cross3Pr_Max95[i] <- quantile((AGB_list[[i]]$AGB_Cross3Pr), 0.975)
      AGB_Results$AGB_Cross2Pr_Mean[i] <- mean(AGB_list[[i]]$AGB_Cross2Pr)
      AGB_Results$AGB_Cross2Pr_Min95[i] <- quantile((AGB_list[[i]]$AGB_Cross2Pr), 0.025)
      AGB_Results$AGB_Cross2Pr_Max95[i] <- quantile((AGB_list[[i]]$AGB_Cross2Pr), 0.975)
      AGB_Results$AGB_CrossMed_Mean[i] <- mean(AGB_list[[i]]$AGB_CrossMed)
      AGB_Results$AGB_CrossMed_Min95[i] <- quantile((AGB_list[[i]]$AGB_CrossMed), 0.025)
      AGB_Results$AGB_CrossMed_Max95[i] <- quantile((AGB_list[[i]]$AGB_CrossMed), 0.975)
    }
    
    # Calculate standardized values -- compared to uncorrected data
      AGB_Results$AGB_MeasTaper_Std <- AGB_Results$AGB_MeasTaper/AGB_Results$AGB_Uncorrected
      AGB_Results$AGB_FullModel_Mean_Std <- AGB_Results$AGB_FullModel_Mean/AGB_Results$AGB_Uncorrected
      AGB_Results$AGB_FullModel_Min95_Std <- AGB_Results$AGB_FullModel_Min95/AGB_Results$AGB_Uncorrected
      AGB_Results$AGB_FullModel_Max95_Std <- AGB_Results$AGB_FullModel_Max95/AGB_Results$AGB_Uncorrected
      AGB_Results$AGB_Cross3Pr_Mean_Std <- AGB_Results$AGB_Cross3Pr_Mean/AGB_Results$AGB_Uncorrected
      AGB_Results$AGB_Cross3Pr_Min95_Std <- AGB_Results$AGB_Cross3Pr_Min95/AGB_Results$AGB_Uncorrected
      AGB_Results$AGB_Cross3Pr_Max95_Std <- AGB_Results$AGB_Cross3Pr_Max95/AGB_Results$AGB_Uncorrected
      AGB_Results$AGB_Cross2Pr_Mean_Std <- AGB_Results$AGB_Cross2Pr_Mean/AGB_Results$AGB_Uncorrected
      AGB_Results$AGB_Cross2Pr_Min95_Std <- AGB_Results$AGB_Cross2Pr_Min95/AGB_Results$AGB_Uncorrected
      AGB_Results$AGB_Cross2Pr_Max95_Std <- AGB_Results$AGB_Cross2Pr_Max95/AGB_Results$AGB_Uncorrected
      AGB_Results$AGB_CrossMed_Mean_Std <- AGB_Results$AGB_CrossMed_Mean/AGB_Results$AGB_Uncorrected
      AGB_Results$AGB_CrossMed_Min95_Std <- AGB_Results$AGB_CrossMed_Min95/AGB_Results$AGB_Uncorrected
      AGB_Results$AGB_CrossMed_Max95_Std <- AGB_Results$AGB_CrossMed_Max95/AGB_Results$AGB_Uncorrected

      # Mean, min, max from applying taper correction
      round(mean(100*abs(AGB_Results$AGB_Uncorrected-AGB_Results$AGB_MeasTaper)/AGB_Results$AGB_Uncorrected),1)
      round(min(100*abs(AGB_Results$AGB_Uncorrected-AGB_Results$AGB_MeasTaper)/AGB_Results$AGB_Uncorrected),1)
      round(max(100*abs(AGB_Results$AGB_Uncorrected-AGB_Results$AGB_MeasTaper)/AGB_Results$AGB_Uncorrected),1)

      # Mean, min, max from applying Model 1
      round(mean(100*abs(AGB_Results$AGB_FullModel_Mean-AGB_Results$AGB_MeasTaper)/AGB_Results$AGB_MeasTaper),1)
      round(min(100*abs(AGB_Results$AGB_FullModel_Mean-AGB_Results$AGB_MeasTaper)/AGB_Results$AGB_MeasTaper),1)
      round(max(100*abs(AGB_Results$AGB_FullModel_Mean-AGB_Results$AGB_MeasTaper)/AGB_Results$AGB_MeasTaper),1)
      
      # Mean, min, max from applying cross-validated model with 3 parameters
      round(mean(100*abs(AGB_Results$AGB_Cross3Pr_Mean-AGB_Results$AGB_MeasTaper)/AGB_Results$AGB_MeasTaper),1)
      round(min(100*abs(AGB_Results$AGB_Cross3Pr_Mean-AGB_Results$AGB_MeasTaper)/AGB_Results$AGB_MeasTaper),1)
      round(max(100*abs(AGB_Results$AGB_Cross3Pr_Mean-AGB_Results$AGB_MeasTaper)/AGB_Results$AGB_MeasTaper),1)

      # Mean, min, max from applying cross-validated model with 2 parameters
      round(mean(100*abs(AGB_Results$AGB_Cross2Pr_Mean-AGB_Results$AGB_MeasTaper)/AGB_Results$AGB_MeasTaper),1)
      round(min(100*abs(AGB_Results$AGB_Cross2Pr_Mean-AGB_Results$AGB_MeasTaper)/AGB_Results$AGB_MeasTaper),1)
      round(max(100*abs(AGB_Results$AGB_Cross2Pr_Mean-AGB_Results$AGB_MeasTaper)/AGB_Results$AGB_MeasTaper),1)

      # Mean, min, max from applying cross-validated model with 2 parameters
      round(mean(100*abs(AGB_Results$AGB_CrossMed_Mean-AGB_Results$AGB_MeasTaper)/AGB_Results$AGB_MeasTaper),1)
      round(min(100*abs(AGB_Results$AGB_CrossMed_Mean-AGB_Results$AGB_MeasTaper)/AGB_Results$AGB_MeasTaper),1)
      round(max(100*abs(AGB_Results$AGB_CrossMed_Mean-AGB_Results$AGB_MeasTaper)/AGB_Results$AGB_MeasTaper),1)

#### Save results ####  
      
  write.csv(AGB_Results, file="ResultsFiles/BiomassEstimates.csv", row.names = F)
      