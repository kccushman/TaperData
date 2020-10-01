## NOTE: if working from the current Git repository, then workflow will function from step 2.
TaperSample <- read.csv("DataFile_TaperParameterSample.csv")
sites <- c("AMA","BCI","BKT","HKK","KCH")
sitesNames <- c("Amacayacu",
                "Barro Colorado",
                "Bukit Timah",
                "Huai Kha Khaeng",
                "Khao Chong")

# Load needed packages ### left for refernce; replaced with package::function notation throughout. 
  # library(splancs)
  # library(MASS)
  # library(lme4)
  # library(MuMIn)

# Load data about sampled trees
  TreeSample <- read.csv("~/Desktop/Taper/Current/TaperCorrection/TreeTaperSample.csv")
  TreeSample$ID <- paste(TreeSample$Site, TreeSample$Tag, sep="-")

##### 1. Calculate area, diameter, and circularity of all-cross sections #####  
  
  # Define function to calculate perimeter of polygon
      calc.perim = function(vertices) {
        dist.sum = sqrt((vertices[1,1]-vertices[length(vertices[,1]),1])^2 + (vertices[1,2]-vertices[length(vertices[,2]),2])^2)
        for (k in 2:length(vertices[,1])) {
          dist.k = sqrt((vertices[k,1]-vertices[k-1,1])^2 + (vertices[k,2]-vertices[k-1,2])^2)
          dist.sum = dist.sum + dist.k
        }
        return(dist.sum)
      }  
    
  # Create data frame structure to store contour results
      contours=data.frame(Tag=NA,
                          Site=NA,
                          ht=NA,
                          perim=NA,
                          area=NA,
                          d=NA,
                          iso.quo=NA)
#### 1.1 Amacayacu ####

    # List all trees for AMA
      AMA.trees <- list.files("~/Desktop/Taper/Current/TaperCorrection//RawTaperData/AMA/contours",
                              full.names = F, no..=F)
      AMA.files <- list.files("~/Desktop/Taper/Current/TaperCorrection//RawTaperData/AMA/contours",
                              full.names = T, no..=F)
      
    # For each tree, read in all cross-section point files
      for(i in 1:length(AMA.trees)){
        # Count numner of cross-sections for this tree
          n=length(list.files(AMA.files[i]))
        # Make a data frame to save results
          tree.data=data.frame(Tag=rep(AMA.trees[i],n),
                               Site=rep("AMA",n),
                               ht=NA,
                               perim=NA,
                               area=NA,
                               d=NA,
                               iso.quo=NA)
        # For each cross-section calculate height, perimeter, area, diameter
          for(j in 1:n){
            # Read in text file to make a table of all points in cross-section. Columns are X,Y,and Z-values of points. 
              vertices=read.table(paste(AMA.files[i],
                                        '/c',
                                        j,
                                        '.txt',
                                        sep=''),header=T)
            # Calculate height of cross-section
              tree.data$ht[j]=mean(vertices[,3])
            # Calculate area of cross-section  
              tree.data$area[j]=splancs::areapl(cbind(vertices[,1],vertices[,2]))
            # Calculate perimeter of cross-section  
              tree.data$perim[j] = calc.perim(vertices)	
            # Find convex hull
              tree.chull <- chull(vertices[,1:2])
              # Calculate area from convex hull
                perim.chull <- calc.perim(vertices[tree.chull,])
              # Estimate diameter from convex hull area
                tree.data$d[j] <- perim.chull/pi
          }
          
        # For all cross-sections, calculate circularity (iso.quo)
          tree.data$iso.quo = (4*pi*tree.data$area)/((tree.data$perim)^2)
          
        # Save info
          contours <- rbind(contours,tree.data)
      }
      contours <- contours[-1,] # get rid of first NA line
      
#### 1.2 Barro Coloardo ####
      # List all trees for BCI
      BCI.trees <- list.files("~/Desktop/Taper/Current/TaperCorrection//RawTaperData/BCI/contours",
                              full.names = F, no..=F)
      BCI.files <- list.files("~/Desktop/Taper/Current/TaperCorrection//RawTaperData/BCI/contours",
                              full.names = T, no..=F)
      
      # For each tree, read in all cross-section point files
      for(i in 1:length(BCI.trees)){
        # Count numner of cross-sections for this tree
        n=length(list.files(BCI.files[i]))
        # Make a data frame to save results
        tree.data=data.frame(Tag=rep(BCI.trees[i],n),
                             Site=rep("BCI",n),
                             ht=NA,
                             perim=NA,
                             area=NA,
                             d=NA,
                             iso.quo=NA)
        # For each cross-section calculate height, perimeter, area, diameter
        for(j in 1:n){
          # Read in text file to make a table of all points in cross-section. Columns are X,Y,and Z-values of points. 
          vertices=read.table(paste(BCI.files[i],
                                    '/c',
                                    j,
                                    '.txt',
                                    sep=''),header=T)
          # Calculate height of cross-section
          tree.data$ht[j]=mean(vertices[,3])
          # Calculate area of cross-section  
          tree.data$area[j]=splancs::areapl(cbind(vertices[,1],vertices[,2]))
          # Calculate perimeter of cross-section  
          tree.data$perim[j] = calc.perim(vertices)	
          # Find convex hull
          tree.chull <- chull(vertices[,1:2])
          # Calculate area from convex hull
          perim.chull <- calc.perim(vertices[tree.chull,])
          # Estimate diameter from convex hull area
          tree.data$d[j] <- perim.chull/pi
        }
        
        # For all cross-sections, calculate circularity (iso.quo)
        tree.data$iso.quo = (4*pi*tree.data$area)/((tree.data$perim)^2)
        
        # Save info
        contours <- rbind(contours,tree.data)
      }    
  
  
#### 1.3 Bukit Timah ####

      # List all trees for BKT
      BKT.trees <- list.files("~/Desktop/Taper/Current/TaperCorrection//RawTaperData/BKT/contours",
                              full.names = F, no..=F)
      BKT.files <- list.files("~/Desktop/Taper/Current/TaperCorrection//RawTaperData/BKT/contours",
                              full.names = T, no..=F)
      
      # For each tree, read in all cross-section point files
      for(i in 1:length(BKT.trees)){
        # Count numner of cross-sections for this tree
        n=length(list.files(BKT.files[i]))
        # Make a data frame to save results
        tree.data=data.frame(Tag=rep(BKT.trees[i],n),
                             Site=rep("BKT",n),
                             ht=NA,
                             perim=NA,
                             area=NA,
                             d=NA,
                             iso.quo=NA)
        # For each cross-section calculate height, perimeter, area, diameter
        for(j in 1:n){
          # Read in text file to make a table of all points in cross-section. Columns are X,Y,and Z-values of points. 
          vertices=read.table(paste(BKT.files[i],
                                    '/c',
                                    j,
                                    '.txt',
                                    sep=''),header=T)
          # Calculate height of cross-section
          tree.data$ht[j]=mean(vertices[,3])
          # Calculate area of cross-section  
          tree.data$area[j]=splancs::areapl(cbind(vertices[,1],vertices[,2]))
          # Calculate perimeter of cross-section  
          tree.data$perim[j] = calc.perim(vertices)	
          # Find convex hull
          tree.chull <- chull(vertices[,1:2])
          # Calculate area from convex hull
          perim.chull <- calc.perim(vertices[tree.chull,])
          # Estimate diameter from convex hull area
          tree.data$d[j] <- perim.chull/pi
        }
        
        # For all cross-sections, calculate circularity (iso.quo)
        tree.data$iso.quo = (4*pi*tree.data$area)/((tree.data$perim)^2)
        
        # Save info
        contours <- rbind(contours,tree.data)
      }    
      
      
#### 1.4 Huai Kha Khaeng ####

      # List all trees for HKK
      HKK.trees <- list.files("~/Desktop/Taper/Current/TaperCorrection//RawTaperData/HKK/contours",
                              full.names = F, no..=F)
      HKK.files <- list.files("~/Desktop/Taper/Current/TaperCorrection//RawTaperData/HKK/contours",
                              full.names = T, no..=F)
      
      # For each tree, read in all cross-section point files
      for(i in 1:length(HKK.trees)){
        # Count numner of cross-sections for this tree
        n=length(list.files(HKK.files[i]))
        # Make a data frame to save results
        tree.data=data.frame(Tag=rep(HKK.trees[i],n),
                             Site=rep("HKK",n),
                             ht=NA,
                             perim=NA,
                             area=NA,
                             d=NA,
                             iso.quo=NA)
        # For each cross-section calculate height, perimeter, area, diameter
        for(j in 1:n){
          # Read in text file to make a table of all points in cross-section. Columns are X,Y,and Z-values of points. 
          vertices=read.table(paste(HKK.files[i],
                                    '/c',
                                    j,
                                    '.txt',
                                    sep=''),header=T)
          # Calculate height of cross-section
          tree.data$ht[j]=mean(vertices[,3])
          # Calculate area of cross-section  
          tree.data$area[j]=splancs::areapl(cbind(vertices[,1],vertices[,2]))
          # Calculate perimeter of cross-section  
          tree.data$perim[j] = calc.perim(vertices)	
          # Find convex hull
          tree.chull <- chull(vertices[,1:2])
          # Calculate area from convex hull
          perim.chull <- calc.perim(vertices[tree.chull,])
          # Estimate diameter from convex hull area
          tree.data$d[j] <- perim.chull/pi
        }
        
        # For all cross-sections, calculate circularity (iso.quo)
        tree.data$iso.quo = (4*pi*tree.data$area)/((tree.data$perim)^2)
        
        # Save info
        contours <- rbind(contours,tree.data)
      }   

#### 1.5 Khao Chong ####

      # List all trees for KCH
      KCH.trees <- list.files("~/Desktop/Taper/Current/TaperCorrection//RawTaperData/KCH/contours",
                              full.names = F, no..=F)
      KCH.files <- list.files("~/Desktop/Taper/Current/TaperCorrection//RawTaperData/KCH/contours",
                              full.names = T, no..=F)
      
      # For each tree, read in all cross-section point files
      for(i in 1:length(KCH.trees)){
        # Count numner of cross-sections for this tree
        n=length(list.files(KCH.files[i]))
        # Make a data frame to save results
        tree.data=data.frame(Tag=rep(KCH.trees[i],n),
                             Site=rep("KCH",n),
                             ht=NA,
                             perim=NA,
                             area=NA,
                             d=NA,
                             iso.quo=NA)
        # For each cross-section calculate height, perimeter, area, diameter
        for(j in 1:n){
          # Read in text file to make a table of all points in cross-section. Columns are X,Y,and Z-values of points. 
          vertices=read.table(paste(KCH.files[i],
                                    '/c',
                                    j,
                                    '.txt',
                                    sep=''),header=T)
          # Calculate height of cross-section
          tree.data$ht[j]=mean(vertices[,3])
          # Calculate area of cross-section  
          tree.data$area[j]=splancs::areapl(cbind(vertices[,1],vertices[,2]))
          # Calculate perimeter of cross-section  
          tree.data$perim[j] = calc.perim(vertices)	
          # Find convex hull
          tree.chull <- chull(vertices[,1:2])
          # Calculate area from convex hull
          perim.chull <- calc.perim(vertices[tree.chull,])
          # Estimate diameter from convex hull area
          tree.data$d[j] <- perim.chull/pi
        }
        
        # For all cross-sections, calculate circularity (iso.quo)
        tree.data$iso.quo = (4*pi*tree.data$area)/((tree.data$perim)^2)
        
        # Save info
        contours <- rbind(contours,tree.data)
      }  

      contours$ID <- paste(contours$Site, contours$Tag, sep="-")
      
      write.csv(contours, file="ContourData.csv", row.names=F)

#### 2. Calculate taper parameter for all sampled trees #####
      
  # Find the name of all trees with cross-sectional measurements
    ContourTrees <- unique(contours$ID)
  
  # Define taper equation
    Eqn1 <- function(h,DBH,b1) {DBH*exp(-b1*(h-1.3))}
    
  # Define negative log likelihood function for optimization
    Lk1 <- function(par,h,d) {
      sigma <- par[1]
      DBH   <- par[2]
      b1    <- par[3]
      res   <- Eqn1(h,DBH,b1)-d
      neglk <- -sum(dnorm(res,sd = sigma,log=T))
      return(neglk)
    }
    
  # Define function to find the best-fit taper parameter by optimizing the likelihood function above  
    Fit.Eqn1= function(par,h,d) {
      return(optim(par,Lk1,h=h,d=d))
    }
  
  # Define columns to save results
    TreeSample$b1.hom <- NA # Taper parameter - measurements above HOM
    TreeSample$b1.iso <- NA # Taper parameter - meaurements above HOM or circular
    TreeSample$range.hom <- NA # Height range - measurements above HOM
    TreeSample$range.iso <- NA # Height range - meaurements above HOM or circular
    TreeSample$iso <- NA # Circularity at height of measurement
  
  # For each tree, loop through to calculate taper    
    for(i in 1:length(ContourTrees)){
      tree.i <- contours[contours$ID==ContourTrees[i],]
      info.i <- TreeSample[TreeSample$ID==ContourTrees[i],]
      
      # Find measurements up to 4 m above the HOM, and find the mean ciruclarity (iso.quo) above the HOM
      tree.hom <- tree.i[tree.i$ht > info.i$HOM & tree.i$ht <= info.i$HOM + 3.6,]
      hom.iso <- mean(tree.hom$iso.quo)
      tree.hom <- tree.hom[!is.na(tree.hom$ht),]
      
      # Find measurements above the HOM or with mean circularity of measurements above the HOM
      tree.iso <- tree.i[(tree.i$ht > info.i$HOM | tree.i$iso.quo >= hom.iso) & tree.i$ht <= info.i$HOM + 3.6,]
      tree.iso <- tree.iso[tree.iso$ht <= min(tree.iso$ht) + 3.6,]
      tree.iso <- tree.iso[!is.na(tree.iso$ht),]
      
      # Store ciruclarity at the height of measurement
      TreeSample[TreeSample$ID==ContourTrees[i],"iso"] <- tree.hom$iso.quo[1]
      
      # Fit taper functions
      results.hom <- ifelse(length(tree.hom$ht)==0,NA,
                            Fit.Eqn1(par=c(sigma=1,DBH=info.i$DBH/1000,b1=0.05),
                                     h=tree.hom$ht, d=tree.hom$d))
      results.iso <- ifelse(length(tree.iso$ht)==0,NA,
                            Fit.Eqn1(par=c(sigma=1,DBH=info.i$DBH/1000,b1=0.05),
                                     h=tree.iso$ht, d=tree.iso$d))
      
      # Save taper parameter results in data frame
      TreeSample[TreeSample$ID==ContourTrees[i],"b1.hom"] <- ifelse(is.na(results.hom), NA, results.hom[[1]][3])
      TreeSample[TreeSample$ID==ContourTrees[i],"b1.iso"] <- ifelse(is.na(results.iso), NA, results.iso[[1]][3])
      
      # Save range of heights used to fit taper parameter
      TreeSample[TreeSample$ID==ContourTrees[i],"range.hom"] <- ifelse(is.na(results.hom),NA,range(tree.hom$ht)[2]-range(tree.hom$ht)[1])
      TreeSample[TreeSample$ID==ContourTrees[i],"range.iso"] <- ifelse(is.na(results.iso),NA,range(tree.iso$ht)[2]-range(tree.iso$ht)[1])
    }
    
      # Convert DBH to cm
      TreeSample$DBH <- TreeSample$DBH/10
    
    # Save only measurements to use in a separate data frame called Taper Sample. These are trees for which
    # a taper parameter was estimated, and which had at least 1.4 m of trunk shape measurements.
    TaperSample <- TreeSample[!is.na(TreeSample$b1.iso),]
    TaperSample <- TaperSample[TaperSample$range.iso >= 1.4,]
    
    # Save only trees with circularity measurement
    CircSample <- TreeSample[!is.na(TreeSample$iso),]
    
   
    
#### 3. Compare models explaining variation in taper parameter ####
    
    # Models WITHOUT family random effect
    
    # Null model
    modelnull <- lme4::lmer(b1.iso~ (1|Site), data = TaperSample)
    # With one variable
    model1a <- lme4::lmer(b1.iso~log(DBH) + (1|Site), data = TaperSample)
    model1b <- lme4::lmer(b1.iso~log(HOM) + (1|Site), data = TaperSample)
    model1c <- lme4::lmer(b1.iso~log(WSG) + (1|Site), data = TaperSample)
    # With two variables
    model2a <- lme4::lmer(b1.iso~log(DBH) + log(HOM) + (1|Site), data = TaperSample)
    model2b <- lme4::lmer(b1.iso~log(DBH) + log(WSG) + (1|Site), data = TaperSample)
    model2c <- lme4::lmer(b1.iso~log(HOM) + log(WSG) + (1|Site), data = TaperSample)
    # With three variables
    model3a <- lme4::lmer(b1.iso~log(DBH) + log(HOM) + log(WSG) + (1|Site), data = TaperSample)
    

    # Make a table comparing models (for main text)    
    ModelCompare <- anova(modelnull, model1a, model1b, model1c, model2a, model2b, model2c, model3a)
    ModelCompare$R2marginal <- c(MuMIn::r.squaredGLMM(modelnull)[1], 
                                 MuMIn::r.squaredGLMM(model1a)[1], 
                                 MuMIn::r.squaredGLMM(model1b)[1], 
                                 MuMIn::r.squaredGLMM(model1c)[1], 
                                 MuMIn::r.squaredGLMM(model2a)[1], 
                                 MuMIn::r.squaredGLMM(model2b)[1], 
                                 MuMIn::r.squaredGLMM(model2c)[1],
                                 MuMIn::r.squaredGLMM(model3a)[1])
    ModelCompare$R2conditional <- c(MuMIn::r.squaredGLMM(modelnull)[2], 
                                    MuMIn::r.squaredGLMM(model1a)[2], 
                                    MuMIn::r.squaredGLMM(model1b)[2], 
                                    MuMIn::r.squaredGLMM(model1c)[2], 
                                    MuMIn::r.squaredGLMM(model2a)[2], 
                                    MuMIn::r.squaredGLMM(model2b)[2], 
                                    MuMIn::r.squaredGLMM(model2c)[2],
                                    MuMIn::r.squaredGLMM(model3a)[2])
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
    
    
#### 4. Save files ####    
    write.csv(ModelResults, file="Table3_TaperMixedModelComparison.csv",
              row.names = F)
    write.csv(TaperSample, file="DataFile_TaperParameterSample.csv",
              row.names = F)
    write.csv(TreeSample, file="DataFile_AllTaperEstimates.csv",
              row.names = F)
    write.csv(CircSample, file="DataFile_CircSample.csv",
              row.names = F)

  # Verify that residuals of best model don't vary by site using a 1-way ANOVA
    ModelResiduals <- data.frame(residual = residuals(model3a),
                                 site = TaperSample[TaperSample$b1.iso>0,"Site"])
    residualANOVA <- lm(residual~site, data=ModelResiduals)
    summary(residualANOVA)

#### Revision edit: Make table with taper by site ####
  
  # Calculate EDBH, basal area, and AGB for each tree
    agb.allometry <- function(E,wsg,dbh) {exp(-1.803-0.976*E+0.976*log(wsg)
                                              + 2.673*log(dbh) - 0.0299*(log(dbh)^2))}
    
    taper.eqn <- function(d,h,b1) {d/(exp(-b1*(h-1.3)))}
    
    TaperSample$EDBH <- taper.eqn(d = TaperSample$DBH,
                                  h = TaperSample$HOM,
                                  b1 = TaperSample$b1.iso)
    
    TaperSample$BA <- pi*(TaperSample$EDBH/2)^2
    
    TaperSample$E <- NA
      TaperSample[TaperSample$Site=="AMA","E"] <- -0.07928769
      TaperSample[TaperSample$Site=="BCI","E"] <- 0.04944549
      TaperSample[TaperSample$Site=="BKT","E"] <- -0.05956875
      TaperSample[TaperSample$Site=="HKK","E"] <- 0.3194663
      TaperSample[TaperSample$Site=="KCH","E"] <- 0.04786947
      
    TaperSample$AGB <- agb.allometry(E = TaperSample$E,
                                     wsg = TaperSample$WSG,
                                     dbh = TaperSample$DBH)
    
    TaperSample$Vol <- TaperSample$AGB/TaperSample$WSG
    
    TaperBySite <- data.frame(Site = c("AMA","BCI","BKT","HKK","KCH"),
                              n = NA,
                              Mean = NA,
                              SD = NA,
                              Mean.BA = NA,
                              Mean.AGB = NA,
                              Mean.Vol = NA)
    
    for(i in 1:length(TaperBySite$Site)){
      TaperBySite[i,"n"] <- length(TaperSample[TaperSample$Site==TaperBySite$Site[i],"b1.iso"])
      TaperBySite[i,"Mean"] <- mean(TaperSample[TaperSample$Site==TaperBySite$Site[i],"b1.iso"])
      TaperBySite[i,"SD"] <- sd(TaperSample[TaperSample$Site==TaperBySite$Site[i],"b1.iso"])
      TaperBySite[i,"Mean.BA"] <- weighted.mean(x = TaperSample[TaperSample$Site==TaperBySite$Site[i],"b1.iso"],
                                                w = TaperSample[TaperSample$Site==TaperBySite$Site[i],"BA"])
      TaperBySite[i,"Mean.AGB"] <- weighted.mean(x = TaperSample[TaperSample$Site==TaperBySite$Site[i],"b1.iso"],
                                                 w = TaperSample[TaperSample$Site==TaperBySite$Site[i],"AGB"])
      TaperBySite[i,"Mean.Vol"] <- weighted.mean(x = TaperSample[TaperSample$Site==TaperBySite$Site[i],"b1.iso"],
                                                 w = TaperSample[TaperSample$Site==TaperBySite$Site[i],"Vol"])
    }
    
    write.csv(TaperBySite, file="TaperVariationTable.csv")
      
#### Revision edit: Make table with circularity by site ####    
# Calculate EDBH, basal area, and AGB for each tree
    agb.allometry <- function(E,wsg,dbh) {exp(-1.803-0.976*E+0.976*log(wsg)
                                              + 2.673*log(dbh) - 0.0299*(log(dbh)^2))}
    
    taper.eqn <- function(d,h,b1) {d/(exp(-b1*(h-1.3)))}
    
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
    
    CircSample$Vol <- CircSample$AGB/CircSample$WSG
    
    CircBySite <- data.frame(Site = c("AMA","BCI","BKT","HKK","KCH"),
                             n = NA,
                              Mean = NA,
                              SD = NA,
                              Mean.BA = NA,
                              Mean.AGB = NA,
                              Mean.Vol = NA)
    
    for(i in 1:length(CircBySite$Site)){
      CircBySite[i,"n"] <- length(CircSample[CircSample$Site==CircBySite$Site[i],"iso"])
      CircBySite[i,"Mean"] <- mean(CircSample[CircSample$Site==CircBySite$Site[i],"iso"])
      CircBySite[i,"SD"] <- sd(CircSample[CircSample$Site==CircBySite$Site[i],"iso"])
      CircBySite[i,"Mean.BA"] <- weighted.mean(x = CircSample[CircSample$Site==CircBySite$Site[i],"iso"],
                                                w = CircSample[CircSample$Site==CircBySite$Site[i],"BA"])
      CircBySite[i,"Mean.AGB"] <- weighted.mean(x = CircSample[CircSample$Site==CircBySite$Site[i],"iso"],
                                                 w = CircSample[CircSample$Site==CircBySite$Site[i],"AGB"])
      CircBySite[i,"Mean.Vol"] <- weighted.mean(x = CircSample[CircSample$Site==CircBySite$Site[i],"iso"],
                                                 w = CircSample[CircSample$Site==CircBySite$Site[i],"Vol"])
    }
    
    write.csv(CircBySite, file="CircVariationTable.csv")
#### Revision edit: Explore alternate taper parameter analysis options ####
    
     # Krusal-Wallis test
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
      
      # one-way ANOVA for family
        familyAnova <- aov(b1.iso~Family, data=TaperSample)
        summary(familyAnova)
        
      # additional variation explained by family
        model3a_fam <- lme4::lmer(b1.iso~log(DBH) + log(HOM) + log(WSG) + (1|Site) + (1|Family), data = TaperSample)
        MuMIn::r.squaredGLMM(model3a_fam)
        
      # Look at distribution of values-- neither normal nor lognormal
        # Taper parameter
        shapiro.test(TaperSample$b1.iso)
        shapiro.test(log(TaperSample$b1.iso))
        shapiro.test(sqrt(TaperSample$b1.iso))
      
        hist(TaperSample$b1.iso,
             breaks = seq(-0.1,0.2,0.01),
             col="black",border="white", xlab = "Taper parameter")
        hist(log(TaperSample$b1.iso),
             breaks = seq(-10,1,0.4),
             col="black",border="white", xlab = "Taper parameter")
        hist(sqrt(TaperSample$b1.iso),
             col="black",border="white", xlab = "Taper parameter")
        
        # DAB -- is lognormal
        shapiro.test(TaperSample$DBH)
        shapiro.test(log(TaperSample$DBH)) 
        
        # HOM -- neither normal nor lognormal
        shapiro.test(TaperSample$HOM)
        shapiro.test(log(TaperSample$HOM))
        
        hist(TaperSample$HOM,
             breaks = seq(1.3,7.9,0.2),
             col="black",border="white", xlab = "HOM")
        
        # WSG -- neither normal nor lognormal
        shapiro.test(TaperSample$WSG)
        shapiro.test(log(TaperSample$WSG))
        
        hist(TaperSample$WSG,
             breaks = seq(0.2,1,0.03),
             col="black",border="white", xlab = "WSG")
        


    
#### Revision edit: Estimate AGB with various models ####
        
  # Define necessary functions
    agb.allometry <- function(E,wsg,dbh) {exp(-1.803-0.976*E+0.976*log(wsg)
                                              + 2.673*log(dbh) - 0.0299*(log(dbh)^2))}
    
    taper.eqn <- function(d,h,b1) {d/(exp(-b1*(h-1.3)))}
    
  # Add E for each site
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
      
    # AGB without site and wood density differences to better sort by size
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
          model.full <- lme4::lmer(b1.iso~log(DBH) + log(HOM) + log(WSG) + (1|Site), data = data.full)
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
          model.3pr <- lme4::lmer(b1.iso~log(DBH) + log(HOM) + log(WSG) + (1|Site), data = data.full)
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
          model.2pr <- lme4::lmer(b1.iso~log(DBH) + log(HOM) + (1|Site), data = data.full)
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
          model.full <- lme4::lmer(b1.iso~log(DBH) + log(HOM) + log(WSG) + (1|Site), data = data.full)
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
          model.3pr <- lme4::lmer(b1.iso~log(DBH) + log(HOM) + log(WSG) + (1|Site), data = data.full)
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
          model.2pr <- lme4::lmer(b1.iso~log(DBH) + log(HOM) + (1|Site), data = data.full)
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
          model.full <- lme4::lmer(b1.iso~log(DBH) + log(HOM) + log(WSG) + (1|Site), data = data.full)
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
          model.3pr <- lme4::lmer(b1.iso~log(DBH) + log(HOM) + log(WSG) + (1|Site), data = data.full)
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
          model.2pr <- lme4::lmer(b1.iso~log(DBH) + log(HOM) + (1|Site), data = data.full)
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
          model.full <- lme4::lmer(b1.iso~log(DBH) + log(HOM) + log(WSG) + (1|Site), data = data.full)
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
          model.3pr <- lme4::lmer(b1.iso~log(DBH) + log(HOM) + log(WSG) + (1|Site), data = data.full)
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
          model.2pr <- lme4::lmer(b1.iso~log(DBH) + log(HOM) + (1|Site), data = data.full)
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
          model.full <- lme4::lmer(b1.iso~log(DBH) + log(HOM) + log(WSG) + (1|Site), data = data.full)
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
          model.3pr <- lme4::lmer(b1.iso~log(DBH) + log(HOM) + log(WSG) + (1|Site), data = data.full)
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
          model.2pr <- lme4::lmer(b1.iso~log(DBH) + log(HOM) + (1|Site), data = data.full)
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
        
  # Summarize results
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
      mean(100*abs(AGB_Results$AGB_Uncorrected-AGB_Results$AGB_MeasTaper)/AGB_Results$AGB_Uncorrected)
      min(100*abs(AGB_Results$AGB_Uncorrected-AGB_Results$AGB_MeasTaper)/AGB_Results$AGB_Uncorrected)
      max(100*abs(AGB_Results$AGB_Uncorrected-AGB_Results$AGB_MeasTaper)/AGB_Results$AGB_Uncorrected)

      # Mean, min, max from applying Model 1
      mean(100*abs(AGB_Results$AGB_FullModel_Mean-AGB_Results$AGB_MeasTaper)/AGB_Results$AGB_MeasTaper)
      min(100*abs(AGB_Results$AGB_FullModel_Mean-AGB_Results$AGB_MeasTaper)/AGB_Results$AGB_MeasTaper)
      max(100*abs(AGB_Results$AGB_FullModel_Mean-AGB_Results$AGB_MeasTaper)/AGB_Results$AGB_MeasTaper)
      
      # Mean, min, max from applying cross-validated model with 3 parameters
      mean(100*abs(AGB_Results$AGB_Cross3Pr_Mean-AGB_Results$AGB_MeasTaper)/AGB_Results$AGB_MeasTaper)
      min(100*abs(AGB_Results$AGB_Cross3Pr_Mean-AGB_Results$AGB_MeasTaper)/AGB_Results$AGB_MeasTaper)
      max(100*abs(AGB_Results$AGB_Cross3Pr_Mean-AGB_Results$AGB_MeasTaper)/AGB_Results$AGB_MeasTaper)

      # Mean, min, max from applying cross-validated model with 2 parameters
      mean(100*abs(AGB_Results$AGB_Cross2Pr_Mean-AGB_Results$AGB_MeasTaper)/AGB_Results$AGB_MeasTaper)
      min(100*abs(AGB_Results$AGB_Cross2Pr_Mean-AGB_Results$AGB_MeasTaper)/AGB_Results$AGB_MeasTaper)
      max(100*abs(AGB_Results$AGB_Cross2Pr_Mean-AGB_Results$AGB_MeasTaper)/AGB_Results$AGB_MeasTaper)

      # Mean, min, max from applying cross-validated model with 2 parameters
      mean(100*abs(AGB_Results$AGB_CrossMed_Mean-AGB_Results$AGB_MeasTaper)/AGB_Results$AGB_MeasTaper)
      min(100*abs(AGB_Results$AGB_CrossMed_Mean-AGB_Results$AGB_MeasTaper)/AGB_Results$AGB_MeasTaper)
      max(100*abs(AGB_Results$AGB_CrossMed_Mean-AGB_Results$AGB_MeasTaper)/AGB_Results$AGB_MeasTaper)

      
      write.csv(AGB_Results, file="DataFile_BiomassEstimates.csv", row.names = F)