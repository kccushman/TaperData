## NOTE: if working from the current Git repository, then workflow will function from step 2.

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
  
    # Save only measurements to use in a separate data frame called Taper Sample. These are trees for which
    # a taper parameter was estimated, and which had at least 1.4 m of trunk shape measurements.
    TaperSample <- TreeSample[!is.na(TreeSample$b1.iso),]
    TaperSample <- TaperSample[TaperSample$range.iso >= 1.4,]
    
    # Convert DBH to cm
    TaperSample$DBH <- TaperSample$DBH/10
    
#### 3. Compare models explaining variation in taper parameter ####
    # Null model
    modelnull <- lmer(log(b1.iso)~ (1|Site) + (1|Family), data = TaperSample)
    # With one variable
    model1a <- lmer(log(b1.iso)~log(DBH) + (1|Site) + (1|Family), data = TaperSample)
    model1b <- lmer(log(b1.iso)~log(HOM) + (1|Site) + (1|Family), data = TaperSample)
    model1c <- lmer(log(b1.iso)~log(WSG) + (1|Site) + (1|Family), data = TaperSample)
    # With two variables
    model2a <- lmer(log(b1.iso)~log(DBH) + log(HOM) + (1|Site) + (1|Family), data = TaperSample)
    model2b <- lmer(log(b1.iso)~log(DBH) + log(WSG) + (1|Site) + (1|Family), data = TaperSample)
    model2c <- lmer(log(b1.iso)~log(HOM) + log(WSG) + (1|Site) + (1|Family), data = TaperSample)
    # With three variables
    model3a <- lmer(log(b1.iso)~log(DBH) + log(HOM) + log(WSG) + (1|Site) + (1|Family), data = TaperSample)
    
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
    ModelCompare$Description <- c(paste("Log(b1) = ",round(summary(modelnull)$coefficients[1],3)," + Site + Family", sep=""),
                                  paste("Log(b1) = ",round(summary(model1a)$coefficients[1,1],3)," + ",round(summary(model1a)$coefficients[2,1],3),"*log(DAB) + Site + Family", sep=""),
                                  paste("Log(b1) = ",round(summary(model1b)$coefficients[1,1],3)," + ",round(summary(model1b)$coefficients[2,1],3),"*log(HOM) + Site + Family", sep=""),
                                  paste("Log(b1) = ",round(summary(model1c)$coefficients[1,1],3)," + ",round(summary(model1c)$coefficients[2,1],3),"*log(WSG) + Site + Family", sep=""),
                                  paste("Log(b1) = ",round(summary(model2a)$coefficients[1,1],3)," + ",round(summary(model2a)$coefficients[2,1],3),"*log(DAB) ",round(summary(model2a)$coefficients[3,1],3),"*log(HOM) ", "+ Site + Family", sep=""),
                                  paste("Log(b1) = ",round(summary(model2b)$coefficients[1,1],3)," + ",round(summary(model2b)$coefficients[2,1],3),"*log(DAB) ",round(summary(model2b)$coefficients[3,1],3),"*log(WSG) ", "+ Site + Family", sep=""),
                                  paste("Log(b1) = ",round(summary(model2c)$coefficients[1,1],3)," + ",round(summary(model2c)$coefficients[2,1],3),"*log(HOM) ",round(summary(model2c)$coefficients[3,1],3),"*log(WSG) ", "+ Site + Family", sep=""),
                                  paste("Log(b1) = ",round(summary(model3a)$coefficients[1,1],3)," + ",round(summary(model3a)$coefficients[2,1],3),"*log(DAB) ",round(summary(model3a)$coefficients[3,1],3),"*log(HOM) ",round(summary(model3a)$coefficients[4,1],3),"*log(WSG) ", "+ Site + Family", sep=""))
    
    ModelResults <- ModelCompare[,c("Description","Df","AIC","R2marginal","R2conditional")]
    ModelResults <- ModelResults[order(ModelResults$AIC),]
    ModelResults$dAIC <- ModelResults$AIC-min(ModelResults$AIC)
    write.csv(ModelResults, file="Table 2_Taper mixed model comparison.csv",
              row.names = F)
    write.csv(TaperSample, file="Data file_taper parameter sample.csv",
              row.names = F)

  # Verify that residuals of best model don't vary by site using a 1-way ANOVA
    ModelResiduals <- data.frame(residual = residuals(model3a),
                                 site = TaperSample[TaperSample$b1.iso>0,"Site"])
    residualANOVA <- lm(residual~site, data=ModelResiduals)
    summary(residualANOVA)