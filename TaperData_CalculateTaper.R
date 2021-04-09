# Use 'groundhog' package to use versions of packages at time of publication
  groundhog::groundhog.library("splancs","2020-12-20")

# Define site abbreviations and names
  sites <- c("AMA","BCI","BKT","HKK","KCH")
  sitesNames <- c("Amacayacu",
                  "Barro Colorado",
                  "Bukit Timah",
                  "Huai Kha Khaeng",
                  "Khao Chong")

# Load data about sampled trees
  TreeSample <- read.csv("DataFiles/DataFile_SampledTreesInfo.csv")
  TreeSample$ID <- paste(TreeSample$Site, TreeSample$Tag, sep="-")

#### 1. Calculate area, diameter, and circularity of all-cross sections ####  
  
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
      AMA.trees <- list.files("DataFiles/contours_AMA/",
                              full.names = F, no..=F)
      AMA.files <- list.files("DataFiles/contours_AMA/",
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
              tree.data$area[j]=areapl(cbind(vertices[,1],vertices[,2]))
            # Calculate perimeter of cross-section  
              tree.data$perim[j] = calc.perim(vertices)	
            # Find convex hull
              tree.chull <- chull(vertices[,1:2])
              # Calculate perimeter of convex hull
                perim.chull <- calc.perim(vertices[tree.chull,])
              # Estimate diameter from convex hull perimeter
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
      BCI.trees <- list.files("DataFiles/contours_BCI/",
                              full.names = F, no..=F)
      BCI.files <- list.files("DataFiles/contours_BCI/",
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
          tree.data$area[j]=areapl(cbind(vertices[,1],vertices[,2]))
          # Calculate perimeter of cross-section  
          tree.data$perim[j] = calc.perim(vertices)	
          # Find convex hull
          tree.chull <- chull(vertices[,1:2])
          # Calculate perimeter of convex hull
          perim.chull <- calc.perim(vertices[tree.chull,])
          # Estimate diameter from convex hull perimeter
          tree.data$d[j] <- perim.chull/pi
        }
        
        # For all cross-sections, calculate circularity (iso.quo)
        tree.data$iso.quo = (4*pi*tree.data$area)/((tree.data$perim)^2)
        
        # Save info
        contours <- rbind(contours,tree.data)
      }    
  
  
#### 1.3 Bukit Timah ####

      # List all trees for BKT
      BKT.trees <- list.files("DataFiles/contours_BKT/",
                              full.names = F, no..=F)
      BKT.files <- list.files("DataFiles/contours_BKT/",
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
          tree.data$area[j]=areapl(cbind(vertices[,1],vertices[,2]))
          # Calculate perimeter of cross-section  
          tree.data$perim[j] = calc.perim(vertices)	
          # Find convex hull
          tree.chull <- chull(vertices[,1:2])
          # Calculate perimeter of convex hull
          perim.chull <- calc.perim(vertices[tree.chull,])
          # Estimate diameter from convex hull perimeter
          tree.data$d[j] <- perim.chull/pi
        }
        
        # For all cross-sections, calculate circularity (iso.quo)
        tree.data$iso.quo = (4*pi*tree.data$area)/((tree.data$perim)^2)
        
        # Save info
        contours <- rbind(contours,tree.data)
      }    
      
      
#### 1.4 Huai Kha Khaeng ####

      # List all trees for HKK
      HKK.trees <- list.files("DataFiles/contours_HKK/",
                              full.names = F, no..=F)
      HKK.files <- list.files("DataFiles/contours_HKK/",
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
          tree.data$area[j]=areapl(cbind(vertices[,1],vertices[,2]))
          # Calculate perimeter of cross-section  
          tree.data$perim[j] = calc.perim(vertices)	
          # Find convex hull
          tree.chull <- chull(vertices[,1:2])
          # Calculate perimeter of convex hull
          perim.chull <- calc.perim(vertices[tree.chull,])
          # Estimate diameter from convex hull perimeter
          tree.data$d[j] <- perim.chull/pi
        }
        
        # For all cross-sections, calculate circularity (iso.quo)
        tree.data$iso.quo = (4*pi*tree.data$area)/((tree.data$perim)^2)
        
        # Save info
        contours <- rbind(contours,tree.data)
      }   

#### 1.5 Khao Chong ####

      # List all trees for KCH
      KCH.trees <- list.files("DataFiles/contours_KCH/",
                              full.names = F, no..=F)
      KCH.files <- list.files("DataFiles/contours_KCH/",
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
          tree.data$area[j]=areapl(cbind(vertices[,1],vertices[,2]))
          # Calculate perimeter of cross-section  
          tree.data$perim[j] = calc.perim(vertices)	
          # Find convex hull
          tree.chull <- chull(vertices[,1:2])
          # Calculate perimeter of convex hull
          perim.chull <- calc.perim(vertices[tree.chull,])
          # Estimate diameter from convex hull perimeter
          tree.data$d[j] <- perim.chull/pi
        }
        
        # For all cross-sections, calculate circularity (iso.quo)
        tree.data$iso.quo = (4*pi*tree.data$area)/((tree.data$perim)^2)
        
        # Save info
        contours <- rbind(contours,tree.data)
      }  

      contours$ID <- paste(contours$Site, contours$Tag, sep="-")

#### 2. Compare candidate taper models #####
      
  # Find the name of all trees with cross-sectional measurements
    ContourTrees <- unique(contours$ID)
      
  ## Define different taper models
  
  eqn.names <- c("Metcalf et al. (2009)",  # 1
				"Kozak et al. (1969)",  # 2
				"Ormerod (1973)",  # 3
				"Forslund (1990)",  # 4
				"Riemer et al. (1995)"  # 5
				)  # References in literature for each equation
				
    # For each model, Eqn defines the taper equation as proposed in the literature
    #  Args:
    #    h: vector of measurement heights (m) (Equations 1-5)
    #    H: measured total tree height (m) (Equations 2-5)
    #    DBH: parameter for diameter at 1.3 m height (cm) (Equations 1-5)
    #    b1: taper parameter (Equations 1-5)
    #    b2: taper parameter (Equations 4,5)
    #    b3: taper parameter (Equation 5)
    #  Returns: 
    #    A vector of calculated diameters for each corresponding measurment height. 
    #
	# For each model, Lk defines the negative log likelihood function for a set of 
	#   measurement heights and diameters given a set of parameter values.
	#  Args:
	#    par: vector of parameter values for each model (Equations 1-5)
	#      sigma: standard deviation between fitted values and expected values
	#      DBH: parameter for diameter at 1.3 m height (cm) (Equations 1-5)
    #      b1: taper parameter (Equations 1-5)
    #      b2: taper parameter (Equations 4,5)
    #      b3: taper parameter (Equation 5)
    #    h: vector of measurement heights (m) (Equations 1-5)
    #    d: vector of measured diameters (cm) (Equations 1-5)
    #    H: measured total tree height (m) (Equations 2-5)
    #  Returns:
    #    The negative log likelihood for the set of measurement heights and diameters given 
    #    the set of parameter values.
    #
	# For each model, Fit.Eqn finds the best-fit paramater values for a set of measurement
	# heights and diamteres by minimizing negative log likelihood
	#  Args:
	#    par: vector of parameter values for each model (Equations 1-5)
	#      sigma: standard deviation between fitted values and expected values
	#      DBH: parameter for diameter at 1.3 m height (cm) (Equations 1-5)
    #      b1: taper parameter (Equations 1-5)
    #      b2: taper parameter (Equations 4,5)
    #      b3: taper parameter (Equation 5)
    #    h: vector of measurement heights (m) (Equations 1-5)
    #    d: vector of measured diameters (cm) (Equations 1-5)
    #    H: measured total tree height (m) (Equations 2-5)
    #  Returns:
    #    The optimization results from minimizing the negative log likelihood function, 
    #    including parameter values ($par) and the negative log likelihood of the resulting 
    #    fit ($value).
    
    					
    #1. Metcalf et al. (2009)
     Eqn1 <- function(h,DBH,b1) {DBH*exp(-b1*(h-1.3))}
     Lk1 <- function(par,h,d) {
       sigma <- par[1]
       DBH   <- par[2]
       b1    <- par[3]
       res   <- Eqn1(h,DBH,b1)-d
       neglk <- -sum(dnorm(res,sd = sigma,log=T))
       return(neglk)
     }
     Fit.Eqn1= function(par,h,d) {
       return(optim(par,Lk1,h=h,d=d))
	 }

    #2. Kozak et al. (1969)
     Eqn2 <- function(h,H,DBH,b1) {DBH*sqrt(b1*(1-2*(h/H)+(h/H)^2))}
     Lk2 <- function(par,h,d,H) {
       sigma <- par[1]
       DBH   <- par[2]
       b1    <- par[3]
       res   <- Eqn2(h,H,DBH,b1)-d
       neglk <- -sum(dnorm(res,sd = sigma,log=T))
       return(neglk)
     }
     Fit.Eqn2= function(par,h,d,H) {
       return(optim(par,Lk2,h=h,d=d,H=H))
     }
		
    #3. Ormerod (1973)
    Eqn3 <- function(h,H,DBH,b1) {DBH*(((H-h)/(H-1.3))^b1)}
    Lk3 <- function(par,h,d,H) {
      sigma	<- par[1]
      DBH   <- par[2]
      b1    <- par[3]
      res <- Eqn3(h,H,DBH,b1)-d
      neglk <- -sum(dnorm(res,sd = sigma,log=T))
      return(neglk)
    }
    Fit.Eqn3= function(par,h,d,H) {
      return(optim(par,Lk3,h=h,d=d,H=H))
    }
		
    #4. Forslund (1990)
    Eqn4 <- function(h,H,DBH,b1,b2) {DBH*(1-(h/H)^b1)^(1/b2)}
    Lk4 <- function(par,h,d,H) {
      sigma	<- par[1]
      DBH   <- par[2]
      b1    <- par[3]
      b2    <- par[4]
      res   <- Eqn4(h,H,DBH,b1,b2)-d
      neglk <- -sum(dnorm(res,sd = sigma,log=T))
      return(neglk)
    }
    Fit.Eqn4= function(par,h,d,H) {
      return(optim(par,Lk4,h=h,d=d,H=H))
    }
	
    #5. Riemer et al. (1995)
    Eqn5 <- function(h,H,DBH,b1,b2,b3) {b1*DBH/(1-exp(b3*(1.3-H)))+(DBH/2-b1*DBH)*(1-1/(1-exp(b2*(1.3-H))))+exp(-b2*h)*(exp(1.3*b2)*(DBH/2-b1*DBH)/(1-exp(b2*(1.3-H))))-exp(b3*h)*(b1*DBH*exp(-b3*H)/(1-exp(b3*(1.3-H))))}
    Lk5 <- function(par,h,d,H) {
      sigma	<- par[1]
      DBH   <- par[2]
      b1    <- par[3]
      b2    <- par[4]
      b3    <- par[5]
      res   <- Eqn5(h,H,DBH,b1,b2,b3)-d
      neglk <- -sum(dnorm(res,sd = sigma,log=T))
      return(neglk)
    }
    Fit.Eqn5= function(par,h,d,H) {
      return(optim(par,Lk5,h=h,d=d,H=H))
    }
    
          # Compare the goodness of fit for each taper model using AICc values
        # CalcAICc calculates the AICc values for the fit of a taper model to a set of 
        # diameter measurements.
	    #  Args:
	    #    fitresult: The return argument from fitting taper models as defined above (i.e. 
	    #    the output from Fit.Eqn1, Fit.Eqn2, Fit.Eqn3, Fit.Eqn4, or Fit.Eqn5)
        #    d: vector of measured diameters (cm)
        #  Returns:
        #    The AICc value (numeric)
        CalcAICc <- function(fitresult,d) {
          AIC <- 2*fitresult$value+2*length(fitresult$par)
          AICc <- AIC+(2*length(fitresult$par)*(length(fitresult$par)+1))/(length(d)-length(fitresult$par)-1)
          return(AICc)
        }

  
  # Define columns to save results
    TreeSample$AIC1 <- NA
    TreeSample$AIC2 <- NA
    TreeSample$AIC3 <- NA
    TreeSample$AIC4 <- NA
    TreeSample$AIC5 <- NA
  
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
      
      # Proceed if there are remaining measurements
      if(!length(tree.iso$ht)==0 & !is.na(info.i$Ht)){
      
        #  Fit each proposed taper function
        # MAKE SURE all diameters are in cm
        results.eqn1 <- Fit.Eqn1(
  	    par=c(sigma=1,DBH=info.i$DBH/10,b1=0.05),
  	    h=tree.iso$ht,
  	    d=tree.iso$d*100)
        results.eqn2 <- Fit.Eqn2(
  	    par=c(sigma=1,DBH=info.i$DBH/10,b1=1),
  	    h=tree.iso$ht,
  	    d=tree.iso$d*100,
  	    H=info.i$Ht)
        results.eqn3 <- Fit.Eqn3(
  	    par=c(sigma=1,DBH=info.i$DBH/10,b1=1),
  	    h=tree.iso$ht,
  	    d=tree.iso$d*100,
  	    H=info.i$Ht)	  
        results.eqn4 <- Fit.Eqn4(
  	    par=c(sigma=1,DBH=info.i$DBH/10,b1=1,b2=0.1),
  	    h=tree.iso$ht,
  	    d=tree.iso$d*100,
  	    H=info.i$Ht)
        results.eqn5 <- Fit.Eqn5(
  	    par=c(sigma=1,DBH=info.i$DBH/10,b1=-1,b2=1,b3=-1),
  	    h=tree.iso$ht,
  	    d=tree.iso$d*100,
  	    H=info.i$Ht)
      
        # Calculate AIC
        AICc=c(CalcAICc(results.eqn1,tree.iso$d),
      	    CalcAICc(results.eqn2,tree.iso$d),
      	    CalcAICc(results.eqn3,tree.iso$d),
      	    CalcAICc(results.eqn4,tree.iso$d),
      	    CalcAICc(results.eqn5,tree.iso$d))
        
        # Calculate whether within 2 of lowest AIC
        AICc <- AICc-min(AICc)
        AICbest <- AICc<=2
        
        # Save results
        TreeSample[TreeSample$ID==ContourTrees[i],"AIC1"] <- AICbest[1]
        TreeSample[TreeSample$ID==ContourTrees[i],"AIC2"] <- AICbest[2]
        TreeSample[TreeSample$ID==ContourTrees[i],"AIC3"] <- AICbest[3]
        TreeSample[TreeSample$ID==ContourTrees[i],"AIC4"] <- AICbest[4]
        TreeSample[TreeSample$ID==ContourTrees[i],"AIC5"] <- AICbest[5]

        # Save range of heights used to fit taper parameter
        TreeSample[TreeSample$ID==ContourTrees[i],"range.iso"] <- range(tree.iso$ht)[2]-range(tree.iso$ht)[1]
      }
    }
    
    # Save only measurements to use in a separate data frame called Taper Sample. These are trees for which
    # a taper parameter was estimated, and which had at least 1.4 m of trunk shape measurements.
    TaperSample <- TreeSample[TreeSample$range.iso >= 1.4,]
    TaperSample <- TaperSample[!is.na(TaperSample$Site),]
    
    # Summarize by site
    bestModels <- data.frame(Site=sites,
                             Model1=NA,
                             Model2=NA,
                             Model3=NA,
                             Model4=NA,
                             Model5=NA)
    for(i in 1:length(sites)){
      bestModels$Model1[i] <- dim(TaperSample[TaperSample$AIC1==T & TaperSample$Site==sites[i],])[1]
      bestModels$Model2[i] <- dim(TaperSample[TaperSample$AIC2==T & TaperSample$Site==sites[i],])[1]
      bestModels$Model3[i] <- dim(TaperSample[TaperSample$AIC3==T & TaperSample$Site==sites[i],])[1]
      bestModels$Model4[i] <- dim(TaperSample[TaperSample$AIC4==T & TaperSample$Site==sites[i],])[1]
      bestModels$Model5[i] <- dim(TaperSample[TaperSample$AIC5==T & TaperSample$Site==sites[i],])[1]
    }
    
    bestModels$Site <- sitesNames
    write.csv(bestModels, file = "ResultsFiles/TableS5_TaperModelComparison.csv", row.names = F)
   
#### 3. Calculate taper parameter (using best model) for all sampled trees #####
    
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
      
      # Find measurements up to 3.6 m above the HOM, and find the mean ciruclarity (iso.quo) above the HOM
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
    
#### 4. Save files ####    
    
    write.csv(TaperSample, file="DataFiles/DataFile_TaperParameterSample.csv",
              row.names = F)
    write.csv(TreeSample, file="DataFiles/DataFile_AllTaperEstimates.csv",
              row.names = F)
    write.csv(CircSample, file="DataFiles/DataFile_CircSample.csv",
              row.names = F)