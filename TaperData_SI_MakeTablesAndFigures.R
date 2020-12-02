###### Load data ######
  # Tree sample and taper parameter information for trees with 3-D models
    TaperSample <- read.csv("DataFile_TaperParameterSample.csv")
    contours <- read.csv("ContourData.csv")
    TreeSample <- read.csv("DataFile_AllTaperEstimates.csv")
  # Data for trees measured with optical dendrometer and DBH tape
    dbh.tape.data <- read.csv("~/Desktop/Taper/Current/TaperCorrection/DBH.tape.data.csv")
    optical.data <- read.csv("~/Desktop/Taper/Current/TaperCorrection/Optical.dendro.data.csv")
  # Plot data
    load("ForestGEO_CensusData.RData")

# Define plot names and colors
  sites <- c("AMA","BCI","BKT","HKK","KCH")
  
  sitesNames <- c("Amacayacu",
                  "Barro Colorado",
                  "Bukit Timah",
                  "Huai Kha Khaeng",
                  "Khao Chong")
  
  site.cols <- data.frame(site=sites,
                          col= c("#efa7a7", # salmon
                                 "#000000", # black
                                 "#158baf", # bright blue"
                                 "#bf0f0f", # cardinal
                                 "#3ead82")) # teal
  
  site.cols$col <- as.character(site.cols$col)
        
###### Table S5:Random intercept values for family and site #####
  
  # Best taper models from model comparison 
  
  # Model 1
    model1 <- lme4::lmer(b1.iso~log(DBH) + log(HOM) + log(WSG) + (1|Site) , data = TaperSample)
  # Model 2
    model2 <- lme4::lmer(b1.iso~log(DBH) + log(HOM) + (1|Site), data = TaperSample)
  # Model 3
    model3 <- lme4::lmer(b1.iso~log(DBH) + log(WSG) + (1|Site), data = TaperSample)
  # Model 4
    model4 <- lme4::lmer(b1.iso~log(DBH) + (1|Site), data = TaperSample)
  # Model 5
    model5 <- lme4::lmer(b1.iso~log(HOM) + log(WSG) + (1|Site), data = TaperSample)
  # Model 6
    model6 <- lme4::lmer(b1.iso~log(HOM) + (1|Site), data = TaperSample)
  # Model 7
    model7 <- lme4::lmer(b1.iso~log(WSG) + (1|Site), data = TaperSample)
  # Model 8
    model8 <- lme4::lmer(b1.iso~ (1|Site), data = TaperSample)
  
    SiteCoefs1 <- data.frame(Model = 1,
                             Group=c(rownames(coef(model1)[[1]])),
                             Intercept=c(coef(model1)[[1]][,1]))
    SiteCoefs2 <- data.frame(Model = 2,
                             Group=c(rownames(coef(model2)[[1]])),
                             Intercept=c(coef(model2)[[1]][,1]))
    SiteCoefs3 <- data.frame(Model = 3,
                             Group=c(rownames(coef(model3)[[1]])),
                             Intercept=c(coef(model3)[[1]][,1]))
    SiteCoefs4 <- data.frame(Model = 4,
                             Group=c(rownames(coef(model4)[[1]])),
                             Intercept=c(coef(model4)[[1]][,1]))
    SiteCoefs5 <- data.frame(Model = 5,
                             Group=c(rownames(coef(model5)[[1]])),
                             Intercept=c(coef(model5)[[1]][,1]))
    SiteCoefs6 <- data.frame(Model = 6,
                             Group=c(rownames(coef(model6)[[1]])),
                             Intercept=c(coef(model6)[[1]][,1]))
    SiteCoefs7 <- data.frame(Model = 7,
                             Group=c(rownames(coef(model7)[[1]])),
                             Intercept=c(coef(model7)[[1]][,1]))
    SiteCoefs8 <- data.frame(Model = 8,
                             Group=c(rownames(coef(model8)[[1]])),
                             Intercept=c(coef(model8)[[1]][,1]))
    SiteCoefs <- rbind(SiteCoefs1,SiteCoefs2,SiteCoefs3,SiteCoefs4,SiteCoefs5,SiteCoefs6,SiteCoefs7,SiteCoefs8)
    
    write.csv(SiteCoefs, row.names = F, file="TableS5_RandomIntercepts.csv")

###### Table S6:Taper values by family #####

  # Are there significant differences among families?
    summary(aov(b1.iso~Family, TaperSample))
    
  # Aggregate - mean Taperularity by family
  FamilyTaper <- aggregate(TaperSample$b1.iso, by=list(TaperSample$Family), FUN="mean")
  # Aggregate - SD Taperularity by family  
  FamilyTaperSD <- aggregate(TaperSample$b1.iso, by=list(TaperSample$Family), FUN="sd")
  # Aggregate - total sample size by family  
  FamilyTaperCount <- aggregate(TaperSample$b1.iso, by=list(TaperSample$Family), FUN="length")
  
  # Put results into data frame
  TaperTable <- data.frame(Family=FamilyTaper$Group.1, 
                          Taperularity=paste(round(FamilyTaper$x,2)," (",round(FamilyTaperSD$x,2),")", sep=""),
                          SampleSize=FamilyTaperCount$x,
                          SampleAMA=NA,
                          SampleBCI=NA,
                          SampleBKT=NA,
                          SampleHKK=NA,
                          SampleKCH=NA)
  # Sort by sample size
  TaperTable <- TaperTable[order(-TaperTable$SampleSize),]
  
  # Find the sample size per plot
  for(i in 1:length(TaperTable$Family)){
    TaperTable[i,"SampleAMA"] <- length(TaperSample[TaperSample$Family==TaperTable$Family[i] & TaperSample$Site=="AMA","Tag"])
    TaperTable[i,"SampleBCI"] <- length(TaperSample[TaperSample$Family==TaperTable$Family[i] & TaperSample$Site=="BCI","Tag"])
    TaperTable[i,"SampleBKT"] <- length(TaperSample[TaperSample$Family==TaperTable$Family[i] & TaperSample$Site=="BKT","Tag"])
    TaperTable[i,"SampleHKK"] <- length(TaperSample[TaperSample$Family==TaperTable$Family[i] & TaperSample$Site=="HKK","Tag"])
    TaperTable[i,"SampleKCH"] <- length(TaperSample[TaperSample$Family==TaperTable$Family[i] & TaperSample$Site=="KCH","Tag"])
  }
  
  write.csv(TaperTable, file="TableS6_TaperByFamily.csv", row.names = F)  
  
###### Table S7:Circularity values by family #####
  # Make a data frame with trees with circularity measurement at the HOM
  CircSample <- TreeSample[!is.na(TreeSample$iso),]
    
  # Are there significant differences among families?
    summary(aov(iso~Family, CircSample))
    
  # Aggregate - mean circularity by family
  FamilyCirc <- aggregate(CircSample$iso, by=list(CircSample$Family), FUN="mean")
  # Aggregate - SD circularity by family  
  FamilyCircSD <- aggregate(CircSample$iso, by=list(CircSample$Family), FUN="sd")
  # Aggregate - total sample size by family  
  FamilyCircCount <- aggregate(CircSample$iso, by=list(CircSample$Family), FUN="length")
  
  # Put results into data frame
  CircTable <- data.frame(Family=FamilyCirc$Group.1, 
                          Circularity=paste(round(FamilyCirc$x,2)," (",round(FamilyCircSD$x,2),")", sep=""),
                          SampleSize=FamilyCircCount$x,
                          SampleAMA=NA,
                          SampleBCI=NA,
                          SampleBKT=NA,
                          SampleHKK=NA,
                          SampleKCH=NA)
  # Sort by sample size
  CircTable <- CircTable[order(-CircTable$SampleSize),]
  
  # Find the sample size per plot
  for(i in 1:length(CircTable$Family)){
    CircTable[i,"SampleAMA"] <- length(CircSample[CircSample$Family==CircTable$Family[i] & CircSample$Site=="AMA","Tag"])
    CircTable[i,"SampleBCI"] <- length(CircSample[CircSample$Family==CircTable$Family[i] & CircSample$Site=="BCI","Tag"])
    CircTable[i,"SampleBKT"] <- length(CircSample[CircSample$Family==CircTable$Family[i] & CircSample$Site=="BKT","Tag"])
    CircTable[i,"SampleHKK"] <- length(CircSample[CircSample$Family==CircTable$Family[i] & CircSample$Site=="HKK","Tag"])
    CircTable[i,"SampleKCH"] <- length(CircSample[CircSample$Family==CircTable$Family[i] & CircSample$Site=="KCH","Tag"])
  }
  
  write.csv(CircTable, file="TableS7_Circularity by family.csv", row.names = F)  
  
###### Table S8: HOM tests by site ######
    HOM.results <- read.csv("Data file_HOMresultsPerPlot.csv")
    load("ForestGEO_CensusData.RData")
    
      # Order census data by StemID and then by HOM (decreasing)
      for(i in 1:length(AMA.cens)){
        AMA.cens[[i]] <- AMA.cens[[i]][order(AMA.cens[[i]]$StemID,-AMA.cens[[i]]$hom),]
      }
      for(i in 1:length(BCI.cens)){
        BCI.cens[[i]] <- BCI.cens[[i]][order(BCI.cens[[i]]$StemID,-BCI.cens[[i]]$hom),]
      }
      for(i in 1:length(BKT.cens)){
        BKT.cens[[i]] <- BKT.cens[[i]][order(BKT.cens[[i]]$StemID,-BKT.cens[[i]]$hom),]
      }
      for(i in 1:length(HKK.cens)){
        HKK.cens[[i]] <- HKK.cens[[i]][order(HKK.cens[[i]]$StemID,-HKK.cens[[i]]$hom),]
      }
      for(i in 1:length(KCH.cens)){
        KCH.cens[[i]] <- KCH.cens[[i]][order(KCH.cens[[i]]$StemID,-KCH.cens[[i]]$hom),]
      }
    # Test for differnces over time within plots

# Test for differnces over time within plots

  ## AMA
  #Kruskal-Wallis test
  AMA.HOM.test <- kruskal.test(list(AMA.cens[[1]][!duplicated(AMA.cens[[1]]$StemID),"hom"],AMA.cens[[2]][!duplicated(AMA.cens[[2]]$StemID),"hom"]))
 
  ## BCI
  #Kruskal-Wallis test
  BCI.HOM.test <- kruskal.test(list(BCI.cens[[2]][!duplicated(BCI.cens[[2]]$StemID),"hom"],
                                    BCI.cens[[3]][!duplicated(BCI.cens[[3]]$StemID),"hom"],BCI.cens[[4]][!duplicated(BCI.cens[[4]]$StemID),"hom"],
                                    BCI.cens[[5]][!duplicated(BCI.cens[[5]]$StemID),"hom"],BCI.cens[[6]][!duplicated(BCI.cens[[6]]$StemID),"hom"]))
  
  ## HKK
  #Kruskal-Wallis test
  HKK.HOM.test <- kruskal.test(list(HKK.cens[[1]][!duplicated(HKK.cens[[1]]$StemID),"hom"],HKK.cens[[2]][!duplicated(HKK.cens[[2]]$StemID),"hom"],
                                    HKK.cens[[3]][!duplicated(HKK.cens[[3]]$StemID),"hom"],HKK.cens[[4]][!duplicated(HKK.cens[[4]]$StemID),"hom"]))
  
  ## KCH
  #Kruskal-Wallis test
  KCH.HOM.test <- kruskal.test(list(KCH.cens[[1]][!duplicated(KCH.cens[[1]]$StemID),"hom"],KCH.cens[[2]][!duplicated(KCH.cens[[2]]$StemID),"hom"],
                                    KCH.cens[[3]][!duplicated(KCH.cens[[3]]$StemID),"hom"]))
  
    HomTestTable <- data.frame(Site=c("Amacayacu", "Barro Colorado", "Huai Kha Khaeng", "Khao Chong"),
                             MinHOM=c(min(HOM.results[HOM.results$Site=="AMA",c("MeanHOM")]),
                                      min(HOM.results[HOM.results$Site=="BCI",c("MeanHOM")]),
                                      min(HOM.results[HOM.results$Site=="HKK",c("MeanHOM")]),
                                      min(HOM.results[HOM.results$Site=="KCH",c("MeanHOM")])),
                             MaxHOM=c(max(HOM.results[HOM.results$Site=="AMA",c("MeanHOM")]),
                                      max(HOM.results[HOM.results$Site=="BCI",c("MeanHOM")]),
                                      max(HOM.results[HOM.results$Site=="HKK",c("MeanHOM")]),
                                      max(HOM.results[HOM.results$Site=="KCH",c("MeanHOM")])),
                             Hval=c(AMA.HOM.test$statistic, BCI.HOM.test$statistic,HKK.HOM.test$statistic,
                                    KCH.HOM.test$statistic),
                             Pval=c(AMA.HOM.test$p.value, BCI.HOM.test$p.value,HKK.HOM.test$p.value,
                                    KCH.HOM.test$p.value))
    write.csv(HomTestTable, file="TableS8_MeasHeightsWithinPlots.csv", row.names = F)

###### Table S9: HOM variation by family #####
  RecentHOM <- rbind(AMA.cens[[2]][!duplicated(AMA.cens[[2]]$StemID),c("Family", "Site","hom","dbh")], BCI.cens[[6]][!duplicated(BCI.cens[[6]]$StemID),c("Family", "Site","hom","dbh")], 
                     BKT.cens[[1]][!duplicated(BKT.cens[[1]]$StemID),c("Family", "Site","hom","dbh")], HKK.cens[[4]][!duplicated(HKK.cens[[4]]$StemID),c("Family", "Site","hom","dbh")],
                     KCH.cens[[3]][!duplicated(KCH.cens[[3]]$StemID),c("Family", "Site","hom","dbh")])
  
  RecentHOM[is.na(RecentHOM$Family),"Family"] <- "Unknown"
  
  summary(aov(hom~Family, data = RecentHOM))
  
  FamilyHOM <- aggregate(RecentHOM$hom, by=list(RecentHOM$Family), FUN="mean")
  FamilyHOMSD <- aggregate(RecentHOM$hom, by=list(RecentHOM$Family), FUN="sd")
  FamilyHOMCount <- aggregate(RecentHOM$hom, by=list(RecentHOM$Family), FUN="length")
  
  # Put results into data frame
  HOMTable <- data.frame(Family=FamilyHOM$Group.1, 
                         HOM=paste(round(FamilyHOM$x,2)," (",round(FamilyHOMSD$x,2),")", sep=""),
                         SampleSize=FamilyHOMCount$x,
                         SampleAMA=NA,
                         SampleBCI=NA,
                         SampleBKT=NA,
                         SampleHKK=NA,
                         SampleKCH=NA)
  
  # Sort by sample size
  HOMTable <- HOMTable[order(-HOMTable$SampleSize),]
  
  # Find the sample size per plot
  for(i in 1:length(HOMTable$Family)){
    HOMTable[i,"SampleAMA"] <- length(RecentHOM[RecentHOM$Family==HOMTable$Family[i] & RecentHOM$Site=="AMA","dbh"])
    HOMTable[i,"SampleBCI"] <- length(RecentHOM[RecentHOM$Family==HOMTable$Family[i] & RecentHOM$Site=="BCI","dbh"])
    HOMTable[i,"SampleBKT"] <- length(RecentHOM[RecentHOM$Family==HOMTable$Family[i] & RecentHOM$Site=="BKT","dbh"])
    HOMTable[i,"SampleHKK"] <- length(RecentHOM[RecentHOM$Family==HOMTable$Family[i] & RecentHOM$Site=="HKK","dbh"])
    HOMTable[i,"SampleKCH"] <- length(RecentHOM[RecentHOM$Family==HOMTable$Family[i] & RecentHOM$Site=="KCH","dbh"])
  }
  
  # Keep families with at least 10 individuals
  HOMTable10 <- HOMTable[HOMTable$SampleSize>=10,]
  
  # Write table to .csv
  write.csv(HOMTable10, file="TableS9_MeasHtByFamily.csv", row.names = F)

###### Figure S1: Taper parameter variation with height #####
  range.lm <- lm(b1.iso~range.iso, data=TaperSample)
  summary(range.lm)
  newx <- seq(min(TaperSample$range.iso), max(TaperSample$range.iso), length.out=100)
  preds <- predict(range.lm, newdata = data.frame(range.iso=newx), 
                   interval = 'confidence')
  
  tiff(width=5, height=5, file="Figure S1_Taper vs measurement height range.tiff",res=300,units="in")
    par(family="serif")
    plot(b1.iso~range.iso, data=TaperSample,
         pch=20,
         xlab=NA,
         ylab=NA)
    mtext("Range of measurement heights (m)", side=1, line=2)
    mtext("Taper parameter", side=2, line=2)
    abline(range.lm, lwd=2)
    lines(newx, preds[ ,3], lty = 'dashed', col = 'red', lwd=2)
    lines(newx, preds[ ,2], lty = 'dashed', col = 'red', lwd=2)
    legend(x=2.5, y=-0.02,
           c("Regression line", "95% CI"),
           lwd=2,
           col=c("black","red"),
           lty=c(1,2),
           bty='n')
  dev.off()
  
  
###### Figure S2: Taper parameter variation including measurements below the HOM ##### 
  IsoCompare <- TaperSample[TaperSample$range.hom>= 1.4 & !(TaperSample$range.hom==TaperSample$range.iso),]
  # Use a Deming regression
    IsoReg=mcreg(IsoCompare$b1.iso,IsoCompare$b1.hom,method.reg="Deming",error.ratio=1,method.ci="analytical")
    getCoefficients(IsoReg)
    
    newx <- seq(min(IsoCompare$b1.iso), max(IsoCompare$b1.iso), length.out=100)
    preds <- MCResultAnalytical.calcResponse(IsoReg, x.levels = newx, alpha=0.05)
    
  tiff(width=5, height=5, file="Figure S2_Taper variation below HOM.tiff",res=300,units="in")
    par(family="serif")
    plot(b1.hom~b1.iso, data=IsoCompare,
         pch=20,
         xlab=NA,
         ylab=NA)
    abline(a=getCoefficients(IsoReg)[1,1],b=getCoefficients(IsoReg)[2,1],lwd=2)
    lines(x=preds[,1],y=preds[,4],lwd=2,lty=2)
    lines(x=preds[,1],y=preds[,5],lwd=2,lty=2)
    abline(a=0,b=1,col="red",lwd=2,lty=2)
    mtext("Taper parameter - with measurements below HOM", side=1, line=2)
    mtext("Taper parameter - without measurements below HOM", side=2, line=2)
    legend(x=0.06, y=0.01,
           c("Regression line", "95% CI","1-1 line"),
           lwd=2,
           col=c("black","black","red"),
           lty=c(1,2,2),
           bty='n')
  dev.off()

###### Figure S3: Optical dendrometer vs. 3-D model histograms #####
  
  # Calculate taper values for trees using optical dendrometer data
  TaperSample$b1.opt <- NA
  
  # define taper equation
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
    
  optical.trees <- unique(as.character(optical.data$tag[optical.data$tag %in% TaperSample$Tag]))
  for(i in 1:length(optical.trees)){
    #select 3-D data for each tree and re-create subset used to fit taper function
    three.d.i <- contours[contours$Tag==optical.trees[i] & contours$Site=="BCI",]
    info.i <- TreeSample[TreeSample$Tag==optical.trees[i] & TreeSample$Site=="BCI",]
    tree.hom <- three.d.i[three.d.i$ht > (info.i$HOM[1]),]
    hom.iso <-  mean(tree.hom$iso.quo)
    tree.iso <-three.d.i[(three.d.i$ht > info.i$HOM | three.d.i$iso.quo >= hom.iso) & three.d.i$ht <= info.i$HOM + 3.6,]
    
    #select optical data for the same range of heights
    opt.i <- optical.data[optical.data$tag==optical.trees[i] & optical.data$ht >= range(tree.iso$ht)[1] & optical.data$ht <= range(tree.iso$ht)[2],]
    
    results.opt <- ifelse(length(opt.i$ht)==0,NA,
                          Fit.Eqn1(
                            par=c(sigma=1,DBH=info.i$DBH[1],b1=0.05),
                            h=opt.i$ht, d=opt.i$diam))
    
    TaperSample[TaperSample$Tag==optical.trees[i],"b1.opt"] <- ifelse(is.na(results.opt), NA,
                                                                      results.opt[[1]][3])
  }
  
  tiff(width=7, height=4, file="Figure S3_Histograms of taper methods.tiff",res=300,units="in")
    par(mfrow=c(1,2), family="serif")
    b1.breaks <- seq(0,0.12,0.02)
    
    histA <- hist(TaperSample[TaperSample$b1.opt>=0 & TaperSample$b1.iso>=0,'b1.opt'],
                  breaks=b1.breaks,
                  xlab = NA,
                  main = 'Optical dendrometer',
                  col="black",border="white")
    densityA <- density(TaperSample[TaperSample$b1.opt>=0 & TaperSample$b1.iso>=0 & !is.na(TaperSample$b1.opt),'b1.opt'])
    lines(x=densityA$x,y=densityA$y/(max(densityA$y))*max(histA$counts),lwd=2, col="red")
    
    histB <- hist(TaperSample[TaperSample$b1.opt>=0 & TaperSample$b1.iso>=0,'b1.iso'],
                  breaks=b1.breaks,
                  xlab = NA,
                  main = '3-D model',
                  col="black",border="white")
    densityB <- density(TaperSample[TaperSample$b1.opt>=0 & TaperSample$b1.iso>=0 & !is.na(TaperSample$b1.opt),'b1.iso'])
    lines(x=densityB$x,y=densityB$y/(max(densityB$y))*max(histB$counts),lwd=2, col="red")
    
    mtext('Taper parameter value', side=1, outer=T, line=-2)
  dev.off()
  
  tapermethod.ttest <- t.test(TaperSample[TaperSample$b1.opt>=0 & TaperSample$b1.iso>=0,'b1.opt'],
                              TaperSample[TaperSample$b1.opt>=0 & TaperSample$b1.iso>=0,'b1.iso'],
                              paired=T)
  
###### Figure S4: Optical dendrometer taper vs. 3-D model taper #####
  
  TaperCompare <- TaperSample[TaperSample$b1.opt>=0 & TaperSample$b1.iso>=0 & !is.na(TaperSample$b1.opt),]
  # Use a Deming regression
  OptReg=mcr::mcreg(x=TaperCompare$b1.iso,y=TaperCompare$b1.opt,method.reg="Deming",error.ratio=1,method.ci="analytical")
  mcr::getCoefficients(OptReg)
  
  newx <- seq(min(TaperCompare$b1.iso), max(TaperCompare$b1.iso), length.out=100)
  preds <- mcr::MCResultAnalytical.calcResponse(OptReg, x.levels = newx, alpha=0.05)
  
  tiff(width=5, height=5, file="FigureS4_Optical dendrometer vs 3-D model.tiff",res=300,units="in")
    par(family="sans")
    plot(b1.opt~b1.iso, data=TaperCompare,
         ylim=c(-0.4,0.4),
         #ylim=range(TaperCompare$b1.opt)+c(-0.015,0),
         pch=20,
         xlab=NA,
         ylab=NA)
     abline(a=mcr::getCoefficients(OptReg)[1,1],b=mcr::getCoefficients(OptReg)[2,1],lwd=2)
     lines(x=preds[,1],y=preds[,4],lwd=2,lty=2)
     lines(x=preds[,1],y=preds[,5],lwd=2,lty=2)
    abline(a=0,b=1,col="red",lwd=2,lty=2)
    mtext("3-D model taper parameter", side=1, line=2)
    mtext("Optical dendrometer taper parameter", side=2, line=2)
     legend(x=0.025, y=-0.1,
            c("Regression line","95% CI","1-1 line"),
            lwd=2,
            col=c("black","black","red"),
            lty=c(1,2,2),
            bty='n')
  dev.off()
  
###### Figure S5: Individual trees, optical dendromter vs. tape vs. 3-D model #####
  dbh.tape.trees <- unique(dbh.tape.data$tag)
  dbh.tape.trees <- dbh.tape.trees[dbh.tape.trees %in% contours[contours$Site=="BCI",'Tag']]
  
  tiff(width=8, height=9, file="Figure S5_Optical dendrometer vs 3-D model individuals.tiff",res=300,units="in")
    par(mfrow=c(3,3),family="serif",mar=c(3,3,1,1),oma=c(2,2,1,1))
    
    for(i in 1:length(dbh.tape.trees)){
      dbh.tape.i <- dbh.tape.data[dbh.tape.data$tag==dbh.tape.trees[i],]
      three.d.i <- contours[contours$Site=="BCI" & contours$Tag==dbh.tape.trees[i],]
        three.d.i$d <- three.d.i$d*100 # convert from meters to cm
      optical.i <- optical.data[optical.data$tag==dbh.tape.trees[i],]
      
      xlims <- range(c(dbh.tape.i$taped,three.d.i$d.area,optical.i$diam),na.rm=T)
      ylims <- range(c(dbh.tape.i$ht,three.d.i$ht,optical.i$ht),na.rm=T)
      
      plot(ht~taped, data=dbh.tape.i, type="n",
           xlab=NA,
           ylab=NA,
           xlim=xlims,ylim=ylims,
           cex.axis=1.5)
      points(ht~taped, data=dbh.tape.i, pch=3, col=adjustcolor(alpha.f=0.8, "black"))
      points(ht~diam, data=optical.i, pch=4, col=adjustcolor(alpha.f=0.8, "red"))
      points(ht~d, data=three.d.i, pch=19, col=adjustcolor(alpha.f=0.6, "blue"))
    }
    par(xpd=T)
    plot(x=0,y=0,type="n",bty="n",
         xlim=c(0,10),ylim=c(0,10),
         xaxt="n",yaxt="n",xlab="",ylab="")
    legend(x=0,y=10,
           c('3-D Model','Diameter tape','Optical dendrometer'),
           pch=c(19,3,4),cex=1.5,
           col=c(adjustcolor(alpha.f=0.6, "blue"),adjustcolor(alpha.f=0.8, "black"),adjustcolor(alpha.f=0.8, "red")), bty='n')
    
    mtext("Diameter measurement (cm)", side=1, line=0.5, outer=T, cex=1.5)
    mtext("Measurement height (m)", side=2, line=0.5, outer=T, cex=1.5)
    
  dev.off()
  
###### Figure S6: Distribution of model residuals #####
  
  # Read in data file for taper parameter sample and define the taper model used in the biomass estimation routine
  model3a <- lme4::lmer(b1.iso~log(DBH) + log(HOM) + log(WSG) + (1|Site), data = TaperSample)

  tiff(width=6, height=4, file="FigureS6_DistributionOfModelResiduals.tiff",res=300,units="in")
  par(mfrow=c(2,2),mar=c(4,3,1,1),oma = c(0,3,0,0), family="sans")
  
  # NO family fixed effects
  plot(x=fitted(model3a),
       y=residuals(model3a),
       xlab=NA, ylab=NA,
       pch=20, col=adjustcolor("black",alpha.f=0.5),
       cex=1.2, cex.axis=1.5)
  mtext("Fitted value", side=1, line=2)
  mtext("Residual value", side=2, outer=T)
  abline(h=0,lty=2)
  text("a", x=0.005, y=0.08, cex=1.4)
  
  plot(x=model3a@frame$`log(DBH)`,
       y=residuals(model3a),
       xlab=NA, ylab= NA,
       pch=20, col=adjustcolor("black",alpha.f=0.5),
       cex=1.2, cex.axis=1.5)
  mtext("log(DAB)", side=1, line=2)
  text("b", x=3, y=0.08, cex=1.4)
  abline(h=0,lty=2)
  
  plot(x=model3a@frame$`log(HOM)`,
       y=residuals(model3a),
       xlab=NA, ylab= NA,
       pch=20, col=adjustcolor("black",alpha.f=0.5),
       cex=1.2, cex.axis=1.5)
  mtext("log(HOM)", side=1, line=2)
  text("c", x=0.4, y=0.08, cex=1.4)
  abline(h=0,lty=2)
  
  plot(x=model3a@frame$`log(WSG)`,
       y=residuals(model3a),
       xlab=NA, ylab= NA,
       pch=20, col=adjustcolor("black",alpha.f=0.5),
       cex=1.2, cex.axis=1.5)
  mtext("log(WSG)", side=1, line=2)
  text("d", x=-1.195, y=0.08, cex=1.4)
  abline(h=0,lty=2)
  
  dev.off()

###### Figure S7: Nonstandard measurement heights for trees >= 30 cm ######
  HOM.results <- read.csv("DataFile_HOMresultsPerPlot30.csv")
  
  axisSize <- 0.9
  cexAx <- 1.1
  textCex <- 1
  ptSize <- 0.8
  pchSig <- c(19,19,19,19,19)
  
  tiff(file="FigureS7_MeasHtChangesOverTime30.tiff",width=3,height=6, units="in", res=300, family = "sans")

    par(mfrow=c(3,1), mar=c(1,6,1,1),oma=c(2,1,0,0), family="sans", xpd=F, las=1)
    
    plot(x=HOM.results$Year, y=HOM.results$Prop*100,
         type='n',
         ylab=NA,ylim=c(0,100),
         xlab=NA,
         xaxt = "n",
         cex.axis=cexAx)
    mtext("% of \nstems", side=2, line=6.8, cex=axisSize, adj=0)
    text(x=1990.5,y=97,"a", cex=textCex+0.3)
    
    for(i in 1:length(sites)){
        points(x=HOM.results[HOM.results$Site==sites[i],"Year"],
               y=HOM.results[HOM.results$Site==sites[i],"Prop"]*100, 
               pch=pchSig[i], col=site.cols[site.cols$site==sites[i],"col"], cex=ptSize)
        lines(x=HOM.results[HOM.results$Site==sites[i],"Year"],
              y=HOM.results[HOM.results$Site==sites[i],"Prop"]*100, col=site.cols[site.cols$site==sites[i],"col"], lty=1)
    }
    
    legend(x=2004, y=105, bty='n',
       sitesNames,
       col=site.cols$col,
       cex=cexAx-0.3,
       pch=pchSig)
    
    plot(x=HOM.results$Year, y=HOM.results$PropBA*100,
         type='n',
         ylab=NA,ylim=c(0,100),
         xlab=NA,
         xaxt = "n",
         cex.axis=cexAx)
    mtext("% of \nbasal \narea", side=2, line=6.8, cex=axisSize, adj=0)
    text(x=1990.5,y=97,"b", cex=textCex+0.3)
    for(i in 1:length(sites)){
        points(x=HOM.results[HOM.results$Site==sites[i],"Year"],
               y=HOM.results[HOM.results$Site==sites[i],"PropBA"]*100, 
               pch=pchSig[i], col=site.cols[site.cols$site==sites[i],"col"], cex=ptSize)
        lines(x=HOM.results[HOM.results$Site==sites[i],"Year"],
              y=HOM.results[HOM.results$Site==sites[i],"PropBA"]*100, col=site.cols[site.cols$site==sites[i],"col"], lty=1)
    }
    
    plot(x=HOM.results$Year,y=HOM.results$MeanHOM,
         type='n',
         ylab=NA,
         xlab=NA,
         cex.axis=cexAx)
    mtext("Census year", side=1, line=2, cex=axisSize)
    mtext("Average \nHOM (m)", side=2, line=6.8, cex=axisSize, adj=0)
    text(x=1990.5,y=3.55,"c", cex=textCex+0.3)
    for(i in 1:length(sites)){
        points(x=HOM.results[HOM.results$Site==sites[i],"Year"],
               y=HOM.results[HOM.results$Site==sites[i],"MeanHOM"], 
               pch=pchSig[i], col=site.cols[site.cols$site==sites[i],"col"], cex=ptSize)
        lines(x=HOM.results[HOM.results$Site==sites[i],"Year"],
              y=HOM.results[HOM.results$Site==sites[i],"MeanHOM"], col=site.cols[site.cols$site==sites[i],"col"], lty=1)
    }
  dev.off()
  
###### Figure S8: Trunk circularity for each tree #####
  contours$Tag <- as.character(contours$Tag) 
    
  tiff(width=8, height=6, file="Figure S8_Trunk circularity for individual trees.tiff",res=300,units="in")
    par(mfrow=c(2,3), family="serif", mar=c(3,3,1,1),oma=c(2,2,1,1))
    circSites <- unique(CircSample$Site)
    for(i in 1:length(circSites)){
      siteTrees <- unique(CircSample[CircSample$Site==circSites[i],"Tag"])
      colorsi <- rainbow(length(siteTrees))
      plot(0,type='n',
           xlim=c(0.01,1),
           ylim=c(-10,10),
           xlab=NA,ylab=NA,
           cex.axis=1.5)
      mtext(sitesNames[i],side=3, cex=1.2)
      for(j in 1:length(siteTrees)){
        treej <- contours[contours$Site==circSites[i] & contours$Tag==siteTrees[j],]
        infoj <- TreeSample[TreeSample$Site==circSites[i] & TreeSample$Tag==siteTrees[j],]
        # Calculate height above HOM for each measurements
          treej$aboveHOM <- treej$ht-infoj$HOM
        # Add circularity by height above HOM to plot
          lines(aboveHOM~iso.quo, data=treej, col=colorsi[j])
      }
      abline(h=0, col="red", lty=2)
    }
    mtext("Circularity", side=1,outer=T, cex=1.5)
    mtext("Distance from measurement height (m)", side=2,outer=T, cex=1.5)
  dev.off()

  
###### Figure S9: HOM variation versus climate and latitude ######
  
    # Vector of MAP in site alphabetical order
    MAP <- c(3200,2600,2500,1500,2800)
  # Vector of dry season length in site alphabetical order
    DSL <- c(0,3,0,5,2.5)
  # Vector of distance from equator
    Lat <- c(3.8, 9.2, 1.3, 15.4, 7.5)
    
  PropStems <- c(0.5512590, 0.5677175, 0.2314316, 0.1071939, 0.6266549)
  PropBA <- c(0.6104313, 0.7050572, 0.2970426, 0.1605881,  0.7186443)
  MeanHOM <- c(2.554515, 3.658446, 1.664472, 1.541162, 3.012815)
  
  cexpt=2
  cexAx = 1.4
  
tiff(width=7, height=6.5, file="FigureS9_MeasHtVersusClimate.tiff",res=300,units="in")
    par(mfrow=c(3,3), family="sans", mar=c(2,2,1,1),oma=c(2,2,1,1))
    
    
  # Proportion of basal area  
  plot(PropBA~MAP, pch =20,
       xlab = NA,
       xaxt="n",
       ylab = NA,
       col=site.cols$col, cex=cexpt, cex.axis = cexAx)
  mtext("Proportion of basal area", side=2, line=2, cex=0.9)
  text("a", x = 1500, y = 0.7, cex = cexAx)
  
      legend(x=1450, y=0.68, bty='n',
       sitesNames,
       col=site.cols$col,
       cex=1.2,
       pch=19)
      
      
  plot(PropBA~DSL, pch =20,
       col=site.cols$col,
       xlab = NA,
       xaxt="n",
       yaxt="n",
       ylab=NA, cex=cexpt, cex.axis = cexAx)
   text("b", x = 0, y = 0.7, cex = cexAx)
   
  plot(PropBA~Lat, pch =20,
       col=site.cols$col,
       xlab = NA,
       xaxt="n",
       yaxt="n",
       ylab=NA, cex=cexpt, cex.axis = cexAx)
   text("c", x = 1.5, y = 0.7, cex = cexAx)

    # Proportion of stems 
  plot(PropStems~MAP, pch =20,
       xlab = NA,
       xaxt="n",
       ylab = NA,
       col=site.cols$col, cex=cexpt, cex.axis = cexAx)
  mtext("Proportion of stems", side=2, line=2, cex=0.9)
  text("d", x = 1500, y = 0.6, cex = cexAx)
   
  plot(PropStems~DSL, pch =20,
       col=site.cols$col,
       xlab = NA,
       xaxt="n",
       yaxt="n",
       ylab=NA, cex=cexpt, cex.axis = cexAx)
  text("e", x = 0, y = 0.6, cex = cexAx) 
  
  plot(PropStems~Lat, pch =20,
       col=site.cols$col,
       xlab = NA,
       xaxt="n",
       yaxt="n",
       ylab=NA, cex=cexpt, cex.axis = cexAx)
  text("f", x = 1.5, y = 0.6, cex = cexAx) 

    # Mean HOM 
  plot(MeanHOM~MAP, pch =20,
       xlab = NA,
       ylab = NA,
       col=site.cols$col, cex=cexpt, cex.axis = cexAx)
  text("g", x = 1500, y = 3.5, cex = cexAx) 
  
  mtext("Mean annual precipitation (mm)", side=1, line=2, cex=0.9)
  mtext("Mean HOM (m)", side=2, line=2, cex=0.9)

  plot(MeanHOM~DSL, pch =20,
       col=site.cols$col,
       xlab = NA,
       yaxt="n",
       ylab=NA, cex=cexpt, cex.axis = cexAx)
  text("h", x = 0, y = 3.5, cex = cexAx)
  mtext("Dry season length (months)", side=1, line=2, cex=0.9)
  plot(MeanHOM~Lat, pch =20,
       col=site.cols$col,
       xlab = NA,
       yaxt="n",
       ylab=NA, cex=cexpt, cex.axis = cexAx)
  text("i", x = 1.5, y = 3.5, cex = cexAx)
  mtext("Latitude (degrees from equator)", side=1, line=2, cex=0.9)
dev.off()