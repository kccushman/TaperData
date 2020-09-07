# Make figures for taper data paper

#### Make vectors of site abbreviates, names, and colors ####

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

#### Figure 2 : Variation of taper parameter with tree qualities ####
  TaperSample <- read.csv("DataFile_TaperParameterSample.csv")

  cexAx <- 1
  axisSize <- 0.8
  textCex <- 1.2
  
  tiff(file="Figure2_TaperParameterVariation.tiff", height=2.3, width=6, units="in", res=300)
    par(mfrow=c(1,3), mar=c(4,2,0,0), xpd=F,family="sans", las=1, oma=c(0,5,1,1))

    plot(log(b1.iso)~log(DBH), data=TaperSample, type='n',
         ylim=range(log(TaperSample$b1.iso),na.rm=T)+c(0,1),
         xlab=NA, ylab=NA,
         cex.axis=cexAx)
    mtext("log DAB (cm)", side=1, line=2.2, cex=axisSize)
    mtext("log taper \nparameter", side=2, line=2, cex=axisSize)
    text(x=3,y=-1.1,"a", cex=textCex+0.1)
    for(i in 1:length(sites)){
      points(log(b1.iso)~log(DBH), data=TaperSample[TaperSample$Site==sites[i],],
             col=site.cols[site.cols$site==sites[i],"col"], pch=20)
    }
    taperDABlm <- lme4::lmer(log(b1.iso)~log(DBH) + (1|Site), data = TaperSample)
    abline(a=summary(taperDABlm)$coefficients[1,1], b=summary(taperDABlm)$coefficients[2,1])
    
        legend(x=2.8, y=-5.7,bty='n',
           sitesNames,
           col=site.cols$col,
           cex=cexAx-0.25,
           pch=19)
        
    plot(log(b1.iso)~log(HOM), data=TaperSample, type='n',
         ylim=range(log(TaperSample$b1.iso),na.rm=T)+c(0,1),
         xlab=NA,ylab=NA,
         yaxt='n',
         cex.axis=cexAx)
  mtext("log HOM (cm)", side=1, line=2.2, cex=axisSize)
  text(x=0.45,y=-1.1,"b", cex=textCex+0.1)
    for(i in  1:length(sites)){
      points(log(b1.iso)~log(HOM), data=TaperSample[TaperSample$Site==sites[i],],
             col=site.cols[site.cols$site==sites[i],"col"], pch=20)
    }
    taperHOMlm <- lme4::lmer(log(b1.iso)~log(HOM) + (1|Site), data = TaperSample)
    abline(a=summary(taperHOMlm)$coefficients[1,1], b=summary(taperHOMlm)$coefficients[2,1])
    
    plot(log(b1.iso)~log(WSG), data=TaperSample, type='n',
         ylim=range(log(TaperSample$b1.iso),na.rm=T)+c(0,1),
         xlab=NA, ylab=NA,yaxt='n',
         cex.axis=cexAx)
    mtext("log WSG (cm)", side=1, line=2.2, cex=axisSize)
    text(x=-1.25,y=-1.1,"c", cex=textCex+0.1)
    for(i in 1:length(sites)){
      points(log(b1.iso)~log(WSG), data=TaperSample[TaperSample$Site==sites[i],],
             col=site.cols[site.cols$site==sites[i],"col"], pch=20)
    }
    taperWSGlm <- lme4::lmer(log(b1.iso)~log(WSG) + (1|Site), data = TaperSample)
    abline(a=summary(taperWSGlm)$coefficients[1,1], b=summary(taperWSGlm)$coefficients[2,1])

  dev.off()  


#### Figure 3: Trunk circularity for each site #####
  # Calculate empirical cumulative distribution functions for each plot
    CircSample <- TreeSample[!is.na(TreeSample$iso),]

    AMA.circCDF <- Hmisc::Ecdf(CircSample[CircSample$Site=='AMA','iso'])
    BCI.circCDF <- Hmisc::Ecdf(CircSample[CircSample$Site=='BCI','iso'])
    BKT.circCDF <- Hmisc::Ecdf(CircSample[CircSample$Site=='BKT','iso'])
    HKK.circCDF <- Hmisc::Ecdf(CircSample[CircSample$Site=='HKK','iso'])
    KCH.circCDF <- Hmisc::Ecdf(CircSample[CircSample$Site=='KCH','iso'])
    
    # Fifth percentile
    AMA.circCDF[[1]][AMA.circCDF[[2]]>0.05][1]
    BCI.circCDF[[1]][BCI.circCDF[[2]]>0.05][1]
    BKT.circCDF[[1]][BKT.circCDF[[2]]>0.05][1]
    HKK.circCDF[[1]][HKK.circCDF[[2]]>0.05][1]
    KCH.circCDF[[1]][KCH.circCDF[[2]]>0.05][1]
    
    # Median values
    AMA.circCDF[[1]][AMA.circCDF[[2]]>=0.5][1]
    BCI.circCDF[[1]][BCI.circCDF[[2]]>=0.5][1]
    BKT.circCDF[[1]][BKT.circCDF[[2]]>=0.5][1]
    HKK.circCDF[[1]][HKK.circCDF[[2]]>=0.5][1]
    KCH.circCDF[[1]][KCH.circCDF[[2]]>=0.5][1]
    
    tiff(width=5, height=4, file="Figure3_TrunkCircularityECDFs.tiff",res=300,units="in")
      
     par(xpd=F,family="sans", las=1, oma=c(1,4,0,0), mar=c(3,4,1,1))
    plot(0,type='n',
           xlim=range(CircSample$iso),
           ylim=c(0.01,1),
           xlab=NA,ylab=NA,
           log='y')
      # Plot ECDF for each site
        lines(AMA.circCDF, col=site.cols[site.cols$site=="AMA","col"], lwd=2)
        lines(BCI.circCDF, col=site.cols[site.cols$site=="BCI","col"], lwd=2)
        lines(BKT.circCDF, col=site.cols[site.cols$site=="BKT","col"], lwd=2)
        lines(HKK.circCDF, col=site.cols[site.cols$site=="HKK","col"], lwd=2)
        lines(KCH.circCDF, col=site.cols[site.cols$site=="KCH","col"], lwd=2)
      abline(h=1,lty=2)
      abline(h=0.5,lty=1)
      # Plot vertical line for mean ciruclarity per site
        abline(v=mean(CircSample[CircSample$Site=='AMA','iso']), col=site.cols[site.cols$site=="AMA","col"])
        abline(v=mean(CircSample[CircSample$Site=='BCI','iso']), col=site.cols[site.cols$site=="BCI","col"])
        abline(v=mean(CircSample[CircSample$Site=='BKT','iso']), col=site.cols[site.cols$site=="BKT","col"])
        abline(v=mean(CircSample[CircSample$Site=='HKK','iso']), col=site.cols[site.cols$site=="HKK","col"])
        abline(v=mean(CircSample[CircSample$Site=='KCH','iso']), col=site.cols[site.cols$site=="KCH","col"])
      legend(x=0.3,y=0.5,
             sitesNames[1:5],
             col=site.cols$col[1:5],
             bty='n',
             lwd=2)
      mtext("Trunk circularity", side=1, line=2)
      mtext("Cumulative\nproportion\nof trees", side=2, line=3.5)
    dev.off()
    
    Circ.Test1 <- kruskal.test(list(CircSample[CircSample$Site=='AMA','iso'],
                                    CircSample[CircSample$Site=='BCI','iso'],
                                    CircSample[CircSample$Site=='BKT','iso'],
                                    CircSample[CircSample$Site=='HKK','iso'],
                                    CircSample[CircSample$Site=='KCH','iso']))
    Circ.Test2 <- kruskal.test(list(CircSample[CircSample$Site=='AMA','iso'],
                                    CircSample[CircSample$Site=='BCI','iso'],
                                    CircSample[CircSample$Site=='HKK','iso'],
                                    CircSample[CircSample$Site=='KCH','iso']))
    
#### Figure 4: Proportion of stems measured at nonstandard heights over time at each plot ####
  HOM.results <- read.csv("Data file_HOMresultsPerPlot.csv")
  
  axisSize <- 0.9
  cexAx <- 1.1
  textCex <- 1
  ptSize <- 0.8
  pchSig <- c(19,19,19,19,19)
  
  tiff(file="Figure4_MeasHtChangesOverTime.tiff",width=3,height=6, units="in", res=300, family = "sans")

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
       pch=19)
    
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
    text(x=1990.5,y=3.14,"c", cex=textCex+0.3)
    for(i in 1:length(sites)){
        points(x=HOM.results[HOM.results$Site==sites[i],"Year"],
               y=HOM.results[HOM.results$Site==sites[i],"MeanHOM"], 
               pch=pchSig[i], col=site.cols[site.cols$site==sites[i],"col"], cex=ptSize)
        lines(x=HOM.results[HOM.results$Site==sites[i],"Year"],
              y=HOM.results[HOM.results$Site==sites[i],"MeanHOM"], col=site.cols[site.cols$site==sites[i],"col"], lty=1)
    }
  dev.off()
  
  