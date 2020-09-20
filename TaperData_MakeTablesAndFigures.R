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

  cexAx <- 1.2
  axisSize <- 1
  textCex <- 1.4
  
  tiff(file="Figure2_TaperParameterVariation.tiff", height=3, width=8, units="in", res=300)
    par(mfrow=c(1,3), mar=c(4,2,0,0), xpd=F,family="sans", las=1, oma=c(0,5,1,1))

    plot(b1.iso~DBH, data=TaperSample, type='n',
         ylim=range(TaperSample$b1.iso,na.rm=T),
         log="x",
         xlab=NA, ylab=NA,
         cex.axis=cexAx)
    mtext("DAB (cm)", side=1, line=2.2, cex=axisSize)
    mtext("taper \nparameter", side=2, line=7, cex=axisSize, adj=0)
    text(x=19,y=0.15,"a", cex=textCex+0.1)
    
    taperDABlm <- lme4::lmer(b1.iso~log(DBH) + (1|Site), data = TaperSample)
    for(i in 1:length(sites)){
      points(b1.iso~DBH, data=TaperSample[TaperSample$Site==sites[i],],
             col=site.cols[site.cols$site==sites[i],"col"], pch=20)
      xrange=seq(from=min(TaperSample[TaperSample$Site==sites[i],"DBH"]),
                 to=max(TaperSample[TaperSample$Site==sites[i],"DBH"]),0.1)
      ypred=predict(taperDABlm,newdata=data.frame(DBH=xrange, Site=sites[i]))
      lines(x=xrange,y=ypred, col=site.cols[site.cols$site==sites[i],"col"],
            lwd=2)
    }

        legend(x=17, y=0,bty='n',
           sitesNames,
           col=site.cols$col,
           cex=0.9,
           pch=19)
        
    plot(b1.iso~HOM, data=TaperSample, type='n',
         ylim=range(TaperSample$b1.iso,na.rm=T),
         xlab=NA,ylab=NA,
         yaxt='n',
         log="x",
         cex.axis=cexAx)
  mtext("HOM (cm)", side=1, line=2.2, cex=axisSize)
  text(x=1.47,y=0.15,"b", cex=textCex+0.1)
    taperHOMlm <- lme4::lmer(b1.iso~log(HOM) + (1|Site), data = TaperSample)
    for(i in 1:length(sites)){
      points(b1.iso~HOM, data=TaperSample[TaperSample$Site==sites[i],],
             col=site.cols[site.cols$site==sites[i],"col"], pch=20)
      xrange=seq(from=min(TaperSample[TaperSample$Site==sites[i],"HOM"]),
                 to=max(TaperSample[TaperSample$Site==sites[i],"HOM"]),0.1)
      ypred=predict(taperHOMlm,newdata=data.frame(HOM=xrange, Site=sites[i]))
      lines(x=xrange,y=ypred, col=site.cols[site.cols$site==sites[i],"col"],
            lwd=2)
    }
    
    plot(b1.iso~WSG, data=TaperSample, type='n',
         ylim=range(TaperSample$b1.iso,na.rm=T),
         xlab=NA,ylab=NA,
         yaxt='n',
         log="x",
         cex.axis=cexAx)
  mtext("WSG (cm)", side=1, line=2.2, cex=axisSize)
  text(x=0.28,y=0.15,"c", cex=textCex+0.1)
    taperWSGlm <- lme4::lmer(b1.iso~log(WSG) + (1|Site), data = TaperSample)
    for(i in 1:length(sites)){
      points(b1.iso~WSG, data=TaperSample[TaperSample$Site==sites[i],],
             col=site.cols[site.cols$site==sites[i],"col"], pch=20)
      xrange=seq(from=min(TaperSample[TaperSample$Site==sites[i],"WSG"]),
                 to=max(TaperSample[TaperSample$Site==sites[i],"WSG"]),0.1)
      ypred=predict(taperWSGlm,newdata=data.frame(WSG=xrange, Site=sites[i]))
      lines(x=xrange,y=ypred, col=site.cols[site.cols$site==sites[i],"col"],
            lwd=2)
   }
    
  dev.off()  


#### Figure 3: CDFs for taper and circularity for each site #####
  # Calculate empirical cumulative distribution functions for each plot
    CircSample <- TreeSample[!is.na(TreeSample$iso),]

    AMA.circCDF <- Hmisc::Ecdf(CircSample[CircSample$Site=='AMA','iso'])
    BCI.circCDF <- Hmisc::Ecdf(CircSample[CircSample$Site=='BCI','iso'])
    BKT.circCDF <- Hmisc::Ecdf(CircSample[CircSample$Site=='BKT','iso'])
    HKK.circCDF <- Hmisc::Ecdf(CircSample[CircSample$Site=='HKK','iso'])
    KCH.circCDF <- Hmisc::Ecdf(CircSample[CircSample$Site=='KCH','iso'])
    
    AMA.taperCDF <- Hmisc::Ecdf(TaperSample[TaperSample$Site=='AMA','b1.iso'])
    BCI.taperCDF <- Hmisc::Ecdf(TaperSample[TaperSample$Site=='BCI','b1.iso'])
    BKT.taperCDF <- Hmisc::Ecdf(TaperSample[TaperSample$Site=='BKT','b1.iso'])
    HKK.taperCDF <- Hmisc::Ecdf(TaperSample[TaperSample$Site=='HKK','b1.iso'])
    KCH.taperCDF <- Hmisc::Ecdf(TaperSample[TaperSample$Site=='KCH','b1.iso'])
    
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
    
    tiff(width=8, height=4, file="Figure3_TrunkTaperCircularityECDFs.tiff",res=300,units="in")
      
     par(xpd=F,family="sans", mfrow=c(1,2), las=1, oma=c(1,4,0,0), mar=c(3,4,1,1))
     
    plot(0,type='n',
           xlim=range(TaperSample$b1.iso),
           ylim=c(0.01,1),
           xlab=NA,ylab=NA)
      # Plot ECDF for each site
        lines(AMA.taperCDF, col=site.cols[site.cols$site=="AMA","col"], lwd=2)
        lines(BCI.taperCDF, col=site.cols[site.cols$site=="BCI","col"], lwd=2)
        lines(BKT.taperCDF, col=site.cols[site.cols$site=="BKT","col"], lwd=2)
        lines(HKK.taperCDF, col=site.cols[site.cols$site=="HKK","col"], lwd=2)
        lines(KCH.taperCDF, col=site.cols[site.cols$site=="KCH","col"], lwd=2)
      abline(h=1,lty=2)
      abline(h=0.5,lty=1)
      text("a",x=-0.06,y=0.95)

      legend(x=0.025,y=0.35,
             sitesNames[1:5],
             col=site.cols$col[1:5],
             bty='n',
             lwd=2, cex=0.8)
      mtext("Tree taper parameter (b)", side=1, line=2)
      mtext("Cumulative\nproportion\nof trees", side=2, line=3.5)
  
      
        plot(0,type='n',
           xlim=range(CircSample$iso),
           ylim=c(0.01,1),
           xlab=NA,ylab=NA)
      # Plot ECDF for each site
        lines(AMA.circCDF, col=site.cols[site.cols$site=="AMA","col"], lwd=2)
        lines(BCI.circCDF, col=site.cols[site.cols$site=="BCI","col"], lwd=2)
        lines(BKT.circCDF, col=site.cols[site.cols$site=="BKT","col"], lwd=2)
        lines(HKK.circCDF, col=site.cols[site.cols$site=="HKK","col"], lwd=2)
        lines(KCH.circCDF, col=site.cols[site.cols$site=="KCH","col"], lwd=2)
      abline(h=1,lty=2)
      abline(h=0.5,lty=1)
      text("b",x=0.32,y=0.95)
      mtext("Trunk circularity", side=1, line=2)
      
    dev.off()
    
    Circ.Test1 <- kruskal.test(list(CircSample[CircSample$Site=='AMA','iso'],
                                    CircSample[CircSample$Site=='BCI','iso'],
                                    CircSample[CircSample$Site=='BKT','iso'],
                                    CircSample[CircSample$Site=='HKK','iso'],
                                    CircSample[CircSample$Site=='KCH','iso']))
    Circ.Test2 <- kruskal.test(list(CircSample[CircSample$Site=='AMA' & CircSample$DBH>=30,'iso'],
                                    CircSample[CircSample$Site=='BCI' & CircSample$DBH>=30,'iso'],
                                    CircSample[CircSample$Site=='BKT' & CircSample$DBH>=30,'iso'],
                                    CircSample[CircSample$Site=='HKK' & CircSample$DBH>=30,'iso'],
                                    CircSample[CircSample$Site=='KCH' & CircSample$DBH>=30,'iso']))
    
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
  
  