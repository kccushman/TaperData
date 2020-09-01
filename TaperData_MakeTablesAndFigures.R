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
  textCex <- 1.2
  tiff(file="Figure2_TaperParameterVariation.tiff", height=2.2, width=8, units="in", res=300)
    layout(matrix(c(1,2,3,4),1,4,byrow = T),
           widths=c(1.5,1.5,1.5,0.8),heights=c(1,1,1,1))
    par(mar=c(4,3,0,0), xpd=F,family="sans", las=1, oma=c(0,5,1,1))

    plot(log(b1.iso)~log(DBH), data=TaperSample, type='n',
         ylim=range(log(TaperSample$b1.iso),na.rm=T)+c(0,1),
         xlab=NA, ylab=NA,
         cex.axis=cexAx)
    mtext("log DAB (cm)", side=1, line=2.2, cex=axisSize)
    mtext("log Taper \nparameter", side=2, line=2, cex=axisSize)
    text(x=3,y=-1.1,"a", cex=textCex+0.1)
    for(i in 1:length(sites)){
      points(log(b1.iso)~log(DBH), data=TaperSample[TaperSample$Site==sites[i],],
             col=site.cols[site.cols$site==sites[i],"col"], pch=20)
    }
    taperDABlm <- lme4::lmer(log(b1.iso)~log(DBH) + (1|Site), data = TaperSample)
    abline(a=summary(taperDABlm)$coefficients[1,1], b=summary(taperDABlm)$coefficients[2,1])
    
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
    
    par(xpd=NA)
    plot(0,type="n",bty="n",
         xlim=c(0,10),ylim=c(0,10),xlab="",ylab="",
         xaxt="n",yaxt="n")
    legend(x=-6.5, y=10,bty='n',
           sitesNames,
           col=site.cols$col,
           cex=cexAx,
           pch=19)
  dev.off()  


#### Figure 3: proportion of stems measured at nonstandard heights over time at each plot ####
  HOM.results <- read.csv("Data file_HOMresultsPerPlot.csv")
  
  axisSize <- 0.9
  cexAx <- 1.1
  textCex <- 1
  ptSize <- 0.8
  
  tiff(file="Figure3_MeasHtChangesOverTime.tiff",width=3,height=6, units="in", res=300, family = "sans")

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
               y=HOM.results[HOM.results$Site==sites[i],"Prop"]*100, pch=19, col=site.cols[site.cols$site==sites[i],"col"], cex=ptSize)
        lines(x=HOM.results[HOM.results$Site==sites[i],"Year"],
              y=HOM.results[HOM.results$Site==sites[i],"Prop"]*100, col=site.cols[site.cols$site==sites[i],"col"], lty=1)
    }
    
    legend(x=2000, y=105, bty='n',
       sitesNames,
       col=site.cols$col,
       cex=cexAx-0.1,
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
               y=HOM.results[HOM.results$Site==sites[i],"PropBA"]*100, pch=19, col=site.cols[site.cols$site==sites[i],"col"], cex=ptSize)
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
    text(x=1990.5,y=3.15,"c", cex=textCex+0.3)
    for(i in 1:length(sites)){
        points(x=HOM.results[HOM.results$Site==sites[i],"Year"],
               y=HOM.results[HOM.results$Site==sites[i],"MeanHOM"], pch=19, col=site.cols[site.cols$site==sites[i],"col"], cex=ptSize)
        lines(x=HOM.results[HOM.results$Site==sites[i],"Year"],
              y=HOM.results[HOM.results$Site==sites[i],"MeanHOM"], col=site.cols[site.cols$site==sites[i],"col"], lty=1)
    }


  dev.off()