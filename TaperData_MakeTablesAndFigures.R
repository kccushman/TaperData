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
  HOM.results <- read.csv("Data file_HOMresultsPerPlot30.csv")
  
  axisSize <- 0.9
  cexAx <- 1.1
  textCex <- 1
  ptSize <- 0.8
  pchSig <- c(1,19,1,1,19)
  
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
  
#### Figure 4: Taper variation with climate ####
  
  # Vector of MAP in site alphabetical order
    MAP <- c(3200,2600,2500,1500,2800)
  # Vector of dry season length in site alphabetical order
    DSL <- c(0,3,0,5,2.5)
  # Vector of distance from equator
    Lat <- c(3.8, 9.2, 1.3, 15.4, 7.5)
    
  # Vector of best model site random effects
    b1 <- c(-0.506965136, -0.371995616, -0.518252059, -0.328540579, -0.458079386)

    labCx = 0.9
    
    tiff(file="Figure4_TaperEffectsVsClimate.tiff", height=3, width=6, units="in", res=300)
      par(mfrow=c(1,2), mar=c(2,1,0,0), oma=c(2,5,1,1),
          xpd=F, family="sans", las=1)
      
      plot(b1 ~ MAP,
           type="n",
           xlim=c(1200,3400),
           ylab=NA, xlab=NA)
      for(i in 1:length(sites)){
        points(x = MAP[i], y = b1[i], pch=20, col=site.cols$col[i])
      }
      text("a", x=1200, y=-0.33)
      mtext("Mean annual precipitation (mm/yr)", side=1, line=2, cex=labCx)
      mtext("Site \nrandom \neffect", side=2, line=6, adj=0, cex=labCx)
      
      plot(b1 ~ DSL,
           type="n", yaxt = "n",
           xlim=c(-0.5,5.5),
           ylab=NA, xlab=NA)
      for(i in 1:length(sites)){
        points(x = DSL[i], y = b1[i], pch=20, col=site.cols$col[i])
      }
      text("b", x=-0.5, y=-0.33)
      mtext("Dry season length (months/yr)", side=1, line=2, cex=labCx)
      
    dev.off()
    
    summary(lm(b1 ~ MAP))
    summary(lm(b1 ~ DSL))
    
#### Figure S: ####    
    

  PropStems <- c(0.5512590, 0.5677175, 0.2314316, 0.1071939, 0.6266549)
  PropBA <- c(0.6104313, 0.7050572, 0.2970426, 0.1605881,  0.7186443)
  MeanHOM <- c(2.554515, 3.658446, 1.664472, 1.541162, 3.012815)
  
  plot(PropBA~MAP, pch =20)
  plot(PropBA~DSL, pch =20)
  plot(PropBA~Lat, pch =20)
  
  summary(lm(PropBA~MAP))  
  summary(lm(PropBA~DSL))  
  summary(lm(PropBA~Lat)) 
  
  plot(PropStems~MAP, pch =20)
  plot(PropStems~DSL, pch =20)
  plot(PropStems~Lat, pch =20)
  
  summary(lm(PropStems~MAP))  
  summary(lm(PropStems~DSL)) 
  summary(lm(PropStems~Lat)) 
  
  plot(MeanHOM~MAP, pch =20)
  plot(MeanHOM~DSL, pch =20)
  plot(MeanHOM~Lat, pch =20)
  
  summary(lm(MeanHOM~MAP))  
  summary(lm(MeanHOM~DSL)) 
  summary(lm(MeanHOM~Lat))
  