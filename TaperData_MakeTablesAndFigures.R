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
  
  model3a <- lme4::lmer(b1.iso~log(DBH) + log(HOM) + log(WSG) + (1|Site), data = TaperSample)

  tiff(file="Figure2_TaperParameterVariation.tiff", height=3, width=8, units="in", res=300)
    par(mfrow=c(1,3), mar=c(4,2,0,0), xpd=F,family="sans", las=1, oma=c(0,5,1,1))

    plot(b1.iso~DBH, data=TaperSample, type='n',
         ylim=range(TaperSample$b1.iso,na.rm=T),
         log="x",
         xlab=NA, ylab=NA,
         cex.axis=cexAx)
    mtext("DAB (cm)", side=1, line=2.2, cex=axisSize)
    par(las=0)
    mtext("Taper parameter (b)", side=2, line=3.5, cex=axisSize)
    par(las=1)
    text(x=19,y=0.15,"a", cex=textCex+0.1)
    
    for(i in 1:length(sites)){
      points(b1.iso~DBH, data=TaperSample[TaperSample$Site==sites[i],],
             col=site.cols[site.cols$site==sites[i],"col"], pch=20)
      
      xrange=seq(from=min(TaperSample[TaperSample$Site==sites[i],"DBH"]),
                 to=max(TaperSample[TaperSample$Site==sites[i],"DBH"]),0.1)
      
      ypred=predict(model3a,newdata=data.frame(DBH=xrange,
                                                  Site=sites[i],
                                                  HOM = mean(TaperSample[TaperSample$Site==sites[i],"HOM"]),
                                                  WSG = mean(TaperSample[TaperSample$Site==sites[i],"WSG"])))
      
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

    for(i in 1:length(sites)){
      points(b1.iso~HOM, data=TaperSample[TaperSample$Site==sites[i],],
             col=site.cols[site.cols$site==sites[i],"col"], pch=20)
      
      xrange=seq(from=min(TaperSample[TaperSample$Site==sites[i],"HOM"]),
                 to=max(TaperSample[TaperSample$Site==sites[i],"HOM"]),0.1)
      
      ypred=predict(model3a,newdata=data.frame(DBH=mean(TaperSample[TaperSample$Site==sites[i],"DBH"]),
                                               Site=sites[i],
                                               HOM = xrange,
                                               WSG = mean(TaperSample[TaperSample$Site==sites[i],"WSG"])))
      
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
      
      ypred=predict(model3a,newdata=data.frame(DBH=mean(TaperSample[TaperSample$Site==sites[i],"DBH"]),
                                               Site=sites[i],
                                               HOM = mean(TaperSample[TaperSample$Site==sites[i],"HOM"]),
                                               WSG = xrange))
      
      lines(x=xrange,y=ypred, col=site.cols[site.cols$site==sites[i],"col"],
            lwd=2)
   }
    
  dev.off()  


#### Figure 3: Violin plots for taper for each site #####
  
    TaperSample <- read.csv("DataFile_TaperParameterSample.csv")
    TaperBySite <- read.csv("TaperVariationTable.csv")

    library(ggplot2)
  
      taperPlot <- ggplot(TaperSample, aes(x=Site, y=b1.iso, fill=Site)) + 
        geom_violin(trim=F) +
        labs(x="", y = "Tree taper parameter (b)") +
        scale_fill_manual(values = as.character(site.cols$col)) +
        scale_y_continuous(breaks=seq(-0.05,0.15,0.05)) +
        theme_bw() +
        theme(panel.border = element_rect(colour = "black"), panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(), axis.line = element_blank(),
              axis.text.x=element_blank(),
              axis.title.x=element_blank(),
              axis.ticks.x=element_blank()) + 
        #geom_point(data=TaperBySite,aes(x=Site, y = Mean), colour="white", shape=17, show.legend=F) + 
        #geom_point(data=TaperBySite,aes(x=Site, y = Mean.BA), colour="white", shape=15, show.legend=F) +
        geom_point(data=TaperBySite,aes(x=Site, y = Mean.AGB), colour="white", shape=16, show.legend=F)
      
    pdf(width=4.5, height=3, file="Figure3_TaperViolins.pdf")
      
     par(xpd=F,family="sans", mfrow=c(1,1), las=1, oma=c(1,4,0,3), mar=c(3,4,1,1))  
     
     taperPlot
     
    dev.off()
    
#### Figure 4: Proportion of stems measured at nonstandard heights over time at each plot ####
  HOM.results <- read.csv("Data file_HOMresultsPerPlot.csv")
  
  axisSize <- 0.9
  cexAx <- 1.1
  textCex <- 1
  ptSize <- 0.8
  pchSig <- c(19,19,19,19,19)
  
  tiff(file="Figure4_MeasHtChangesOverTime.tiff",width=3,height=6, units="in", res=300, family = "sans")

    par(mfrow=c(3,1), mar=c(1,3,1,1),oma=c(2,1,0,0), family="sans", xpd=F, las=1)
    
    plot(x=HOM.results$Year, y=HOM.results$Prop*100,
         type='n',
         ylab=NA,ylim=c(0,100),
         xlab=NA,
         xaxt = "n",
         cex.axis=cexAx)
    par(las=0)
    mtext("Stems (%)", side=2, line=2.5, cex=axisSize)
    par(las=1)
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
    par(las=0)
    mtext("Basal area (%)", side=2, line=2.5, cex=axisSize)
    par(las=1)
    text(x=1990.5,y=97,"b", cex=textCex+0.3)
    for(i in 1:length(sites)){
        points(x=HOM.results[HOM.results$Site==sites[i],"Year"],
               y=HOM.results[HOM.results$Site==sites[i],"PropBA"]*100, 
               pch=pchSig[i], col=site.cols[site.cols$site==sites[i],"col"], cex=ptSize)
        lines(x=HOM.results[HOM.results$Site==sites[i],"Year"],
              y=HOM.results[HOM.results$Site==sites[i],"PropBA"]*100, col=site.cols[site.cols$site==sites[i],"col"], lty=1)
    }
    
    plot(x=HOM.results$Year,y=HOM.results$MeanHOM.AGB,
         type='n',
         ylab=NA,
         xlab=NA,
         cex.axis=cexAx)
    mtext("Census year", side=1, line=2, cex=axisSize)
    par(las=0)
    mtext("Average HOM (m)", side=2, line=2.5, cex=axisSize)
    par(las=1)
    text(x=1990.5,y=3.14,"c", cex=textCex+0.3)
    for(i in 1:length(sites)){
        points(x=HOM.results[HOM.results$Site==sites[i],"Year"],
               y=HOM.results[HOM.results$Site==sites[i],"MeanHOM.AGB"], 
               pch=pchSig[i], col=site.cols[site.cols$site==sites[i],"col"], cex=ptSize)
        lines(x=HOM.results[HOM.results$Site==sites[i],"Year"],
              y=HOM.results[HOM.results$Site==sites[i],"MeanHOM.AGB"], col=site.cols[site.cols$site==sites[i],"col"], lty=1)
    }
  dev.off()
  
  
#### Figure 5: Violin plots for circularity for each site #####

    CircSample <- read.csv("DataFile_CircSample.csv")
    
    CircBySite <- read.csv("CircVariationTable.csv")

    library(ggplot2)

      circPlot <- ggplot(CircSample, aes(x=Site, y=iso, fill=Site)) + 
        geom_violin(trim=F) +
        labs(x="", y = "Tree circularity") +
        scale_fill_manual(values = as.character(site.cols$col)) +
        theme_bw() +
        theme(panel.border = element_rect(colour = "black"), panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(), axis.line = element_blank(),
              axis.text.x=element_blank(),
              axis.title.x=element_blank(),
              axis.ticks.x=element_blank()) + 
        #geom_point(data=CircBySite,aes(x=Site, y = Mean), colour="white", shape=17, show.legend=F) + 
        #geom_point(data=CircBySite,aes(x=Site, y = Mean.BA), colour="white", shape=15, show.legend=F) +
        geom_point(data=CircBySite,aes(x=Site, y = Mean.AGB), colour="white", shape=16, show.legend=F)
      

    
    pdf(width=4.5, height=3, file="Figure5_CircularityViolins.pdf")
      
     par(xpd=F,family="sans", mfrow=c(1,1), las=1, oma=c(1,4,0,3), mar=c(3,4,1,1))  
     
     circPlot
     
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
    
#### Figure 6: Compare biomass estimates from different taper models ####
  AGB_Results <- read.csv("DataFile_BiomassEstimates.csv")
  
  tiff(file="Figure6_BiomassEstimates.tiff",width=6,height=4, units="in", res=400, family = "sans")

    par(mfrow=c(1,1), mar=c(3,3,1,1),oma=c(1,1,1,1), family="sans", xpd=F, las=1)
        
      plot(x=1:5,
           y=AGB_Results$AGB_MeasTaper_Std,
           ylab="AGB relative to uncorrected",
           xlim=c(0.9,5.5),
           ylim=c(1,1.35),
           pch=19,
           col="#253494",
           xaxt="na",
           xlab = NA)
      # Add full taper model
      arrows(x0=1:5+0.1, x1=1:5+0.1,
             y0=AGB_Results$AGB_FullModel_Min95_Std, y1=AGB_Results$AGB_FullModel_Max95_Std,
             length=0,
             angle=90,
             col="#2c7fb8",
             lwd=2)
      points(x=1:5+0.1, y=AGB_Results$AGB_FullModel_Mean_Std, pch=19, col="#2c7fb8")
      # 3 parameter cross-validation
      arrows(x0=1:5+0.2, x1=1:5+0.2,
             y0=AGB_Results$AGB_Cross3Pr_Min95_Std, y1=AGB_Results$AGB_Cross3Pr_Max95_Std,
             length=0,
             angle=90,
             col="#41b6c4",
             lwd=2)
      points(x=1:5+0.2, y=AGB_Results$AGB_Cross3Pr_Mean_Std, pch=19, col="#41b6c4")
      # 2 parameter cross-validation
      arrows(x0=1:5+0.3, x1=1:5+0.3,
             y0=AGB_Results$AGB_Cross2Pr_Min95_Std, y1=AGB_Results$AGB_Cross2Pr_Max95_Std,
             length=0,
             angle=90,
             col="#a1dab4",
             lwd=2)
      points(x=1:5+0.3, y=AGB_Results$AGB_Cross2Pr_Mean_Std, pch=19, col="#a1dab4")
      # Add "median" taper value
      arrows(x0=1:5+0.4, x1=1:5+0.4,
             y0=AGB_Results$AGB_CrossMed_Min95_Std, y1=AGB_Results$AGB_CrossMed_Max95_Std,
             length=0,
             angle=90,
             col="grey",
             lwd=2)
      points(x=1:5+0.4, y=AGB_Results$AGB_CrossMed_Mean_Std, pch=19, col="grey")
      
      abline(h=1,lty=2)
      
    legend(x=1,y=1.37,
           c("Measured taper",
             "Taper modeled from DAB, HOM, and WSG, all sites (Model 1, Table 3)",
             "Taper modeled from DAB, HOM and WSG, other sites",
             "Taper modeled from DAB and HOM, other sites",
             "Biomass-weighted median taper, other sites"),
           y.intersp = 0.75,
           cex = 0.8,
           pch=19,
           col=c("#253494","#2c7fb8","#41b6c4","#a1dab4","grey"),
           bty="n")
    
    text(x=1.2,y=1.01,sitesNames[1], cex=0.8)
    text(x=2.2,y=1.03,sitesNames[2], cex=0.8)
    text(x=3.2,y=1.01,sitesNames[3], cex=0.8)
    text(x=4.2,y=1.03,sitesNames[4], cex=0.8)
    text(x=5.2,y=1.01,sitesNames[5], cex=0.8)
    par(las=0)
    mtext("Relative AGB", side=2,line = 3)
    
  dev.off()