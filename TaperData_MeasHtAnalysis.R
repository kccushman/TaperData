# Load species list of figs so the stranglers can be omitted
  FigList <- read.csv("~/Desktop/Taper/Current/RawCensusData_CTFS/all.fig.sp.csv")
  StranglerFigs <- FigList[FigList$Strangler=="Yes","Mnemonic"]

## Name columns to include for all censuses (EDIT??)
  taper.cols <- c('TreeID','StemID','Tag','StemTag','Mnemonic','Genus','Family','WSG',
                  'dbh','hom','ExactDate','date','DFstatus','status','codes','notes',
                  'Site','Census')

#### 1. Format data to include for each plot ####

  ## Load and format species lists
  sp.list <- read.csv("~/Desktop/Taper/Current/RawCensusData_CTFS/Species list.csv")
    AMA.sp.list <- sp.list[sp.list$plot=='Amacayacu',]
    BKT.sp.list <- sp.list[sp.list$plot=='Bukit Timah Big Trees',]
    HKK.sp.list <- sp.list[sp.list$plot=='Huai Kha Khaeng',]
  ## Load species list for Barro Colorado
  BCI.sp.list <- read.csv("~/Desktop/Taper/Current/RawCensusData_CTFS/BCI.species.list.csv")
    BCI.sp.list[BCI.sp.list$mnemonic=='protte','subspecies'] <- ""
    BCI.sp.list[BCI.sp.list$mnemonic=='quaras','subspecies'] <- ""
    BCI.sp.list$LatinAll <- paste(BCI.sp.list$genus, BCI.sp.list$species, BCI.sp.list$subspecies)
  ## Load species list for Khao Chong
  KCH.sp.list <- read.csv("~/Desktop/Taper/Current/RawCensusData_CTFS/KCH.species.table.csv")
  ## Combine into one data frame
  AMA.sp.list$site <- 'AMA'
  BCI.sp.list$site <- 'BCI'
  BKT.sp.list$site <- 'BKT'
  HKK.sp.list$site <- 'HKK'
  KCH.sp.list$site <- 'KCH'
  taper.sp.list <- rbind(AMA.sp.list[,c('mnemonic','genus','species','family','idlevel','site')],
                         BCI.sp.list[,c('mnemonic','genus','species','family','idlevel','site')],
                         BKT.sp.list[,c('mnemonic','genus','species','family','idlevel','site')],
                         HKK.sp.list[,c('mnemonic','genus','species','family','idlevel','site')],
                         KCH.sp.list[,c('mnemonic','genus','species','family','idlevel','site')])
  colnames(taper.sp.list) <- c('Mnemonic','Genus','Species','Family','idlevel','Site')
  
#### 1.1 Amacayacu ####
  load("~/Desktop/Taper/Current/RawCensusData_CTFS/amacayacu.stem1.RData")
  load("~/Desktop/Taper/Current/RawCensusData_CTFS/AMAcen2.RData")
    AMA.cen1 <- amacayacu.stem1
    AMA.cen2 <- AMAcen2
  
  # Rename columns
  AMA.cen1 <- AMA.cen1[,c('treeID','stemID','tag','StemTag','sp','dbh','hom','ExactDate','date','DFstatus','status')]
    colnames(AMA.cen1) <-c('TreeID','StemID','Tag','StemTag','Mnemonic','dbh','hom','ExactDate','date','DFstatus','status')
    AMA.cen1$codes <- NA
    AMA.cen1$notes <- NA
    AMA.cen1$Site <- 'AMA'
    AMA.cen1$Census <- 1

  AMA.cen2<- AMA.cen2[c('StemTag','sp','dbh.2..cm.','hom..cm..2','date.2','C.digos','Observaciones.remedicion')]
    colnames(AMA.cen2) <- c('StemTag','Mnemonic','dbh','hom','date','codes','notes')
    AMA.cen2$TreeID <- NA
    AMA.cen2$StemID <- NA
    AMA.cen2$Tag <- NA
    AMA.cen2$ExactDate <- NA
    AMA.cen2$DFstatus <- NA
    AMA.cen2$Site <- 'AMA'
    AMA.cen2$Census <- 2
    AMA.cen2$StemTag <- gsub("_" , "", AMA.cen2$StemTag)
    AMA.cen2$status <- "A"
      AMA.cen2[grep('D',AMA.cen2$codes),'status'] <- "D"

  # Make sure that HOM is in m and DBH is in cm
    AMA.cen1$hom <- AMA.cen1$hom/100
    AMA.cen1$dbh <- AMA.cen1$dbh/10
    AMA.cen2$hom <- AMA.cen2$hom/100

  # Keep stems that are alive and greater than 10 cm (or greater than 10 cm in census1)
    AMA.cen1 <- AMA.cen1[AMA.cen1$DFstatus=='alive',]
    AMA.cen1 <- AMA.cen1[AMA.cen1$status=='A',]
    AMA.cen1 <- AMA.cen1[AMA.cen1$dbh>=10 & !is.na(AMA.cen1$dbh),]

    AMA.cen2 <- AMA.cen2[AMA.cen2$status=='A',]
    AMA.cen2 <-  AMA.cen2[(AMA.cen2$dbh>=10 | AMA.cen2$StemTag %in% AMA.cen1$StemTag) & !is.na(AMA.cen2$dbh),]

  # Uniformly capitalize species mnemonic
  AMA.cen1$Mnemonic <- R.utils::capitalize(AMA.cen1$Mnemonic)
  AMA.cen2$Mnemonic <- R.utils::capitalize(as.character(AMA.cen2$Mnemonic))
  
  # Remove strangler figs
  AMA.cen1 <- AMA.cen1[!(AMA.cen1$Mnemonic %in% StranglerFigs), ]
  AMA.cen2 <- AMA.cen2[!(AMA.cen2$Mnemonic %in% StranglerFigs), ]
  
  # "StemID" should be unique "StemTag" value
  AMA.cen1$StemID <- AMA.cen1$StemTag
  AMA.cen2$StemID <- AMA.cen2$StemTag
  
  AMA.cens <- list(AMA.cen1, AMA.cen2)


#### 1.2 Barro Colorado ####

    # BCI.cen2 <- read.csv('~/Desktop/Taper/Current/RawCensusData_CTFS/BCI.census2.csv') # omit 1985 census now
    BCI.cen3 <- read.csv('~/Desktop/Taper/Current/RawCensusData_CTFS/BCI.census3.csv')  #1990
    BCI.cen4 <- read.csv('~/Desktop/Taper/Current/RawCensusData_CTFS/BCI.census4.csv')  #1995
    BCI.cen5 <- read.csv('~/Desktop/Taper/Current/RawCensusData_CTFS/BCI.census5.csv')  #2000
    BCI.cen6 <- read.csv('~/Desktop/Taper/Current/RawCensusData_CTFS/BCI.census6.csv')  #2005
    BCI.cen7 <- read.csv('~/Desktop/Taper/Current/RawCensusData_CTFS/BCI.census7.csv')  #2010
    BCI.cen8 <- read.csv('~/Desktop/Taper/Current/RawCensusData_CTFS/BCI.census8.csv')  #2015
    #BCI.cens <- list(BCI.cen2,BCI.cen3,BCI.cen4,BCI.cen5,BCI.cen6,BCI.cen7, BCI.cen8) 
    BCI.cens <- list(BCI.cen3,BCI.cen4,BCI.cen5,BCI.cen6,BCI.cen7, BCI.cen8) 
    nCensus <- length(BCI.cens)

  for (i in 1:length(BCI.cens)){
    # Only keep stems with a diameter measurement
      BCI.cens[[i]] <- BCI.cens[[i]][!is.na(BCI.cens[[i]]$DBH),]
    
    # Merge with species list to get species mnemonic
    BCI.cens[[i]]$LatinAll <- paste(BCI.cens[[i]]$Latin,BCI.cens[[i]]$SubSpecies)
    BCI.cens[[i]] <- merge(BCI.cens[[i]],BCI.sp.list[c('LatinAll','mnemonic','genus','family')], by="LatinAll", all.x=T, all.y=F)
    

    # Keep correct columns
    BCI.cens[[i]] <- BCI.cens[[i]][,c('TreeID','StemID','Tag','StemTag','mnemonic','genus','family',
                                      'DBH','HOM','Date','Status','Codes')]
    colnames(BCI.cens[[i]]) <- c('TreeID','StemID','Tag','StemTag','Mnemonic','Genus','Family',
                                      'dbh','hom','date','status','codes')
    BCI.cens[[i]]$notes <- NA
    BCI.cens[[i]]$Site <- "BCI"
    BCI.cens[[i]]$Census <- i
    BCI.cens[[i]]$ExactDate <- NA
    BCI.cens[[i]]$DFstatus <- NA
      
    # Make sure that HOM is in m and DBH is in cm    
    BCI.cens[[i]]$dbh <- BCI.cens[[i]]$dbh/10
    # Keep stems that are alive
    BCI.cens[[i]] <- BCI.cens[[i]][BCI.cens[[i]]$status=='alive',]
    
    # Remove strangler figs
    BCI.cens[[i]] <- BCI.cens[[i]][!(BCI.cens[[i]]$Mnemonic %in% StranglerFigs), ]
    
    # Remove duplicate measurements where the new HOM is unknown (one entry in census 3 and one in census 7)
    BCI.cens[[i]] <- BCI.cens[[i]][!is.na(BCI.cens[[i]]$hom),]
  }

  # Keep stems that are greater than 10 cm DBH, or were greater than 10 cm DBH in a previous census
  BCI.bigstems <- c()
  for (i in 1:length(BCI.cens)){
    BCI.bigstems <- unique(c(BCI.bigstems,BCI.cens[[i]][BCI.cens[[i]]$dbh >= 10 & !is.na(BCI.cens[[i]]$dbh),'StemID']))
    BCI.cens[[i]] <- BCI.cens[[i]][(BCI.cens[[i]]$dbh >= 10 | BCI.cens[[i]]$StemID %in% BCI.bigstems) & !is.na(BCI.cens[[i]]$dbh),]
  }
  
#### 1.3 Bukit Timah ####
  
  load('~/Desktop/Taper/Current/RawCensusData_CTFS/bukittimahbigtrees.full1.rdata')
  BKT.cen1 <- bukittimahbigtrees.full1
    # Keep correct columns
    BKT.cen1 <- BKT.cen1[,c('treeID','stemID','tag','StemTag','sp',
                            'dbh','hom','ExactDate','date','DFstatus','status','codes')]
    colnames(BKT.cen1) <- c('TreeID','StemID','Tag','StemTag','Mnemonic',
                            'dbh','hom','ExactDate','date','DFstatus','status','codes')
    BKT.cen1$notes <- NA
    BKT.cen1$Site <- "BKT"
    BKT.cen1$Census <- 1
  # Make sure that HOM is in m and DBH is in cm  
  BKT.cen1$dbh <-  BKT.cen1$dbh/10
  #Keep only trees that are alive and greater than 10 cm dbh
  BKT.cen1 <- BKT.cen1[BKT.cen1$DFstatus=='alive',]
  BKT.cen1 <- BKT.cen1[BKT.cen1$dbh >=10 & !is.na(BKT.cen1$dbh),]
  # Remove strangler figs
  BKT.cen1 <- BKT.cen1[!(BKT.cen1$Mnemonic %in% StranglerFigs), ]
  # Round HOM to 1 decimal place
  BKT.cen1$hom <- round(BKT.cen1$hom,1)
  

#### 1.4 Huai Kha Khaeng ####
  load('~/Desktop/Taper/Current/RawCensusData_CTFS/hkk.stem1.RData')
  load('~/Desktop/Taper/Current/RawCensusData_CTFS/hkk.stem2.RData')
  load('~/Desktop/Taper/Current/RawCensusData_CTFS/hkk.stem3.RData')
  load('~/Desktop/Taper/Current/RawCensusData_CTFS/hkk.stem4.RData')
  HKK.cen1 <- hkk.stem1
  HKK.cen2 <- hkk.stem2
  HKK.cen3 <- hkk.stem3
  HKK.cen4 <- hkk.stem4
  #merge to get real height data for census 4
  HOM.HKK <- read.csv('~/Desktop/Taper/Current/RawCensusData_CTFS/HKK.census.data.csv')
  HOM.HKK$POM09 <- HOM.HKK$POM09/100
  HKK.cen4 <- merge(HKK.cen4,
    HOM.HKK[,c('Tag','POM09')],
    by.x='tag',by.y='Tag',all.x=T,all.y=F)
  HKK.cen4$hom <- ifelse(is.na(HKK.cen4$POM09),1.3,HKK.cen4$POM09)
  HKK.cens <- list(HKK.cen1,HKK.cen2,HKK.cen3,HKK.cen4)

  for (i in 1:length(HKK.cens)){
    # Keep correct columns
    HKK.cens[[i]] <- HKK.cens[[i]][,c('treeID','stemID','tag','StemTag','sp',
                                      'dbh','hom','ExactDate','date','DFstatus','status','codes')]
    colnames(HKK.cens[[i]]) <- c('TreeID','StemID','Tag','StemTag','Mnemonic',
                                      'dbh','hom','ExactDate','date','DFstatus','status','codes')
    HKK.cens[[i]]$notes <- NA
    HKK.cens[[i]]$Site <- "HKK"
    HKK.cens[[i]]$Census <- i
    # Make sure that HOM is in m and DBH is in cm    
    HKK.cens[[i]]$dbh <- HKK.cens[[i]]$dbh/10
    
    # Change precision of HOM to one decimal place (two significant figures)
    HKK.cens[[i]]$hom <- round(HKK.cens[[i]]$hom,1)
    
    # Keep stems that are alive
    HKK.cens[[i]] <- HKK.cens[[i]][HKK.cens[[i]]$DFstatus=='alive',]
    # Assume anything without a HOM is measured at 1.3 m
    HKK.cens[[i]][is.na(HKK.cens[[i]]$hom),'hom'] <- 1.3
    
    # Remove strangler figs
    HKK.cens[[i]] <- HKK.cens[[i]][!(HKK.cens[[i]]$Mnemonic %in% StranglerFigs), ]
  }

  # Keep stems that are greater than 10 cm DBH, or were greater than 10 cm DBH in a previous census
  HKK.bigstems <- c()
  for (i in 1:length(HKK.cens)){
    HKK.bigstems <- unique(c(HKK.bigstems,HKK.cens[[i]][HKK.cens[[i]]$dbh >= 10 & !is.na(HKK.cens[[i]]$dbh),'StemID']))
    HKK.cens[[i]] <- HKK.cens[[i]][(HKK.cens[[i]]$dbh >= 10 | HKK.cens[[i]]$StemID %in% HKK.bigstems) & !is.na(HKK.cens[[i]]$dbh),]
  }


#### 1.5 Khao Chong ####
  KCH.cens1and2 <- read.csv('~/Desktop/Taper/Current/RawCensusData_CTFS/ViewFullTable.csv',as.is=T)
    KCH.cens1and2$HOM <- as.numeric(KCH.cens1and2$HOM)
    KCH.cens1and2$DBH <- as.numeric(KCH.cens1and2$DBH)
  load('~/Desktop/Taper/Current/RawCensusData_CTFS/khaochong.stem3.RData')

  # Census 1
  KCH.cen1 <- KCH.cens1and2[KCH.cens1and2$CensusID==1,]
    # Keep correct columns
    KCH.cen1 <- KCH.cen1[,c('TreeID','StemID','Tag','StemTag','Mnemonic',
                            'DBH','HOM','ExactDate','Status','ListOfTSM')]
    colnames(KCH.cen1) <- c('TreeID','StemID','Tag','StemTag','Mnemonic',
                            'dbh','hom','ExactDate','status','codes')
    KCH.cen1$date <- NA
    KCH.cen1$DFstatus <- NA
    KCH.cen1$notes <- NA
    KCH.cen1$Site <- 'KCH'
    KCH.cen1$Census <- 1
    # Make sure that HOM is in m and DBH is in cm  
    KCH.cen1$dbh <- KCH.cen1$dbh/10
    # Keep stems that are alive
    KCH.cen1 <- KCH.cen1[KCH.cen1$status=='alive',]
    # Remove strangler figs
    KCH.cen1 <- KCH.cen1[!(KCH.cen1$Mnemonic %in% StranglerFigs), ]

  # Census 2
  KCH.cen2 <- KCH.cens1and2[KCH.cens1and2$CensusID==2,]
    # Keep correct columns
    KCH.cen2 <- KCH.cen2[,c('TreeID','StemID','Tag','StemTag','Mnemonic',
                            'DBH','HOM','ExactDate','Status','ListOfTSM')]
    colnames(KCH.cen2) <- c('TreeID','StemID','Tag','StemTag','Mnemonic',
                            'dbh','hom','ExactDate','status','codes')
    KCH.cen2$date <- NA
    KCH.cen2$DFstatus <- NA
    KCH.cen2$notes <- NA
    KCH.cen2$Site <- 'KCH'
    KCH.cen2$Census <- 2
    # Make sure that HOM is in m and DBH is in cm  
    KCH.cen2$dbh <- KCH.cen2$dbh/10
    # Keep stems that are alive
    KCH.cen2 <- KCH.cen2[KCH.cen2$status=='alive',]
    # Remove strangler figs
    KCH.cen2 <- KCH.cen2[!(KCH.cen2$Mnemonic %in% StranglerFigs), ]

    # Census 3
    KCH.cen3 <- khaochong.stem3
    # Keep correct columns
    KCH.cen3 <- KCH.cen3[,c('treeID','stemID','tag','StemTag','sp',
                                      'dbh','hom','ExactDate','date','DFstatus','status','codes')]
    colnames(KCH.cen3) <- c('TreeID','StemID','Tag','StemTag','Mnemonic',
                                      'dbh','hom','ExactDate','date','DFstatus','status','codes')
    KCH.cen3$notes <- NA
    KCH.cen3$Site <- "KCH"
    KCH.cen3$Census <- 3
    # Make sure that HOM is in m and DBH is in cm    
    KCH.cen3$dbh <- KCH.cen3$dbh/10
    # Keep stems that are alive
    KCH.cen3 <- KCH.cen3[KCH.cen3$DFstatus=='alive',]
    # Remove strangler figs
    KCH.cen3 <- KCH.cen3[!(KCH.cen3$Mnemonic %in% StranglerFigs), ]

    # Keep stems that are (or were) greater than 10 cm
    KCH.cen1 <- KCH.cen1[KCH.cen1$dbh >=10 & !is.na(KCH.cen1$dbh),] 
      KCH.bigstems <- unique(KCH.cen1$StemTag)
    KCH.cen2 <- KCH.cen2[(KCH.cen2$dbh >=10 | KCH.cen2$StemTag %in% KCH.bigstems) & !is.na(KCH.cen2$dbh),]
      KCH.bigstems <- unique(c(KCH.bigstems,KCH.cen2$StemTag))
    KCH.cen3 <- KCH.cen3[(KCH.cen3$dbh >=10 | KCH.cen3$StemTag %in% KCH.bigstems) & !is.na(KCH.cen3$dbh),]
    
    # Change StemID to unique StemTag value
    KCH.cen1$StemID <- KCH.cen1$StemTag
    KCH.cen2$StemID <- KCH.cen2$StemTag
    KCH.cen3$StemID <- KCH.cen3$StemTag
    
    # Change presicion of hom to two significant figures (one decimal place)
    KCH.cen1$hom <- round(KCH.cen1$hom,1)
    KCH.cen2$hom <- round(KCH.cen2$hom,1)
    KCH.cen3$hom <- round(KCH.cen3$hom,1)

  KCH.cens <- list(KCH.cen1,KCH.cen2,KCH.cen3)

#### 2. For all census data, merge with information on family and WSG and calculate AGB ####
  
  # Family: taper.sp.list
  # WSG
  load('~/Desktop/Taper/Current/RawCensusData_CTFS/wsg.ctfs2.Rdata')

#### 2.1 Amacayacu ####
  # Species list
  AMA.sp <- droplevels(taper.sp.list[taper.sp.list$Site=='AMA',c('Mnemonic','Genus','Family')])
    AMA.sp$Mnemonic <- R.utils::capitalize(as.character(AMA.sp$Mnemonic))
  # Wood density file
  AMA.wd <- droplevels(wsg.ctfs2[wsg.ctfs2$site=="Amacayacu",])
    AMA.wd.sp <- aggregate.data.frame(AMA.wd$wsg,by=list(AMA.wd$sp),FUN='mean',na.rm=T)
     colnames(AMA.wd.sp) <- c('Mnemonic','WSG')
     AMA.wd.sp$Mnemonic <- R.utils::capitalize(as.character(AMA.wd.sp$Mnemonic))
    AMA.wd.genus <- aggregate.data.frame(AMA.wd$wsg,by=list(AMA.wd$genus),FUN='mean',na.rm=T)
      colnames(AMA.wd.genus) <- c('Genus','g.wsg')
    AMA.wd.family <- aggregate.data.frame(AMA.wd$wsg,by=list(AMA.wd$fam),FUN='mean',na.rm=T)
      colnames(AMA.wd.family) <- c('Family','f.wsg')
  for(i in 1:length(AMA.cens)){
    AMA.cens[[i]] <- merge(AMA.cens[[i]],AMA.sp,by='Mnemonic',all.x=T,all.y=F)
    AMA.cens[[i]] <- merge(AMA.cens[[i]],AMA.wd.sp,by='Mnemonic',all.x=T,all.y=F)
    AMA.cens[[i]] <- merge(AMA.cens[[i]],AMA.wd.genus,by='Genus',all.x=T,all.y=F)
    AMA.cens[[i]] <- merge(AMA.cens[[i]],AMA.wd.family,by='Family',all.x=T,all.y=F)
    AMA.cens[[i]]$WSG <- ifelse(is.na(AMA.cens[[i]]$WSG),AMA.cens[[i]]$g.wsg,AMA.cens[[i]]$WSG)
    AMA.cens[[i]]$WSG <- ifelse(is.na(AMA.cens[[i]]$WSG),AMA.cens[[i]]$f.wsg,AMA.cens[[i]]$WSG)
    AMA.cens[[i]]$WSG <- ifelse(is.na(AMA.cens[[i]]$WSG),
                                weighted.mean(AMA.cens[[i]]$WSG, w=AMA.cens[[i]]$dbh, na.rm=T),
                                AMA.cens[[i]]$WSG)
    AMA.cens[[i]] <- AMA.cens[[i]][,taper.cols]
    
  }

#### 2.2 Barro Colorado Island ####
  # Species list
  BCI.sp <- droplevels(taper.sp.list[taper.sp.list$Site=='BCI',c('Mnemonic','Genus','Family')])
  # Wood density file
  BCI.wd <- droplevels(wsg.ctfs2[wsg.ctfs2$site=="bci",])
    BCI.wd.sp <- aggregate.data.frame(BCI.wd$wsg,by=list(BCI.wd$sp),FUN='mean',na.rm=T)
     colnames(BCI.wd.sp) <- c('Mnemonic','WSG')
    BCI.wd.genus <- aggregate.data.frame(BCI.wd$wsg,by=list(BCI.wd$genus),FUN='mean',na.rm=T)
      colnames(BCI.wd.genus) <- c('Genus','g.wsg')
    BCI.wd.family <- aggregate.data.frame(BCI.wd$wsg,by=list(BCI.wd$fam),FUN='mean',na.rm=T)
      colnames(BCI.wd.family) <- c('Family','f.wsg')
  for(i in 1:length(BCI.cens)){
    BCI.cens[[i]] <- merge(BCI.cens[[i]],BCI.wd.sp,by='Mnemonic',all.x=T,all.y=F)
    BCI.cens[[i]] <- merge(BCI.cens[[i]],BCI.wd.genus,by='Genus',all.x=T,all.y=F)
    BCI.cens[[i]] <- merge(BCI.cens[[i]],BCI.wd.family,by='Family',all.x=T,all.y=F)
    BCI.cens[[i]]$WSG <- ifelse(is.na(BCI.cens[[i]]$WSG),BCI.cens[[i]]$g.wsg,BCI.cens[[i]]$WSG)
    BCI.cens[[i]]$WSG <- ifelse(is.na(BCI.cens[[i]]$WSG),BCI.cens[[i]]$f.wsg,BCI.cens[[i]]$WSG)
    BCI.cens[[i]]$WSG <- ifelse(is.na(BCI.cens[[i]]$WSG),
                                weighted.mean(BCI.cens[[i]]$WSG, w=BCI.cens[[i]]$dbh, na.rm=T),
                                BCI.cens[[i]]$WSG)
    BCI.cens[[i]] <- BCI.cens[[i]][,taper.cols]
    
  }

#### 2.3 Bukit Timah ####
  # Species list
  BKT.sp <- droplevels(taper.sp.list[taper.sp.list$Site=='BKT',c('Mnemonic','Genus','Family')])
  # Wood density file
  BKT.wd <- droplevels(wsg.ctfs2[wsg.ctfs2$site=="Bukit Timah Big Tree",])
    BKT.wd.sp <- aggregate.data.frame(BKT.wd$wsg,by=list(BKT.wd$sp),FUN='mean',na.rm=T)
     colnames(BKT.wd.sp) <- c('Mnemonic','WSG')
    BKT.wd.genus <- aggregate.data.frame(BKT.wd$wsg,by=list(BKT.wd$genus),FUN='mean',na.rm=T)
      colnames(BKT.wd.genus) <- c('Genus','g.wsg')
    BKT.wd.family <- aggregate.data.frame(BKT.wd$wsg,by=list(BKT.wd$fam),FUN='mean',na.rm=T)
      colnames(BKT.wd.family) <- c('Family','f.wsg')

    BKT.cen1 <- merge(BKT.cen1,BKT.sp,by='Mnemonic',all.x=T,all.y=F)
    BKT.cen1 <- merge(BKT.cen1,BKT.wd.sp,by='Mnemonic',all.x=T,all.y=F)
    BKT.cen1 <- merge(BKT.cen1,BKT.wd.genus,by='Genus',all.x=T,all.y=F)
    BKT.cen1 <- merge(BKT.cen1,BKT.wd.family,by='Family',all.x=T,all.y=F)
    BKT.cen1$WSG <- ifelse(is.na(BKT.cen1$WSG),BKT.cen1$g.wsg,BKT.cen1$WSG)
    BKT.cen1$WSG <- ifelse(is.na(BKT.cen1$WSG),BKT.cen1$f.wsg,BKT.cen1$WSG)
    BKT.cen1$WSG <- ifelse(is.na(BKT.cen1$WSG),
                                weighted.mean(BKT.cen1$WSG, w=BKT.cen1$dbh, na.rm=T),
                                BKT.cen1$WSG)
    BKT.cen1 <- BKT.cen1[,taper.cols]
    BKT.cens <- list(BKT.cen1)

#### 2.4 Huai Kha Khaeng ####
  # Species list
  HKK.sp <- droplevels(taper.sp.list[taper.sp.list$Site=='HKK',c('Mnemonic','Genus','Family')])
  # Wood density file
  HKK.wd <- droplevels(wsg.ctfs2[wsg.ctfs2$site=="Huai Kha Khaeng",])
    HKK.wd.sp <- aggregate.data.frame(HKK.wd$wsg,by=list(HKK.wd$sp),FUN='mean',na.rm=T)
     colnames(HKK.wd.sp) <- c('Mnemonic','WSG')
    HKK.wd.genus <- aggregate.data.frame(HKK.wd$wsg,by=list(HKK.wd$genus),FUN='mean',na.rm=T)
      colnames(HKK.wd.genus) <- c('Genus','g.wsg')
    HKK.wd.family <- aggregate.data.frame(HKK.wd$wsg,by=list(HKK.wd$fam),FUN='mean',na.rm=T)
      colnames(HKK.wd.family) <- c('Family','f.wsg')
  for(i in 1:length(HKK.cens)){
    HKK.cens[[i]] <- merge(HKK.cens[[i]],HKK.sp,by='Mnemonic',all.x=T,all.y=F)
    HKK.cens[[i]] <- merge(HKK.cens[[i]],HKK.wd.sp,by='Mnemonic',all.x=T,all.y=F)
    HKK.cens[[i]] <- merge(HKK.cens[[i]],HKK.wd.genus,by='Genus',all.x=T,all.y=F)
    HKK.cens[[i]] <- merge(HKK.cens[[i]],HKK.wd.family,by='Family',all.x=T,all.y=F)
    HKK.cens[[i]]$WSG <- ifelse(is.na(HKK.cens[[i]]$WSG),HKK.cens[[i]]$g.wsg,HKK.cens[[i]]$WSG)
    HKK.cens[[i]]$WSG <- ifelse(is.na(HKK.cens[[i]]$WSG),HKK.cens[[i]]$f.wsg,HKK.cens[[i]]$WSG)
    HKK.cens[[i]]$WSG <- ifelse(is.na(HKK.cens[[i]]$WSG),
                                weighted.mean(HKK.cens[[i]]$WSG, w=HKK.cens[[i]]$dbh, na.rm=T),
                                HKK.cens[[i]]$WSG)
    HKK.cens[[i]] <- HKK.cens[[i]][,taper.cols]
  }

#### 2.5 Khao Chong ####
  # Species list
  KCH.sp <- droplevels(taper.sp.list[taper.sp.list$Site=='KCH',c('Mnemonic','Genus','Family')])
  # Wood density file
  KCH.wd <- droplevels(wsg.ctfs2[wsg.ctfs2$site=="Khao Chong",])
    KCH.wd.sp <- aggregate.data.frame(KCH.wd$wsg,by=list(KCH.wd$sp),FUN='mean',na.rm=T)
     colnames(KCH.wd.sp) <- c('Mnemonic','WSG')
    KCH.wd.genus <- aggregate.data.frame(KCH.wd$wsg,by=list(KCH.wd$genus),FUN='mean',na.rm=T)
      colnames(KCH.wd.genus) <- c('Genus','g.wsg')
    KCH.wd.family <- aggregate.data.frame(KCH.wd$wsg,by=list(KCH.wd$fam),FUN='mean',na.rm=T)
      colnames(KCH.wd.family) <- c('Family','f.wsg')
  for(i in 1:length(KCH.cens)){
    KCH.cens[[i]] <- merge(KCH.cens[[i]],KCH.sp,by='Mnemonic',all.x=T,all.y=F)
    KCH.cens[[i]] <- merge(KCH.cens[[i]],KCH.wd.sp,by='Mnemonic',all.x=T,all.y=F)
    KCH.cens[[i]] <- merge(KCH.cens[[i]],KCH.wd.genus,by='Genus',all.x=T,all.y=F)
    KCH.cens[[i]] <- merge(KCH.cens[[i]],KCH.wd.family,by='Family',all.x=T,all.y=F)
    KCH.cens[[i]]$WSG <- ifelse(is.na(KCH.cens[[i]]$WSG),KCH.cens[[i]]$g.wsg,KCH.cens[[i]]$WSG)
    KCH.cens[[i]]$WSG <- ifelse(is.na(KCH.cens[[i]]$WSG),KCH.cens[[i]]$f.wsg,KCH.cens[[i]]$WSG)
    KCH.cens[[i]]$WSG <- ifelse(is.na(KCH.cens[[i]]$WSG),
                                weighted.mean(KCH.cens[[i]]$WSG, w=KCH.cens[[i]]$dbh, na.rm=T),
                                KCH.cens[[i]]$WSG)
    KCH.cens[[i]] <- KCH.cens[[i]][,taper.cols]
  }
      
      
#### 2.6 Add E value to all entries for Chave et al. 2014 AGB allometry ####
    
    # NOTE: This may not be necessary if AGB anlyses are totally dropped    
      
    load("~/Desktop/Taper/Current/RawCensusData_CTFS/Chave.Evals.RData")
  
    # Define function (agb.allometry) to estimate aboveground biomass from wood specific 
    # gravity (wsg, in units of g/(cm^3)), tree diameter (dbh, in units of cm) from Chave et al. 2014.     
    agb.allometry <- function(E,wsg,dbh) {exp(-1.803-0.976*E+0.976*log(wsg) + 2.673*log(dbh) - 0.0299*(log(dbh)^2))}
    
    for(i in 1:length(AMA.cens)){
      AMA.cens[[i]]$E <- E.AMA
      AMA.cens[[i]]$AGB <- agb.allometry(AMA.cens[[i]]$E, AMA.cens[[i]]$WSG, AMA.cens[[i]]$dbh)
    }
    for(i in 1:length(BCI.cens)){
      BCI.cens[[i]]$E <- E.BCI
      BCI.cens[[i]]$AGB <- agb.allometry(BCI.cens[[i]]$E, BCI.cens[[i]]$WSG, BCI.cens[[i]]$dbh)
    }
    for(i in 1:length(BKT.cens)){
      BKT.cens[[i]]$E <- E.BKT
      BKT.cens[[i]]$AGB <- agb.allometry(BKT.cens[[i]]$E, BKT.cens[[i]]$WSG, BKT.cens[[i]]$dbh)
    }
    for(i in 1:length(HKK.cens)){
      HKK.cens[[i]]$E <- E.HKK
      HKK.cens[[i]]$AGB <- agb.allometry(HKK.cens[[i]]$E, HKK.cens[[i]]$WSG, HKK.cens[[i]]$dbh)
    }
    for(i in 1:length(KCH.cens)){
      KCH.cens[[i]]$E <- E.KCH
      KCH.cens[[i]]$AGB <- agb.allometry(KCH.cens[[i]]$E, KCH.cens[[i]]$WSG, KCH.cens[[i]]$dbh)
    }
    
#### 3. Make a table of stats per census for each plot ####

  # Define function to estimate the proportion of biomass measured at nonstandard heights over time
  
  # Input data required for AGB change estimation function
  # 1. CensusData = list of census data. Each element of the list is a data frame with data for a separate census. 
      # The following fields are used in the function and should be included in each census data frame:
          # "StemID" = unique identifier for each stem in the census.
          # "hom" = height of diameter measurement (m).
          # "dbh" = diameter measurement at hom in (cm). 
      # There should be no more than one measurement for each StemID and hom for each census. There can be multiple measurements if one StemID was measured at 
      # multiple heights (more than one hom). If there are multiple measurement heights, the uncorrected routine and taper correction routine uses the tallest
      # height.
  
  # Output data are a data frame:
      # 1. propDBH = data frame with the following columns:
          # "Census" = census number
          # "Prop" = proportion of stems with HOM above or below 1.3 m (+/- height threshold)
          # "PropBA" = proportion of basal area with HOM above or below 1.3 m (+/- height threshold)
    NonStdProp <- function(CensusData){
      
      # Count the number of censuses
      nCensus <- length(CensusData)
      
      # Create data frames to store results
        propDBH <- data.frame(Census=1:nCensus,
                             Prop=NA,
                             PropBA=NA)
  
      # Calculate the proportion of stems and basal area measured at non-standard heights
      for (i in 1:nCensus){
        # Only keep one (highest) diameter measurement per unique stem
          # Order by StemID (increasing), then measurement height (decreasing)
          CensusData[[i]] <- CensusData[[i]][order(CensusData[[i]]$StemID, -as.numeric(CensusData[[i]]$hom)),]
          # Only keep one measurement per StemID
          CensusData[[i]] <- CensusData[[i]][!duplicated(CensusData[[i]]$StemID),]
        
        # Proportion of stems
        propDBH$Prop[i] <- length(CensusData[[i]][!(CensusData[[i]]$hom==1.3),"hom"])/length(CensusData[[i]][,"hom"])
        
        # Proportion of basal area
        CensusData[[i]]$BasalArea <- pi*((CensusData[[i]]$dbh/2)^2)
        propDBH$PropBA[i] <- sum(CensusData[[i]][!(CensusData[[i]]$hom==1.3),"BasalArea"])/sum(CensusData[[i]][,"BasalArea"])
        
      }
      
      return(propDBH)
    }
      
    # Use function to calculate proportion of stems and basal area measured at non-standard heights for each plot
      HOM.AMA <- NonStdProp(CensusData=AMA.cens); HOM.AMA$Site <- "AMA"
      HOM.BCI <- NonStdProp(CensusData=BCI.cens); HOM.BCI$Site <- "BCI"
      HOM.BKT <- NonStdProp(CensusData=BKT.cens); HOM.BKT$Site <- "BKT"
      HOM.HKK <- NonStdProp(CensusData=HKK.cens); HOM.HKK$Site <- "HKK"
      HOM.KCH <- NonStdProp(CensusData=KCH.cens); HOM.KCH$Site <- "KCH"
      
    # Calculate the mean HOM for each plot and year
      HOM.AMA$MeanHOM <- NA
      for(i in 1:length(AMA.cens)){
        HOM.AMA$MeanHOM[i] <- mean(x = AMA.cens[[i]][!duplicated(AMA.cens[[i]]$StemID),"hom"])
      }
      HOM.BCI$MeanHOM <- NA
      for(i in 1:length(BCI.cens)){
        HOM.BCI$MeanHOM[i] <- mean(x = BCI.cens[[i]][!duplicated(BCI.cens[[i]]$StemID),"hom"])
      }
      HOM.BKT$MeanHOM <- NA
      for(i in 1:length(BKT.cens)){
        HOM.BKT$MeanHOM[i] <- mean(x = BKT.cens[[i]][!duplicated(BKT.cens[[i]]$StemID),"hom"])
      }
      HOM.HKK$MeanHOM <- NA
      for(i in 1:length(HKK.cens)){
        HOM.HKK$MeanHOM[i] <- mean(x = HKK.cens[[i]][!duplicated(HKK.cens[[i]]$StemID),"hom"])
      }
      HOM.KCH$MeanHOM <- NA
      for(i in 1:length(KCH.cens)){
        HOM.KCH$MeanHOM[i] <- mean(x = KCH.cens[[i]][!duplicated(KCH.cens[[i]]$StemID),"hom"])
      }
      
      # Calculate the basal area-weighted mean HOM for each plot and year
      HOM.AMA$MeanHOM.BA <- NA
      for(i in 1:length(AMA.cens)){
        AMA.cens[[i]]$BA <- pi*(AMA.cens[[i]]$dbh/2)^2
        HOM.AMA$MeanHOM.BA[i] <- weighted.mean(x = AMA.cens[[i]][!duplicated(AMA.cens[[i]]$StemID),"hom"],
                                            w = AMA.cens[[i]][!duplicated(AMA.cens[[i]]$StemID),"BA"])
      }
      HOM.BCI$MeanHOM.BA <- NA
      for(i in 1:length(BCI.cens)){
        BCI.cens[[i]]$BA <- pi*(BCI.cens[[i]]$dbh/2)^2
        HOM.BCI$MeanHOM.BA[i] <- weighted.mean(x = BCI.cens[[i]][!duplicated(BCI.cens[[i]]$StemID),"hom"],
                                            w = BCI.cens[[i]][!duplicated(BCI.cens[[i]]$StemID),"BA"])
      }
      HOM.BKT$MeanHOM.BA <- NA
      for(i in 1:length(BKT.cens)){
        BKT.cens[[i]]$BA <- pi*(BKT.cens[[i]]$dbh/2)^2
        HOM.BKT$MeanHOM.BA[i] <- weighted.mean(x = BKT.cens[[i]][!duplicated(BKT.cens[[i]]$StemID),"hom"],
                                            w = BKT.cens[[i]][!duplicated(BKT.cens[[i]]$StemID),"BA"])
      }
      HOM.HKK$MeanHOM.BA <- NA
      for(i in 1:length(HKK.cens)){
        HKK.cens[[i]]$BA <- pi*(HKK.cens[[i]]$dbh/2)^2
        HOM.HKK$MeanHOM.BA[i] <- weighted.mean(x = HKK.cens[[i]][!duplicated(HKK.cens[[i]]$StemID),"hom"],
                                            w = HKK.cens[[i]][!duplicated(HKK.cens[[i]]$StemID),"BA"])
      }
      HOM.KCH$MeanHOM.BA <- NA
      for(i in 1:length(KCH.cens)){
        KCH.cens[[i]]$BA <- pi*(KCH.cens[[i]]$dbh/2)^2
        HOM.KCH$MeanHOM.BA[i] <- weighted.mean(x = KCH.cens[[i]][!duplicated(KCH.cens[[i]]$StemID),"hom"],
                                            w = KCH.cens[[i]][!duplicated(KCH.cens[[i]]$StemID),"BA"])
      }
      
    # Calculate the biomass-weighted mean HOM for each plot and year
      HOM.AMA$MeanHOM.AGB <- NA
      for(i in 1:length(AMA.cens)){
        HOM.AMA$MeanHOM.AGB[i] <- weighted.mean(x = AMA.cens[[i]][!duplicated(AMA.cens[[i]]$StemID),"hom"],
                                            w = AMA.cens[[i]][!duplicated(AMA.cens[[i]]$StemID),"AGB"])
      }
      HOM.BCI$MeanHOM.AGB <- NA
      for(i in 1:length(BCI.cens)){
        HOM.BCI$MeanHOM.AGB[i] <- weighted.mean(x = BCI.cens[[i]][!duplicated(BCI.cens[[i]]$StemID),"hom"],
                                            w = BCI.cens[[i]][!duplicated(BCI.cens[[i]]$StemID),"AGB"])
      }
      HOM.BKT$MeanHOM.AGB <- NA
      for(i in 1:length(BKT.cens)){
        HOM.BKT$MeanHOM.AGB[i] <- weighted.mean(x = BKT.cens[[i]][!duplicated(BKT.cens[[i]]$StemID),"hom"],
                                            w = BKT.cens[[i]][!duplicated(BKT.cens[[i]]$StemID),"AGB"])
      }
      HOM.HKK$MeanHOM.AGB <- NA
      for(i in 1:length(HKK.cens)){
        HOM.HKK$MeanHOM.AGB[i] <- weighted.mean(x = HKK.cens[[i]][!duplicated(HKK.cens[[i]]$StemID),"hom"],
                                            w = HKK.cens[[i]][!duplicated(HKK.cens[[i]]$StemID),"AGB"])
      }
      HOM.KCH$MeanHOM.AGB <- NA
      for(i in 1:length(KCH.cens)){
        HOM.KCH$MeanHOM.AGB[i] <- weighted.mean(x = KCH.cens[[i]][!duplicated(KCH.cens[[i]]$StemID),"hom"],
                                            w = KCH.cens[[i]][!duplicated(KCH.cens[[i]]$StemID),"AGB"])
      }
      
    # Combine into a single data frame
      HOM.results <- rbind(HOM.AMA,HOM.BCI,HOM.BKT,HOM.HKK,HOM.KCH)
      
     # Add actual census years for each site
      HOM.results$Year <- NA
        # Amacayacu
        HOM.results[HOM.results$Site=="AMA","Year"] <- c(2007,2013)
        # Barro Colorado
        HOM.results[HOM.results$Site=="BCI","Year"] <- c(1990,1995,2000,2005,2010,2015)
        # Bukit Timah
        HOM.results[HOM.results$Site=="BKT","Year"] <- c(2006)
        # Huai Kha Khaeng
        HOM.results[HOM.results$Site=="HKK","Year"] <- c(1994,1999,2004,2009)
        # Khao Chong
        HOM.results[HOM.results$Site=="KCH","Year"] <- c(2000,2005,2010)

      # Write a .csv file with these results to make a figure
        write.csv(HOM.results, file = "Data file_HOMresultsPerPlot.csv", row.names = F)

#### 4. Repeat HOM analysis but only using stems greater than 30 cm DBH (comparable across plots) ####
  # Make new lists of data greater than 30 cm
      AMA.cens30 <- AMA.cens
      for(i in 1:length(AMA.cens)){
        AMA.cens30[[i]] <- AMA.cens[[i]][AMA.cens[[i]]$dbh >= 30,]
      }
      BCI.cens30 <- BCI.cens
      for(i in 1:length(BCI.cens)){
        BCI.cens30[[i]] <- BCI.cens[[i]][BCI.cens[[i]]$dbh >= 30,]
      }
      BKT.cens30 <- BKT.cens
      for(i in 1:length(BKT.cens)){
        BKT.cens30[[i]] <- BKT.cens[[i]][BKT.cens[[i]]$dbh >= 30,]
      }
      HKK.cens30 <- HKK.cens
      for(i in 1:length(HKK.cens)){
        HKK.cens30[[i]] <- HKK.cens[[i]][HKK.cens[[i]]$dbh >= 30,]
      }
      KCH.cens30 <- KCH.cens
      for(i in 1:length(KCH.cens)){
        KCH.cens30[[i]] <- KCH.cens[[i]][KCH.cens[[i]]$dbh >= 30,]
      }
       
    # Use function to calculate proportion of stems and basal area measured at non-standard heights for each plot
      HOM30.AMA <- NonStdProp(CensusData=AMA.cens30); HOM30.AMA$Site <- "AMA"
      HOM30.BCI <- NonStdProp(CensusData=BCI.cens30); HOM30.BCI$Site <- "BCI"
      HOM30.BKT <- NonStdProp(CensusData=BKT.cens30); HOM30.BKT$Site <- "BKT"
      HOM30.HKK <- NonStdProp(CensusData=HKK.cens30); HOM30.HKK$Site <- "HKK"
      HOM30.KCH <- NonStdProp(CensusData=KCH.cens30); HOM30.KCH$Site <- "KCH"
      
    # Calculate the biomass-weighted mean HOM for each plot and year
      HOM30.AMA$MeanHOM <- NA
      for(i in 1:length(AMA.cens30)){
        HOM30.AMA$MeanHOM[i] <- weighted.mean(x = AMA.cens30[[i]]$hom, w = AMA.cens30[[i]]$AGB)
      }
      HOM30.BCI$MeanHOM <- NA
      for(i in 1:length(BCI.cens30)){
        HOM30.BCI$MeanHOM[i] <- weighted.mean(x = BCI.cens30[[i]]$hom, w = BCI.cens30[[i]]$AGB)
      }
      HOM30.BKT$MeanHOM <- NA
      for(i in 1:length(BKT.cens30)){
        HOM30.BKT$MeanHOM[i] <- weighted.mean(x = BKT.cens30[[i]]$hom, w = BKT.cens30[[i]]$AGB)
      }
      HOM30.HKK$MeanHOM <- NA
      for(i in 1:length(HKK.cens30)){
        HOM30.HKK$MeanHOM[i] <- weighted.mean(x = HKK.cens30[[i]]$hom, w = HKK.cens30[[i]]$AGB)
      }
      HOM30.KCH$MeanHOM <- NA
      for(i in 1:length(KCH.cens30)){
        HOM30.KCH$MeanHOM[i] <- weighted.mean(x = KCH.cens30[[i]]$hom, w = KCH.cens30[[i]]$AGB)
      }
      
    # Combine into a single data frame
      HOM30.results <- rbind(HOM30.AMA,HOM30.BCI,HOM30.BKT,HOM30.HKK,HOM30.KCH)
      
     # Add actual census years for each site
      HOM30.results$Year <- NA
        # Amacayacu
        HOM30.results[HOM30.results$Site=="AMA","Year"] <- c(2007,2013)
        # Barro Colorado
        HOM30.results[HOM30.results$Site=="BCI","Year"] <- c(1990,1995,2000,2005,2010,2015)
        # Bukit Timah
        HOM30.results[HOM30.results$Site=="BKT","Year"] <- c(2006)
        # Huai Kha Khaeng
        HOM30.results[HOM30.results$Site=="HKK","Year"] <- c(1994,1999,2004,2009)
        # Khao Chong
        HOM30.results[HOM30.results$Site=="KCH","Year"] <- c(2000,2005,2010)
        
      # Write a .csv file with these results to make a figure
        write.csv(HOM30.results, file = "Data file_HOMresultsPerPlot30.csv", row.names = F)
#### 5. Statistical tests for HOM variation among plots and over time ####

  # Kruskal - Wallis test for differences among distributions of measurement heights among plots

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
  
  # Using the most recent census from each plot
  
    # With Bukit Timah (different minimum tree censused)            
    HOM.list = list(AMA.cens[[2]][!duplicated(AMA.cens[[2]]$StemID),"hom"], BCI.cens[[6]][!duplicated(BCI.cens[[6]]$StemID),"hom"], 
                    BKT.cens[[1]][!duplicated(BKT.cens[[1]]$StemID),"hom"], HKK.cens[[4]][!duplicated(HKK.cens[[4]]$StemID),"hom"],
                    KCH.cens[[3]][!duplicated(KCH.cens[[3]]$StemID),"hom"])
    KW.HOM.test = kruskal.test(HOM.list)
    
    # Without Bukit Timah
    HOM.list2 = list(AMA.cens[[2]][!duplicated(AMA.cens[[2]]$StemID),"hom"], BCI.cens[[6]][!duplicated(BCI.cens[[6]]$StemID),"hom"], 
                     HKK.cens[[4]][!duplicated(HKK.cens[[4]]$StemID),"hom"],KCH.cens[[3]][!duplicated(KCH.cens[[3]]$StemID),"hom"])
    KW.HOM.test2 = kruskal.test(HOM.list2)
  

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
  
#### 6. Statistical tests for HOM variation among plots and over time only using stems greater than 30 cm DBH ####

  # Kruskal - Wallis test for differences among distributions of measurement heights among plots

  # Order census data by StemID and then by HOM (decreasing)
      for(i in 1:length(AMA.cens30)){
        AMA.cens30[[i]] <- AMA.cens30[[i]][order(AMA.cens30[[i]]$StemID,-AMA.cens30[[i]]$hom),]
      }
      for(i in 1:length(BCI.cens30)){
        BCI.cens30[[i]] <- BCI.cens30[[i]][order(BCI.cens30[[i]]$StemID,-BCI.cens30[[i]]$hom),]
      }
      for(i in 1:length(BKT.cens30)){
        BKT.cens30[[i]] <- BKT.cens30[[i]][order(BKT.cens30[[i]]$StemID,-BKT.cens30[[i]]$hom),]
      }
      for(i in 1:length(HKK.cens30)){
        HKK.cens30[[i]] <- HKK.cens30[[i]][order(HKK.cens30[[i]]$StemID,-HKK.cens30[[i]]$hom),]
      }
      for(i in 1:length(KCH.cens30)){
        KCH.cens30[[i]] <- KCH.cens30[[i]][order(KCH.cens30[[i]]$StemID,-KCH.cens30[[i]]$hom),]
      }
  
  # Using the most recent census from each plot
  
    # With Bukit Timah (different minimum tree censused)            
    HOM.list = list(AMA.cens30[[2]][!duplicated(AMA.cens30[[2]]$StemID),"hom"], BCI.cens30[[6]][!duplicated(BCI.cens30[[6]]$StemID),"hom"], 
                    BKT.cens30[[1]][!duplicated(BKT.cens30[[1]]$StemID),"hom"], HKK.cens30[[4]][!duplicated(HKK.cens30[[4]]$StemID),"hom"],
                    KCH.cens30[[3]][!duplicated(KCH.cens30[[3]]$StemID),"hom"])
    KW.HOM.test = kruskal.test(HOM.list)


# Test for differnces over time within plots

  ## AMA
  #Kruskal-Wallis test
  AMA.HOM.test <- kruskal.test(list(AMA.cens30[[1]][!duplicated(AMA.cens30[[1]]$StemID),"hom"],AMA.cens30[[2]][!duplicated(AMA.cens30[[2]]$StemID),"hom"]))
 
  ## BCI
  #Kruskal-Wallis test
  BCI.HOM.test <- kruskal.test(list(BCI.cens30[[2]][!duplicated(BCI.cens30[[2]]$StemID),"hom"],
                                    BCI.cens30[[3]][!duplicated(BCI.cens30[[3]]$StemID),"hom"],BCI.cens30[[4]][!duplicated(BCI.cens30[[4]]$StemID),"hom"],
                                    BCI.cens30[[5]][!duplicated(BCI.cens30[[5]]$StemID),"hom"],BCI.cens30[[6]][!duplicated(BCI.cens30[[6]]$StemID),"hom"]))
  
  ## HKK
  #Kruskal-Wallis test
  HKK.HOM.test <- kruskal.test(list(HKK.cens30[[1]][!duplicated(HKK.cens30[[1]]$StemID),"hom"],HKK.cens30[[2]][!duplicated(HKK.cens30[[2]]$StemID),"hom"],
                                    HKK.cens30[[3]][!duplicated(HKK.cens30[[3]]$StemID),"hom"],HKK.cens30[[4]][!duplicated(HKK.cens30[[4]]$StemID),"hom"]))
  
  ## KCH
  #Kruskal-Wallis test
  KCH.HOM.test <- kruskal.test(list(KCH.cens30[[1]][!duplicated(KCH.cens30[[1]]$StemID),"hom"],KCH.cens30[[2]][!duplicated(KCH.cens30[[2]]$StemID),"hom"],
                                    KCH.cens30[[3]][!duplicated(KCH.cens30[[3]]$StemID),"hom"]))
  
#### 7. Save formatted census data ####
save(AMA.cens,BCI.cens,BKT.cens,HKK.cens,KCH.cens,file="ForestGEO_CensusData.RData")

#### Revision edit: Look at first interval for BCI ####
  # 1990-1995
  HOMchange <- BCI.cens[[1]]
  HOMchange <- HOMchange[!duplicated(HOMchange$StemID),]
  HOMchange <- merge(x = HOMchange, y = BCI.cens[[2]][!duplicated(BCI.cens[[2]]$StemID),c("StemID","dbh","hom")],
                     by = "StemID", all.x = T, all.y = T)
  
  HOMchange$dHOM <- HOMchange$hom.y-HOMchange$hom.x
  hist(HOMchange[HOMchange$dHOM>0,"dHOM"], breaks=seq(-12,6,0.2), col="black",border="white",
       main="1990 - 1995 HOM increases",xlab="HOM change (m)",
       xlim=c(0,6),ylim=c(0,160))
  weighted.mean(x = HOMchange$dHOM, w = HOMchange$AGB, na.rm=T)

  #1995-2000
  HOMchange <- BCI.cens[[2]]
  HOMchange <- HOMchange[!duplicated(HOMchange$StemID),]
  HOMchange <- merge(x = HOMchange, y = BCI.cens[[3]][!duplicated(BCI.cens[[3]]$StemID),c("StemID","dbh","hom")],
                     by = "StemID", all.x = T, all.y = T)
  
  HOMchange$dHOM <- HOMchange$hom.y-HOMchange$hom.x
  hist(HOMchange[HOMchange$dHOM>0,"dHOM"], breaks=seq(-12,6,0.2), col="black",border="white",
       main="1995 - 2000 HOM increases",xlab="HOM change (m)",
       xlim=c(0,6),ylim=c(0,160))
  HOMchange[HOMchange$dHOM>3 & !is.na(HOMchange$dHOM),]
  weighted.mean(x = HOMchange$dHOM, w = HOMchange$AGB, na.rm=T)
  
  #2000-2005
  HOMchange <- BCI.cens[[3]]
  HOMchange <- HOMchange[!duplicated(HOMchange$StemID),]
  HOMchange <- merge(x = HOMchange, y = BCI.cens[[4]][!duplicated(BCI.cens[[4]]$StemID),c("StemID","dbh","hom")],
                     by = "StemID", all.x = T, all.y = T)
  
  HOMchange$dHOM <- HOMchange$hom.y-HOMchange$hom.x
  hist(HOMchange[HOMchange$dHOM>0,"dHOM"], breaks=seq(-12,6,0.2), col="black",border="white",
       main="2000 - 2005 HOM increases",xlab="HOM change (m)",
       xlim=c(0,6),ylim=c(0,160))
  HOMchange[HOMchange$dHOM>3 & !is.na(HOMchange$dHOM),]
  weighted.mean(x = HOMchange$dHOM, w = HOMchange$AGB, na.rm=T)
  
  weighted.mean(x = BCI.cens[[2]][!duplicated(BCI.cens[[2]]$StemID),"hom"],
                w = BCI.cens[[2]][!duplicated(BCI.cens[[2]]$StemID),"AGB"])
  