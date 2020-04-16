########################################
# CSA Niches Analysis
# v0.1 
#
# Christine Lamanna
# c.lamanna@cgiar.org
# Date Created: 18 September 2019
# Date Edited: 25 September 2019
#########################################

#1.) Libraries and Tools Needed
#Depends on outputs from Basic Analysis.R script


require(ggplot2)
require(viridis)
require(maps)
require(raster)
require(dismo)
require(BiodiversityR)
require(SDMTools)
require(miceadds)
require(Hmisc)
require(lme4)

#2.) Input Data
  load.Rdata(filename=list.files(pattern="Prepared.RData",recursive=T,full.names=T),"Meta.Table2")
  Aggegrate.By<-c("PrName","Out.Ind","Product")
  DATA<-Meta.Table2[grepl("Solo",Pr.Class) & grepl("Product Yield",Out.Ind)]

#3.) Choose Subset of Data to Analyze
  #PracticexCrop combinations with at least 30 sites
  Niche.List <- ANALYSED.DATA$PrName[ANALYSED.DATA$Product=="Maize"&ANALYSED.DATA$Sites>=30] #11 practices
  Niche.List <- Niche.List[2:11]
  Niche.List <- Niche.List[order(Niche.List)]

  #Extract Data for Analysis
  Niche.Data <- Meta.Table2[Meta.Table2$PrName %in% Niche.List & Meta.Table2$Product=="Maize" & Meta.Table2$Out.Ind=="Product Yield",]
  Niche.Data <- Niche.Data[,c("Code","Latitude","Longitude","Site.Type","Site.ID","ID","M.Year","PrName","Product","Rep","yi")]
  
  #Visualize Available Data
  africa_map <- map_data("world", region = c("Angola","Botswana","Lesotho","Madagascar","Malawi","Mauritius",
                      "Mozambique","Namibia","Seychelles","South Africa","Swaziland", "Tanzania","Zambia", "Zimbabwe",
                      "Democratic Republic of the Congo", "Republic of Congo","South Sudan","Eritrea",
                      "Ethiopia","Guinea-Bissau","Egypt","Libya","Sudan","Chad","Tunisia","Algeria","Mali","Mauritania",
                      "Morocco","Western Sahara","Senegal","Gambia","Guinea","Ivory Coast", "Liberia", "Sierra Leone", "Ghana",
                      "Burkina Faso","Niger","Benin","Togo","Nigeria","Cameroon","Central African Republic", "Gabon",
                      "Equatorial Guinea","Uganda","Rwanda","Burundi","Kenya","Somalia","Djibouti","Lake Albert","Lake Malawi",
                      "Lake Tanganyika","Lake Victoria","Lake Kariba","Cape Verde","Comoros","Sao Tome and Principe"))

#Get High Resolution Country maps
uganda <- raster::getData('GADM', country='UGA', level=1)
uganda <- map_data("world", region=c("Uganda", "Lake Victoria"))
  
  MapThis <- as.data.frame(Niche.Data[Niche.Data$PrName==Niche.List[8],c("Longitude","Latitude","yi")])

  ggplot(africa_map, aes(x=long, y=lat, group=group))+
    geom_polygon(colour="black",fill="lightgrey")+
    coord_map("polyconic", ylim=c(-35,35), xlim=c(-25,55))+
    scale_color_viridis(option="magma")+
    geom_point(data=MapThis, aes(x=Longitude, y=Latitude, color=yi, group=NULL), size=1)+
    ggtitle("Data For Reduced Tillage")

ggplot(uganda, aes(x=long, y=lat, group=group))+
  geom_polygon(colour="black",fill="lightgrey")
  
#4) Aggregate Data up to Site Level
  ROUND<-4
  Aggegrate.By<-c("Site.ID","PrName")

  #4.1) Identify outliers
  Out.Calc<-function(Vals){
  return((Vals < quantile(Vals)[2] - 3 *  IQR(Vals)  | Vals > quantile(Vals)[4] + 3 * IQR(Vals)))
  }

  Outliers<-unlist(Niche.Data[,R:=1:nrow(Niche.Data)
                      ][,list(Outliers=list(R[Out.Calc(yi)])), by=Aggegrate.By
                        ][,Outliers])

  # 4.2) Calculate Weighted mean effect size by practice and site combination
  Weight.Group<-unique(c("Code",Aggegrate.By))

  Niche.Data.Sites<-Niche.Data[!Outliers
                    ][,N.Obs.Study:=.N,by=Weight.Group # Recalculate Weightings by study within obervations grouping 
                      ][,Weight.Study:=(Rep^2/(2*Rep))/N.Obs.Study # Recalculate Weightings
                        ][,list(Observations=.N,
                                Studies=length(unique(Code)),
                                Sites=length(unique(Site.ID)),
                                Lat=mean(Latitude),
                                Long=mean(Longitude),
                                RR=round(weighted.mean(yi,Weight.Study,na.rm=T),ROUND),
                                RR.sd=suppressWarnings(round(abs(wtd.var(yi,Weight.Study,na.rm=T))^0.5,ROUND))
                        ),
                        by=Aggegrate.By  # this 
                        ][,pc:=round(100*exp(RR)-100,ROUND-1)
                          ][,pc.sd:=round(100*exp(RR.sd)-100,ROUND-1)]
  Niche.Data.Sites$Label <- paste(Niche.Data.Sites$Site.ID, Niche.Data.Sites$PrName)


#5.) Import Spatial Data
  #5.1a) Cropland Area
  #Import Cropland Area, this will be the mask for all other layers
  library(maptools)
  data(wrld_simpl)
  cropland <- raster('./Data/Spatial Data/cpland_prob--ssa.tif/cpland_prob--SSA.tif')
  plot(wrld_simpl, add=TRUE)
  points(Niche.Data[Niche.Data$PrName=="Inorganic Fertilizer",c("Longitude","Latitude")], col='blue')
  
  #5.1b) Climate Data
  #started with AFRICLIM:
  #tbio - worldclim baseline, 5' resolution
  #mbio - CHIRPS baseline, 5' resolution
  cfiles <- list.files(path="./Data/Spatial Data/res5min/climate/",pattern=".tif", full.names = T)
  climate <- stack(cfiles)
  
  #5.1c) Physical data
  
  #5.1d) Soil data
  sfiles <- list.files(path="./Data/Spatial Data/res5min/soils/", pattern=".tif", full.names=T)
  soils <- stack(sfiles)
  
  #5.1e) Social data
  pfiles <- list.files(path="./Data/Spatial Data/res5min/social/", pattern=".tif", full.names=T)
  market <- raster("./Data/Spatial Data/TravelTimeToMarket_SSA_GeoTiff/traveltimetomarket_ssa_100k.tif")
  maize <- raster("./Data/Spatial Data/maiz_p--SSA.tif/maiz_p--SSA.tif")
  plot(market)
  plot(maize)
  social <- stack(pfiles)
  social <- crop(social, cropland)


  #5.2) Visualize Spatial Data
  #plot(predictors)
  #plot(predictors, 21)
  #plot(wrld_simpl, add=TRUE)
  #points(Niche.Data[Niche.Data$PrName=="Inorganic Fertilizer",c("Longitude","Latitude")], col='blue')

  #5.3) Extract Spatial Data and check for correlations
  #library(corrplot)
  #presvals <- extract(predictors, Niche.Data.Sites[,c("Long","Lat")])
  #presvals <- as.data.frame(presvals)
  #row.names(presvals) <- Niche.Data.Sites$Label
  #presvals <- na.omit(presvals)
  #M<-cor(presvals)
  #corrplot(M, method="color")
  
  #5.4) Mask Spatial Data to Cropland in SSA
  cropland <- raster('./Data/Spatial Data/cpland_prob--ssa.tif/cpland_prob--SSA.tif')
  plot(wrld_simpl, add=TRUE)
  points(Niche.Data[Niche.Data$PrName=="Inorganic Fertilizer",c("Longitude","Latitude")], col='blue')
  MAT <- raster("./Data/Spatial Data/tbio_wc150s/bio1_wc150s.tif")
  MATr <- resample(MAT,cropland)
  extent(MAT)
  extent(cropland)
  e <- intersect(extent(MAT), extent(cropland)) #find minimum extent
  MATcr <- crop(MATr,e)
  croplandc <- crop(cropland,e)
  masked1 <- raster::mask(MATcr,croplandc)
  masked2 <- raster::mask(masked1,croplandc,maskvalue=0)
  

  #match extent and resolution of all layers
  croplandR <- resample(cropland,climate)
  soilsR <- resample(soils,climate)
  
  #crop for SSA
  e <- extent(-25.5,60.5,-40,30)
  croplandR <- crop(croplandR,e)
  soilR <- crop(soilsR,e)
  climateC <- crop(climate,e)
  
  #Mask using croplands
  climateM <- raster::mask(climateC,croplandR) #masks on NA, to SSA countries
  soilsM <- raster::mask(soilR,croplandR)
  climateM <- raster::mask(climateM,croplandR, maskvalue=0) #masks on 0 to eliminate non-crop areas
  soilsM <- raster::mask(soilsM,croplandR, maskvalue=0)
  
  #5.6) Create stack of predictors
  predictors <- stack(climateM,soilsM)
  names(predictors) <- c("MAT","MAP","RFS","ISO","MTWM","LLDS","PET","CLAY","ROOT","SOC","PH","SAND")
  
  
  #5.7) Remove observations from North Africa (n=1), Cabo Verde (n=1), and NAs (n=3)
  Niche.Data.Sites <- Niche.Data.Sites[Niche.Data.Sites$Lat < 20,]
  Niche.Data.Sites <- Niche.Data.Sites[Niche.Data.Sites$Long > -20,]
  Niche.Data.Sites <- na.omit(Niche.Data.Sites)

  Niche.Data <- Niche.Data[Niche.Data$Lat < 20,]
  Niche.Data <- Niche.Data[Niche.Data$Long > -20,]
  Niche.Data <- na.omit(Niche.Data)

  
#6) Create Presence & Absence Datasets
  #visualizing response ratios to determine bins for abundance analysis
  hist(Niche.Data.Sites$RR)
  summary(Niche.Data.Sites$RR)
  length(Niche.Data.Sites$RR[Niche.Data.Sites$RR < 0.3])/length(Niche.Data.Sites$RR)
  length(Niche.Data.Sites$RR[Niche.Data.Sites$RR>0.75 & Niche.Data.Sites$RR<1])/length(Niche.Data.Sites$RR)
  hist(Niche.Data.Sites$RR[Niche.Data.Sites$PrName=="Inorganic Fertilizer"])
 
  hist(Niche.Data$yi)
  summary(Niche.Data$yi)
  length(Niche.Data$yi[Niche.Data$yi < 0.3])/length(Niche.Data$yi)
  length(Niche.Data$yi[Niche.Data$yi> 0.75 & Niche.Data$yi<1])/length(Niche.Data$yi)
  hist(Niche.Data$yi[Niche.Data$PrName=="Inorganic Fertilizer"])
  
  
   
  #6.1a)
  #Create Two binary Suitability variables (Presence/Absence)
  #Threshold for Suitability is a positive RR, RR >0
  Niche.Data.Sites$Suit0 <- TRUE
  Niche.Data.Sites$Suit0[Niche.Data.Sites$RR<0] <- FALSE
  #Threshold for Suitability is RR>0.15, at least 16% increase in yield
  Niche.Data.Sites$Suit15 <- TRUE
  Niche.Data.Sites$Suit15[Niche.Data.Sites$RR<0.15] <- FALSE
 
  #Threshold for Suitability is a positive RR, RR >0
  Niche.Data$Suit0 <- TRUE
  Niche.Data$Suit0[Niche.Data$yi<0] <- FALSE
  #Threshold for Suitability is RR>0.15, at least 16% increase in yield
  Niche.Data$Suit15 <- TRUE
  Niche.Data$Suit15[Niche.Data$yi<0.15] <- FALSE
  
   
  #6.1b)
  #alternatively assign each RR to an "abudance" class and binary present/absent
  # Create a named vector that reflects the class midpoints 
  a<-c("-0.125", "0.125", "0.375","0.625","1")
  b<-as.list(a)
  names(b)<-c("0", "1", "2","3","4")
  
  #Abundance is a factor, representing the effect size bin the observation belongs to
  #using 5 classes with relatively even distribution across all the data
  #Present is a binary factor, 1 = positive effect, 0 = negative effect
  for (i in (1:length(Niche.Data.Sites$Site.ID))){
    ifelse (Niche.Data.Sites$RR[i] < 0, Niche.Data.Sites$Abundance[i] <- "0",
            ifelse(Niche.Data.Sites$RR[i] < 0.25, Niche.Data.Sites$Abundance[i] <- "1",
                   ifelse(Niche.Data.Sites$RR[i] < 0.5, Niche.Data.Sites$Abundance[i] <- "2",
                          ifelse(Niche.Data.Sites$RR[i] <0.75, Niche.Data.Sites$Abundance[i] <- "3",
                                 ifelse(Niche.Data.Sites$RR[i]>0.75, Niche.Data.Sites$Abundance[i] <- "4")))))

  }
  Niche.Data.Sites$Abundance <- as.ordered(Niche.Data.Sites$Abundance)
  Niche.Data.Sites$Presence <- 1
  Niche.Data.Sites$Presence[Niche.Data.Sites$Abundance=="0"] <- 0
  Niche.Data.Sites$Presence <- as.factor(Niche.Data.Sites$Presence)
  
  # A vector of weights based on their individual abundance classification
    d<-ifelse(Niche.Data.Sites$Abundance=="0", 1/(unname(table(Niche.Data.Sites$Abundance)["0"])/length(Niche.Data.Sites$Abundance)),
            ifelse(Niche.Data.Sites$Abundance=="1", 1/(unname(table(Niche.Data.Sites$Abundance)["1"])/length(Niche.Data.Sites$Abundance)),
                   ifelse(Niche.Data.Sites$Abundance=="2", 1/(unname(table(Niche.Data.Sites$Abundance)["2"])/length(Niche.Data.Sites$Abundance)),
                          ifelse(Niche.Data.Sites$Abundance=="3", 1/(unname(table(Niche.Data.Sites$Abundance)["3"])/length(Niche.Data.Sites$Abundance)),
                                 ifelse(Niche.Data.Sites$Abundance=="4", 1/(unname(table(Niche.Data.Sites$Abundance)["4"])/length(Niche.Data.Sites$Abundance)),0)))))

#7) Extract Predictor Data for Suitable & Unsuitable Points    
    #7.0) Download extracted data for site/year combinations
    load(filename="POWER.Annual.RData") #actual temperature data
    load(filename="CHIRPS Annual.RData") #actual rainfall data
    Soils.Data <- read.csv("SoilGrids Data.csv",header=T) #actual soil data (higher accuracy)
    
    #7.1) Split M.Year in Niche.Data and assign to first year, ignoring seasons for now
    x <- strsplit(Niche.Data$M.Year,"[.]")
    a <- unlist(lapply(x, FUN = function(y){
      y[1]
    }))
    Niche.Data$Year <- a
    
    #7.2) Append climate and soil data
    Niche.Data$Key <- paste(Niche.Data$ID,"-",Niche.Data$Year)
    POWER.Annual$Key <- paste(POWER.Annual$Site.Key,"-",POWER.Annual$Year)
    CHIRPS.Annual$Key <- paste(CHIRPS.Annual$Site.Key,"-",CHIRPS.Annual$YEAR)
    
    
    Niche.Data$MAP <- CHIRPS.Annual$Total.Rain[match(Niche.Data$Key,CHIRPS.Annual$Key)] #actual annual precip to match with MAP
    Niche.Data$MAT <- POWER.Annual$Temp.Mean.Mean[match(Niche.Data$Key,POWER.Annual$Key)]*10 #actual annual temp to match with MAT AFRICLIM
    Niche.Data$MTWM <- POWER.Annual$Temp.Max[match(Niche.Data$Key,POWER.Annual$Key)]*10 #maximum temp to match MTWM
    Niche.Data$SOC <- Soils.Data$ORCDRC_M_sl2_250m[match(Niche.Data$ID,Soils.Data$Site.Key)] #accurate SOC to match SOC
    Niche.Data$PH <- Soils.Data$PHIHOX_M_sl2_250m[match(Niche.Data$ID,Soils.Data$Site.Key)] #accurate PH to match PH
    Niche.Data$CLAY <- Soils.Data$CLYPPT_M_sl2_250m[match(Niche.Data$ID,Soils.Data$Site.Key)] #accurate CLAY % to match CLAY
    Niche.Data$SAND <- Soils.Data$SNDPPT_M_sl2_250m[match(Niche.Data$ID,Soils.Data$Site.Key)] #accurate SAND % to match SAND
    Niche.Data$RFS <- extract(predictors[["RFS"]],Niche.Data[,c("Longitude","Latitude")])
    Niche.Data$ISO <- extract(predictors[["ISO"]],Niche.Data[,c("Longitude","Latitude")])
    Niche.Data$LLDS <- extract(predictors[["LLDS"]],Niche.Data[,c("Longitude","Latitude")])
    Niche.Data$PET <- extract(predictors[["PET"]],Niche.Data[,c("Longitude","Latitude")])
    Niche.Data$ROOT <- extract(predictors[["ROOT"]],Niche.Data[,c("Longitude","Latitude")])
    
    
#8) Develop Suitability Models using SDM Methods 
    
    #8.0) Run SDM with SDMtools/dismo and random background points
      #8.0a) Create random points on the continent
#      bg <- as.data.frame(randomPoints(predictors, 500))
#      colnames(bg) <- c("Long","Lat")
#      plot(predictors[[1]])
#      plot(wrld_simpl, add=TRUE)
#      points(Niche.Data[Niche.Data$PrName=="Tree Management",c("Longitude","Latitude")], col='blue')
#      points(bg, col='black', size=1)
    
#      #8.0b) Extract predictors for presense & absense values
#      presvals <- extract(predictors,Niche.Data.Sites[Niche.Data.Sites$PrName=="Tree Management",c("Long","Lat")])
#      absvals <- extract(predictors, bg)
#      pb <- c(rep(1,nrow(presvals)),rep(0,nrow(absvals)))
#      sdmdata<- data.frame(cbind(pb, rbind(presvals,absvals)))
    
#      #8.0c) Run SDMs
#      m.glm = glm(pb ~ ., data=sdmdata)
    
#      #8.0d) Predict & Plot Suitabilty
#      p <- predict(predictors, m.glm) #takes ~10seconds at 10' resolution
#      plot(p)
#      plot(wrld_simpl, add=T)

    #8.1) with Absence data & dismo methods
      #8.1a) extract predictors for suitable and unsuitable points
      suitable <- Niche.Data[Niche.Data$PrName=="Reduced Tillage" & Niche.Data$Suit15==TRUE,
                             c("MAT","MAP","RFS","ISO","MTWM","LLDS","PET","CLAY","ROOT","SOC","PH","SAND")]
      unsuitable <- Niche.Data[Niche.Data$PrName=="Reduced Tillage" & Niche.Data$Suit15==FALSE,
                               c("MAT","MAP","RFS","ISO","MTWM","LLDS","PET","CLAY","ROOT","SOC","PH","SAND")]
      st <- c(rep(1,nrow(suitable)),rep(0,nrow(unsuitable)))
      sdmdata.rt<- data.frame(cbind(st, rbind(suitable,unsuitable)))
      sdmdata.rt<-na.omit(sdmdata.rt)
      
      #8.1b) run the model
      #binomial regression
      m.glm.if3 = glm(st ~ ., family=binomial(link = "logit"), data=sdmdata.if)
      summary(m.glm.if3)
      
      m.glm.rt3 = glm(st ~ ., family=binomial(link = "logit"), data=sdmdata.rt)
      summary(m.glm.rt3)
      m.glm.af3 = glm(st ~ ., family=binomial(link = "logit"), data=sdmdata.af)
      summary(m.glm.af3)
      
      m.glm.ic3 = glm(st ~ ., family=binomial(link = "logit"), data=sdmdata.ic)
      summary(m.glm.ic3)
      
      m.glm.cr3 = glm(st ~ ., family=binomial(link = "logit"), data=sdmdata.cr)
      summary(m.glm.cr3)
      
      m.glm.cres3 = glm(st ~ ., family=binomial(link = "logit"), data=sdmdata.cres)
      summary(m.glm.cres3)
      
      m.glm.wh3 = glm(st ~ ., family=binomial(link = "logit"), data=sdmdata.wh)
      summary(m.glm.wh3)
      
      #8.1c) Binomial Regression with Random Effects
      suitable <- Niche.Data[Niche.Data$PrName=="Water Harvesting" & Niche.Data$Suit15==TRUE,
                             c("Code","ID","MAT","MAP","RFS","ISO","MTWM","LLDS","PET","CLAY","ROOT","SOC","PH","SAND")]
      unsuitable <- Niche.Data[Niche.Data$PrName=="Water Harvesting" & Niche.Data$Suit15==FALSE,
                               c("Code","ID","MAT","MAP","RFS","ISO","MTWM","LLDS","PET","CLAY","ROOT","SOC","PH","SAND")]
      st <- c(rep(1,nrow(suitable)),rep(0,nrow(unsuitable)))
      sdmdata.whr<- data.frame(cbind(st, rbind(suitable,unsuitable)))
      sdmdata.whr<-na.omit(sdmdata.whr)
      
      require(lmerTest)
      require(MuMIn)
      require(lme4)
      control = lmerControl(optimizer = "optim", calc.derivs = F)
      
      m.glmer.if <- glmer(st ~ MAT + MAP * CLAY + RFS + ISO + MTWM + LLDS + PET + CLAY + ROOT + SOC + PH + SAND + (1|ID/Code),
                           binomial, data=sdmdata.ifr)
      r.squaredGLMM(m.glmer.if)
      anova(m.glmer.if,test="LR")
      
      m.glmer.rt <- glmer(st ~ MAT + MAP * CLAY + RFS + ISO + MTWM + LLDS + PET + CLAY + ROOT + SOC + PH + SAND + (1|ID/Code),
                           binomial, data=sdmdata.rtr)
      
      m.glmer.res <- glmer(st ~ MAT + MAP * CLAY + RFS + ISO + MTWM + LLDS + PET + CLAY + ROOT + SOC + PH + SAND + (1|ID/Code),
                          binomial, data=sdmdata.resr)

      m.glmer.rot <- glmer(st ~ MAT + MAP * CLAY + RFS + ISO + MTWM + LLDS + PET + CLAY + ROOT + SOC + PH + SAND + (1|ID/Code),
                           binomial, data=sdmdata.rotr)
      
      m.glmer.gm <- glmer(st ~ MAT + MAP * CLAY + RFS + ISO + MTWM + LLDS + PET + CLAY + ROOT + SOC + PH + SAND + (1|ID/Code),
                           binomial, data=sdmdata.gmr)
      
      m.glmer.ic <- glmer(st ~ MAT + MAP * CLAY + RFS + ISO + MTWM + LLDS + PET + CLAY + ROOT + SOC + PH + SAND + (1|ID/Code),
                           binomial, data=sdmdata.icr)
      
      m.glmer.ml <- glmer(st ~ MAT + MAP * CLAY + RFS + ISO + MTWM + LLDS + PET + CLAY + ROOT + SOC + PH + SAND + (1|ID/Code),
                           binomial, data=sdmdata.mlr)
      
      m.glmer.of <- glmer(st ~ MAT + MAP * CLAY + RFS + ISO + MTWM + LLDS + PET + CLAY + ROOT + SOC + PH + SAND + (1|ID/Code),
                           binomial, data=sdmdata.ofr)
      
      m.glmer.af <- glmer(st ~ MAT + MAP * CLAY + RFS + ISO + MTWM + LLDS + PET + CLAY + ROOT + SOC + PH + SAND + (1|ID/Code),
                           binomial, data=sdmdata.tmr)
      
      m.glmer.wh <- glmer(st ~ MAT + MAP * CLAY + RFS + ISO + MTWM + LLDS + PET + CLAY + ROOT + SOC + PH + SAND + (1|ID/Code),
                           binomial, data=sdmdata.whr)
      
            
#      anova(m.glmer.if4,m.glmer.if4a,test="LR")
#      options(na.action='na.fail' )
#      dredge<-MuMIn::dredge(m.glmer.if4a )
      
      #8.1c) Predict & Plot Suitabilty
      p.if <- predict(predictors, m.glmer.if, re.form=NA) 
      p.rt <- predict(predictors, m.glmer.rt, re.form=NA)
      p.res <- predict(predictors, m.glmer.res, re.form=NA)
      p.rot <- predict(predictors, m.glmer.rot, re.form=NA)
      p.gm <- predict(predictors, m.glmer.gm, re.form=NA)
      p.ic <- predict(predictors, m.glmer.ic, re.form=NA)
      p.ml <- predict(predictors, m.glmer.ml, re.form=NA)
      p.of <- predict(predictors, m.glmer.of, re.form=NA)
      p.af <- predict(predictors, m.glmer.af, re.form=NA)
      p.wh <- predict(predictors, m.glmer.wh, re.form=NA)
      
      
      colr <- colorRampPalette(brewer.pal(11, 'RdYlBu'))
      plot(p.if, col=viridis(100), main="Inorganic Fertilizer Suitability")
      plot(p.rt, col=viridis(100),main="Reduced Tillage Suitability")
      plot(p.res, col=viridis(100),main="Crop Residue Suitability")
      plot(p.rot, col=viridis(100),main="Crop Rotation Suitability")
      plot(p.gm, col=viridis(100),main="Green Manure Suitability")
      plot(p.ic, col=viridis(100),main="Intercropping Suitability")
      plot(p.ml, col=viridis(100),main="Mulching Suitability")
      plot(p.of, col=viridis(100),main="Organic Fertilizer Suitability")
      plot(p.af, col=viridis(100),main="Agroforestry Suitability")
      plot(p.wh, col=viridis(100),main="Water Harvesting Suitability")
      plot(wrld_simpl, add=T)
      
      s <- stack(p.res, p.rot, p.gm, p.if, p.ic, p.ml,p.of,p.rt,p.af,p.wh)
      names(s) <- Niche.List
#      s[s>2.94] <- 2.94 
#      s[s<-2.94] <- -2.94
      s <- exp(s)/(1+exp(s)) #Change odds to chance suitable

      
      p <- gplot(s) + geom_tile(aes(fill = value)) +
        facet_wrap(~ variable) +
        scale_fill_divergent(low='tomato', 
                             mid = 'lemonchiffon', 
                             high = 'darkcyan', 
                             midpoint=0.5, 
                             na.value="white"
                             ) +
        coord_equal()+
        labs(fill="Chance Suitable")+
        geom_polygon(data=africa_map, 
                     aes(x=long, y=lat, group=group), 
                     color="grey20", size=0.25, fill=NA)+
        theme(
          axis.line=element_blank(),axis.text.x=element_blank(),
          axis.text.y=element_blank(),axis.ticks=element_blank(),
          axis.title.x=element_blank(),axis.title.y=element_blank(),
          panel.background=element_blank(),panel.border=element_blank(),
          panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
          plot.background=element_blank(),
          strip.background = element_rect(color="white", fill="white"),
          strip.text = element_text(size=12, face = "bold")
          )
    p
      
    ggsave("CurrentChanceSuitability-All.png", device="png")


    
    
    #8.2) with CART random forest + abundance classes
    library(party)
    
      #8.2a) extract raster values for each observation
      names(predictors)<- c("MAT","MAP","RFS","ISOTHERM","MTWM","LLDS","PET","CLAY","ROOT","SOC","PH","SAND")
      rastervals <- extract(predictors,Niche.Data.Sites[,c("Long","Lat")])
      Niche.Data.Sites.Ext <- cbind(Niche.Data.Sites,rastervals,d)
    
      #8.2b) run conditional random forest
      training.data <- as.data.frame(Niche.Data.Sites.Ext[Niche.Data.Sites.Ext$PrName=="Intercropping",])
      abundance.cforest<-cforest(Abundance~.,data=training.data[,c("MAT","MAP","RFS","ISOTHERM","MTWM","LLDS","PET","CLAY","SAND","ROOT","SOC","PH","Abundance")],
                               weights=training.data$d,scores=b,controls=cforest_unbiased(ntree=1000, mtry=5))
    
      #8.2c)predict suitability using random forest model
      testing.data <- as.data.frame(Niche.Data.Sites.Ext[Niche.Data.Sites.Ext$PrName=="Inorganic Fertilizer",]) 
      # Use model to make probabilistic predictions
      probability.pred<-predict(abundance.cforest, newdata=testing.data, type="prob") 
      probability.pred <- data.frame(matrix(unlist(probability.pred), nrow=dim(testing.data)[1], byrow=T))
    
      # Probability of occurrence is the sum of the non-zero predictions
      g<-dim(probability.pred)[2]
      ifelse(g<=2, probability.pred$presence<-(probability.pred[,2]),
           ifelse(g>2, probability.pred$presence<-rowSums(probability.pred[,c(2:g)], 0)))
      # These values of probability of occurence can then be used to calculate AUC and calibration curves
    
      # Make predictions on the scale of the response variable, in this case the categorical abundance scale
      response.pred<-predict(abundance.cforest, newdata=testing.data, type="response") 
      response.pred<-as.numeric(as.character(response.pred))
      
    #8.3) Run SDM with BiodiversityR ensemble methods
      require(BiodiversityR)
#      forSDM <- c("Inorganic Fertilizer","Intercropping","Tree Management","Reduced Tillage")
      forSDM <- c("Inorganic Fertilizer")
      all.suit <- Niche.Data.Sites[Niche.Data.Sites$Suit15==TRUE & Niche.Data.Sites$PrName %in% forSDM,c("PrName","Long","Lat")]
      all.unsuit <- Niche.Data.Sites[Niche.Data.Sites$Suit15==FALSE & Niche.Data.Sites$PrName %in% forSDM,c("PrName","Long","Lat")]
      all.suit <- na.omit(all.suit)
      all.unsuit <- na.omit(all.unsuit)
      esm1.if <- ensemble.batch(x=predictors, xn=predictors, species.presence = all.suit, species.absence = all.unsuit,
                     BIOCLIM=0, DOMAIN=0, MAHAL=0, MAHAL01=0,
                     MAXENT=0)
      #for inorganic fertilizer it's only successfully running GLM with correct input data...

      #plot esemble models
      pres <- all.suit[,c("Long","Lat")]
      abs <- all.unsuit[,c("Long","Lat")]
      plot1 <- ensemble.plot(RASTER.species.name="Inorganic Fertilizer", RASTER.stack.name=predictors, p=pres, a=abs)
      
  
  
  
#9.) Create a Raster Stack of Future Climate Layers
    #9.1) Future Climate Data
      #AFRICLIM:
      #tbio - worldclim baseline, 5' resolution, RCP8.5, 2085
      #mbio - CHIRPS baseline, 5' resolution, RCP8.5, 2085
      fcfiles <- list.files(path="./Data/Spatial Data/Future Climate Data/RCP45/for model/",pattern=".tif", full.names = T)
      rcp45 <- stack(fcfiles)
      fcfiles <- list.files(path="./Data/Spatial Data/Future Climate Data/RCP85/for model/",pattern=".tif", full.names = T)
      rcp85 <- stack(fcfiles)
      
      #crop for SSA
      e <- extent(-25.5,60.5,-40,30)
      rcp45C <- crop(rcp45,e)
      rcp85C <- crop(rcp85,e)
      
      #Mask using croplands
      rcp45M <- raster::mask(rcp45C,croplandR) #masks on NA, to SSA countries
      rcp45M <- raster::mask(rcp45M,croplandR, maskvalue=0) #masks on 0 to eliminate non-crop areas
      rcp85M <- raster::mask(rcp85C,croplandR) #masks on NA, to SSA countries
      rcp85M <- raster::mask(rcp85M,croplandR, maskvalue=0) #masks on 0 to eliminate non-crop areas
      
      #5.2) Create stack of predictors for future scenarios
      #using current soils data and future climate data
      rcp45all <- stack(rcp45M,soilsM)
      names(rcp45all) <- c("MAT","MAP","RFS","ISO","MTWM","LLDS","PET","CLAY","ROOT","SOC","PH","SAND")
      rcp85all <- stack(rcp85M,soilsM)
      names(rcp85all) <- c("MAT","MAP","RFS","ISO","MTWM","LLDS","PET","CLAY","ROOT","SOC","PH","SAND")
      
      #5.3) Run models on future climate stack
      #5.3c) rcp4.5
      p.if45 <- predict(rcp45all, m.glmer.if, re.form=NA) 
      p.rt45 <- predict(rcp45all, m.glmer.rt, re.form=NA)
      p.res45 <- predict(rcp45all, m.glmer.res, re.form=NA)
      p.rot45 <- predict(rcp45all, m.glmer.rot, re.form=NA)
      p.gm45 <- predict(rcp45all, m.glmer.gm, re.form=NA)
      p.ic45 <- predict(rcp45all, m.glmer.ic, re.form=NA)
      p.ml45 <- predict(rcp45all, m.glmer.ml, re.form=NA)
      p.of45 <- predict(rcp45all, m.glmer.of, re.form=NA)
      p.af45 <- predict(rcp45all, m.glmer.af, re.form=NA)
      p.wh45 <- predict(rcp45all, m.glmer.wh, re.form=NA)
      s45 <- stack(p.res45, p.rot45, p.gm45, p.if45, p.ic45, p.ml45,p.of45,p.rt45,p.af45,p.wh45)
      names(s45) <- Niche.List
      s45 <- exp(s45)/(1+exp(s45)) #Change odds to chance suitable
      s45d <- s45-s #Change in suitability in 2085 vs. present, RCP4.5
      
      #5.3d) rcp85
      p.if85 <- predict(rcp85all, m.glmer.if, re.form=NA) 
      p.rt85 <- predict(rcp85all, m.glmer.rt, re.form=NA)
      p.res85 <- predict(rcp85all, m.glmer.res, re.form=NA)
      p.rot85 <- predict(rcp85all, m.glmer.rot, re.form=NA)
      p.gm85 <- predict(rcp85all, m.glmer.gm, re.form=NA)
      p.ic85 <- predict(rcp85all, m.glmer.ic, re.form=NA)
      p.ml85 <- predict(rcp85all, m.glmer.ml, re.form=NA)
      p.of85 <- predict(rcp85all, m.glmer.of, re.form=NA)
      p.af85 <- predict(rcp85all, m.glmer.af, re.form=NA)
      p.wh85 <- predict(rcp85all, m.glmer.wh, re.form=NA)
      s85 <- stack(p.res85, p.rot85, p.gm85, p.if85, p.ic85, p.ml85,p.of85,p.rt85,p.af85,p.wh85)
      names(s85) <- Niche.List
      s85 <- exp(s85)/(1+exp(s85)) #Change odds to chance suitable
      s85d <- s85-s #Change in suitability in 2085 vs. present, RCP4.5
      
      
      
            
      #5.3e) Plot Change in Suitability under RCP 8.5
      p85 <- gplot(s85) + geom_tile(aes(fill = value)) +
        facet_wrap(~ variable) +
        scale_fill_divergent(low='tomato', 
                             mid = 'lemonchiffon', 
                             high = 'darkcyan', 
                             midpoint=0, 
                             na.value="white"
        ) +
        coord_equal()+
        labs(fill="Change in Suitablility")+
        labs(title="Change in Suitability 2085 vs. 2000, RCP 8.5")+
        geom_polygon(data=africa_map, 
                     aes(x=long, y=lat, group=group), 
                     color="grey20", size=0.25, fill=NA)+
        theme(
          axis.line=element_blank(),axis.text.x=element_blank(),
          axis.text.y=element_blank(),axis.ticks=element_blank(),
          axis.title.x=element_blank(),axis.title.y=element_blank(),
          panel.background=element_blank(),panel.border=element_blank(),
          panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
          plot.background=element_blank(),
          strip.background = element_rect(color="white", fill="white"),
          strip.text = element_text(size=12, face = "bold")
        )
      p85
      
      ggsave("ChangeInSuitability-All-RCP85.png", device="png")
      
      #5.3e) Plot Change in Suitability under RCP 8.5
      p45 <- gplot(s45d) + geom_tile(aes(fill = value)) +
        facet_wrap(~ variable) +
        scale_fill_divergent(low='tomato', 
                             mid = 'lemonchiffon', 
                             high = 'darkcyan', 
                             midpoint=0, 
                             na.value="white"
        ) +
        coord_equal()+
        labs(fill="Change in Suitablility")+
        labs(title="Change in Suitability 2085 vs. 2000, RCP 4.5")+
        geom_polygon(data=africa_map, 
                     aes(x=long, y=lat, group=group), 
                     color="grey20", size=0.25, fill=NA)+
        theme(
          axis.line=element_blank(),axis.text.x=element_blank(),
          axis.text.y=element_blank(),axis.ticks=element_blank(),
          axis.title.x=element_blank(),axis.title.y=element_blank(),
          panel.background=element_blank(),panel.border=element_blank(),
          panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
          plot.background=element_blank(),
          strip.background = element_rect(color="white", fill="white"),
          strip.text = element_text(size=12, face = "bold")
        )
      p45
      
      ggsave("ChangeInSuitability-All-RCP45.png", device="png")
      
      
      
      
      
      
      
      
#10.) map for Uganda
      uganda <- map_data("world", region=c("Uganda", "Lake Victoria"))
      uganda <- getData("GADM", country='UGA', level=1)
      UGA <- gadm
      p.res <- predict(predictors, m.glmer.res, re.form=NA)
      s.res <- exp(p.res)/(1+exp(p.res))
      
      test <- gplot(s.res) + geom_tile(aes(fill = value)) +
        facet_wrap(~ variable) +
        scale_fill_divergent(low='tomato', 
                             mid = 'lemonchiffon', 
                             high = 'darkcyan', 
                             midpoint=0, 
                             na.value="white"
        ) +
        coord_map("polyconic", ylim=c(-5,20), xlim=c(15,35))+
        labs(fill="Suitablility")+
        labs(title="Suitability of Crop Residue")+
        geom_polygon(data=uganda, 
                     aes(x=long, y=lat, group=group), 
                     color="grey20", size=0.25, fill=NA)+
        theme(
          axis.line=element_blank(),axis.text.x=element_blank(),
          axis.text.y=element_blank(),axis.ticks=element_blank(),
          axis.title.x=element_blank(),axis.title.y=element_blank(),
          panel.background=element_blank(),panel.border=element_blank(),
          panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
          plot.background=element_blank(),
          strip.background = element_rect(color="white", fill="white"),
          strip.text = element_text(size=12, face = "bold")
        )
      
      
      plot(s.res,
           xlim=c(28,36), ylim=c(-2,5))
      plot(UGA, axes=TRUE, col=NA, add=TRUE)
      
  

