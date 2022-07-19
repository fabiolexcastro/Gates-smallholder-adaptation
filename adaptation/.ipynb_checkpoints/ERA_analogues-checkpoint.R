#JRV - ERA analyses for adaptation atlas

#load libraries & create functions ----
require(analogues)
require(rgdal)
require(raster)
require(maptools); data(wrld_simpl)
require(data.table)
require(Hmisc)
require(tidyverse)
require(ggplot2)
require(metR)

require(rworldmap)
require(doSNOW)
require(miceadds)
require(parallel)
require(foreach)

#analysis version
vr <- 4
# Set number of cores for parallel processing ====
Cores<-15
# Run full or streamlined analysis?
DoLite<-T


#set directories - MAKE SURE YOU START CONSOLE FROM ANALOGUES FOLDER-----
wd <- "~/work/ONECGIAR/Atlas_MVP/adaptation_options"
if(!dir.exists(wd)){dir.create(wd)}
cimdir <- paste(wd,"/ERA_analogues",sep="")
if(!dir.exists(cimdir)){dir.create(cimdir)}

#get Africa shapefile -----
sh_ctry <- readOGR("~/work/ONECGIAR/Data/Africa_shp/African_continet.shp")
sh_ctry <- spTransform(sh_ctry, crs("+proj=longlat +ellps=WGS84 +no_defs"))
sh_xt <- extent(sh_ctry)
sh_xt@xmin <- sh_xt@xmin-1; sh_xt@ymin <- sh_xt@ymin-1; sh_xt@xmax <- sh_xt@xmax+1; sh_xt@ymax <- sh_xt@ymax+1

#load Soils Data -----
if(DoLite){
  load(paste(cimdir,"/input_data/soils2.rda",sep=""))
}else{
  load(paste(cimdir,"/input_data/soils.rda",sep=""))
}

# Create a Mask- ---
load(paste(cimdir,"/input_data/wc_prec.rda",sep=""))
msk <- wc_prec[[1]]; rm(wc_prec)
msk[which(!is.na(msk[]))] <- 1

#read in ERA data =====
all_data <- read.csv(paste(cimdir,"/ERA_Test_Data.csv",sep=""))
all_data <- droplevels(all_data)

#aggregate data up to site level =====
rnd <- 4
agg_by <- c("Site.ID","PrName","Out.SubInd","Product.Simple","Product.Type")

#determine outliers =====
outlier_calc <- function(x){
  return((x < quantile(x)[2] - 3 *  IQR(x)  | x > quantile(x)[4] + 3 * IQR(x)))
}

all_data <- as.data.table(all_data)
outliers <- unlist(all_data[,R:=1:nrow(all_data)][,list(outliers=list(R[outlier_calc(yi)])), by=agg_by][,outliers])

#Calculate Weighted mean effect size by practice and site combination  =====
wt_grp <- unique(c("Code",agg_by))
data_sites <- all_data[!outliers
                       ][,N.Obs.Study:=.N,by=wt_grp # Recalculate Weightings by study within obervations grouping 
                         ][,Weight.Study:=(Rep^2/(2*Rep))/N.Obs.Study # Recalculate Weightings
                           ][,list(Observations=.N,
                                   Studies=length(unique(Code)),
                                   Country=unique(Country),
                                   AEZ16=unique(AEZ16simple),
                                   Lat=mean(Latitude),
                                   Lon=mean(Longitude),
                                   RR=round(weighted.mean(yi,Weight.Study,na.rm=T),rnd),
                                   RR.sd=suppressWarnings(round(abs(wtd.var(yi,Weight.Study,na.rm=T))^0.5,rnd))),
                             by=agg_by  # this 
                             ][,pc:=round(100*exp(RR)-100,rnd-1)
                               ][,pc.sd:=round(100*exp(RR.sd)-100,rnd-1)]
data_sites$Label <- paste(data_sites$Site.ID, data_sites$PrName)
data_sites<-data_sites[!(is.na(Lat) | is.na(Lon))]

#back to data.frame
data_sites[,PrName:=as.character(PrName)]
data_sites <- as.data.frame(data_sites)
data_sites$Npracs<-unlist(lapply(strsplit(data_sites$PrName,"-"),length))


#Create Scenarios x Years x Tresholds Loop ####
Scenarios <- c("rcp4.5","rcp8.5")
Years <- c(2030,2050)
Thresholds <- c(0.0,0.15,0.27,0.41)
Vars <- expand.grid(Years=Years, Scenarios=Scenarios, Threshold=Thresholds)
Vars$Scenarios <- as.character(Vars$Scenarios)
Vars <- rbind(Vars,expand.grid(Years=NA, Scenarios="baseline", Threshold=Thresholds))

#options
DoNeg<-F # Produce negative suitability?
DoLite<-T # To save time average the soil rasters across depths (rather than using multiple depths) and cut out min class and quantiles from run_points function

ExcludeProducts<-c("Fodder Legume","Fodder Tree","Okra","Olive","Other Bean","Pepper","Pumpkin","Tomato","Watermelon","African Yam Bean","Amaranth",
                   "Amaranth Grain","Apple","Cabbage","Capsicum","Carrot & Parsnip","Chickpea","Chili","Cucumber","Date","Eggplant","Fava Bean","Firewood",
                   "Fodder Tree","Fonio","Garlic","Grape","Grapefruit & Pomelo","Jatropha","Jujube","Lablab","Melon","Napier Grass","Onion","Other Leafy Green",
                   "Other Spice","Other Veg","Peas","Spinach","Turnip","Zucchini")
MinSites<-1
MaxPracs<-1 # Max number of practices to consider

for(k in 1:nrow(Vars)){
    Scenario<-Vars$Scenarios[k]
    Year<-Vars$Years[k]
    Variable<-Vars$Vars[k]
    Threshold<-Vars$Threshold[k]
    print(paste0("Running: Scenario = ",Vars$Scenarios[k]," | Year = ",Vars$Years[k]," | Threshold = ",Vars$Threshold[k]))
    
    ############################################################
    #load and resample historical and future hazards data -----
    # dry days mean (hzrs_m): Load Baseline ####
    File<-paste0(cimdir,"/input_data/hzrs_m.rda")
    if(!file.exists(File)){
        hzrs_m <- stack(paste(hzdir,"/dry_days_hist/monthly/dry_days_1981_2019_",1:12,"_mean.tif",sep="")) %>%
                raster::readAll(.) %>%
                raster::crop(., msk) %>%
                resample (.,msk) %>%
                readAll()
        save(hzrs_m, file=File,compress="xz",compression_level=9)
    } else {
        load(File)
    }

    # hzrs_m: Load Future Data #####
    File.Future<-paste0(cimdir,"/input_data/hzrs_m_",Year,"_",Scenario,".rda")
    if (Scenario != 'baseline') {
       if(!file.exists(File.Future)){
           if(Year==2030){
               hzrs_m_fut <- stack(paste(hzdir,"/dry_days_future/monthly/dry_days_",Scenario,"_",Year,"/dry_days_2020_2049_",1:12,"_mean.tif",sep=""))
           } else {
               hzrs_m_fut <- stack(paste(hzdir,"/dry_days_future/monthly/dry_days_",Scenario,"_",Year,"/dry_days_2040_2069_",1:12,"_mean.tif",sep=""))
           }
           hzrs_m_fut<-hzrs_m_fut %>%
                       raster::readAll(.) %>%
                       raster::crop(., msk) %>%
                       raster::resample(msk, method = 'bilinear') %>%
                       raster::stack(.) %>%
                       raster::readAll()
           save(hzrs_m_fut, file=File.Future,compress="xz",compression_level=9)
       } else {
           load(File.Future)
       }
    }

    # dry days cv (hzrs_v): Load Baseline ####
    File<-paste0(cimdir,"/input_data/hzrs_v.rda")
    if(!file.exists(File)){
        hzrs_v <- stack(paste(hzdir,"/dry_days_hist/monthly/dry_days_1981_2019_",1:12,"_cv.tif",sep="")) %>%
                    raster::readAll(.) %>%
                    raster::crop(., msk) %>%
                    raster::resample(.,msk) %>%
                    raster::readAll()
        save(hzrs_v, file=File,compress="xz",compression_level=9)
    } else {
        load(File)
    }
    
    # dry days cv (hzrs_v): Load Future #####
    File.Future<-paste0(cimdir,"/input_data/hzrs_v_",Year,"_",Scenario,".rda")
    if (Scenario != 'baseline') {
       if(!file.exists(File.Future)){
           if(Year==2030){
               hzrs_v_fut <- stack(paste(hzdir,"dry_days_future/monthly/dry_days_",Scenario,"_",Year,"/dry_days_2020_2049_",1:12,"_cv.tif",sep=""))
           } else {
               hzrs_v_fut <- stack(paste(hzdir,"dry_days_future/monthly/dry_days_",Scenario,"_",Year,"/dry_days_2040_2069_",1:12,"_cv.tif",sep=""))
           }

           hzrs_v_fut<-hzrs_v_fut %>%
                       raster::readAll(.) %>%
                       raster::crop(., msk) %>%
                       raster::resample(msk, method = 'bilinear') %>%
                       raster::stack(.) %>%
                       raster::readAll()
           save(hzrs_v_fut, file=File.Future,compress="xz",compression_level=9)
       } else {
           load(File.Future)
       }
    }
    
    #ppt driest quarter data =====
    #ppt driest quarter: Load Baseline ####
    File<-paste0(cimdir,"/input_data/pdq_rs.rda")
    if(!file.exists(File)){
        Data.File<-paste(hzdir,"/ppt_driest_quarter/pdq_hist.tif",sep="")
        pdq_rs <- raster(Data.File) %>%
                raster::readAll(.) %>%
                raster::resample(.,msk) %>%
                raster::readAll()
        save(pdq_rs, file=File,compress="xz",compression_level=9)
    } else {
        load(File)
    }
  
    #ppt driest quarter: Load Future ####
    File<-paste0(cimdir,"/input_data/pdq_rs_",Year,"_",Scenario,".rda")
    if (Scenario != 'baseline') {
       if(!file.exists(File)){
           Data.File<-paste0(hzdir,"/ppt_driest_quarter/pdq_fut_mmm_",Year,"_",Scenario,".tif")
           pdq_rs_fut <- raster(Data.File) %>%
                       raster::readAll(.) %>%
                       raster::resample(.,msk) %>%
                       raster::readAll()
           save(pdq_rs_fut, file=File,compress="xz",compression_level=9)
       } else {
           load(File)
       }
    }
    
    #aridity index data =====
    #aridity index data: Load Baseline ####
    File<-paste0(cimdir,"/input_data/ai_rs.rda")
    if(!file.exists(File)){
        Data.File<-paste(hzdir,"aridity/aridity_thornthwaite_hist.tif",sep="")
        ai_rs <- raster(Data.File) %>%
                raster::readAll(.) %>%
                raster::resample (msk, method = 'bilinear') %>%
                raster::readAll()
        save(ai_rs, file=File,compress="xz",compression_level=9)
    } else {
        load(File)
    }
      
    #aridity index data: Load Future ####
    File<-paste0(cimdir,"/input_data/ai_rs.rda_",Year,"_",Scenario,".rda")
    if (Scenario != 'baseline') {
       if(!file.exists(File)){
           Data.File<-paste0(hzdir,"aridity/aridity_thornthwaite_fut_mmm_",Year,"_",Scenario,".tif")
           ai_rs_fut <- raster(Data.File) %>%
                       raster::readAll(.) %>%
                       raster::resample(msk, method = 'bilinear') %>%
                       raster::readAll()
           save(ai_rs_fut, file=File,compress="xz",compression_level=9)
       } else {
           load(File)
       }
    }
    
    #chirps cv data =====
    #chirps cv data: Load Baseline ####
    File<-paste0(cimdir,"/input_data/chcv_rs.rda")
    if(!file.exists(File)){
        chcv_rs <- raster(paste(hzdir,"/chirps_cv/cvr_africa.tif",sep="")) %>%
                    raster::readAll(.) %>%
                    raster::crop(., msk) %>%
                    raster::resample(msk, method = 'bilinear') %>%
                    raster::readAll()
        save(chcv_rs, file=File,compress="xz",compression_level=9)
    } else {
        load(File)
    }
      
    #chirps cv data: Load Future #####
    File<-paste0(cimdir,"/input_data/chcv_",Year,"_",Scenario,".rda")
    if (Scenario != 'baseline') {
       if(!file.exists(File)){
           chcv_rs_fut <- raster(paste(hzdir,"/chirps_cv/cv_",Scenario,"_",Year,".tif",sep="")) %>%
                       raster::readAll(.) %>%
                       raster::resample (msk, method = 'bilinear') %>%
                       raster::readAll()
           save(chcv_rs_fut, file=File,compress="xz",compression_level=9)
       } else {
           load(File)
       }
    }
    
    #load monthly climate data -----
    #climate: load baseline ####
    load(paste(cimdir,"/input_data/wc_prec.rda",sep="")); wc_prec <- stack(wc_prec)
    load(paste(cimdir,"/input_data/wc_tmean.rda",sep="")); wc_tmean <- stack(wc_tmean)
    
    #climate: load future ####
    File<-paste0(cimdir,"/input_data/prec_",gsub("[.]","_",Scenario),"_",Year,"_","Resamp.rda")
    if (Scenario != 'baseline') {
       if(!file.exists(File)){
           wc_prec_fut<-load.Rdata2(filename=paste0("prec_",gsub("[.]","_",Scenario),"_",Year,"_","2.5min.RData"),path=cimdir) 
           wc_prec_fut<-resample (wc_prec_fut,msk, method = 'bilinear')
           wc_prec_fut<-round(wc_prec_fut,0)
           wc_prec_fut<-stack(wc_prec_fut)
           if(!inMemory(wc_prec_fut)){wc_prec_fut<-readAll(wc_prec_fut)}
           save(wc_prec_fut,file=File)
       } else {
           wc_prec_fut<-load.Rdata2(filename=paste0("prec_",gsub("[.]","_",Scenario),"_",Year,"_","Resamp.rda"),path=paste0(cimdir,"/input_data"))
       }
    }
    
    File<-paste0(cimdir,"/input_data/tmean_",gsub("[.]","_",Scenario),"_",Year,"_","Resamp.rda")
    if (Scenario != 'baseline') {
       if(!file.exists(File)){
           wc_tmean_fut<-load.Rdata2(filename=paste0("tmean_",gsub("[.]","_",Scenario),"_",Year,"_","2.5min.RData"),path=cimdir) 
           wc_tmean_fut<-resample (wc_tmean_fut,msk, method = 'bilinear')
           wc_tmean_fut<-round(wc_tmean_fut)
           wc_tmean_fut<-stack(wc_tmean_fut)
           if(!inMemory(wc_tmean_fut)){wc_tmean_fut<-readAll(wc_tmean_fut)}
           save(wc_tmean_fut,file=File)
       } else {
           wc_tmean_fut<-load.Rdata2(filename=paste0("tmean_",gsub("[.]","_",Scenario),"_",Year,"_","Resamp.rda"),path=paste0(cimdir,"/input_data")) 
       }
    }
   
    ############################################################
    ############################################################
    #select practice and organize practice data -----
    
    #calculate similarity (both temp and precip for positive sites) - Function -----
    #adjust so that for each site
    #1. similarity for mean climate
    #2. similarity for each soil variable
    #3. similarity for hazard layer (mean and var)
    #4. combine all in a sensible layer and save individual layers and output as .RData
    
    #Subset Era Data ====
    Y<-data.table(data_sites)
    
    #mulch-reduced tillage interaction
    M.RT<-Y[Product.Type=="Plant Product" &  Out.SubInd == "Crop Yield" & RR>=Threshold &!PrName=="" & !Product.Simple=="" & PrName == "Mulch-Reduced Tillage" & !Product.Simple %in% ExcludeProducts,
            list(N.Sites=length(unique(Site.ID)),N.Countries=length(unique(Country)),N.AEZ16=length(unique(AEZ16))),
            by=c("PrName","Product.Simple","Out.SubInd")][N.Sites >= MinSites]
    
    #individual practices
    Y<-Y[Product.Type=="Plant Product" &  Out.SubInd == "Crop Yield" & RR>=Threshold &!PrName=="" & !Product.Simple=="" & Npracs<=MaxPracs & !Product.Simple %in% ExcludeProducts,
         list(N.Sites=length(unique(Site.ID)),N.Countries=length(unique(Country)),N.AEZ16=length(unique(AEZ16))),
         by=c("PrName","Product.Simple","Out.SubInd")
         ][N.Sites >= MinSites]
    
    #combine both tables
    Y<-rbind(Y,M.RT)
    Y[,Product.Simple:=as.character(Product.Simple)]
    
    # Points <=10 - Run practices in parallel (gets clogged up when some cores have a practice with large numbers of sites) =====
    run_points <- function(i, pr_df, data_sites, cimdir, Threshold, Year, Scenario, vr, etype='pos', DoLite=F) {
        #get practice and product name
        prname<-pr_df[i,PrName]
        Product<-pr_df[i,Product.Simple]
        
        #Restructure data and classifiy points-----
        pdata <- data_sites[which(data_sites$PrName == prname & data_sites$Product.Simple == Product & data_sites$Out.SubInd == "Crop Yield"),]
        pdata <- pdata[which(!is.na(pdata$Lat)),]
        pdata <- pdata[which(!is.na(pdata$Lon)),]
        pdata <- droplevels(pdata)
        pdata <- pdata[,c("Lon","Lat","RR")]
        
        #cleanup any data points outside African continent
        pdata$isvalid <- raster::extract(msk, pdata[,c("Lon","Lat")])
        pdata <- pdata[which(!is.na(pdata$isvalid)),]
        pdata$isvalid <- NULL
        
        #positive yield impact  =====
        in_data <- pdata[which(pdata$RR >= Threshold),]
        in_data$ID <- paste("s",1:nrow(in_data),sep="")
        
        #output directory  =====
        pr_odir <- paste(cimdir,"/T",Threshold,"/",Product,"/",Year,"/",gsub("[.]","_",Scenario),"/",gsub(" ", "_", prname, fixed=T),"_v",vr,sep="")
        if (!file.exists(pr_odir)) {dir.create(pr_odir,recursive=T)}
        
        #verbose what i'm running
        cat("processing practice=", prname, i, "of", nrow(pr_df), "/ crop=", Product, "/ n=", nrow(in_data),"\n")
        
        #loop through points for practice
        allstk <- c()
        for (j in 1:nrow(in_data)) {
          #cat("processing site=",j,"\n")
          if (!file.exists(paste(pr_odir,"/out_",etype,"_",in_data$ID[j],".RData",sep=""))) {
            #1. mean climate ####
            #cat("...mean climate\n")
            par1 <- createParameters(x=in_data$Lon[j], y=in_data$Lat[j], vars=c("prec","tmean"),weights=c(0.75,0.25),
                                     ndivisions=c(12,12),growing.season=c(1,12),
                                     rotation="prec",threshold=1,env.data.ref=list(wc_prec,wc_tmean), 
                                     env.data.targ=if(is.na(Year)){list(wc_prec,wc_tmean)}else{list(wc_prec_fut,wc_tmean_fut)},outfile=pr_odir,
                                     fname=NA,writefile=F)
            sim1 <- calc_similarity(par1)
            if (!inMemory(sim1)) {sim1 <- readAll(sim1)}
            
            #2. each soil variable, then combine in one by taking the median ####
            #cat("...soils ORCDRC, PHIHOX, CLYPPT, SNDPPT\n")
            simsol <- c()
            for (svar in c("ORCDRC","PHIHOX","CLYPPT","SNDPPT")) {
                #print(svar)
                if(!DoLite){soil_data <- stack(soilstk[[svar]])} else {soil_data <- soilstk2[[svar]]}
                parx <- createParameters(x=in_data$Lon[j], y=in_data$Lat[j], vars=c(svar),weights=1,
                                         ndivisions=c(3),growing.season=c(1,3),rotation="none",
                                         threshold=1,env.data.ref=list(soil_data), 
                                         env.data.targ=list(soil_data),outfile=pr_odir,
                                         fname=NA,writefile=F)
                simsol <- c(simsol,calc_similarity(parx))
            }
            sim2 <- stack(simsol) %>% min(., na.rm=TRUE)
            if (!inMemory(sim2)) {sim2 <- readAll(sim2)} 
            rm(soil_data)
            gc()
            
            #3. similarity hazard layer (mean and variability) ####
            #cat("...hazard layers\n")
            #a. mean
            #cat("   ...mean\n")
            par3a <- createParameters(x=in_data$Lon[j], y=in_data$Lat[j], vars=c("tmean"),weights=1,
                                      ndivisions=c(12),growing.season=c(1,12),
                                      rotation="tmean",threshold=1,env.data.ref=list(stack(hzrs_m)), 
                                      env.data.targ=if(is.na(Year)){list(hzrs_m)}else{list(hzrs_m_fut)},outfile=pr_odir,
                                      fname=NA,writefile=F)
            sim3a <- calc_similarity(par3a)
            if (!inMemory(sim3a)) {sim3a <- readAll(sim3a)}
            
            #b. variability
            #cat("   ...variability\n")
            par3b <- createParameters(x=in_data$Lon[j], y=in_data$Lat[j], vars=c("tmean"),weights=1,
                                      ndivisions=c(12),growing.season=c(1,12),
                                      rotation="tmean",threshold=1,env.data.ref=list(stack(hzrs_v)), 
                                      env.data.targ=if(is.na(Year)){list(stack(hzrs_v))}else{list(stack(hzrs_v_fut))},outfile=pr_odir,
                                      fname=NA,writefile=F)
            sim3b <- calc_similarity(par3b)
            if (!inMemory(sim3b)) {sim3b <- readAll(sim3b)}
            
            #c. ppt driest quarter
            #cat("   ...prec driest quarter\n")
            par3c <- createParameters(x=in_data$Lon[j], y=in_data$Lat[j], vars=c("bio_1"),
                                      weights=c(1),ndivisions=c(1),growing.season=NA,
                                      rotation="none",threshold=1,env.data.ref=list(pdq_rs), 
                                      env.data.targ=if(is.na(Year)){list(pdq_rs)}else{list(pdq_rs_fut)},outfile=pr_odir,
                                      fname=NA,writefile=F)
            sim3c <- calc_similarity(par3c)
            if (!inMemory(sim3c)) {sim3c <- readAll(sim3c)}
            
            #d. aridity index
            #cat("   ...aridity index\n")
            par3d <- createParameters(x=in_data$Lon[j], y=in_data$Lat[j], vars=c("bio_1"),
                                      weights=c(1),ndivisions=c(1),growing.season=NA,
                                      rotation="none",threshold=1,env.data.ref=list(ai_rs), 
                                      env.data.targ=if(is.na(Year)){list(ai_rs)}else{list(ai_rs_fut)},outfile=pr_odir,
                                      fname=NA,writefile=F)
            sim3d <- calc_similarity(par3d)
            if (!inMemory(sim3d)) {sim3d <- readAll(sim3d)}
            
            #e. chirps cv
            #cat("   ...chirps c.v.\n")
            par3e <- createParameters(x=in_data$Lon[j], y=in_data$Lat[j], vars=c("bio_1"),
                                      weights=c(1),ndivisions=c(1),growing.season=NA,
                                      rotation="none",threshold=1,env.data.ref=list(chcv_rs), 
                                      env.data.targ=if(is.na(Year)){list(chcv_rs)}else{list(chcv_rs_fut)},outfile=pr_odir,
                                      fname=NA,writefile=F)
            sim3e <- calc_similarity(par3e)
            if (!inMemory(sim3e)) {sim3d <- readAll(sim3e)}
            
            #calculate low limit for hazards
            #cat("   ...aggregating\n")
            sim3 <- stack(c(sim3a, sim3b, sim3c, sim3d, sim3e)) %>% min(., na.rm=TRUE)
            if (!inMemory(sim3)) {sim3 <- readAll(sim3)}
            
            #4. combine all in a sensible layer and save individual layers and output as .RData ####
            #soil weight 50%, and climate indicators are divided equally over the remaining %
            #sim_i <- sim1 * 0.25 + sim2 * 0.5 + sim3 * 0.25
            #cat("...minimum similarity across all factors\n")
            sim_i <- stack(c(sim1, sim2, sim3)) %>%
              min(., na.rm=TRUE)
            if (!inMemory(sim_i)) {sim_i <- readAll(sim_i)}
            
            #5. put results in stack ####
            simstk <- stack(sim1, sim2, sim3, sim_i)
            names(simstk) <- c("meanclim","soil","hazard","overall")
            if (!inMemory(simstk)) {simstk <- readAll(simstk)}
            
            #6. save results  ####
            #cat("...saving\n")
            save(list=c("simstk"), 
                 file=paste(pr_odir,"/out_",etype,"_",in_data$ID[j],".RData",sep=""),
                 compress="xz",compression_level=9)
            
            #cleanup ####
            #cat("...clean-up\n")
            rm(list=c("sim1","sim2","sim3a","sim3b","sim3c","sim3d","sim3e","sim3"))
            gc(full=TRUE, verbose=FALSE)
          } else {
            #cat("...already processed, hence loading\n")
            load(paste(pr_odir,"/out_",etype,"_",in_data$ID[i],".RData",sep=""))
          }
          allstk <- c(allstk, simstk[["overall"]])
          rm(simstk)
        }
        
        #maximum values  - see original code for more information
        if(nlayers(stack(allstk))>1){
          maxsim <- max(stack(allstk), na.rm=TRUE) %>%
                    raster::writeRaster(., paste(pr_odir,"/max_similarity_",etype,".tif",sep=""),
                                        overwrite=T)
        }else{
          maxsim <- allstk[[1]] %>%
                    raster::writeRaster(., paste(pr_odir,"/max_similarity_",etype,".tif",sep=""),
                                        overwrite=T)
        }
        
        #verbose what i just did
        cat("completed practice=",prname," ",i,"of",nrow(pr_df),"/n=",nrow(in_data),"\n")
        
        #return object
        return(maxsim)
    }
    
    #practice list with only two fields
    X<-Y[,c("PrName","Product.Simple")]
    
    #run in parallel if 1 practice or more
    if(nrow(X)>0){
        xres <- parallel::mclapply(1:nrow(X), run_points, pr_df = X, data_sites, cimdir, 
                                   Threshold, Year, Scenario, vr, etype = "pos", DoLite,
                                   mc.cores = Cores, mc.preschedule = FALSE)
    }

    #clean-up
    rm(X,Y)
    gc()
}


