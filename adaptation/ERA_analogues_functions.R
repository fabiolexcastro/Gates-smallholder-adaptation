#load data
load_ERA_data <- function(cimdir, DoLite, Year, Scenario) {
    #load soil data
    #load Soils Data -----
    if (DoLite) {
        load(paste(cimdir,"/input_data/soils2.rda",sep=""), envir=.GlobalEnv)
    } else {
        load(paste(cimdir,"/input_data/soils.rda",sep=""), envir=.GlobalEnv)
    }
    
    #load and resample historical and future hazards data -----
    # dry days mean (hzrs_m): Load Baseline ####
    File<-paste0(cimdir,"/input_data/hzrs_m.rda")
    load(File, envir=.GlobalEnv)
    
    # hzrs_m: Load Future Data #####
    if (!is.na(Year)) {
        File.Future<-paste0(cimdir,"/input_data/hzrs_m_",Year,"_",Scenario,".rda")
        load(File.Future, envir=.GlobalEnv)
    }

    # dry days cv (hzrs_v): Load Baseline ####
    File<-paste0(cimdir,"/input_data/hzrs_v.rda")
    load(File, envir=.GlobalEnv)
    
    # dry days cv (hzrs_v): Load Future #####
    if (!is.na(Year)) {
        File.Future<-paste0(cimdir,"/input_data/hzrs_v_",Year,"_",Scenario,".rda")
        load(File.Future, envir=.GlobalEnv)
    }
    
    #ppt driest quarter data =====
    #ppt driest quarter: Load Baseline ####
    File<-paste0(cimdir,"/input_data/pdq_rs.rda")
    load(File, envir=.GlobalEnv)
  
    #ppt driest quarter: Load Future ####
    if (!is.na(Year)) {
        File<-paste0(cimdir,"/input_data/pdq_rs_",Year,"_",Scenario,".rda")
        load(File, envir=.GlobalEnv)
    }
  
    #aridity index data =====
    #aridity index data: Load Baseline ####
    File<-paste0(cimdir,"/input_data/ai_rs.rda")
    load(File, envir=.GlobalEnv)
      
    #aridity index data: Load Future ####
    if (!is.na(Year)) {
        File<-paste0(cimdir,"/input_data/ai_rs.rda_",Year,"_",Scenario,".rda")
        load(File, envir=.GlobalEnv)
    }
    
    #chirps cv data =====
    #chirps cv data: Load Baseline ####
    File<-paste0(cimdir,"/input_data/chcv_rs.rda")
    load(File, envir=.GlobalEnv)
      
    #chirps cv data: Load Future #####
    if (!is.na(Year)) {
        File<-paste0(cimdir,"/input_data/chcv_",Year,"_",Scenario,".rda")
        load(File, envir=.GlobalEnv)
    }
    
    #load monthly climate data -----
    #climate: load baseline ####
    load(paste(cimdir,"/input_data/wc_prec.rda",sep=""), envir=.GlobalEnv); wc_prec <<- stack(wc_prec)
    load(paste(cimdir,"/input_data/wc_tmean.rda",sep=""), envir=.GlobalEnv); wc_tmean <<- stack(wc_tmean)
    
    #climate: load future ####
    #precipitation
    if (!is.na(Year)) {
        wc_prec_fut <<- load.Rdata2(filename=paste0("prec_",gsub("[.]","_",Scenario),"_",Year,"_","Resamp.rda"), path=paste(cimdir,"/input_data",sep=""))
    }
    
    #temperature
    if (!is.na(Year)) {
        wc_tmean_fut <<- load.Rdata2(filename=paste0("tmean_",gsub("[.]","_",Scenario),"_",Year,"_","Resamp.rda"), path=paste(cimdir,"/input_data",sep=""))
    }
}

# Points <=10 - Run practices in parallel (gets clogged up when some cores have a practice with large numbers of sites) =====
run_points <- function(pr_i, pr_df, data_sites, cimdir, Threshold, Year, Scenario, vr, etype='pos', DoLite=FALSE) {
    #get practice and product name
    prname<-pr_df[pr_i,PrName]
    Product<-pr_df[pr_i,Product.Simple]
    
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
    cat("processing practice=", prname, pr_i, "of", nrow(pr_df), "/ crop=", Product, "/ n=", nrow(in_data),"\n")
    
    #loop through points for practice
    allstk <- c()
    for (i in 1:nrow(in_data)) {
      #cat("processing site=",i,"\n")
      if (!file.exists(paste(pr_odir,"/out_",etype,"_",in_data$ID[i],".RData",sep=""))) {
        #1. mean climate ####
        #cat("...mean climate\n")
        par1 <- createParameters(x=in_data$Lon[i], y=in_data$Lat[i], vars=c("prec","tmean"),weights=c(0.75,0.25),
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
            if(!DoLite){
                soil_data <- stack(soilstk[[svar]])
                ndiv <- 3; gseas <- c(1,3)
            } else {
                soil_data <- soilstk2[[svar]]
                ndiv <- 1; gseas <- NA
            }
            parx <- createParameters(x=in_data$Lon[i], y=in_data$Lat[i], vars=c(svar),weights=1,
                                     ndivisions=ndiv,growing.season=gseas,rotation="none",
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
        par3a <- createParameters(x=in_data$Lon[i], y=in_data$Lat[i], vars=c("tmean"),weights=1,
                                  ndivisions=c(12),growing.season=c(1,12),
                                  rotation="tmean",threshold=1,env.data.ref=list(stack(hzrs_m)), 
                                  env.data.targ=if(is.na(Year)){list(hzrs_m)}else{list(hzrs_m_fut)},outfile=pr_odir,
                                  fname=NA,writefile=F)
        sim3a <- calc_similarity(par3a)
        if (!inMemory(sim3a)) {sim3a <- readAll(sim3a)}
        
        #b. variability
        #cat("   ...variability\n")
        par3b <- createParameters(x=in_data$Lon[i], y=in_data$Lat[i], vars=c("tmean"),weights=1,
                                  ndivisions=c(12),growing.season=c(1,12),
                                  rotation="tmean",threshold=1,env.data.ref=list(stack(hzrs_v)), 
                                  env.data.targ=if(is.na(Year)){list(stack(hzrs_v))}else{list(stack(hzrs_v_fut))},outfile=pr_odir,
                                  fname=NA,writefile=F)
        sim3b <- calc_similarity(par3b)
        if (!inMemory(sim3b)) {sim3b <- readAll(sim3b)}
        
        #c. ppt driest quarter
        #cat("   ...prec driest quarter\n")
        par3c <- createParameters(x=in_data$Lon[i], y=in_data$Lat[i], vars=c("bio_1"),
                                  weights=c(1),ndivisions=c(1),growing.season=NA,
                                  rotation="none",threshold=1,env.data.ref=list(pdq_rs), 
                                  env.data.targ=if(is.na(Year)){list(pdq_rs)}else{list(pdq_rs_fut)},outfile=pr_odir,
                                  fname=NA,writefile=F)
        sim3c <- calc_similarity(par3c)
        if (!inMemory(sim3c)) {sim3c <- readAll(sim3c)}
        
        #d. aridity index
        #cat("   ...aridity index\n")
        par3d <- createParameters(x=in_data$Lon[i], y=in_data$Lat[i], vars=c("bio_1"),
                                  weights=c(1),ndivisions=c(1),growing.season=NA,
                                  rotation="none",threshold=1,env.data.ref=list(ai_rs), 
                                  env.data.targ=if(is.na(Year)){list(ai_rs)}else{list(ai_rs_fut)},outfile=pr_odir,
                                  fname=NA,writefile=F)
        sim3d <- calc_similarity(par3d)
        if (!inMemory(sim3d)) {sim3d <- readAll(sim3d)}
        
        #e. chirps cv
        #cat("   ...chirps c.v.\n")
        par3e <- createParameters(x=in_data$Lon[i], y=in_data$Lat[i], vars=c("bio_1"),
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
             file=paste(pr_odir,"/out_",etype,"_",in_data$ID[i],".RData",sep=""),
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
    cat("completed practice=",prname," ",pr_i,"of",nrow(pr_df),"/n=",nrow(in_data),"\n")
    
    #return object
    return(maxsim)
}

