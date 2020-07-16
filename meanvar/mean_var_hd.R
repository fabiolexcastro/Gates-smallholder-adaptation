#JRV - calculate mean and variability for number of hot days
#July 2020

#working directory
wd <- "~/work"
hdir <- paste(wd,"/hazard_layers",sep="")
caldir <- paste(wd,"/crop_calendar",sep="")

#libraries
library(raster)
library(tidyverse)
library(tibble)

#layer name
lname <- "heat_stress_days"
period <- "hist"

#years
yi <- 1983
yf <- 2016

#given classes
m <- c(0, 1, 1,  
       1, 5, 2,
       5, 10, 3,
       10, 50, 4)
rclmat <- matrix(m, ncol=3, byrow=TRUE)

#list of crops
crop_list <- list.files(paste(hdir,"/",lname,"_",period,sep=""))

#list 
cal_df <- data.frame(
           crop=c("arabica_coffee","banana","barley_s1","barley_s2","bean","cassava","chickpea","cocoa",
                  "coconut","cotton","cowpea","generic","groundnut","lentil","maiz_s1","maiz_s2","oilpalm",
                  "pearl_millet","pigeonpea","plantain","potato","rapeseed","rice_s1","rice_s2","robusta_coffee",
                  "sesameseed","small_millet","sorghum_s1","sorghum_s2","sugar_cane","sugarbeet","sunflower",
                  "sweet_potato","teas","tobacco","wheat_s1","wheat_s2","yams"),
           calendar=c("perennial","perennial","Barley","Barley.Winter","Pulses","Cassava","Pulses","perennial",
                      "perennial","Cotton","Pulses","perennial","Groundnuts","Pulses","Maize","Maize.2","perennial",
                      "Millet","Pulses","perennial","Potatoes","Rapeseed","Rice","Rice.2","perennial","Pulses","Millet",
                      "Sorghum","Sorghum.2","perennial","Sugarbeets","Sunflower","Sweet.Potatoes","perennial","Maize",
                      "Wheat","Wheat.Winter","Yams"))

#loop crops
for (crop_i in crop_list) {
    #crop_i <- crop_list[1]
    cat("processing crop=",crop_i,"\n")
    
    #calendar name
    calname <- paste(cal_df$calendar[which(cal_df$crop == crop_i)])
    
    #output name
    obdir <- paste(hdir,"/",lname,"_",period,"/",crop_i,sep="")
    obname <- paste("heat_crop_",crop_i,"_",yi,"_",yf,"_Africa_", sep="")

    #load all yearly layers, compute long-term mean, c.v. and 95th percentile
    rstk <- stack(paste(obdir,"/heat_crop_",crop_i,"_",yi:yf,"_Africa.tif",sep=""))
    
    #calculate monthly average (12 months in a year if perennial, else from crop calendar)
    if (calname == "perennial") {
        rstk <- rstk / 12
    } else {
        #work out calendar
        cal_plant <- raster(paste(caldir,"/",calname,".crop.calendar.fill.nc",sep=""), varname="plant") %>%
                raster::crop(., rstk) %>%
                raster::resample(., rstk, method="ngb")
        cal_harv <- raster(paste(caldir,"/",calname,".crop.calendar.fill.nc",sep=""), varname="harvest") %>%
                crop(., rstk) %>%
                resample(., rstk, method="ngb")
        cal_dur <- cal_harv - cal_plant #season length
        cal_dur[which(cal_dur[] < 0)] <- cal_dur[which(cal_dur[] < 0)] + 365 #correct where harv before plant
        cal_dur <- cal_dur / 30 #season length in months
        
        #days per month
        rstk <- rstk / cal_dur
    }
    #x <- terra::classify(rsmean, rclmat, include.lowest=TRUE)

    #calculate statistics
    #mean
    rsmean <- mean(rstk, na.rm=TRUE)
    rsmean <- raster::writeRaster(rsmean, 
                                 filename=file.path(path.expand(obdir), paste(obname,"mean.tif",sep="")), 
                                 overwrite=TRUE)

    #coefficient of variation
    rscv <- raster::calc(rstk, fun=function(x) {sd(x,na.rm=TRUE)/mean(x,na.rm=TRUE)})
    rscv <- raster::writeRaster(rscv, 
                               filename=file.path(path.expand(obdir), paste(obname,"cv.tif",sep="")), 
                               overwrite=TRUE)

    #median (2.5 in 5 years)
    rsmedian <- raster::calc(rstk, fun=function(x) {quantile(x, probs=0.5, na.rm=TRUE)})
    rsmedian <- raster::writeRaster(rsmedian, 
                                    filename=file.path(path.expand(obdir), paste(obname,"p50.tif",sep="")),
                                    overwrite=TRUE)

    rm(list=c("rsmean","rscv","rsmedian")); gc(T)

    #write probability rasters (upper and lower)
    #2 in 5 years (probability of 40%)
    #1 in 5 years (20%)
    #1 in 10 years (10%)
    #1 in 20 years (5%)
    #1 in 50 years (2%)
    #1 in 100 years (1%)
    for (pval in c(0.4, 0.2, 0.1, 0.05, 0.02, 0.01)) {
        #pval <- 0.4
        cat("processing and writing pval=",format(pval, nsmall=2),"\n")
        fn_u <- file.path(path.expand(obdir), 
                          paste(obname,"p",formatC(pval*100, width = 2, format = "d", flag = "0"),"u.tif",sep=""))
        if (!file.exists(fn_u)) {
            rs_u <- terra::calc(rstk, fun=function(x) {quantile(x, probs=pval, na.rm=T)},
                               filename=fn_u,
                               overwrite=TRUE,
                               wopt=list(gdal=c("COMPRESS=LZW")))
        }

        fn_l <- file.path(path.expand(obdir), 
                          paste(obname,"p",formatC(pval*100, width = 2, format = "d", flag = "0"),"l.tif",sep=""))
        if (!file.exists(fn_l)) {
            rs_l <- terra::app(rstk, fun=function(x) {quantile(x, probs=(1-pval), na.rm=T)},
                                   filename=fn_l,
                                   overwrite=TRUE,
                                   wopt=list(gdal=c("COMPRESS=LZW")))
        }
        rm(list=c("rs_u","rs_l","fn_u","fn_l"))
    }

    #clean up
    rm(rstk); gc()
}

#tar.bz2 everything
setwd(hdir)
system(paste("tar -cjvf ", lname, "_", period, ".tar.bz2 ", lname, "_", period, sep=""))

