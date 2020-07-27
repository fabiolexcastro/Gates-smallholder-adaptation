#JRV - calculate mean and variability indicators for number of dry days
#July 2020

#working directory
wd <- "~/work"
hdir <- paste(wd,"/hazard_layers",sep="")

#libraries
library(raster)
library(tidyverse)

#layer name
lname <- "dry_days"
period <- "50s"
rcp <- "rcp45"

#years
if (period == "30s") {
    yi <- 2020
    yf <- 2049
} else {
    yi <- 2040
    yf <- 2069
}

#Africa mask
msk <- raster(paste(hdir,"/chirps_cv/cvr_africa.tif",sep=""))

#output base name
obdir <- paste(hdir,"/",lname,"_future/yearly/",lname,"_",rcp,"_",period,sep="")
obname <- paste(lname,"_",rcp,"_",period,"_",yi,"_",yf,"_", sep="")

#given classes
m <- c(0, 15, 1,  
       15, 20, 2,
       20, 25, 3,
       25, 35, 4)
rclmat <- matrix(m, ncol=3, byrow=TRUE)

#rename files (dryDays to dry_days)
flist <- list.files(obdir, pattern="\\.tif")
setwd(obdir)
for (fl in flist) {
    if (file.exists(fl)) {
        nfl <- gsub("dryDays", lname, fl)
        system(paste("mv ", fl, " ", nfl, sep=""))
    }
}
setwd("~")

#load data
rstk <- stack(paste(obdir, "/", lname, "_", rcp, "_", period, "_", yi:yf, ".tif", sep="")) %>%
        raster::resample(., msk, method="ngb")

#calculate monthly average
rstk <- rstk / 12

#calculate statistics
#mean
rsmean <- mean(rstk, na.rm=TRUE)
rsmean <- raster::writeRaster(rsmean, 
                             filename=file.path(path.expand(obdir), paste(obname,"mean.tif",sep="")), 
                             overwrite=TRUE)

#coefficient of variation
if (!file.exists(file.path(path.expand(obdir), paste(obname,"cv.tif",sep="")))) {
    rscv <- raster::calc(rstk, fun=function(x) {sd(x,na.rm=TRUE) / mean(x,na.rm=TRUE) * 100},
                         filename=file.path(path.expand(obdir), paste(obname,"cv.tif",sep="")), 
                         overwrite=TRUE)
}

#median (2.5 in 5 years)
if (!file.exists(file.path(path.expand(obdir), paste(obname,"p50.tif",sep="")))) {
    rsmedian <- raster::calc(rstk, fun=function(x) {quantile(x, probs=0.5, na.rm=TRUE)},
                             filename=file.path(path.expand(obdir), paste(obname,"p50.tif",sep="")),
                             overwrite=TRUE)
}

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
        rs_u <- raster::calc(rstk, fun=function(x) {quantile(x, probs=pval, na.rm=T)},
                           filename=fn_u,
                           overwrite=TRUE)
    }
    
    fn_l <- file.path(path.expand(obdir), 
                      paste(obname,"p",formatC(pval*100, width = 2, format = "d", flag = "0"),"l.tif",sep=""))
    if (!file.exists(fn_l)) {
        rs_l <- raster::calc(rstk, fun=function(x) {quantile(x, probs=(1-pval), na.rm=T)},
                               filename=fn_l,
                               overwrite=TRUE)
    }
    rm(list=c("rs_u","rs_l","fn_u","fn_l"))
}

#clean up
rm(rstk); gc()



################################################################
################################################################
#now calculate for each month
#output base name
obdir <- paste(hdir,"/",lname,"_future/monthly/",lname,"_",rcp,"_",period,sep="")
obname <- paste(lname,"_",yi,"_",yf,"_", sep="")

#rename files (dryDays to dry_days)
flist <- list.files(obdir, pattern="dryDays_")
setwd(obdir)
for (fl in flist) {
    if (file.exists(fl)) {
        nfl <- gsub("dryDays", lname, fl)
        system(paste("mv ", fl, " ", nfl, sep=""))
    }
}
setwd("~")

#function to parallelise
run_month <- function(i, obdir, obname, lname, yi, yf, msk) {
    #i <- 1
    cat("processing month=",i,"\n")
    
    #load rasters for all years
    rstk <- stack(paste(obdir,"/",lname,"_",yi:yf,"_",i,".tif",sep="")) %>%
            raster::resample(., msk, method="ngb")
    
    #calculate statistics
    #mean
    if (!file.exists(file.path(path.expand(obdir), paste(obname,i,"_mean.tif",sep="")))) {
        rsmean <- mean(rstk, na.rm=TRUE)
        rsmean <- raster::writeRaster(rsmean, 
                                     filename=file.path(path.expand(obdir), paste(obname,i,"_mean.tif",sep="")), 
                                     overwrite=TRUE)
        rm(rsmean)
    }
    
    #coefficient of variation
    if (!file.exists(file.path(path.expand(obdir), paste(obname,i,"_cv.tif",sep="")))) {
        rsmean <- raster(paste(file.path(path.expand(obdir), paste(obname,i,"_mean.tif",sep=""))))
        rsstd <- raster::calc(rstk, fun=sd, na.rm=TRUE)
        rscv <- rsstd / rsmean * 100
        rscv <- raster::writeRaster(rscv, 
                                   filename=file.path(path.expand(obdir), paste(obname,i,"_cv.tif",sep="")), 
                                   overwrite=TRUE)
        rm(list=c("rsstd","rsmean","rscv"))
    }
    
    #median (2.5 in 5 years)
    if (!file.exists(file.path(path.expand(obdir), paste(obname,i,"_p50.tif",sep="")))) {
        rsmedian <- raster::calc(rstk, fun=median, na.rm=TRUE, 
                               filename=file.path(path.expand(obdir), paste(obname,i,"_p50.tif",sep="")),
                               overwrite=TRUE)
        rm(rsmedian)
    }
    
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
                          paste(obname,i,"_p",formatC(pval*100, width = 2, format = "d", flag = "0"),"u.tif",sep=""))
        if (!file.exists(fn_u)) {
            rs_u <- raster::calc(rstk, fun=function(x) {quantile(x, probs=pval, na.rm=T)},
                               filename=fn_u,
                               overwrite=TRUE)
        }
        
        fn_l <- file.path(path.expand(obdir), 
                          paste(obname,i,"_p",formatC(pval*100, width = 2, format = "d", flag = "0"),"l.tif",sep=""))
        if (!file.exists(fn_l)) {
            rs_l <- raster::calc(rstk, fun=function(x) {quantile(x, probs=(1-pval), na.rm=T)},
                               filename=fn_l,
                               overwrite=TRUE)
        }
        rm(list=c("rs_u","rs_l","fn_u","fn_l"))
    }
    
    #clean up
    rm(rstk); gc()
}

#parallelization
require(doSNOW)
require(parallel)
require(foreach)

cl <- makeCluster(6)
registerDoSNOW(cl)

foreach(i = 1:12, .packages = c('raster', 'rgdal', 'sp', 'tidyverse'), .verbose = TRUE) %dopar% {
  print(i)
  run_month(i, obdir, obname, lname, yi, yf, msk)
}

stopCluster(cl)


#tar.bz2 everything
setwd(hdir)
system(paste("tar -cjvf ", lname, "_future.tar.bz2 ", lname, "_future", sep=""))

