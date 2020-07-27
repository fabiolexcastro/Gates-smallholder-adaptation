#JRV - calculate mean and variability indicators for max. cons. dry days future climate
#June 2020

#working directory
wd <- "~/work"
hdir <- paste(wd,"/hazard_layers",sep="")

#libraries
library(raster)
library(tidyverse)

#layer name
lname <- "max_cons_dry_days"
period <- "50s"
rcp <- "rcp85"

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
obdir <- paste(hdir,"/",lname,"_future/",lname,"_",rcp,"_",period,sep="")
obname <- paste("drySsnYear_",yi,"_",yf,"_", sep="")

#rename files (dryDays to dry_days)
flist <- list.files(obdir, pattern="consecutiveDryDays")
setwd(obdir)
for (fl in flist) {
    if (file.exists(fl)) {
        nfl <- gsub("consecutiveDryDays", "drySsnYear", fl)
        system(paste("mv ", fl, " ", nfl, sep=""))
    }
}
setwd("~")

#load all yearly layers, compute long-term mean, c.v. and 95th percentile
rstk <- stack(paste(hdir,"/",lname,"_future/",lname,"_",rcp,"_",period,"/drySsnYear_",yi:yf,".tif",sep="")) %>%
        raster::resample(., msk, method="ngb")

#given classes
m <- c(0, 7, 1,  
       7, 21, 2,
       21, 42, 3,
       42, 100, 4)
rclmat <- matrix(m, ncol=3, byrow=TRUE)

#calculate seasonal average (4 seasons in a year)
rstk <- rstk / 4

#calculate statistics
#mean
rsmean <- mean(rstk, na.rm=TRUE)
rsmean <- raster::writeRaster(rsmean, 
                             filename=file.path(path.expand(obdir), paste(obname,"mean.tif",sep="")), 
                             overwrite=TRUE)

#coefficient of variation
rsstd <- raster::calc(rstk, fun=sd, na.rm=TRUE)
rscv <- rsstd / rsmean * 100
rscv <- raster::writeRaster(rscv, 
                           filename=file.path(path.expand(obdir), paste(obname,"cv.tif",sep="")), 
                           overwrite=TRUE)

#median (2.5 in 5 years)
rsmedian <- raster::calc(rstk, fun=function(x) {quantile(x, probs=0.5, na.rm=T)},
                            filename=file.path(path.expand(obdir), paste(obname,"p50.tif",sep="")),
                            overwrite=TRUE)

rm(list=c("rsmean","rsstd","rscv","rsmedian")); gc(T)

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

#tar.bz2 everything
setwd(hdir)
system(paste("tar -cjvf ", lname, "_future.tar.bz2 ", lname, "_future", sep=""))
