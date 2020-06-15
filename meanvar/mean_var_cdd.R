#JRV - calculate mean and variability indicators for max. cons. dry days
#June 2020

#working directory
wd <- "~/work"
hdir <- paste(wd,"/hazard_layers",sep="")

#libraries
library(terra)

#layer name
lname <- "max_cons_dry_days"
period <- "hist"

#years
yi <- 1981
yf <- 2019

#output base name
obdir <- paste(hdir,"/",lname,"_",period,sep="")
obname <- paste("drySsnYear_",yi,"_",yf,"_", sep="")

#load all yearly layers, compute long-term mean, c.v. and 95th percentile
rstk <- rast(paste(hdir,"/",lname,"_",period,"/drySsnYear_",yi:yf,".tif",sep=""))

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
rsmean <- terra::writeRaster(rsmean, 
                             filename=file.path(path.expand(obdir), paste(obname,"mean.tif",sep="")), 
                             overwrite=TRUE,
                             wopt=list(gdal=c("COMPRESS=LZW")))

#coefficient of variation
rsstd <- terra::app(rstk, fun=sd, na.rm=TRUE)
rscv <- rsstd / rsmean * 100
rscv <- terra::writeRaster(rscv, 
                           filename=file.path(path.expand(obdir), paste(obname,"cv.tif",sep="")), 
                           overwrite=TRUE,
                           wopt=list(gdal=c("COMPRESS=LZW")))

#median (2.5 in 5 years)
rsmedian <- terra::quantile(rstk, probs=0.5, na.rm=TRUE, 
                            filename=file.path(path.expand(obdir), paste(obname,"p50.tif",sep="")),
                            overwrite=TRUE,
                            wopt=list(gdal=c("COMPRESS=LZW")))

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
        rs_u <- terra::app(rstk, fun=function(x) {quantile(x, probs=pval, na.rm=T)},
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

#tar.bz2 everything
setwd(hdir)
system(paste("tar -cjvf ", lname, "_", period, ".tar.bz2 ", lname, "_", period, sep=""))
