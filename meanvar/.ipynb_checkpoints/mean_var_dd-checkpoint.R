#JRV - calculate mean and variability indicators for number of dry days
#April 2020

#working directory
wd <- "~/work"
hdir <- paste(wd,"/hazard_layers",sep="")

#libraries
library(raster)
library(tidyverse)

#layer name
lname <- "dry_days"
period <- "hist"

#years
yi <- 1981
yf <- 2019

#Africa mask
msk <- raster(paste(hdir,"/chirps_cv/cvr_africa.tif",sep=""))

#output base name
obdir <- paste(hdir,"/",lname,"_",period,"/yearly",sep="")
obname <- paste(lname,"_",yi,"_",yf,"_", sep="")

#given classes
m <- c(0, 15, 1,  
       15, 20, 2,
       20, 25, 3,
       25, 35, 4)
rclmat <- matrix(m, ncol=3, byrow=TRUE)

#replace .asc by .tif files
if (file.exists(paste(obdir,"/",lname,"_",yi,".asc",sep=""))) {
    #load all yearly layers, compute long-term mean, c.v. and 95th percentile
    rstk <- stack(paste(obdir,"/",lname,"_",yi:yf,".asc",sep=""))

    #write all these rasters as GTiff
    for (i in 1:nlyr(rstk)) {
        raster::writeRaster(rstk[[i]], 
                           filename=file.path(path.expand(obdir), paste(names(rstk)[i],".tif",sep="")),
                           overwrite=TRUE)
    }

    #reload stack from GTiff files
    rstk <- stack(paste(obdir,"/",lname,"_",yi:yf,".tif",sep=""))

    #remove all .asc files
    setwd(obdir)
    system("rm -f *.asc")
    setwd("~")
} else {
    rstk <- stack(paste(obdir,"/",lname,"_",yi:yf,".tif",sep=""))
}

#calculate monthly average
rstk <- rstk / 12
#rstk <- crop(rstk, msk)

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



################################################################
################################################################
#now calculate for each month
#output base name
obdir <- paste(hdir,"/",lname,"_",period,"/monthly",sep="")
obname <- paste(lname,"_",yi,"_",yf,"_", sep="")

for (i in 1:12) {
    #i <- 1
    cat("processing month=",i,"\n")
    
    #replace .asc by .tif files
    if (file.exists(paste(obdir,"/",lname,"_",yi,"_",i,".asc",sep=""))) {
        #load all yearly layers, compute long-term mean, c.v. and 95th percentile
        rstk <- stack(paste(obdir,"/",lname,"_",yi:yf,"_",i,".asc",sep=""))
    
        #write all these rasters as GTiff
        for (j in 1:nlyr(rstk)) {
            raster::writeRaster(rstk[[j]], 
                               filename=file.path(path.expand(obdir), paste(names(rstk)[j],".tif",sep="")),
                               overwrite=TRUE)
        }
    
        #reload stack from GTiff files
        rstk <- stack(paste(obdir,"/",lname,"_",yi:yf,"_",i,".tif",sep=""))
    
        #remove all .asc files
        setwd(obdir)
        system(paste("rm -f *_",i,".asc",sep=""))
        setwd("~")
    } else {
        rstk <- stack(paste(obdir,"/",lname,"_",yi:yf,"_",i,".tif",sep=""))
    }
    
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

#tar.bz2 everything
setwd(hdir)
system(paste("tar -cjvf ", lname, "_", period, ".tar.bz2 ", lname, "_", period, sep=""))
