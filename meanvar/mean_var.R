#JRV - calculate mean and variability indicators for given layer

#working directory
wd <- "~/work"
hdir <- paste(wd,"/hazard_layers",sep="")

#libraries
library(raster)

#layer name
lname <- "dry_days"

#years
yi <- 1981
yf <- 2019

#load all yearly layers, compute long-term mean, c.v. and 95th percentile
rstk <- stack(paste(hdir,"/",lname,"/yearly/",lname,"_",yi:yf,".asc",sep=""))
rstk <- readAll(rstk)

#function to calculate statistics
calc_stats <- function(x) {
    #calculate statistics
    rsmean <- calc(x, fun=mean, na.rm=T)
    rsstd <- calc(x, fun=sd, na.rm=T)
    rscv <- rsmean/rsstd * 100
    q95fun <- function(x){quantile(x, probs = .95, na.rm=TRUE)}
    q05fun <- function(x){quantile(x, probs = .05, na.rm=TRUE)}
    rs95 <- calc(x, fun=q95fun, forceapply=T)
    rs05 <- calc(x, fun=q05fun, forceapply=T)
    return(stack(rsmean,rscv,rs95,rs05))
}

#calculate statistics
rsst <- calc_stats(rstk)

#write rasters
writeRaster(rsst[[1]], paste(hdir,"/",lname,"/yearly/",lname,"_",yi,"_",yf,"_mean.tif",sep=""))
writeRaster(rsst[[2]], paste(hdir,"/",lname,"/yearly/",lname,"_",yi,"_",yf,"_cv.tif",sep=""))
writeRaster(rsst[[3]], paste(hdir,"/",lname,"/yearly/",lname,"_",yi,"_",yf,"_95p.tif",sep=""))
writeRaster(rsst[[4]], paste(hdir,"/",lname,"/yearly/",lname,"_",yi,"_",yf,"_05p.tif",sep=""))

#clean up
rm(rstk); rm(rsst); gc()

#now calculate for each month
for (i in 1:12) {
    cat("processing month=",i,"\n")
    
    #load all yearly layers, compute long-term mean, c.v. and 95th percentile
    rstk <- stack(paste(hdir,"/",lname,"/monthly/",lname,"_",yi:yf,"_",i,".asc",sep=""))
    rstk <- readAll(rstk)
    
    #calculate statistics
    rsst <- calc_stats(rstk)
    
    #write rasters
    writeRaster(rsst[[1]], paste(hdir,"/",lname,"/monthly/",lname,"_",yi,"_",yf,"_",i,"_mean.tif",sep=""))
    writeRaster(rsst[[2]], paste(hdir,"/",lname,"/monthly/",lname,"_",yi,"_",yf,"_",i,"_cv.tif",sep=""))
    writeRaster(rsst[[3]], paste(hdir,"/",lname,"/monthly/",lname,"_",yi,"_",yf,"_",i,"_95p.tif",sep=""))
    writeRaster(rsst[[4]], paste(hdir,"/",lname,"/monthly/",lname,"_",yi,"_",yf,"_",i,"_05p.tif",sep=""))
    
    #clean up
    rm(rstk); rm(rsst); gc()
}

