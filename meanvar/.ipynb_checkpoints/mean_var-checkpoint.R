#JRV - calculate mean and variability indicators for given layer
#April 2020

#working directory
wd <- "~/work"
hdir <- paste(wd,"/hazard_layers",sep="")

#libraries
library(raster)

#layer name
lname <- "dry_days"
period <- "hist"

#years
yi <- 1981
yf <- 2019

#load all yearly layers, compute long-term mean, c.v. and 95th percentile
rstk <- stack(paste(hdir,"/",lname,"_",period,"/yearly/",lname,"_",yi:yf,".asc",sep=""))
rstk <- readAll(rstk)

#function to calculate statistics
calc_stats <- function(x) {
    #calculate statistics
    rsmean <- calc(x, fun=mean, na.rm=T)
    rsstd <- calc(x, fun=sd, na.rm=T)
    rscv <- rsstd/rsmean * 100
    
    #median (2.5 in 5 years)
    rsmedian <- calc(x, fun=median, na.rm=T)
    
    #2 in 5 years (probability of 40%)
    q40fun <- function(x){quantile(x, probs = .40, na.rm=TRUE)}
    q60fun <- function(x){quantile(x, probs = .60, na.rm=TRUE)}
    rs40 <- calc(x, fun=q40fun, forceapply=T)
    rs60 <- calc(x, fun=q60fun, forceapply=T)
    
    #1 in 5 years (20%)
    q20fun <- function(x){quantile(x, probs = .20, na.rm=TRUE)}
    q80fun <- function(x){quantile(x, probs = .80, na.rm=TRUE)}
    rs20 <- calc(x, fun=q20fun, forceapply=T)
    rs80 <- calc(x, fun=q80fun, forceapply=T)
    
    #1 in 10 years (10%)
    q10fun <- function(x){quantile(x, probs = .10, na.rm=TRUE)}
    q90fun <- function(x){quantile(x, probs = .90, na.rm=TRUE)}
    rs10 <- calc(x, fun=q10fun, forceapply=T)
    rs90 <- calc(x, fun=q90fun, forceapply=T)
    
    #1 in 20 years (5%)
    q95fun <- function(x){quantile(x, probs = .95, na.rm=TRUE)}
    q05fun <- function(x){quantile(x, probs = .05, na.rm=TRUE)}
    rs05 <- calc(x, fun=q05fun, forceapply=T)
    rs95 <- calc(x, fun=q95fun, forceapply=T)
    
    #1 in 50 years (2%)
    q02fun <- function(x){quantile(x, probs = .02, na.rm=TRUE)}
    q98fun <- function(x){quantile(x, probs = .98, na.rm=TRUE)}
    rs02 <- calc(x, fun=q02fun, forceapply=T)
    rs98 <- calc(x, fun=q98fun, forceapply=T)
    
    #1 in 100 years (1%)
    q01fun <- function(x){quantile(x, probs = .01, na.rm=TRUE)}
    q99fun <- function(x){quantile(x, probs = .99, na.rm=TRUE)}
    rs01 <- calc(x, fun=q01fun, forceapply=T)
    rs99 <- calc(x, fun=q99fun, forceapply=T)
    
    #return object
    outstk <- stack(rsmean, rscv, rsmedian, rs40, rs60, rs20, rs80, rs10, rs90, rs05, rs95, rs02, rs98, rs01, rs99)
    names(outstk) <- c("mean","cv","median","40p","60p","20p","80p","10p","90p","05p","95p","02p","98p","01p","99p")
    return(outstk)
}

#calculate statistics
rsst <- calc_stats(rstk)

#write rasters
for (k in 1:nlayers(39)) {
    kname <- names(rsst)[k]
    writeRaster(rsst[[k]], paste(hdir,"/",lname,"_",period,"/yearly/",lname,"_",yi,"_",yf,"_",kname,".asc",sep=""), overwrite=T)
}

#clean up
rm(rstk); rm(rsst); gc()

#now calculate for each month
for (i in 1:12) {
    cat("processing month=",i,"\n")
    
    #load all yearly layers, compute long-term mean, c.v. and 95th percentile
    rstk <- stack(paste(hdir,"/",lname,"_",period,"/monthly/",lname,"_",yi:yf,"_",i,".asc",sep=""))
    rstk <- readAll(rstk)
    
    #calculate statistics
    rsst <- calc_stats(rstk)
    
    #write rasters
    for (k in 1:nlayers(39)) {
        kname <- names(rsst)[k]
        writeRaster(rsst[[k]], paste(hdir,"/",lname,"_",period,"/monthly/",lname,"_",yi,"_",yf,"_",i,"_",kname,".asc",sep=""))
    }
    
    #clean up
    rm(rstk); rm(rsst); gc(T)
}

#double check results for given month
rstk <- stack(paste(hdir,"/",lname,"_",period,"/monthly/",lname,"_",yi,"_",yf,"_",1,"_",c("mean","cv","95p","05p"),".asc",sep=""))
plot(rstk)

#tar.bz2 everything
setwd(hdir)
system(paste("tar -cjvf ", lname, "_", period, ".tar.bz2 ", lname, "_", period), sep="")
