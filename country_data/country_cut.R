#JRV - calculate mean and variability indicators for number of dry days
#June 2020

#working directory
wd <- "~/work"
hdir <- paste(wd,"/hazard_layers",sep="")
odir <- paste(hdir, "/country_data", sep="")
if (!file.exists(odir)) {dir.create(odir)}

#libraries
library(raster)
library(rgdal)
library(tidyverse)

#country
ctry_name <- "Ethiopia"

#output directory for country
out_cdir <- paste(odir, "/", ctry_name, sep="")
if (!file.exists(out_cdir)) {dir.create(out_cdir)}

#layer name
lnames <- c("dry_days_hist/yearly", "hi", "chirps_cv", "heat_stress_flips", "max_cons_dry_days_hist",
           "ppt_driest_month", "ppt_driest_quarter", "aridity", "thi")

#load shapefile of Africa
sh_ctry <- readOGR(paste(wd,"/Africa_shp/African_continet.shp",sep=""))
sh_ctry <- sh_ctry[which(sh_ctry$ADM0_NAME %in% ctry_name),]
sh_xt <- extent(sh_ctry)
sh_xt@xmin <- sh_xt@xmin-1; sh_xt@ymin <- sh_xt@ymin-1; sh_xt@xmax <- sh_xt@xmax+1; sh_xt@ymax <- sh_xt@ymax+1

#load and crop historical layer
msk1 <- raster(paste(hdir, "/aridity/aridity_thornthwaite_hist.tif", sep="")) %>%
            setValues(., 1, index=which(!is.na(.[]))) %>%
            crop(., sh_xt) %>%
            mask(., sh_ctry)
msk2 <- raster(paste(hdir, "/dry_days_hist/yearly/dry_days_1981.tif",sep="")) %>%
            setValues(., 1, index=which(!is.na(.[]))) %>%
            crop(., sh_xt) %>%
            mask(., sh_ctry)
msk3 <- raster(paste(hdir, "/hi/hi_av_2000.asc", sep="")) %>%
        setValues(., 1, index=which(!is.na(.[]))) %>%
        crop(., sh_xt) %>%
        mask(., sh_ctry)

#loop through layers
for (lname in lnames) {
    #lname <- lnames[1]
    cat("processing layer=", lname, "\n")
    
    #layer output name
    loname <- gsub("\\/", "_", lname)
    
    #pattern
    if (lname %in% c("hi","thi")) {
        patt <- ".asc"
    } else if (lname %in% c("aridity")) {
        patt <- "aridity_thornthwaite_"
    } else {
        patt <- ".tif"
    }
    
    #mask
    if (lname %in% c("hi","thi")) {
        msk <- msk3
    } else if (lname %in% c("dry_days_hist/yearly", "max_cons_dry_days_hist", "chirps_cv")) {
        msk <- msk2
    } else {
        msk <- msk1
    }
    
    #output layer directory
    out_lydir <- paste(out_cdir, "/", loname, sep="")
    if (!file.exists(out_lydir)) {dir.create(out_lydir)}
    
    #load rasters
    rstk <- stack(list.files(paste(hdir, "/", lname, "/",sep=""), pattern=patt, full.names=TRUE)) %>%
            crop(., msk) %>% mask(., msk)
    
    #write rasters
    for (i in 1:nlayers(rstk)) {
        if (nlayers(rstk) == 1) {rsname <- lname} else {rsname <- names(rstk)[i]}
        x <- writeRaster(rstk[[i]], 
                         paste(out_lydir, "/", rsname, ".tif", sep=""), 
                         overwrite=TRUE)
        rm(x)
    }
}

#tar.bz2 everything
setwd(odir)
system(paste("tar -cjvf ", ctry_name, ".tar.bz2 ", ctry_name, sep=""))
