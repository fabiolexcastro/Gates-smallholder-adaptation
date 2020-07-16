#JRV - Fix extent of number of hot days
#July 2020

#working directory
wd <- "~/work"
hdir <- paste(wd,"/hazard_layers",sep="")

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
yrs <- yi:yf

#load mask
msk <- raster(paste(hdir,"/heat_stress_flips/heat_stress_hotspots_already_rcp45_2030s.tif",sep=""))
msk[which(!is.na(msk[]))] <- 1
msk_vls <- raster::rasterToPoints(msk) %>%
            as_tibble() %>% 
            setNames(c('x', 'y', 'value')) %>%
            dplyr::select(., x, y)

#list of crops
crop_list <- list.files(paste(hdir,"/",lname,"_",period,"_v1",sep=""))

#loop crops
for (crop_i in crop_list) {
    #crop_i <- crop_list[1]
    cat("processing crop=",crop_i,"\n")
    
    #output name
    obdir <- paste(hdir,"/",lname,"_",period,"/",crop_i,sep="")
    if (!file.exists(obdir)) {dir.create(obdir, recursive=TRUE)}
    
    #load all yearly layers, compute long-term mean, c.v. and 95th percentile
    rstk <- stack(paste(hdir,"/",lname,"_",period,"_v1/",crop_i,"/heat_crop_",crop_i,"_",yi:yf,"_Africa.tif",sep=""))

    #calculate monthly average (12 months in a year)
    for (i in 1:nlayers(rstk)) {
        if (!file.exists(paste(obdir,"/heat_crop_",crop_i,"_",yrs[i],"_Africa.tif",sep=""))) {
            msk_vls <- dplyr::mutate(msk_vls, value=raster::extract(rstk[[i]], msk_vls[,c('x','y')]))
            out_rs <- raster(msk)
            out_rs[cellFromXY(msk, as.data.frame(msk_vls[,c('x','y')]))]  <- msk_vls$value
            out_rs <- writeRaster(out_rs, 
                                  paste(obdir,"/heat_crop_",crop_i,"_",yrs[i],"_Africa.tif",sep=""),
                                  overwrite=TRUE)
            msk_vls <- dplyr::select(msk_vls, x, y)
        }
    }
}
