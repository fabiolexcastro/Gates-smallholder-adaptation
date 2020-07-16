#JRV - calculate mean and variability indicators for number of dry days
#April 2020

#working directory
wd <- "~/work"
hdir <- paste(wd,"/hazard_layers",sep="")

#libraries
library(raster)
library(tidyverse)

#layer name
lname <- "heat_stress_flips"

#output base name
lydir <- paste(hdir,"/",lname,sep="")

#list of crops
crop_list <- list.files(lydir, pattern="heat_stress_thr_")
crop_list <- unique(unlist(lapply(strsplit(crop_list, split="_", fixed=TRUE), 
                                  FUN=function(x) {
                                        if (length(x)==7) {
                                            return(paste(x[4],"_",x[5],sep=""))
                                        } else {
                                            return(x[4])
                                        }
                                      })))

#mask
msk <- stack(list.files(lydir, pattern="\\.tif", full.names=TRUE)[1])
msk[which(!is.na(msk[]))] <- 1

#crop area
carea <- stack(stack(list.files(lydir, pattern="heat_stress_thr_", full.names=TRUE)))
carea <- mean(carea, na.rm=T)
carea[which(carea[] == 0)] <- NA
carea[which(!is.na(carea[]))] <- 1

#scenario and period list
sce_list <- expand.grid(rcp=c("rcp45","rcp85"),period=c("2050s","2030s"))

#loop through crops to create summary outputs
rsum1 <- rsum2 <- rsum3 <- c()
for (i in 1:nrow(sce_list)) {
    #i <- 1
    sce <- paste(sce_list$rcp[i])
    per <- paste(sce_list$period[i])
    
    cat("processing scenario=", sce, "/ period=", per, "\n")
    
    for (cname in crop_list) {
        #cname <- crop_list[1]
        cat("crop=", cname, "\n")
        
        #load raster
        #0: no crop
        #1: temperature flip
        #2: safe levels
        #3: both above threshold
        rs <- raster(paste(lydir,"/heat_stress_thr_",cname,"_",sce,"_",per,".tif",sep=""))
        
        #summary 1: count of crops that suffer a temperature flip
        rs1 <- raster(rs)
        rs1[which(!is.na(rs[]))] <- 0
        rs1[which(rs[] == 1)] <- 1
        rsum1 <- c(rsum1, rs1)
        
        #summary 2: count of crops that experience both periods above threshold
        rs2 <- raster(rs)
        rs2[which(!is.na(rs[]))] <- 0
        rs2[which(rs[] == 3)] <- 1
        rsum2 <- c(rsum2, rs2)
        
        #summary 3: count of crops that experience either period above threshold
        rs3 <- raster(rs)
        rs3[which(!is.na(rs[]))] <- 0
        rs3[which(rs[] == 3 | rs[] == 1)] <- 1
        rsum3 <- c(rsum3, rs3)
    }
    
    #stack and compute sum
    rsum1 <- stack(rsum1) %>%
                sum(., na.rm=T) %>%
                mask(., carea) %>%
                writeRaster(., filename=paste(lydir, "/heat_stress_hotspots_flipping_",sce,"_",per,".tif",sep=""), overwrite=TRUE)
    rsum2 <- stack(rsum2) %>%
                sum(., na.rm=TRUE) %>%
                mask(., carea) %>%
                writeRaster(., filename=paste(lydir, "/heat_stress_hotspots_already_",sce,"_",per,".tif",sep=""), overwrite=TRUE)
    rsum3 <- stack(rsum3) %>%
                sum(., na.rm=TRUE) %>%
                mask(., carea) %>%
                writeRaster(., filename=paste(lydir, "/heat_stress_hotspots_combined_",sce,"_",per,".tif",sep=""), overwrite=TRUE)
}

#tar.bz2 everything
setwd(hdir)
system(paste("rm -f ", lname, ".tar.bz2", sep=""))
system(paste("tar -cjvf ", lname, ".tar.bz2 ", lname, sep=""))

