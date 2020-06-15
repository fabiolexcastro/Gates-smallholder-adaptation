#JRV - calculate mean and variability indicators for max. cons. dry days
#June 2020

#load libraries
library(raster)
library(analogues)

#working directory
wd <- "~/work"
hdir <- paste(wd,"/hazard_layers",sep="")

#layer directory
adir <- paste(hdir,"/ppt_driest_quarter",sep="")
if (!file.exists(adir)) {dir.create(adir)}

#given classes
m <- c(0, 40, 1,  
       40, 60, 2,
       60, 80, 3,
       80, 100, 4)
rclmat <- matrix(m, ncol=3, byrow=TRUE)

#get tmax, tmin, and prec from WCLv1.4
if (!file.exists(paste(adir,"/pdq_hist.tif",sep=""))) {
    pdq_rs <- getData('worldclim', var='bio', res=2.5, path=adir)
    pdq_rs <- pdq_rs[["bio17"]]
    pdq_rs <- writeRaster(pdq_rs, paste(adir,"/pdq_hist.tif",sep=""))
} else {
    pdq_rs <- raster(paste(adir,"/pdq_hist.tif",sep=""))
}

#clean-up folder
if (file.exists(paste(adir,"/wc2-5",sep=""))) {unlink(paste(adir,"/wc2-5",sep=""), recursive=TRUE, force=TRUE)}

#clean up all objects except srad
rm(list=c("pdq_rs"))
gc(T)

##### do the future scenarios
#MOHC_HADGEM2_ES (29), CESM1_CAM5 (6), GFDL_CM3 (13), MPI ESM_LR (30), MIROC_MIROC5 (27)
#for RCP4.5, RCP8.5, 2030 and 2050
for (sce in c(4.5, 8.5)) {
    for (yr in c(2030, 2050)) {
        for (k in c(6,13,27,29,30)) {
            #sce <- 4.5; yr <- 2030; k <- 6
            cat("processing rcp",sce,"/ year=",yr,"/ model=",k,"\n")
            
            #get CMIP5 bioclim data from CCAFS-Climate
            if (!file.exists(paste(adir,"/pdq_fut_m",k,"_",yr,"_rcp",sce,".tif",sep=""))) {
                pdq_rs <- getCMIP5(var="bio", rcp=sce, model=k, year=yr, res=2.5, path=adir)
                pdq_rs <- pdq_rs[["bio_17"]]
                pdq_rs <- writeRaster(pdq_rs, paste(adir,"/pdq_fut_m",k,"_",yr,"_rcp",sce,".tif",sep=""))
            } else {
                pdq_rs <- raster(paste(adir,"/pdq_fut_m",k,"_",yr,"_rcp",sce,".tif",sep=""))
            }
            
            #delete ascii files
            if (file.exists(paste(adir,"/cmip5",sep=""))) {unlink(paste(adir,"/cmip5",sep=""), recursive=TRUE, force=TRUE)}
            
            #clean up all objects except srad
            rm(list=c("pdq_rs"))
            gc(T)
        }
        
        #calculate muti-model mean of aridity index
        cat("calculate multi-model mean\n")
        if (!file.exists(paste(adir,"/pdq_fut_mmm_",yr,"_rcp",sce,".tif",sep=""))) {
            aindex <- stack(paste(adir,"/pdq_fut_m",c(6,13,27,29,30),"_",yr,"_rcp",sce,".tif",sep=""))
            aindex <- mean(aindex, na.rm=T)
            aindex <- writeRaster(aindex, paste(adir,"/pdq_fut_mmm_",yr,"_rcp",sce,".tif",sep=""))
            rm(aindex); gc(T)
        }
    }
}

setwd(hdir)
system("tar -cjvf ppt_driest_quarter.tar.bz2 ppt_driest_quarter")

