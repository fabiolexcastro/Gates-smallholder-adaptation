#JRV - calculate mean and variability indicators for max. cons. dry days
#June 2020

#load libraries
library(raster)
library(envirem)
library(analogues)

#assignNames
assignNames(solrad='srad.##', tmin='tmin.##', tmax='tmax.##', precip='prec.##', tmean='tmean.##')

#working directory
wd <- "~/work"
hdir <- paste(wd,"/hazard_layers",sep="")

#layer directory
adir <- paste(hdir,"/aridity",sep="")
if (!file.exists(adir)) {dir.create(adir)}

#given classes
m <- c(0, 40, 1,  
       40, 60, 2,
       60, 80, 3,
       80, 100, 4)
rclmat <- matrix(m, ncol=3, byrow=TRUE)

#get tmax, tmin, and prec from WCLv1.4
if (!file.exists(paste(adir,"/tmin_hist.tif",sep=""))) {
    tmin <- getData('worldclim', var='tmin', res=2.5, path=adir)
    tmin <- tmin * 0.1
    tmin <- writeRaster(tmin, paste(adir,"/tmin_hist.tif",sep=""))
} else {
    tmin <- stack(paste(adir,"/tmin_hist.tif",sep=""))
}

if (!file.exists(paste(adir,"/tmax_hist.tif",sep=""))) {
    tmax <- getData('worldclim', var='tmax', res=2.5, path=adir)
    tmax <- tmax * 0.1
    tmax <- writeRaster(tmax, paste(adir,"/tmax_hist.tif",sep=""))
} else {
    tmax <- stack(paste(adir,"/tmax_hist.tif",sep=""))
}

if (!file.exists(paste(adir,"/tmean_hist.tif",sep=""))) {
    tmean <- getData('worldclim', var='tmean', res=2.5, path=adir)
    tmean <- tmean * 0.1
    tmean <- writeRaster(tmean, paste(adir,"/tmean_hist.tif",sep=""))
} else {
    tmean <- stack(paste(adir,"/tmean_hist.tif",sep=""))
}

if (!file.exists(paste(adir,"/prec_hist.tif",sep=""))) {
    prec <- getData('worldclim', var='prec', res=2.5, path=adir)
    prec <- writeRaster(prec, paste(adir,"/prec_hist.tif",sep=""))
} else {
    prec <- stack(paste(adir,"/prec_hist.tif",sep=""))
}

if (file.exists(paste(adir,"/wc2-5",sep=""))) {unlink(paste(adir,"/wc2-5",sep=""), recursive=TRUE, force=TRUE)}

#get srad from WCLv2.0, in kJ m-2 day-1
if (!file.exists(paste(adir,"/srad.tif",sep=""))) {
    setwd(adir)
    system("unrar x ET_SolRad.rar")
    srad <- stack(paste(srdir,"/ET_SolRad/et_solrad_",1:12,sep=""))
    names(srad) <- paste("srad.",1:12,sep="")
    srad <- resample(srad, prec, method='bilinear')
    srad <- writeRaster(srad, paste(adir,"/srad.tif",sep=""))
    unlink(paste(adir,"/ET_SolRad",sep=""), recursive=TRUE, force=TRUE)
} else {
    srad <- stack(paste(adir,"/srad.tif",sep=""))
}

#diurnal t. range
if (!file.exists(paste(adir,"/dtr_hist.tif",sep=""))) {
    dtr <- c()
    for (i in 1:12) {
        dtr <- c(dtr, (tmax[[i]] - tmin[[i]]))
    }
    dtr <- stack(dtr)
    names(dtr) <- paste("dtr.",1:12,sep="")
    dtr <- writeRaster(dtr, paste(adir,"/dtr_hist.tif",sep=""))
} else {
    dtr <- stack(paste(adir,"/dtr_hist.tif",sep=""))
}

#calculate PET
if (!file.exists(paste(adir,"/pet_hist.tif",sep=""))) {
    pet <- monthlyPET(tmean, srad, dtr)
    pet <- writeRaster(pet, paste(adir,"/pet_hist.tif",sep=""))
} else {
    pet <- stack(paste(adir,"/pet_hist.tif",sep=""))
}

#calculate Thornthwaite aridity index
#Thornthwaite aridity index = 100 * (d / n)
#d = sum of monthly differences between prec and PET for months where precip < PET
#n = sum of monthly PET for those months
if (!file.exists(paste(adir,"/aridity_thornthwaite_hist.tif",sep=""))) {
    aindex <- aridityIndexThornthwaite(prec, pet)
    aindex <- writeRaster(aindex, paste(adir,"/aridity_thornthwaite_hist.tif",sep=""))
}

#clean up all objects except srad
rm(list=c("aindex","pet","dtr","tmin","tmax","tmean","prec"))
gc(T)

##### do the future scenarios
#MOHC_HADGEM2_ES (29), CESM1_CAM5 (6), GFDL_CM3 (13), MPI ESM_LR (30), MIROC_MIROC5 (27)
#for RCP4.5, RCP8.5, 2030 and 2050
for (sce in c(4.5, 8.5)) {
    for (yr in c(2030, 2050)) {
        for (k in c(6,13,27,29,30)) {
            #sce <- 4.5; yr <- 2030; k <- 6
            cat("processing rcp",sce,"/ year=",yr,"/ model=",k,"\n")
            
            #get all data
            if (!file.exists(paste(adir,"/tmin_fut_m",k,"_",yr,"_rcp",sce,".tif",sep=""))) {
                tmin <- getCMIP5(var="tmin", rcp=sce, model=k, year=yr, res=2.5, path=adir)
                tmin <- tmin * 0.1
                tmin <- writeRaster(tmin, paste(adir,"/tmin_fut_m",k,"_",yr,"_rcp",sce,".tif",sep=""))
                names(tmin) <- paste("tmin.",1:12,sep="")
            } else {
                tmin <- stack(paste(adir,"/tmin_fut_m",k,"_",yr,"_rcp",sce,".tif",sep=""))
            }
            
            if (!file.exists(paste(adir,"/tmax_fut_m",k,"_",yr,"_rcp",sce,".tif",sep=""))) {
                tmax <- getCMIP5(var="tmax", rcp=sce, model=k, year=yr, res=2.5, path=adir)
                tmax <- tmax * 0.1
                tmax <- writeRaster(tmax, paste(adir,"/tmax_fut_m",k,"_",yr,"_rcp",sce,".tif",sep=""))
                names(tmax) <- paste("tmax.",1:12,sep="")
            } else {
                tmax <- stack(paste(adir,"/tmax_fut_m",k,"_",yr,"_rcp",sce,".tif",sep=""))
            }
            
            if (!file.exists(paste(adir,"/prec_fut_m",k,"_",yr,"_rcp",sce,".tif",sep=""))) {
                prec <- getCMIP5(var="prec", rcp=sce, model=k, year=yr, res=2.5, path=adir)
                prec <- writeRaster(prec, paste(adir,"/prec_fut_m",k,"_",yr,"_rcp",sce,".tif",sep=""))
                names(prec) <- paste("prec.",1:12,sep="")
            } else {
                prec <- stack(paste(adir,"/prec_fut_m",k,"_",yr,"_rcp",sce,".tif",sep=""))
            }
            
            if (!file.exists(paste(adir,"/dtr_fut_m",k,"_",yr,"_rcp",sce,".tif",sep=""))) {
                dtr <- c(); for (i in 1:12) {dtr <- c(dtr, (tmax[[i]] - tmin[[i]]))}
                dtr <- stack(dtr)
                names(dtr) <- paste("dtr.",1:12,sep="")
                dtr <- writeRaster(dtr, paste(adir,"/dtr_fut_m",k,"_",yr,"_rcp",sce,".tif",sep=""))
            } else {
                dtr <- stack(paste(adir,"/dtr_fut_m",k,"_",yr,"_rcp",sce,".tif",sep=""))
            }
            
            if (!file.exists(paste(adir,"/tmean_fut_m",k,"_",yr,"_rcp",sce,".tif",sep=""))) {
                tmean <- c()
                for (i in 1:12) {tmean <- c(tmean, (tmax[[i]] * 0.5 + tmin[[i]] * 0.5))}
                tmean <- stack(tmean)
                names(tmean) <- paste("tmean.",1:12,sep="")
                tmean <- writeRaster(tmean, paste(adir,"/tmean_fut_m",k,"_",yr,"_rcp",sce,".tif",sep=""))
            } else {
                tmean <- stack(paste(adir,"/tmean_fut_m",k,"_",yr,"_rcp",sce,".tif",sep=""))
            }
            
            #delete ascii files
            if (file.exists(paste(adir,"/cmip5",sep=""))) {unlink(paste(adir,"/cmip5",sep=""), recursive=TRUE, force=TRUE)}
            
            #calculate PET
            if (!file.exists(paste(adir,"/pet_fut_m",k,"_",yr,"_rcp",sce,".tif",sep=""))) {
                pet <- monthlyPET(tmean, srad, dtr)
                pet <- writeRaster(pet, paste(adir,"/pet_fut_m",k,"_",yr,"_rcp",sce,".tif",sep=""))
            } else {
                pet <- stack(paste(adir,"/pet_fut_m",k,"_",yr,"_rcp",sce,".tif",sep=""))
            }

            #calculate Thornthwaite aridity index
            if (!file.exists(paste(adir,"/aridity_thornthwaite_fut_m",k,"_",yr,"_rcp",sce,".tif",sep=""))) {
                aindex <- aridityIndexThornthwaite(prec, pet)
                aindex <- writeRaster(aindex, paste(adir,"/aridity_thornthwaite_fut_m",k,"_",yr,"_rcp",sce,".tif",sep=""))
            }
            
            #clean up all objects except srad
            rm(list=c("aindex","pet","dtr","tmin","tmax","tmean","prec"))
            gc(T)
        }
        
        #calculate muti-model mean of aridity index
        cat("calculate multi-model mean\n")
        if (!file.exists(paste(adir,"/aridity_thornthwaite_fut_mmm_",yr,"_rcp",sce,".tif",sep=""))) {
            aindex <- stack(paste(adir,"/aridity_thornthwaite_fut_m",c(6,13,27,29,30),"_",yr,"_rcp",sce,".tif",sep=""))
            aindex <- mean(aindex, na.rm=T)
            aindex <- writeRaster(aindex, paste(adir,"/aridity_thornthwaite_fut_mmm_",yr,"_rcp",sce,".tif",sep=""))
            rm(aindex); gc(T)
        }
    }
}

setwd(adir)
system("tar -cjvf aridity_thornthwaite.tar.bz2 aridity_thornthwaite_*.tif")
system("mv aridity_thornthwaite.tar.bz2 ./../.")
