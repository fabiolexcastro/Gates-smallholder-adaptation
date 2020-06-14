#JRV, April 2020
#ecocrop analysis for cassava, test for A. Ghosh

#libraries
library(raster)
library(analogues)
library(rgdal)

#working directory
wd <- "~/work"

#source EcoCrop function
source(paste(wd,"/Gates-smallholder-adaptation/crop_suitability/EcoCrop-model_WCl.R",sep=""))

#hazard dir
hdir <- paste(wd,"/hazard_layers/ecocrop",sep="")
if (!file.exists(hdir)) {dir.create(hdir)}

#load global shapefile
shp <- readOGR(paste(wd,"/world_shp/ne_50m_admin_0_countries.shp",sep=""))

#get all worldclim data
for (tvar in c("prec","tmin","tmean","tmax")) {
  #tvar <- "prec"
  if (!file.exists(paste(hdir,"/wc_",tvar,".rda",sep=""))) {
    wc_data <- getData('worldclim', var=tvar, res=5, path=hdir)
    wc_data <- readAll(wc_data)
    assign(paste("wc_",tvar,sep=""),wc_data)
    rm(wc_data)
    save(list=c(paste("wc_",tvar,sep="")), file=paste(hdir,"/wc_",tvar,".rda",sep=""), compress="xz",compression_level=9)
    unlink(paste(hdir,"/wc5",sep=""),recursive=T,force=T)
  }
}

#run ecocrop
#write climate files first
cldir <- paste(hdir,"/wc_data",sep="")
if (!file.exists(cldir)) {dir.create(cldir)}

#loop variables to write rasters
for (tvar in c("prec","tmin","tmean","tmax")) {
  cat("writing wcl var=",tvar,"\n")
  load(paste(hdir,"/wc_",tvar,".rda",sep=""))
  x <- get(paste("wc_",tvar,sep=""))
  rm(list=c(paste("wc_",tvar,sep="")))
  for (j in 1:12) {
    if (!file.exists(paste(cldir,"/",tvar,"_",j,".tif",sep=""))) {
      writeRaster(x[[j]], paste(cldir,"/",tvar,"_",j,".tif",sep=""))
    }
  }
  rm(x)
}

#run ecocrop (Ceballos et al. parameters)
if (!file.exists(paste(wd,"/wcl_suit/cassava_suit.tif",sep=""))) {
  sx <- suitCalc(climPath=cldir,Gmin=240,Gmax=240,
                 Tkmp=0,Tmin=150,Topmin=220,Topmax=320,Tmax=450,
                 Rmin=300,Ropmin=800,Ropmax=2200,Rmax=2800,
                 outfolder=paste(hdir,"/wc_suit",sep=""),ext=".tif",cropname='cassava')
}

#load data
wc_sx <- raster(paste(hdir,"/wc_suit/cassava_suit.tif",sep=""))
wc_sx[which(wc_sx[] == 0)] <- NA
plot(wc_sx, breaks=seq(0,100,by=20), col=rev(terrain.colors(5)))

#delete monthly data
#unlink(paste(wd,"/wc_data",sep=""),recursive=T,force=T)
#unlink(paste(wd,"/c5_data",sep=""),recursive=T,force=T)


