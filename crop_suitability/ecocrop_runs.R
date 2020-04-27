#JRV, Feb 2020
#ecocrop analysis for bean

#libraries
library(raster); library(analogues)

#working directory
wd <- "~/Documents/CIMMYT-frijol-para-Mexico/ecocrop_beans"

#get Mex shapefile
if (!file.exists(paste(wd,"/gadm36_MEX_0_sp.rds",sep=""))) {
  shp <- getData('GADM', country='MEX', level=0, path=wd)
} else {
  shp <- readRDS(file=paste(wd,"/gadm36_MEX_0_sp.rds",sep=""))
}

#extent
xt <- extent(shp)
xt@xmin <- xt@xmin - 2; xt@ymin <- xt@ymin - 2
xt@xmax <- xt@xmax + 2; xt@ymax <- xt@ymax + 2

#get, crop, mask and save global worldclim data
wc_prec <- getData('worldclim', var='prec', res=2.5, path=wd)
wc_prec <- crop(wc_prec, xt)

#create mask
if (!file.exists(paste(wd,"/mask.rda",sep=""))) {
  msk <- rasterize(shp, wc_prec)
  msk <- readAll(x)
  save(list=c("msk"), file=paste(wd,"/mask.rda",sep=""), compress="xz",compression_level=9)
} else {
  load(paste(wd,"/mask.rda",sep=""))
}

if (!file.exists(paste(wd,"/wc_prec.rda",sep=""))) {
  wc_prec <- mask(wc_prec, msk)
  wc_prec <- readAll(wc_prec)
  save(list=c("wc_prec"), file=paste(wd,"/wc_prec.rda",sep=""), compress="xz",compression_level=9)
} else {
  load(file=paste(wd,"/wc_prec.rda",sep=""))
}

#delete folder
unlink(paste(wd,"/wc2-5",sep=""),recursive=T,force=T)

for (tvar in c("tmin","tmean","tmax")) {
  #tvar <- "tmin"
  if (!file.exists(paste(wd,"/wc_",tvar,".rda",sep=""))) {
    wc_data <- getData('worldclim', var=tvar, res=2.5, path=wd)
    wc_data <- crop(wc_data, xt)
    wc_data <- mask(wc_data, msk)
    wc_data <- readAll(wc_data)
    assign(paste("wc_",tvar,sep=""),wc_data)
    rm(wc_data)
    save(list=c(paste("wc_",tvar,sep="")), file=paste(wd,"/wc_",tvar,".rda",sep=""), compress="xz",compression_level=9)
    unlink(paste(wd,"/wc2-5",sep=""),recursive=T,force=T)
  }
}

#climate change models
modlist <- cmip5_table$id[which(cmip5_table$rcp8_5==1)]
modlist <- modlist[which(!modlist %in% c(7,22))]
for (i in modlist) {
  for (tvar in c("prec","tmin","tmean","tmax")) {
    cat("model i=",i,"/ variable=",tvar,"\n")
    if (!file.exists(paste(wd,"/c5_",tvar,"_m",i,".rda",sep=""))) {
      x <- getCMIP5(var=tvar,model=i,res=2.5,rcp=8.5,year=2030,path=wd)
      x <- crop(x, msk)
      x <- mask(x, msk)
      x <- readAll(x)
      assign(paste("c5_",tvar,"_m",i,sep=""), x)
      rm(x)
      save(list=c(paste("c5_",tvar,"_m",i,sep="")), file=paste(wd,"/c5_",tvar,"_m",i,".rda",sep=""), compress="xz",compression_level=9)
      unlink(paste(wd,"/cmip5",sep=""),recursive=T,force=T)
      rm(list=c(paste("c5_",tvar,"_m",i,sep="")))
    }
  }
}

#calculate average GCM
for (tvar in c("prec","tmin","tmean","tmax")) {
  if (!file.exists(paste(wd,"/c5_",tvar,".rda",sep=""))) {
    rsvar <- c()
    for (j in 1:12) {
      cat("averaging var=",tvar," / month=",j,"\n")
      rsj <- c() #list of all GCM rasters
      for (i in modlist) {
        load(paste(wd,"/c5_",tvar,"_m",i,".rda",sep=""))
        x <- get(paste("c5_",tvar,"_m",i,sep=""))
        rm(list=c(paste("c5_",tvar,"_m",i,sep="")))
        rsj <- c(rsj,x[[j]])
        rm(x)
      }
      rsj <- stack(rsj)
      rsj_m <- calc(rsj, mean, na.rm=T)
      rsvar <- c(rsvar, rsj_m)
      rm(rsj); rm(rsj_m)
    }
    rsvar <- stack(rsvar)
    rsvar <- readAll(rsvar)
    assign(paste("c5_",tvar,sep=""), rsvar)
    rm(rsvar)
    save(list=c(paste("c5_",tvar,sep="")), file=paste(wd,"/c5_",tvar,".rda",sep=""), compress="xz",compression_level=9)
    rm(list=c(paste("c5_",tvar,sep="")))
  }
}

#run ecocrop
#write climate files first
source(paste(wd,"/EcoCrop-model_WCl.R",sep=""))
for (cty in c("wc","c5")) {
  cldir <- paste(wd,"/",cty,"_data",sep="")
  if (!file.exists(cldir)) {dir.create(cldir)}
  
  #loop variables to write rasters
  for (tvar in c("prec","tmin","tmean","tmax")) {
    cat("writing type=",cty,"/ var=",tvar,"\n")
    load(paste(wd,"/",cty,"_",tvar,".rda",sep=""))
    x <- get(paste(cty,"_",tvar,sep=""))
    rm(list=c(paste(cty,"_",tvar,sep="")))
    #if (tvar != "prec") {x <- x * 0.1}
    for (j in 1:12) {
      if (!file.exists(paste(cldir,"/",tvar,"_",j,".tif",sep=""))) {
        writeRaster(x[[j]], paste(cldir,"/",tvar,"_",j,".tif",sep=""))
      }
    }
    rm(x)
  }
  
  #run ecocrop (Beebe et al. parameters)
  if (!file.exists(paste(wd,"/",cty,"_suit/drybean1_suit.tif",sep=""))) {
    sx <- suitCalc(climPath=cldir,Gmin=90,Gmax=120,
                   Tkmp=0,Tmin=135,Topmin=175,Topmax=231,Tmax=256,
                   Rmin=200,Ropmin=363,Ropmax=450,Rmax=710,
                   outfolder=paste(wd,"/",cty,"_suit",sep=""),ext=".tif",cropname='drybean1')
  }
  
  #run ecocrop (JRV modified parameters to reduce Rmin and Ropmin as test for increased drought tolerance)
  if (!file.exists(paste(wd,"/",cty,"_suit/drybean2_suit.tif",sep=""))) {
    sx <- suitCalc(climPath=cldir,Gmin=90,Gmax=120,
                   Tkmp=0,Tmin=135,Topmin=175,Topmax=231,Tmax=256,
                   Rmin=100,Ropmin=300,Ropmax=450,Rmax=710,
                   outfolder=paste(wd,"/",cty,"_suit",sep=""),ext=".tif",cropname='drybean2')
  }
}

#drybean1
wc_sx <- raster(paste(wd,"/wc_suit/drybean1_suit.tif",sep=""))
c5_sx <- raster(paste(wd,"/c5_suit/drybean1_suit.tif",sep=""))
chg_sx1 <- c5_sx-wc_sx
chg_sx1[which(wc_sx[] == 0 & c5_sx[] == 0)] <- NA
plot(chg_sx1,breaks=c(-100,-50,-5,5,50,100),col=rev(terrain.colors(6)))

#drybean2
wc_sx <- raster(paste(wd,"/wc_suit/drybean2_suit.tif",sep=""))
c5_sx <- raster(paste(wd,"/c5_suit/drybean2_suit.tif",sep=""))
chg_sx2 <- c5_sx-wc_sx
chg_sx2[which(wc_sx[] == 0 & c5_sx[] == 0)] <- NA
plot(chg_sx2,breaks=c(-100,-50,-5,5,50,100),col=rev(terrain.colors(6)))

setwd(wd)
system("zip -r ecocrop_beans_mx.zip c5_suit wc_suit")

#delete monthly data
unlink(paste(wd,"/wc_data",sep=""),recursive=T,force=T)
unlink(paste(wd,"/c5_data",sep=""),recursive=T,force=T)


