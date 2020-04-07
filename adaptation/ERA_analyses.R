#JRV - analogues analysis for bean sites

#working directory
wd <- "~/Documents/Gates-Ag-Adapt-Atlas"
cimdir <- paste(wd,"/ERA_comparison",sep="")

#load library
library(analogues)
library(rgdal)
library(raster)
library(maptools); data(wrld_simpl)

#read ERA data file
all_data <- read.csv(paste(cimdir,"/ERA_Test_Data.csv",sep=""))
all_data <- all_data[which(all_data$Product.Simple == "Maize"),]
all_data <- all_data[which(all_data$PrName %in% c("Mulch","Reduced Tillage","Green Manure","Water Harvesting","Agroforestry Pruning")),]
all_data <- droplevels(all_data)

#get Africa shapefile
sh_ctry <- readOGR(paste("~/Documents/analogues/avisa_sites/Bean_Corridor_Data/Africa_Bean_corridor/African_continet.shp",sep=""))
sh_ctry <- spTransform(sh_ctry, crs("+proj=longlat +ellps=WGS84 +no_defs"))
sh_xt <- extent(sh_ctry)
sh_xt@xmin <- sh_xt@xmin-1; sh_xt@ymin <- sh_xt@ymin-1; sh_xt@xmax <- sh_xt@xmax+1; sh_xt@ymax <- sh_xt@ymax+1

#get worldclim climate data
for (tvar in c("prec","tmin","tmean","tmax")) {
  #tvar <- "prec"
  cat("downloading",tvar,"\n")
  if (!file.exists(paste(cimdir,"/wc_",tvar,".rda",sep=""))) {
    wc_data <- getData('worldclim', var=tvar, res=2.5, path=cimdir)
    wc_data <- crop(wc_data, sh_xt)
    wc_data <- readAll(wc_data)
    assign(paste("wc_",tvar,sep=""), wc_data)
    rm(wc_data)
    save(list=c(paste("wc_",tvar,sep="")), file=paste(cimdir,"/wc_",tvar,".rda",sep=""), compress="xz",compression_level=9)
    unlink(paste(cimdir,"/wc2-5",sep=""),recursive=T,force=T)
  }
}

#create a mask
load(paste(cimdir,"/wc_prec.rda",sep=""))
msk <- wc_prec[[1]]; rm(wc_prec)
msk[which(!is.na(msk[]))] <- 1

#select practice and organize data
#Practice (n) as follows
#   1. Reduced Tillage (47)
#   2. Mulch (42)
#   3. Agroforestry Pruning (38)
#   4. Green Manure (37)
#   5. Water Harvesting (35)

prname <- "Mulch"
pdata <- all_data[which(all_data$PrName == prname),]
pdata <- pdata[which(!is.na(pdata$Latitude)),]
pdata <- pdata[which(!is.na(pdata$Longitude)),]
pdata <- droplevels(pdata)
pdata <- pdata[,c("Longitude","Latitude","yi")]
pdata <- stats::aggregate(pdata[,c("yi")], list(Lon=pdata$Longitude,Lat=pdata$Latitude), mean, na.rm=T)

#cleanup any data points outside African continent
pdata$isvalid <- raster::extract(msk, pdata[,c("Lon","Lat")])
pdata <- pdata[which(!is.na(pdata$isvalid)),]
pdata$isvalid <- NULL

#positive and negative yield impact
pdata_p <- pdata[which(pdata$x >= 0.15),]
pdata_p$ID <- paste("s",1:nrow(pdata_p),sep="")
pdata_n <- pdata[which(pdata$x < 0.15),]
pdata_n$ID <- paste("s",1:nrow(pdata_n),sep="")

#plot
par(mar=c(5,5,1,1))
plot(pdata$Lon, pdata$Lat, pch=20, col="white")
plot(msk,add=T,legend=F)
plot(sh_ctry, add=T)
points(pdata_p$Lon, pdata_p$Lat, pch=20, col="blue", cex=1.5)
points(pdata_n$Lon, pdata_n$Lat, pch=20, col="red", cex=1.5)

pr_odir <- paste(cimdir,"/",gsub(" ","_",prname,fixed=T),sep="")
if (!file.exists(pr_odir)) {dir.create(pr_odir,recursive=T)}

#load climate data for calculation
load(paste(cimdir,"/wc_prec.rda",sep="")); wc_prec <- stack(wc_prec)
load(paste(cimdir,"/wc_tmean.rda",sep="")); wc_tmean <- stack(wc_tmean)

#calculate similarity (both temp and precip for positive sites)
for (i in 1:nrow(pdata_p)) {
  #i <- 1
  cat("processing site=",i,"\n")
  if (!file.exists(paste(pr_odir,"/out_pos_",pdata_p$ID[i],".tif",sep=""))) {
    params <- createParameters(x=pdata_p$Lon[i], y=pdata_p$Lat[i], vars=c("prec","tmean"),weights=c(0.70,0.30),
                               ndivisions=c(12,12),growing.season=c(1,12),
                               rotation="prec",threshold=1,env.data.ref=list(wc_prec,wc_tmean), 
                               env.data.targ=list(wc_prec,wc_tmean),outfile=pr_odir,
                               fname=paste("pos_",pdata_p$ID[i],sep=""),writefile=T)
    sim_out <- calc_similarity(params)
  }
}

#calculate maximum similarity across sites
for (i in 1:nrow(pdata_p)) {
  #i <- 1
  cat("...processing site=",i,"\n")
  rs <- raster(paste(pr_odir,"/out_pos_",pdata_p$ID[i],".tif",sep=""))
  if (i == 1) {
    maxsim <- rs
  } else {
    maxsim <- calc(stack(maxsim, rs), fun=max)
  }
}

#write max similarity rasters
writeRaster(maxsim, paste(pr_odir,"/max_similarity_pos.tif",sep=""),overwrite=T)

#cut similarity to growing areas
aharv <- raster(paste(wd,"/spam2010v1r1_global_harv_area.geotiff/spam2010V1r1_global_H_MAIZ_A.tif",sep=""))
aharv <- crop(aharv, maxsim)
aharv <- resample(aharv, maxsim)
aharv[which(aharv[] < 0.01)] <- NA
aharv[which(!is.na(aharv[]))] <- 1
maxsim_m <- mask(maxsim, aharv)
writeRaster(maxsim_m, paste(pr_odir,"/max_similarity_neg_maizeareas.tif",sep=""),overwrite=T)

#plot
par(mar=c(5,5,1,1))
plot(pdata$Lon, pdata$Lat, pch=20, col="white")
plot(maxsim_m,add=T,legend=F,zlim=c(0,1),breaks=seq(0,1,by=0.2),col=rev(terrain.colors(6)))
plot(sh_ctry, add=T)
points(pdata_p$Lon, pdata_p$Lat, pch=20, col="blue", cex=1)


###
#calculate similarity (both temp and precip for negative sites)
for (i in 1:nrow(pdata_n)) {
  #i <- 1
  cat("processing site=",i,"\n")
  if (!file.exists(paste(pr_odir,"/out_neg_",pdata_n$ID[i],".tif",sep=""))) {
    params <- createParameters(x=pdata_n$Lon[i], y=pdata_n$Lat[i], vars=c("prec","tmean"),weights=c(0.70,0.30),
                               ndivisions=c(12,12),growing.season=c(1,12),
                               rotation="prec",threshold=1,env.data.ref=list(wc_prec,wc_tmean), 
                               env.data.targ=list(wc_prec,wc_tmean),outfile=pr_odir,
                               fname=paste("neg_",pdata_n$ID[i],sep=""),writefile=T)
    sim_out <- calc_similarity(params)
  }
}

#calculate maximum similarity across sites
for (i in 1:nrow(pdata_n)) {
  #i <- 1
  cat("...processing site=",i,"\n")
  rs <- raster(paste(pr_odir,"/out_neg_",pdata_n$ID[i],".tif",sep=""))
  if (i == 1) {
    maxsim <- rs
  } else {
    maxsim <- calc(stack(maxsim, rs), fun=max)
  }
}

#write max similarity rasters
writeRaster(maxsim, paste(pr_odir,"/max_similarity_neg.tif",sep=""),overwrite=T)

#cut similarity to growing areas
aharv <- raster(paste(wd,"/spam2010v1r1_global_harv_area.geotiff/spam2010V1r1_global_H_MAIZ_A.tif",sep=""))
aharv <- crop(aharv, maxsim)
aharv <- resample(aharv, maxsim)
aharv[which(aharv[] < 0.01)] <- NA
aharv[which(!is.na(aharv[]))] <- 1
maxsim_m <- mask(maxsim, aharv)
writeRaster(maxsim_m, paste(pr_odir,"/max_similarity_neg_maizeareas.tif",sep=""),overwrite=T)

#plot
par(mar=c(5,5,1,1))
plot(pdata$Lon, pdata$Lat, pch=20, col="white")
plot(maxsim_m,add=T,legend=F,zlim=c(0,1),breaks=seq(0,1,by=0.2),col=rev(terrain.colors(6)))
plot(sh_ctry, add=T)
points(pdata_n$Lon, pdata_n$Lat, pch=20, col="red", cex=1)



