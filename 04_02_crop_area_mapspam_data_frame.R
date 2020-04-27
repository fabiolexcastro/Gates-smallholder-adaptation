g=gc()
rm(list=ls())
crop <- "BEAN"
continent <- "South America"

shp  <- shapefile("//dapadfs/workspace_cluster_13/GATES/data/shp/base/continents.shp")
shp_oce  <- shp[which(shp@data$CONTINENT == continent),]
root <- "//dapadfs/workspace_cluster_13/GATES/" 

base <- raster("//dapadfs/workspace_cluster_13/GATES/raster/indicators/dryDays/monthly/dry_days_1981_1.asc")
cat(paste0("cultivos Mapspam------>",crop,"\n"))
mapspam  <- raster(paste0("//dapadfs/workspace_cluster_13/GATES/raster/mapSPAM_BM/spam2010V1r1_global_H_",crop,"_A.tif"))
mapspam  <- raster::resample(mapspam, base)
oce   <- crop(mapspam, shp_oce)
oce   <- mask(oce, shp_oce)
oce[oce > 0 ] <- 1
oce[oce <= 0] <- NA
map_oce  <- as.data.frame(rasterToPoints(oce))
map_oce <- na.omit(map_oce)
cell <- cellFromXY(base, map_oce[,1:2])
map_oce <- cbind(cell, map_oce)
colnames(map_oce) <- c("cell","lon","lat","area")

saveRDS(map_oce,paste0(root,"raster/mapSPAM_BM/data.frame/",crop,"_",continent,".rds"))


