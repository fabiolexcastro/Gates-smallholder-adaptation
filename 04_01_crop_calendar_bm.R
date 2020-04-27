#Tests indicadores de temperatura "BEAN"
g=gc(); rm(list=ls())
require(raster);require(ggplot2); require(lubridate);require(tidyverse);require(dplyr)
base <- raster("//dapadfs/workspace_cluster_13/GATES/raster/indicators/dryDays/monthly/dry_days_1981_1.asc")
shp  <- shapefile("//dapadfs/workspace_cluster_13/GATES/data/shp/base/continents.shp")
root <- "//dapadfs/workspace_cluster_13/GATES/" 

continent <- "Africa"

heat_stress <- function(crop, continent){
  shp_oce  <- shp[which(shp@data$CONTINENT == continent),]
  
  ### Bases de datos de clima
  tmax     <- readRDS(paste0(root,"/rds/chirts/Tmax_",continent,"_chirts_1983.rds"))
  tmin     <- readRDS(paste0(root,"/rds/chirts/Tmin_",continent,"_chirts_1983.rds"))
  tmean <- (tmax[,3:ncol(tmax)] + tmin[,3:ncol(tmin)]) *0.5
  tmean <- cbind(tmax[,1:2], tmean)
  cellID <- cellFromXY(base, tmean[,1:2])
  tmean <- cbind(cellID, tmean)
  rm(tmax,tmin, cellID)
  
  
  ### Bases de datos de cultivos 
  mapspam  <- raster("//dapadfs/workspace_cluster_13/GATES/raster/mapSPAM/spam2010V1r1_global_H_BARL_A.tif")
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
  rm(mapspam, oce, cell)
  
  ### Filtro de base de clima con base mapspam 
  tmean_o <- dplyr::filter(tmean, cellID %in% map_oce$cell)
  
  
}

### Filtro de base de clima con base mapspam 
tmean_o <- dplyr::filter(tmean, cellID %in% map_oce$cell)

### Ahora si a cargar la base de calendario 

planting <- raster("//dapadfs/workspace_cluster_13/GATES/netCDF/crop_calendar/crops/filled/Barley.crop.calendar.fill.nc", varname = 'plant.start')
harvest  <- raster("//dapadfs/workspace_cluster_13/GATES/netCDF/crop_calendar/crops/filled/Barley.crop.calendar.fill.nc",varname = 'harvest.end')

tmean_o$Planting <- raster::extract(x = planting, y = tmean_o[,c("x", "y")])
tmean_o$Harvest  <- raster::extract(x = harvest,  y = tmean_o[,c("x", "y")]) 
tmean_o$Duration <- ifelse(test = tmean_o$Planting < tmean_o$Harvest, yes = "One year", no = "Two years")
tmean_o <- na.omit(tmean_o)

library(doSNOW) ; library(foreach) ; library(parallel); library(doParallel)
cores<- detectCores();cl<- makeCluster(cores-10);registerDoParallel(cl) 

system.time(indexes_been <- foreach(i=1:nrow(tmean_o)) %dopar%{ 
  mylist <- list()
  # Parameters
  duration <- tmean_o$Duration[i]
  start    <- tmean_o$Planting[i]
  end      <- tmean_o$Harvest[i]
  
  # Just one pixel
  time.serie <- tmean_o[i, 1:(ncol(tmean_o)-3)]
  nms <- gsub("y_", "",colnames(time.serie[,-c(1:3)]))
  nms <- gsub('\\_', '-', nms)
  nms <- c('cellID',"x", "y", nms)
  colnames(time.serie) <- nms
  
  if(duration== "One year"){
    suppressMessages(library(tidyverse))
    suppressMessages(library(compiler))
    
    X <- time.serie
    X <- X %>% gather(key = Date, value = Value, -(cellID:y))
    X$Year <- lubridate::year(as.Date(X$Date))
    X$Yday <- lubridate::yday(as.Date(X$Date))
    X <- X%>% dplyr::filter(Yday >= start & Yday <= end)
    heat <- sum(X$Value > 23.1)
    results<- data.frame(cellID= unique(X$cellID),lon= unique(X$x),lat= unique(X$y), value = heat, variable= "Heat_Stress")
    mylist[[i]] <- results
  }else{
    ###  
    X <- time.serie
    X <- X %>% gather(key = Date, value = Value, -(cellID:y))
    X$Year <- lubridate::year(as.Date(X$Date))
    X$Yday <- lubridate::yday(as.Date(X$Date))
    X <- X %>% filter(Yday %in% c(start:365, 1:end))
    
    heat <- sum(X$Value > 23.1)
    results<- data.frame(cellID= unique(X$cellID),lon= unique(X$x),lat= unique(X$y), value = heat, variable= "Heat_Stress")
    mylist[[i]] <- results
  }})

stopCluster(cl)
df <- do.call(rbind, indexes_been)
system.time(r <- rasterFromXYZ(df[,c(2:4)]))

