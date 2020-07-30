#JRV - Fix extent of number of hot days
#July 2020
#libraries
g=gc(); rm(list=ls())
library(raster); library(tidyverse);library(tibble)
dim <- function(rcp, periodo){
path <- "//dapadfs/workspace_cluster_13/GATES/raster/indicators/heat_stress/future/"  
#load mask
cat("LOAD MASK", "\n")

msk <- raster("//dapadfs/workspace_cluster_13/GATES/workspace/change_temperature/tif/future/rcp45_2020_2049/tmax_2020_01_1.tif")
msk[which(!is.na(msk[]))] <- 1
msk_vls <- raster::rasterToPoints(msk) %>%
  as_tibble() %>% 
  setNames(c('x', 'y', 'value')) %>%
  dplyr::select(., x, y)

#list of crops
crops <- list.files(paste0(path,rcp,"_",periodo,"/"),pattern = ".tif")
crop_list <- list.files(paste0(path,rcp,"_",periodo,"/"),pattern = ".tif", full.names = T)
lapply(1 :length(crop_list), function(i){
     cat(paste0("Procesando archivo:: " ,i, " de ", length(crop_list), "\n"))
     file <- raster(crop_list[i])
     msk_vls <- msk_vls %>% mutate(value = raster::extract(file, msk_vls[1:2]))
     out_rs <- raster(msk)
     out_rs[cellFromXY(msk, as.data.frame(msk_vls[,c('x','y')]))]  <- msk_vls$value
     out_rs <- writeRaster(out_rs,paste0(path,rcp,"_",periodo,"/correccion_2/",crops[i]))
 })
}

dim(rcp= "rcp85", periodo="2020_2049")



