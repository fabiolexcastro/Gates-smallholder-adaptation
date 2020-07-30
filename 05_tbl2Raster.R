

# 
require(pacman)
pacman::p_load(raster, rgdal, rgeos, stringr, velox, sf, tidyverse)
rm(list = ls())

fls <- list.files('../rds/chirps/dryDays', full.names = TRUE)
yrs <- 2010:2019

fle <- grep('2019', fls, value = T)


# 
rst2tbl <- function(x){
  # x <- yrs[1]
  fle <- grep(x, fls, value = T)
  lapply(1:length(fle), function(k){
    print(fle[k])
    fl <- readRDS(fle[k])  
    ly <- rasterFromXYZ(fl)
    writeRaster(ly, paste0('../raster/indicators/dryDays/', gsub('.rds', '.tif', basename(fle[k]))), overwrite = TRUE)
  })
  print('Done!') 
}

map(.x = yrs, .f = rst2tbl)

