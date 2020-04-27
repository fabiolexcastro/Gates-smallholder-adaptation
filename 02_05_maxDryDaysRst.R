
# Load libraries ------------------------------------------------------------
require(pacman)
pacman::p_load(raster, rgdal, rgeos, stringr, terra, sf, tidyverse, glue, gtools, hablar, readxl)
rm(list = ls())

# Functions to use --------------------------------------------------------
tbl2rst <- function(cn, yr){
  # yr <- yrs[1]
  # cn <- 'North America'
  print(yr)
  po <- cnt %>%
    filter(continent == cn) %>% 
    pull(1)
  tb <- tbl[[po]] %>% 
    filter(year == yr)
  rs <- rasterFromXYZ(tb[,c(3,4,2)])
  writeRaster(rs, paste0('../raster/indicators/drySeasonYear/drySsnYear_', yr, '_', cn, '.asc'), overwrite = T)
  print('Done...!')
}
mosaicRaster <- function(yr){
  
  fl <- grep(yr, fls, value = T) 
  ly <- map(.x = fl, .f = raster)
  
  print(length(ly))
  
  ms1 <- mosaic(x = ly[[1]], y = ly[[2]], fun = 'mean')
  ms2 <- mosaic(x = ly[[3]], y = ly[[4]], fun = 'mean')
  ms3 <- mosaic(x = ms1, y = ms2, fun = 'mean')
  ms4 <- mosaic(x = ms3, y = ly[[5]], fun = 'mean')
  
  print('To write')
  writeRaster(ms4, paste0('../raster/indicators/drySeasonYear/mosaic/drySsnYear_', yr, '.tif'), overwrite = TRUE)
  print('Done!!!')
}


# Load Data ---------------------------------------------------------------
tbl <- list.files('../rds/chirps/consDaysYearly', full.names = T) %>% 
  map(.x = ., .f = readRDS)
cnt <- data.frame(position = 1:5, continent = c('Africa', 'Asia', 'North America', 'Oceania', 'South America'))
yrs <- as.character(1981:2019)

lapply(1:length(yrs), function(k) tbl2rst(cn = 'Africa', yr = yrs[k]))
lapply(1:length(yrs), function(k) tbl2rst(cn = 'Oceania', yr = yrs[k]))
lapply(1:length(yrs), function(k) tbl2rst(cn = 'North America', yr = yrs[k]))
lapply(1:length(yrs), function(k) tbl2rst(cn = 'Asia', yr = yrs[k]))
lapply(1:length(yrs), function(k) tbl2rst(cn = 'South America', yr = yrs[k]))

fls <- list.files('../raster/indicators/drySeasonYear', full.names = T, pattern = '.asc')
map(.x = yrs, .f = mosaicRaster)
