
#
#library(pacman)
#pacman::p_load(tidyverse, stringr, raster, velox, foreach, doSNOW, glue)
library(tidyverse)
library(stringr)
library(raster)
library(velox)
library(glue)
setwd('/dapadfs/workspace_cluster_13/GATES/codes')
rm(list = ls())

# Functions
rst2tbl <- function(x, y){
  
  # Proof
  # x <- yrs[15]
  # y <- cnt[5]
  cutRaster <- function(pos, msk){
    cut <- raster::crop(r[[pos]], msk) %>% 
      raster::mask(., msk)
    cut[which(cut[] < 0)] <- NA
    return(cut)
  }
  
  # Starting
  r <- grep(x, fls, value = T) %>% stack() 
  y <- shp[shp@data$CONTINENT %in% y,]
  
  # Cutting
  r_cut <- lapply(1:nlayers(r), function(i){
    print(paste0('Processing -----> ', i))
    cutRaster(pos = i, msk = y)
  })
  r_cut <- stack(r_cut) 
  
  # Raster to Table
  t <- rasterToPoints(r_cut, spatial = FALSE) %>% 
    as_tibble()
  colnames(t) <- c('x', 'y', paste0('y_', str_sub(names(r), start = 13, end = nchar(names(r)))))
  colnames(t) <- gsub('\\.', '_', colnames(t))
  saveRDS(object = t, file = glue('../rds/chirps/{y@data$CONTINENT}_chirps_{x}.rds'))
  print('Done!') 
}

# Load
path <- '/dapadfs/data_cluster_4/observed/gridded_products/chirps/daily'
fls <- list.files(path, full.names = T, pattern = '.tif$')
shp <- shapefile('/dapadfs/workspace_cluster_13/GATES/data/shp/base/continents.shp')
yrs <- str_sub(string = fls, start = 77, end = 80) %>% unique() %>% as.numeric()

# Continents
cnt <- unique(shp$CONTINENT)
#lapply(1:length(yrs), function(k){print(yrs[k]); rst2tbl(x = yrs[k], y = cnt[3])})
lapply(6:length(yrs), function(k){print(yrs[k]); rst2tbl(x = yrs[k], y = cnt[5])})
