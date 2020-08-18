

# Load libraries ----------------------------------------------------------
require(pacman)
pacman::p_load(raster, rgdal, rgeos, stringr, tidyverse, gtools, foreach, doSNOW, parallel)

g <- gc(reset = TRUE)
rm(list = ls())
options(scipen = 999)

# Functions to use --------------------------------------------------------
add_change <- function(yar, rcp, prd, mnt){
  
  # Proof
  yar <- 2004
  rcp <- 'rcp45'
  prd <- '2020_2049'
  mnt <- 6
  
  print('History data')
  chr <- grep(yar, chrp, value = TRUE)
  chr <- grep(paste0('_', ifelse(mnt < 10, paste0('0', mnt), mnt), '_'), chr, value = TRUE)
  
  print('Future data')
  prc <- grep(paste0(rcp, '_', prd), prcn, value = TRUE)
  prc <- grep(paste0('_', mnt, '_'), prc, value = TRUE)
  prc <- raster(prc)
  prc <- resample(x = prc, y = mask, method = 'ngb')
  
  yar <- date %>% 
    filter(year_1 == yar & period == prd) %>% 
    pull(year_2)
  
  for(i in 26:length(chr)){
    
    print(chr[i])
    
    ch <- grep(paste0('_', ifelse(i < 10, paste0('0', i), i), '.tif'), chr, value = TRUE) 
    ch <- raster(ch)
    st <- stack(ch, prc)
    tb <- rasterToPoints(st)
    tb <- as_tibble(tb)
    colnames(tb) <- c('x', 'y', 'base', 'change')
    tb <- tb %>% 
      mutate(future = round(base + round(base * change /  100, 1), 1)) %>% 
      mutate(future = ifelse(future < 0, 0, future))
    rs <- rasterFromXYZ(tb[,c(1, 2, 5)])
    writeRaster(x = rs,
                filename = paste0('../tif/future/', rcp, '_', prd, '/ppt_', yar, '_', mnt, '_', i, '.tif'),
                overwrite = TRUE)
    
  }
  
}

# Load data ---------------------------------------------------------------
root <- '//dapadfs/workspace_cluster_13/GATES'
shpf <- shapefile(paste0(root, '/data/shp/base/africa_continent.shp'))

# CHIRPS Data
chrp <- paste0(root, '/workspace/dryDays/tif/chirps/continents') %>% 
  list.files(., full.names = TRUE, pattern = '.tif$') %>% 
  grep('Africa', ., value = TRUE)
year <- 1981:2019
mask <- raster(chrp[1])
mask <- mask * 0

# Climate change data
prcn <- list.files('../tif/difference_wc/porc', full.names = TRUE, pattern = '.tif$')
rcps <- c('rcp45', 'rcp85')
prds <- c('2020_2049', '2040_2069')

# Dates
date <- data.frame(period = c(rep(c('2020_2049'), length(1981:2010)), rep(c('2040_2069'), length(1981:2010))),
                   year_1 = c(1981:2010, 1981:2010),
                   year_2 = c(2020:2049, 2040:2069))

# Apply the function ------------------------------------------------------

##
cl <- makeCluster(10)
registerDoSNOW(cl)

foreach(j = 1:length(year), .verbose = TRUE) %dopar% {
  
  require(pacman)
  pacman::p_load(raster, rgdal, rgeos, stringr, tidyverse, gtools, foreach, doSNOW, parallel)
  
  foreach(k = 1:12) %do% {
    
    add_change(yar = year[j], rcp = 'rcp85', prd = '2020_2049', mnt = k)
    
  }
  
}

stopCluster(cl)


cl <- makeCluster(10)
registerDoSNOW(cl)

foreach(j = 1:length(year), .verbose = TRUE) %dopar% {
  
  require(pacman)
  pacman::p_load(raster, rgdal, rgeos, stringr, tidyverse, gtools, foreach, doSNOW, parallel)
  
  foreach(k = 1:12) %do% {
    
    add_change(yar = year[j], rcp = 'rcp85', prd = '2040_2069', mnt = k)
    
  }
  
}

stopCluster(cl)

cl <- makeCluster(10)
registerDoSNOW(cl)

foreach(j = 2:length(year), .verbose = TRUE) %dopar% {
  
  require(pacman)
  pacman::p_load(raster, rgdal, rgeos, stringr, tidyverse, gtools, foreach, doSNOW, parallel)
  
  foreach(k = 1:12) %do% {
    
    add_change(yar = year[j], rcp = 'rcp45', prd = '2020_2049', mnt = k)
    
  }
  
}

stopCluster(cl)




