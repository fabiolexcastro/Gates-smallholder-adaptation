

# Load libraries ----------------------------------------------------------
require(pacman)
pacman::p_load(raster, rgdal, rgeos, stringr, sf, tidyverse, gtools,
               foreach, doSNOW, parallel)

g <- gc(reset = TRUE)
rm(list = ls())

# Functions to use --------------------------------------------------------
add_change <- function(year, rcpe, prdo, moth){
  
  year <- 1981
  rcpe <- 'rcp45'
  prdo <- '2020_2049'
  moth <- 5
  
  print('Future data')
  chng <- grep(paste0(rcpe, '_', prdo), prc, value = TRUE)
  chng <- grep(paste0('_', moth, '_'), chng, value = TRUE)
  chng <- raster(chng)
  
  print('History data')
  base <- grep(year, chr, value = TRUE)
  moth <- ifelse(moth < 10, paste0('0', moth), moth)
  base <- grep(paste0('_', moth, '_'), base, value = TRUE)
  
  print('Path base')
  path <- paste0('../tif/future/', paste0(rcpe, '_', prdo))
  year <- dates %>% filter(year_1 == year) %>% filter(period == prdo) %>% pull(year_2)
  
  print('Table to raster')
  lapply(1:length(base), function(k){
    
    print(k)
    k <- ifelse(k < 10, paste0('0', k), k)
    print(k)
    r <- grep(paste0('_', k, '.tif'), base, value = TRUE)
    r <- raster(r)
    
    mask <- r * 0
    chng <- raster::resample(chng, mask, method = 'ngb')
    
    r <- stack(r, chng)
    r <- rasterToPoints(r)
    r <- as_tibble(r)
    names(r) <- c('x', 'y', 'base', 'change')
    r <- r %>% 
      mutate(future = base * change,
             future = future / 100,
             future = round(future, 1),
             future = base + future)
    r <- rasterFromXYZ(r[,c(1, 2, 5)])
    writeRaster(x = r,
                filename = paste0(path, '/ppt_', year, '_', moth, '_', k, '.tif'),
                overwrite = TRUE)
    
  })
 
  print('Done!')
  
}

# Load data ---------------------------------------------------------------
root <- '//dapadfs/workspace_cluster_13/GATES'
shp <- shapefile(paste0(root, '/data/shp/base/africa_continent.shp'))
chr <- list.files(paste0(root, '/workspace/dryDays/tif/chirps/continents'), full.names = TRUE, pattern = '.tif$') %>% 
  grep('Africa', ., value = TRUE)
prc <- list.files('../tif/difference_wc/porc', full.names = T, pattern = '.tif$')
rcp <- c('rcp45', 'rcp85')
prd <- c('2020_2049', '2040_2069')
yrs <- 1981:2019

# Dates data frame
dates <- data.frame(
  period = c(rep(c('2020_2049'), length(1981:2010)), rep(c('2040_2069'), length(1981:2010))),
  year_1 = c(1981:2010, 1981:2010),
  year_2 = c(2020:2049, 2040:2069)
)

# Adding the change -------------------------------------------------------

for(j in 5:12){
  
  add_change(year = 1980, rcpe = 'rcp45', prdo = '2020_2049', moth = j)
  
}


cl <- makeCluster(4)
registerDoSNOW(cl)

foreach(j = 3:12, .verbose = TRUE) %dopar% {
  require(pacman)
  pacman::p_load(raster, rgdal, rgeos, stringr, sf, tidyverse, gtools)
  add_change(year = 1980, rcpe = 'rcp45', prdo = '2020_2049', moth = j)
}

stopCluster(cl)


