
require(raster)
require(tidyverse)
rm(list = ls())

# Functions
myMean <- function(x){
  # x <- 'mnt1_'
  fl <- grep(x, fls, value = T) 
  st <- stack(fl)
  av <- mean(st)
  print('Done!')
  return(av)
}

# Load data
fls <- list.files('../raster/indicators/dryDays', full.names = TRUE, pattern = '.tif$') %>% 
  grep('Africa', ., value = T)
mns <- paste0('mnt', 1:12, '_')
rsl <- map(.x = mns, .f = myMean)
Map('writeRaster', x = rsl, filename =  paste0('../raster/indicators/dryDays/meanMonthly/Africa_countDryDays_', 1:12, '.tif'), overwrite = FALSE)
  
