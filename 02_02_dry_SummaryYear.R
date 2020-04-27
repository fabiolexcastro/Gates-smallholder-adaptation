
print('To load the libraries')
require(raster)
require(rgdal)
require(rgeos)
require(stringr)
require(velox)
require(sf)
require(tidyverse)
require(glue)
require(gtools)
rm(list = ls())


print('Functions to use')
mySum <- function(yr){
  print(yr)
  st <- grep(yr, fls, value = T) %>% 
    mixedsort %>% 
    stack()
  df <- st %>%
    rasterToPoints() %>%
    as_tibble()
  d2 <- df %>%
    mutate(id = 1:nrow(.)) %>% 
    gather(var, value, -x, -y, -id) %>% 
    group_by(id, x, y) %>% 
    dplyr::summarise(value = sum(value)) %>% 
    ungroup()
  rs <- rasterFromXYZ(d2[,2:4])
  writeRaster(rs, paste0('../raster/indicators/dryDays/yearly/dry_days_', yr, '.asc'), overwrite = T)
  print('Done')
} 


print('Load data')
fls <- list.files('../raster/indicators/dryDays/monthly', full.names = T, pattern = '.asc$')
map(.x = 1981:2019, .f = mySum)
