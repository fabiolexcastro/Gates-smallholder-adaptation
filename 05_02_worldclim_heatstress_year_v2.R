
# Load libraries -----------------------------------------------------
require(pacman)
pacman::p_load(raster, rgdal, gtools, sf, tidyverse, terra, ncdf4, magrittr, foreach, parallel, doSNOW)
g <- gc(reset = TRUE)
rm(list= ls())


# Functions to use --------------------------------------------------------

my_function <- function(climate, crp, map, rcp, yrs, thr){
  
  # climate <- crn_r45_30s
  # crp <- 'banana'
  # map <- 'BANA'
  # thr <- 27
  # rcp <- 'rcp45'
  # yrs <- '2030s'
  
  print('To start... ---> Extracting the areas for the crop')
  map_rst <- grep(map, map_spm, value = TRUE) %>% raster() 
  map_rst <- raster::resample(x = map_rst, y = msk, method = 'bilinear')
  map_rst <- raster::crop(map_rst, shp) %>% raster::mask(., shp)
  map_tbl <- as_tibble(as.data.frame(rasterToPoints(map_rst))) %>% setNames(c('x', 'y', 'value'))
  
  cell <- cellFromXY(map_rst, xy = as.data.frame(map_tbl[,1:2]))
  map_tbl <- as_tibble(cbind(cell = cell, map_tbl))
  map_tbl <- map_tbl %>% mutate(value = ifelse(value <= 0, 0, 1))
  
  cell_sub <- filter(map_tbl, value == 1) %>% pull(1)
  clim_sub <- climate %>% filter(cll %in% cell_sub)
  
  print('Climate sub')
  clim_sub %<>% 
    gather(var, value, -cll, -x, -y) %>% 
    mutate(month = parse_number(var)) %>% 
    drop_na %>% 
    mutate(period = str_sub(var, 1, 3)) %>% 
    dplyr::select(-var) %>% 
    group_by(cll, x, y, period) %>% 
    dplyr::summarise(value = mean(value)) %>% 
    ungroup() %>% 
    spread(period, value)%>% 
    mutate(class = 
             case_when(
               crn < thr & ftr >= thr ~ 'Temperature flips',
               crn <= thr & ftr <= thr ~ 'Safe levels',
               crn >= thr & ftr < thr ~ 'Safe levels',
               crn >= thr & ftr >= thr ~ 'Both above threshold',
               TRUE ~ "others"
             )
    ) %>% 
    inner_join(., lbl, by = c('class' = 'label'))
  
  rst <- msk
  rst[pull(clim_sub, cll)] <- pull(clim_sub, value)
  
  print('To write the raster final')
  writeRaster(rst,
              paste0('../raster/indicators/heat_stress_wc/heat_stress_thr_', crp, '_', rcp, '_', yrs, '.tif'),
              overwrite = TRUE)
  print('Done!')
  
 }

# Load data ---------------------------------------------------------------

# Base files
shp <- shapefile('../data/shp/base/continents.shp')
shp <- shp[shp@data$CONTINENT == 'Africa',]
msk <- raster::raster('../input/worldclim/future/rcp45_30s/tmean_1.tif') * 0
msk <- raster::crop(msk, shp) %>% raster::mask(., shp)

# Load climate
crn_r45_30s <- readRDS(file = '../rds/worldclim/current_future_africar45_30s.rds')
crn_r45_50s <- readRDS(file = '../rds/worldclim/current_future_africar45_50s.rds')
crn_r85_30s <- readRDS(file = '../rds/worldclim/current_future_africar85_30s.rds')
crn_r85_50s <- readRDS(file = '../rds/worldclim/current_future_africar85_50s.rds')

# MapSPAM data
map_spm <- list.files('../raster/mapSPAM', full.names = TRUE, pattern = '.tif$') 
map_spm <- grep('_H_', map_spm, value = TRUE) %>% grep('_A.tif', ., value = TRUE)
map_lst <- basename(map_spm) %>% str_split(., pattern = '_') %>% map(.x = ., .f = function(x) x[[4]]) %>% unlist()
map_abb <- read.csv('../tbl/abb_mapspam.csv')[,1:3]
map_abb <- map_abb %>% mutate(name = toupper(name))

# Label of the raster
lbl <- data.frame(value = 1:3, label = c('Temperature flips', 'Safe levels', 'Both above threshold'))
lbl <- lbl %>% mutate(label = factor(label, levels = c('Temperature flips', 'Safe levels', 'Both above threshold')))
jln <- read_csv('../tbl/julian.csv')

# Apply the function ------------------------------------------------------
climate_list <- list(crn_r45_30s, crn_r45_50s, crn_r85_30s, crn_r85_50s)
rcps <- c('rcp45', 'rcp45', 'rcp85', 'rcp85')
years <- c('2030s', '2050s', '2030s', '2050s')

# climate <- crn_r45_30s
# crp <- 'banana'
# map <- 'BANA'
# thr <- 27
# rcp <- 'rcp45'
# yrs <- '2030s'


cl <- makeCluster(4)
registerDoSNOW(cl)
foreach(i = 1:4, .packages = c('raster', 'rgdal', 'rgeos', 'stringr', 'tidyverse'), .verbose = TRUE) %dopar% {
  print(i)
  my_function(climate = climate_list[[i]], 
              crp = 'banana',
              map = 'BANA',
              rcp = rcps[i],
              yrs = years[i],
              thr = 27)
}
stopCluster(cl)

map_abb

cl <- makeCluster(4)
registerDoSNOW(cl)
foreach(i = 1:4, .packages = c('raster', 'rgdal', 'rgeos', 'stringr', 'tidyverse'), .verbose = TRUE) %dopar% {
  print(i)
  my_function(climate = climate_list[[i]], 
              crp = 'cocoa',
              map = 'COCO',
              rcp = rcps[i],
              yrs = years[i],
              thr = 32)
}
stopCluster(cl)

map_abb

cl <- makeCluster(4)
registerDoSNOW(cl)
foreach(i = 1:4, .packages = c('raster', 'rgdal', 'rgeos', 'stringr', 'tidyverse'), .verbose = TRUE) %dopar% {
  print(i)
  my_function(climate = climate_list[[i]], 
              crp = 'coconut',
              map = 'CNUT',
              rcp = rcps[i],
              yrs = years[i],
              thr = 34)
}
stopCluster(cl)


cl <- makeCluster(4)
registerDoSNOW(cl)
foreach(i = 1:4, .packages = c('raster', 'rgdal', 'rgeos', 'stringr', 'tidyverse'), .verbose = TRUE) %dopar% {
  print(i)
  my_function(climate = climate_list[[i]], 
              crp = 'robusta',
              map = 'RCOF',
              rcp = rcps[i],
              yrs = years[i],
              thr = 30)
}
stopCluster(cl)

cl <- makeCluster(4)
registerDoSNOW(cl)
foreach(i = 1:4, .packages = c('raster', 'rgdal', 'rgeos', 'stringr', 'tidyverse'), .verbose = TRUE) %dopar% {
  print(i)
  my_function(climate = climate_list[[i]], 
              crp = 'arabica',
              map = 'ACOF',
              rcp = rcps[i],
              yrs = years[i],
              thr = 25)
}
stopCluster(cl)

cl <- makeCluster(4)
registerDoSNOW(cl)
foreach(i = 1:4, .packages = c('raster', 'rgdal', 'rgeos', 'stringr', 'tidyverse'), .verbose = TRUE) %dopar% {
  print(i)
  my_function(climate = climate_list[[i]], 
              crp = 'oilpalm',
              map = 'OILP',
              rcp = rcps[i],
              yrs = years[i],
              thr = 32)
}
stopCluster(cl)

cl <- makeCluster(4)
registerDoSNOW(cl)
foreach(i = 1:4, .packages = c('raster', 'rgdal', 'rgeos', 'stringr', 'tidyverse'), .verbose = TRUE) %dopar% {
  print(i)
  my_function(climate = climate_list[[i]], 
              crp = 'plantain',
              map = 'PLNT',
              rcp = rcps[i],
              yrs = years[i],
              thr = 27)
}
stopCluster(cl)

cl <- makeCluster(4)
registerDoSNOW(cl)
foreach(i = 1:4, .packages = c('raster', 'rgdal', 'rgeos', 'stringr', 'tidyverse'), .verbose = TRUE) %dopar% {
  print(i)
  my_function(climate = climate_list[[i]], 
              crp = 'sugarcane',
              map = 'SUGC',
              rcp = rcps[i],
              yrs = years[i],
              thr = 35)
}
stopCluster(cl)

cl <- makeCluster(4)
registerDoSNOW(cl)
foreach(i = 1:4, .packages = c('raster', 'rgdal', 'rgeos', 'stringr', 'tidyverse'), .verbose = TRUE) %dopar% {
  print(i)
  my_function(climate = climate_list[[i]], 
              crp = 'tea',
              map = 'TEAS',
              rcp = rcps[i],
              yrs = years[i],
              thr = 30)
}
stopCluster(cl)

