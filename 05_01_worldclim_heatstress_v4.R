
# Load libraries -----------------------------------------------------
require(pacman)
pacman::p_load(raster, rgdal, gtools, sf, tidyverse, magrittr, terra, ncdf4)
g <- gc(reset = TRUE)
rm(list= ls())

# Functions to use ---------------------------------------------------
extract_julian <- function(x){
  x <- rasterToPoints(x) %>% as_tibble() %>% setNames(c('x', 'y', 'julian')) %>% inner_join(., jln, by = c('julian', 'julian'))
  x <- rasterFromXYZ(x[,c(1, 2, 4)])
  return(x)
  print('Done')
}
red_ftr <- function(x){
  rsl <- mixedsort(x) %>% grep('tmean', ., value = TRUE) %>% stack()
  rsl <- rsl / 10
}
myFunction <- function(x, st, ed){
  
  if(st > ed){
    y <- c(st:12, 1:ed)
    x %>% 
      filter(month %in% y) 
  } else {
    x %>% 
      filter(month %in% st:ed)
  }
  
} 
my_function <- function(climate, cln, map, rcp, yrs, thr){ 
  
  # Proof
  climate <- crn_r85_50s
  cln <- 'Maize'
  map <- 'MAIZ'
  rcp <- 'rcp85'
  yrs <- '2050s'
  thr <- 28.6
  
  print('To start... ---> Extracting the areas for the crop')
  map_rst <- grep(map, map_spm, value = TRUE) %>% raster() 
  map_rst <- raster::resample(x = map_rst, y = msk, method = 'bilinear')
  map_tbl <- as_tibble(as.data.frame(rasterToPoints(map_rst))) %>% setNames(c('x', 'y', 'value'))
  
  cell <- cellFromXY(map_rst, xy = as.data.frame(map_tbl[,1:2]))
  map_tbl <- as_tibble(cbind(cell = cell, map_tbl))
  map_tbl <- map_tbl %>% mutate(value = ifelse(value <= 0, 0, 1))
  
  cell_sub <- filter(map_tbl, value == 1) %>% pull(1)
  clim_sub <- climate %>% filter(cll %in% cell_sub)
  
  print('Extracting the calendar zones')
  cln_rst <- grep(paste0(cln, '.crop'), crp_cln, value = TRUE) 
  cln_str <- raster::brick(cln_rst, varname = 'plant')[[1]] * 1
  cln_str <- raster::resample(x = cln_str, y = msk, method = 'bilinear')
  cln_str <- round(cln_str, 0)
  cln_end <- raster::brick(cln_rst, varname = 'harvest')[[1]] * 1
  cln_end <- raster::resample(x = cln_end, y = msk, method = 'bilinear')
  cln_end <- round(cln_end, 0)
  cln_stk <- raster::stack(cln_str, cln_end)
  names(cln_stk) <- c('start', 'end')
  
  cln_tbl <- as_tibble(as.data.frame(rasterToPoints(cln_stk))) %>% setNames(c('x', 'y', 'start', 'end'))
  cell <- cellFromXY(cln_stk, xy = as.data.frame(cln_tbl)[,1:2])
  cln_tbl <- as_tibble(cbind(cell = cell, cln_tbl))
  cln_tbl <- cln_tbl %>% filter(cell %in% cell_sub)
  cln_tbl <- inner_join(cln_tbl, jln, by = c('start' = 'julian')) %>% rename(month_start = month)
  cln_tbl <- inner_join(cln_tbl, jln, by = c('end' = 'julian')) %>% rename(month_end = month)
  cln_tbl <- cln_tbl %>% dplyr::select(cell, x, y, start = month_start, end = month_end)
  
  print('Checking the row length')
  nrow(map_tbl %>% filter(value == 1))
  nrow(clim_sub)
  nrow(cln_tbl)
  
  print('Tidy the climate table')
  clim_sub %<>% 
    gather(var, value, -cll, -x, -y) %>% 
    mutate(month = parse_number(var)) %>% 
    drop_na %>% 
    mutate(period = str_sub(var, 1, 3)) %>% 
    dplyr::select(-var)
  
  print('To make the summarise')
  tbl <- inner_join(clim_sub, cln_tbl[,c(1, 4, 5)], by = c('cll' = 'cell'))
  rsl <- tbl %>% 
    nest(-cll, -x, -y, -start, -end, -period) %>% 
    mutate(output = pmap(list(data, start, end), .f = myFunction)) %>% 
    dplyr::select(-data)
  rs2 <- rsl %>% 
    unnest(output) %>% 
    spread(period, value) 
  rs2 <- rs2 %>% 
    group_by(cll, x, y) %>% 
    dplyr::summarise(crn = mean(crn), ftr = mean(ftr)) %>% 
    ungroup() %>% 
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
  
  print('Table to raster!')
  
  rst <- msk
  rst[pull(rs2, 1)] <- pull(rs2, value)
  
  print('To write the raster final')
  writeRaster(rst,
              paste0('../raster/indicators/heat_stress_wc/heat_stress_thr_', cln, '_', rcp, '_', yrs, '.tif'),
              overwrite = TRUE)
  print('Done!')
  
}

# Prepare data -----------------------------------------------------------------------------------------------------------

# Base shape
shp <- shapefile('../data/shp/base/continents.shp')

# Base mask
msk <- raster::raster('../input/worldclim/future/rcp45_30s/tmean_1.tif') * 0

# Load climate
crn_r45_30s <- readRDS(file = '../rds/worldclim/current_future_r45_30s.rds')
crn_r45_50s <- readRDS(file = '../rds/worldclim/current_future_r45_50s.rds')
crn_r85_30s <- readRDS(file = '../rds/worldclim/current_future_r85_30s.rds')
crn_r85_50s <- readRDS(file = '../rds/worldclim/current_future_r85_50s.rds')

# Crop calendar data
crp_cln <- list.files('../input/crop_calendar', full.names = TRUE, pattern = '.nc$')
crp_lst <- basename(crp_cln)
crp_lst <- str_sub(crp_lst, start = 1, end = nchar(crp_lst) - 22)

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
for(i in 1:4){
  my_function(cln = 'Cassava', 
              map = 'CASS', 
              rcp = rcps[i], 
              yrs = year[i], 
              ftr_rst = clim[i], 
              thr = 32)    
} 
