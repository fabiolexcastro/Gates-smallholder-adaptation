
# Load libraries -----------------------------------------------------
require(pacman)
pacman::p_load(raster, rgdal, gtools, sf, tidyverse, terra, ncdf4)
g <- gc(reset = TRUE)
rm(list= ls())

# Functions to use ---------------------------------------------------
extract_julian <- function(x){
  x <- rasterToPoints(x) %>% as_tibble() %>% setNames(c('x', 'y', 'julian')) %>% inner_join(., jln, by = c('julian', 'julian'))
  x <- rasterFromXYZ(x[,c(1, 2, 4)])
  return(x)
  print('Done')
}
red_ftr <- function(x){mixedsort(x) %>% grep('tmean', ., value = TRUE) %>% stack()}
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
my_function <- function(cln, map, rcp, yrs, ftr_rst, thr){ 
  
  # # Proof
  # cln <- 'Cassava'
  # map <- 'CASS'
  # rcp <- 'rcp45'
  # yrs <- '2030s'
  # ftr_rst <- r45_30s
  # thr <- 32
  
  print('To start... ---> Extracting the areas for the crop')
  map_rst <- grep(map, map_spm, value = TRUE) %>% raster() 
  map_rst[which(map_rst[] <= 0)] <- NA
  map_pol <- rasterToPolygons(map_rst, dissolve = TRUE)
  map_pol@data$gid <- 1
  map_pol <- aggregate(map_pol, 'gid')
  
  print('Extracting the climate zones')
  crn_rst <- bse_crn %>% stack %>% raster::crop(., map_pol) %>% raster::mask(., map_pol)
  ftr_rst <- ftr_rst %>% stack %>% raster::crop(., map_pol) %>% raster::mask(., map_pol)
  ftr_rst <- ftr_rst / 10
  names(crn_rst) <- c(paste0('crn_tavg_', 1:12))
  names(ftr_rst) <- c(paste0('ftr_tavg_', 1:12))
  crn_ftr <- stack(crn_rst, ftr_rst)
  
  print('Extracting the calendar zones')
  cln_rst <- grep(cln, crp_cln, value = TRUE) 
  cln_str <- raster::brick(cln_rst, varname = 'plant')[[1]] * 1
  cln_str <- raster::crop(cln_str, map_pol) %>% raster::mask(., map_pol) %>% round() %>% extract_julian(.)
  cln_end <- raster::brick(cln_rst, varname = 'harvest')[[1]] * 1
  cln_end <- raster::crop(cln_end, map_pol) %>% raster::mask(., map_pol) %>% round() %>% extract_julian(.)
  
  print('Calendar zones')
  cln_stk <- stack(cln_str, cln_end)
  names(cln_stk) <- c('start', 'end')
  cln_stk <- raster::crop(cln_stk, map_pol) %>% raster::mask(., map_pol) 
  cln_tbl <- cln_stk %>% rasterToPoints(.) %>% as_tibble %>% mutate(id = 1:nrow(.))
  
  print('To extract the temperature values')
  tbl <- raster::extract(crn_ftr, cln_tbl[,1:2])
  tbl <- cbind(cln_tbl, tbl) %>% 
    as_tibble %>% 
    gather(var, value, -id, -x, -y, -start, -end) %>% 
    mutate(month = parse_number(var)) %>% 
    drop_na() %>% 
    mutate(period = str_sub(var, 1, 3)) %>% 
    dplyr::select(-var)
  
  print('To make the summary')
  rsl <- tbl %>% 
    nest(-id, -x, -y,-start, -end, -period) %>% 
    mutate(output = pmap(list(data, start, end), .f = myFunction)) %>% 
    dplyr::select(-data)
  rsl <- rsl %>% 
    unnest(output) %>% 
    spread(period, value) 
  rsl <- rsl %>% 
    group_by(id, x, y) %>% 
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
  rst <- rasterFromXYZ(rsl[,c(2, 3, 7)])
  
  print('To write the raster final')
  writeRaster(rst,
              paste0('../raster/indicators/heat_stress_wc/heat_stress_thr_', cln, '_', rcp, '_', yrs, '.tif'),
              overwrite = TRUE)
  print('Done!')
  
}

# Prepare data -------------------------------------------------------
bse_crn <- list.files('../input/worldclim/current', full.names = TRUE) %>% mixedsort() %>% stack()
bse_crn <- rast(bse_crn)

# Future data
dir_ftr <- list.files('../input/worldclim/future', full.names = TRUE)
fls_ftr <- map(.x = dir_ftr, function(x) list.files(x, full.names = TRUE, pattern = '.tif$'))
r45_30s <- red_ftr(fls_ftr[[1]]); r45_50s <- red_ftr(fls_ftr[[2]])
r85_30s <- red_ftr(fls_ftr[[3]]); r85_50s <- red_ftr(fls_ftr[[4]])        

# Crop calendar data
crp_cln <- list.files('../input/crop_calendar', full.names = TRUE, pattern = '.nc$')
crp_lst <- basename(crp_cln)
crp_lst <- str_sub(crp_lst, start = 1, end = nchar(crp_lst) - 22)

# MapSPAM data
map_spm <- list.files('../input/mapspam', full.names = TRUE, pattern = '.tif$') 
map_lst <- basename(map_spm) %>% str_split(., pattern = '_') %>% map(.x = ., .f = function(x) x[[4]]) %>% unlist()
map_abb <- read.csv('../tbl/abb_mapspam.csv')[,1:3]

# Label of the raster
lbl <- data.frame(value = 1:3, label = c('Temperature flips', 'Safe levels', 'Both above threshold'))
lbl <- lbl %>% mutate(label = factor(label, levels = c('Temperature flips', 'Safe levels', 'Both above threshold')))
jln <- read_csv('../tbl/julian.csv')


# Apply the function ------------------------------------------------------
rcps <- c('rcp45', 'rcp45', 'rcp85', 'rcp85')
year <- c('2030s', '2050s', '2030s', '2040s')
clim <- list(r45_30s, r45_50s, r85_30s, r85_50s)

for(i in 1:4){
  my_function(cln = 'Rice.2', map = 'RICE', rcp = rcps[i], yrs = year[i], ftr_rst = clim[i], thr = 30)    
}

for(i in 1:4){
  my_function(cln = 'Maize', map = 'MAIZ', rcp = rcps[i], yrs = year[i], ftr_rst = clim[i], thr = 28.6)    
}
