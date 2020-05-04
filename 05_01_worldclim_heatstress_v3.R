
# Load libraries ----------------------------------------------------------
require(pacman)
pacman::p_load(raster, rgdal, rgeos, stringr, gtools, tidyverse, lubridate, oce, terra, epitools, forcats)
rm(list = ls())
source('05_01_functions.R')

# Climate -----------------------------------------------------------------

# Some main files -----------
cnt <- shapefile('../data/shp/base/continents_1.shp')
cnts <- c('Africa', 'Asia', 'North America', 'Oceania', 'South America', 'Europe')
jln <- read_csv('../tbl/julian.csv')
tbl_thr <- read_csv('../tbl/thresholds_temperature.csv')

# Current -------------------
wcl_crn <- list.files('//dapadfs/data_cluster_4/observed/gridded_products/worldclim/Global_2_5min_v2', full.names = T) %>% 
  grep('tavg', ., value = T) %>% 
  mixedsort() %>% 
  stack()
wcl_crn <- rast(wcl_crn)

# Future -------------------
wcl_ftr <- list(read_future(x = 'rcp45', y = '2020_2049'), 
                read_future(x = 'rcp45', y = '2040_2069'), 
                read_future(x = 'rcp85', y = '2020_2049'), 
                read_future(x = 'rcp85', y = '2040_2069'))

# Crop calendar ------------
crp_cln <- list.files('../netCDF/crop_calendar/crops/filled', full.names = T, pattern = '.nc$')
crp_lst <- str_sub(crp_cln, 38, nchar(crp_cln) - 22)

# MapSPAM ------------------
fls_map <- list.files('../raster/mapSPAM', full.names = T, pattern = '.tif')
abb_map <- read_csv('../tbl/abb_mapspam.csv')
names(abb_map) <- c('crop', 'full_name', 'name')
crops <- toupper(abb_map$name)
# map(.x = crops[32:length(crops)], .f = zones_crop)

fls_map <- list.files('../data/shp/map_spam', full.names = T, pattern = '.shp$')
crops <- sort(crops)

lbl <- data.frame(
  value = 1:3, 
  label = c('Temperature flips', 'Safe levels', 'Both above threshold')
)  %>% 
  mutate(label = factor(label, levels = c('Temperature flips', 'Safe levels', 'Both above threshold')))


# Extracting the months from crop calendar --------------------------------
  # Proof
  crop <- 'BARL'
  crn <- wcl_crn
  ftr <- rast(wcl_ftr[[2]])
  thr <- tbl_thr %>% filter(Crop == 'Barley') %>% pull(2)
  cln <- 'Barley.Winter'
  rcp <- 'rcp45'
  yr <- '2050s'

extract_crop <- function(crop, crn, thr, cln, ftr, rcp, yr){
  
  print('Zones')
  zne_crp <- grep(crop, fls_map, value = T) %>% vect()
  crn_zne <- terra::crop(crn, zne_crp) %>% terra::mask(., zne_crp) %>% raster::stack()
  ftr_zne <- terra::crop(ftr, zne_crp) %>% terra::mask(., zne_crp) %>% raster::stack()
  
  names(crn_zne) <- c(paste0('crn_tavg_', 1:12))
  names(ftr_zne) <- c(paste0('ftr_tavg_', 1:12))
  
  print('Calendar')=
  cln <- grep(cln, crp_cln, value =T) 
  cln_str <- extract_julian(round(raster::brick(cln, varname = 'plant')[[1]] * 1))
  cln_end <- extract_julian(round(raster::brick(cln, varname = 'harvest')[[1]] * 1))
  
  print('Calendar stack')
  cln_stk <- rast(stack(cln_str, cln_end))
  names(cln_stk) <- c('start', 'end')
  cln_stk <- terra::crop(cln_stk, zne_crp) %>% terra::mask(., zne_crp)
  cln_stk <- stack(cln_stk)
  cln_tbl <- cln_stk %>% rasterToPoints(.) %>% as_tibble() %>% mutate(id = 1:nrow(.))
  
  print('To extract the temperature values')
  crn_ftr <- stack(crn_zne, ftr_zne)
  tbl <- raster::extract(crn_ftr, cln_tbl[,1:2])
  tb2 <- cbind(cln_tbl, tbl) %>% 
    as_tibble(.) %>% 
    gather(var, value, -id, -x, -y, -start, -end) %>% 
    mutate(month = parse_number(var))
  tb2 <- drop_na(tb2)
  tb2 <- tb2 %>% mutate(period = str_sub(var, 1, 3)) %>%  dplyr::select(-var)
  
  print('To make the summary')
  rsl <- tb2 %>% 
    nest(-id, -x, -y, -start, -end, -period) %>% 
    mutate(output = pmap(list(data, start, end), .f = myFunction)) %>% 
    dplyr::select(-data) %>% 
    unnest(output) %>% 
    spread(period, value) %>% 
    drop_na() %>% 
    group_by(id, x, y) %>% 
    dplyr::summarise(crn = mean(crn),
                     ftr = mean(ftr)) %>% 
    ungroup()
  rsl <- rsl %>% 
    mutate(class =   
             case_when(
               crn < thr & ftr >= thr ~ 'Temperature flips',
               crn <= thr & ftr <= thr ~ 'Safe levels',
               crn >= thr & ftr < thr ~ 'Safe levels',
               crn >= thr & ftr >= thr ~ 'Both above threshold',
               TRUE ~ "others"
             )
    )
  
  rsl <- inner_join(rsl, lbl, by = c('class' = 'label'))
  rst <- rasterFromXYZ(rsl[,c(2, 3, 7)])
  writeRaster(rst, paste0('../raster/indicators/heat_stress_wc/heat_stress_thr_', 'Barley.Winter', '_', rcp, '_', yr, '.tif'), overwrite = T)
  print('Done!!!')
}


# Apply the function

climate <- list(r45_30s, r45_50s, r85_30s, r85_50s)
climate_name <- c('rcp45', 'rcp45', 'rcp85', 'rcp85')
years <- c('2030s', '2050s', '2030s', '2050s')

for(i in 1:4){
  extract_crop(
    crop = 'BARL',
    crn = wcl_crn,
    ftr = climate[[i]],
    thr = tbl_thr %>% filter(Crop == 'Barley') %>% pull(2),
    cln = 'Barley.Winter',
    rcp = climate_name[i],
    yr = years[i]
  )
  
}

