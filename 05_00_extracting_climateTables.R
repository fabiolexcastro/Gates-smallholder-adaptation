

# Load libraries ------------------------------------------------------
require(pacman)
pacman::p_load(raster, rgdal, gtools, sf, tidyverse, terra, ncdf4)
g <- gc(reset = TRUE)
rm(list= ls())
source('05_01_functions_v5.R')
options(scipen = 999)

# Main function -----------------------------------------------------------
create_table <- function(crn, ftr, nme){
  # crn <- bse_crn
  # ftr <- r45_30s
  # nme <- 'r45_30s'
  
  print('To start...!')
  stk <- addLayer(crn, ftr)
  vls <- rasterToPoints(stk, spatial = FALSE)
  print('Cell from XY')
  cll <- cellFromXY(stk[[1]], xy = vls[,1:2])
  vls <- cbind(cll, vls)
  vls <- as.data.frame(vls)
  vl2 <- drop_na(vls)
  vl2 <- as_tibble(vl2)
  saveRDS(object = vl2, file = paste0('../rds/worldclim/current_future_', nme, '.rds'))
  print('Done!')

}

# Prepare data --------------------------------------------------------

# Current
bse_crn <- list.files('../input/worldclim/current', full.names = TRUE) %>% mixedsort() %>% stack() %>% setNames(c(paste0('crn_tavg_', 1:12)))
shp <- shapefile('../data/shp/base/continents.shp')

# Future data
fls_ftr <- list.files('../input/worldclim/future', full.names = TRUE) %>% 
  map(.x = ., function(x) list.files(x, full.names = TRUE, pattern = '.tif$'))
r45_30s <- red_ftr(fls_ftr[[1]]); r45_50s <- red_ftr(fls_ftr[[2]])
r85_30s <- red_ftr(fls_ftr[[3]]); r85_50s <- red_ftr(fls_ftr[[4]])      
bse_crn <- raster::crop(bse_crn, r45_30s[[1]])

create_table(crn = bse_crn, ftr = r45_50s, nme = 'r45_50s')
create_table(crn = bse_crn, ftr = r85_30s, nme = 'r85_30s')
create_table(crn = bse_crn, ftr = r85_50s, nme = 'r85_50s')



