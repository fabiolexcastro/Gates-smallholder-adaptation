
# Load libraries ----------------------------------------------------------
require(pacman)
pacman::p_load(raster, rgdal, glue, rgeos, terra, sf, tidyverse, gtools)
rm(list = ls())

#

make_ensemble <- function(rc, yr){
  
  # Proof
  rc <- rcp[2]
  yr <- yrs[2]
  
  print(paste0('To start ---- > ', rc, ' ', yr))
  
  rst <- lapply(1:length(gcm), function(k){
    print(gcm[k])
    stk <- glue('{pth_ftr}/{rc}/global_2_5min/{gcm[k]}/r1i1p1/{yr}') %>% 
      list.files(., full.names = TRUE) %>% 
      grep(paste0('prec_', 1:12, collapse = '|'), ., value = TRUE) %>% 
      mixedsort() %>% 
      stack(.)
  })
  
  avg <- lapply(1:12, function(k){
    print(k)
    rsl <- mean(rst[[1]][[k]], rst[[2]][[k]], rst[[3]][[k]], rst[[4]][[k]])
  })
  avg <- stack(avg)
  names(avg) <- paste0('prec_', 1:12)
  Map('writeRaster', x = unstack(avg), filename = paste0('../data/rst/climate/ftr/ensemble/prec_', 1:12, '_', rc, '_', yr, '.tif'), overwrite = TRUE)
  print('Done...!')
  
  
}

# Load data ---------------------------------------------------------------
pth_crn <- '//DAPADFS/data_cluster_4/observed/gridded_products/worldclim/Global_2_5min_v2'
pth_ftr <- '//DAPADFS/data_cluster_2/gcm/cmip5/downscaled'
gcm <- c('MOHC_HADGEM2_ES', 'CESM1_CAM5', 'GFDL_CM3', 'MPI_ESM_LR', 'MIROC_MIROC5')

# Precipitation current
ppt_crn <- list.files(pth_crn, full.names = TRUE, pattern = '.tif$') %>% 
  mixedsort() %>% 
  grep('prec_', ., value = TRUE) %>% 
  rast()

# Precipitation future
rcp <- c('rcp45', 'rcp85')
yrs <- c('2020_2049', '2040_2069')


