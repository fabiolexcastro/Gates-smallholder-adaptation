


zones_crop <- function(crop){
  print(crop)
  rst <- grep(crop, fls_map, value = T) %>% 
    grep('_A.tif', ., value = T) %>% 
    raster()
  rst[which(rst[] == 0)] <- NA
  rst[which(rst[] > 0)] <- 1
  pol <- rasterToPolygons(rst, dissolve = T)
  shapefile(pol, paste0('../data/shp/map_spam/mapSPAM_', crop, '.shp'))
  print('Done!')
}

extract_julian <- function(x){
  x <- rasterToPoints(x) %>% as_tibble() %>% setNames(c('x', 'y', 'julian')) %>% inner_join(., jln, by = c('julian', 'julian'))
  x <- rasterFromXYZ(x[,c(1, 2, 4)])
  return(x)
  print('Done')
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
ensembling <- function(rc, yr){
  
  # Proof
  # rc <- 'rcp45'
  # yr <- '2020_2049'
  
  pth <- list.files(paste0(wcl_ftr, '/', rc, '/global_2_5min'), full.names = T)
  gcm <- tolower(c('MOHC_HADGEM2_ES', 'CESM1_CAM5', 'GFDL_CM3', 'MPI_ESM_LR', 'MIROC_MIROC5'))
  pth <- grep(paste0(gcm, collapse = '|'), pth, value = T)
  vrs <- paste0('tmean_', 1:12, '$')
  
  fls <- lapply(1:12, function(k){
    print(vrs[k])
    fl <- list.files(paste0(pth, '/r1i1p1/', yr), full.names = T) %>% 
      grep(vrs[k], ., value = T) %>% 
      mixedsort() 
    st <- stack(fl)
    st <- mean(st)
    writeRaster(st, paste0('../data/rst/climate/ftr/', rc, '/', yr, '/', unique(basename(fl)), '.tif'), overwrite = T)
    print('Done')
    return(st)
  })
  
  print('Done...!')
  
}
