

# Load libraries ----------------------------------------------------------
require(pacman)
pacman::p_load(raster, rgdal, rgeos, stringr, sf, tidyverse)
rm(list = ls())

# Functions to use ---------------------------------------------------------
make_mosaic <- function(yr){
  
  print(yr)
  fle <- grep(yr, fls, value = T)
  rst <- map(.x = fle, .f = raster)
  
  ms1 <- raster::mosaic(rst[[1]], rst[[2]], fun = 'mean')
  ms2 <- raster::mosaic(rst[[3]], rst[[4]], fun = 'mean')
  ms3 <- raster::mosaic(ms1, rst[[5]], fun = 'mean')
  ms4 <- raster::mosaic(ms2, ms3, fun = 'mean')
  ms5 <- raster::mosaic(ms3, ms4, fun = 'mean')
  # ms6 <- raster::mosaic(ms4, ms5, fun = 'mean')
  
  # ms1 <- raster::mosaic(rst[[1]], rst[[2]], fun = 'mean')
  # ms2 <- raster::mosaic(rst[[3]], rst[[4]], fun = 'mean')
  # ms3 <- raster::mosaic(rst[[5]], rst[[6]], fun = 'mean')
  # ms4 <- raster::mosaic(ms1, ms2, fun = 'mean')
  # ms5 <- raster::mosaic(ms3, ms4, fun = 'mean')
  # ms6 <- raster::mosaic(ms4, ms5, fun = 'mean')
  
  # rm(ms1, ms2, ms3, ms4, ms5)
  
  # writeRaster(ms6, paste0('../raster/indicators/heat_stress/mosaic/heat_barley_s2_', yr, '.tif'), overwrite = T)
  writeRaster(ms5, paste0('../raster/indicators/heat_stress/mosaic/heat_cowpea_', yr, '.tif'), overwrite = T)
  print('Done!!!')
  
}
make_raster <- function(crp){
  
  crp <- 'cowpea' 
  stk <- grep(crp, fls, value = T)
  stk <- grep(paste0(1983:2016, collapse = '|'), stk, value = T)
  stk <- stk[-grep('_s2_', stk, value = F)]
  stk <- stack(stk)
  avg <- mean(stk)

  cfv <- calc(stk, sd) / avg
  cfv <- cfv * 100
  
  q95fun <- function(x){quantile(x, probs = .95, na.rm=TRUE)}
  p95 <- calc(stk, fun=q95fun, forceapply=T)
  
  print('To write the raster')
  writeRaster(avg,  paste0('../raster/indicators/heat_stress/mosaic/heat_crop_', crp, '_mean.tif'), overwrite = T)
  writeRaster(cfv,  paste0('../raster/indicators/heat_stress/mosaic/heat_crop_', crp, '_cv.tif'), overwrite = T)
  writeRaster(p95, paste0('../raster/indicators/heat_stress/mosaic/heat_crop_', crp, '_p95.tif'), overwrite = T)
  
  print('To calculate the raster stack')
  stk <- stack(avg, cfv, p95)
  names(stk) <- c('mean', 'CV', 'percentil95')
  vls <- rasterToPoints(stk) %>% 
    as_tibble() %>% 
    mutate(id = 1:nrow(.)) %>% 
    gather(var, value, -id, -x, -y)
  
  g1 <- ggplot(data = filter(vls, var == 'mean')) +
    geom_tile(aes(x = x, y =  y, fill = value)) +
    scale_fill_gradientn(colours = RColorBrewer::brewer.pal(n = 8, name = "YlOrRd"), 
                         na.value = 'white') +
    theme_bw() +
    geom_polygon(data = shp, aes(x = long, y = lat, group = group), color = 'grey', fill = 'NA') +
    coord_equal() +
    labs(title = '', fill = 'Days',  x = 'Longitude', y = 'Latitude') +
    theme(legend.position = 'bottom',
          plot.title = element_text(hjust = 0.5),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.key.width = unit(5, 'line')) +
    guides(shape = guide_legend(override.aes = list(size = 10)))
  
  cfv_tbl <- rasterToPoints(cfv) %>% 
    as_tibble() %>% 
    mutate(id = 1:nrow(.)) %>% 
    gather(var, value, -id, -x, -y)
  cfv_tbl <- cfv_tbl %>% mutate(value = ifelse(value > 100, 101, value))
  
  g2 <- ggplot(data = cfv_tbl) +
    geom_tile(aes(x = x, y =  y, fill = value)) +
    scale_fill_gradientn(colours = RColorBrewer::brewer.pal(n = 8, name = "YlOrRd"), 
                         na.value = 'white',
                         labels = c(0, 25, 50, 75, '>100')) +
    theme_bw() +
    geom_polygon(data = shp, aes(x = long, y = lat, group = group), color = 'grey', fill = 'NA') +
    coord_equal() +
    labs(title = '', fill = 'CV',  x = 'Longitude', y = 'Latitude') +
    theme(legend.position = 'bottom',
          plot.title = element_text(hjust = 0.5),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.key.width = unit(5, 'line')) +
    guides(shape = guide_legend(override.aes = list(size = 10)))
  
  g3 <- ggplot(data = filter(vls, var == 'percentil95')) +
    geom_tile(aes(x = x, y =  y, fill = value)) +
    scale_fill_gradientn(colours = RColorBrewer::brewer.pal(n = 8, name = "YlOrRd"), 
                         na.value = 'white') +
    theme_bw() +
    geom_polygon(data = shp, aes(x = long, y = lat, group = group), color = 'grey', fill = 'NA') +
    coord_equal() +
    labs(title = '', fill = 'CV',  x = 'Longitude', y = 'Latitude') +
    theme(legend.position = 'bottom',
          plot.title = element_text(hjust = 0.5),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.key.width = unit(5, 'line')) +
    guides(shape = guide_legend(override.aes = list(size = 10)))
  
  library(ggpubr)
  gg <- ggarrange(g1, g2, g3, ncol = 1, nrow = 3, labels = c('Mean', 'CV', 'Percentil95'))
  ggsave(plot = gg, filename = paste0('../png/maps/mean_cv_p95_', crp, '.png'), units = 'in', width = 12, height = 25, dpi = 300)
  print('DOne!')
  
}

# Load data ---------------------------------------------------------------
fls <- list.files('../raster/indicators/heat_stress', full.names = T, pattern = '.tif$')
fls <- grep('cowpea', fls, value = T)
fls <- fls[grep('s2.tif', fls, value = F)]
yrs <- 1983:2016
shp <- shapefile('../data/shp/base/continents_1.shp')

# Apply the function ------------------------------------------------------
map(.x = yrs, .f = make_mosaic)

# Mean --------------------------------------------------------------------
fls <- list.files('../raster/indicators/heat_stress/mosaic', full.names = T, pattern = '.tif$')

make_raster(crp = 'cowpea')

