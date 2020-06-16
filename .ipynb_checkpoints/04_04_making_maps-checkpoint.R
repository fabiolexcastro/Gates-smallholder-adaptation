

# Load libraries ----------------------------------------------------------
require(pacman)
pacman::p_load(raster, rgdal, rgeos, stringr, velox, sf, tidyverse, 
               maptools)
rm(list = ls())

# Load data ---------------------------------------------------------------
pth <- '../raster/indicators/heat_stress/mosaic'
fls <- list.files(pth, full.names = T, pattern = '.tif$')
shp <- shapefile('../data/shp/base/continents_1.shp')

make_map <- function(crp){
  crp <- 'cowpea'
  fle <- grep(crp, fls, value = T)
  fle <- grep(paste0(1983:2016, collapse = '|'), fle, value = T)
  fle <- fle[-grep('_s2_', fle, value = F)]
  stk <- stack(fle)

  vls <- rasterToPoints(stk) %>% as_tibble() 
  vls <- vls %>% gather(var, value, -x, -y)
  vls <- vls %>% mutate(year = parse_number(var))
  # vls <- vls %>% mutate(year = str_sub(var, 13, nchar(var)))
  vls <- vls %>% mutate(year = factor(year, levels = 1983:2016))
  
}

prd <- c(1983:1994, 1995:2006, 2007:2016)
prd[1]

mapping <- function(x){
  x <- vls %>% filter(year %in% 1995:2006)
  g <- ggplot(x)  +
    geom_tile(aes(x = x, y =  y, fill = value)) +
    facet_wrap(~ year) +
    scale_fill_gradientn(colours = RColorBrewer::brewer.pal(n = 8, name = "YlOrRd"), 
                         na.value = 'white', limits = c(1, 60)) +
    theme_bw() +
    geom_polygon(data = shp, aes(x = long, y = lat, group = group), color = 'grey', fill = 'NA') +
    # scale_x_continuous(breaks = c(-78, -77, -76, -75)) +
    coord_equal() +
    # coord_equal(xlim = extent(vll)[1:2], ylim = extent(vll)[3:4]) +
    labs(title = 'Heat stress Cowpea', fill = 'days',  x = 'Longitude', y = 'Latitude') +
    theme(legend.position = 'bottom',
          plot.title = element_text(hjust = 0.5),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.key.width = unit(5, 'line')) +
    guides(shape = guide_legend(override.aes = list(size = 10)))
  ggsave(plot = g, filename = paste0('../png/maps/heat_stress_cowpea_95_06.png'), units = 'in', width = 15, height = 10, dpi = 300)
  print('Done!!!')
}