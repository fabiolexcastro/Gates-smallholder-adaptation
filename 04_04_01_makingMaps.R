
# Load libraries ----------------------------------------------------------
require(pacman)
pacman::p_load(raster, rgdal, rgeos, stringr, sf, tidyverse, 
               maptools, foreach, parallel, doSNOW, RColorBrewer)
rm(list = ls())

# Functions to use --------------------------------------------------------
extract_table <- function(cr){
  # cr <- 'cotton'
  print(cr)
  fle <- grep(cr, fls, value = T)
  fle <- grep(paste0(1983:2016, collapse = '|'), fle, value = T)
  fle <- fle[1:34]
  stk <- stack(fle)
  vls <- rasterToPoints(stk) %>% as_tibble() 
  vls <- vls %>% gather(var, value, -x, -y)
  # vls <- vls %>% mutate(year = str_sub(var, start = nchar(var) - 3, end = nchar(var)))
  # vls <- vls %>% mutate(year = parse_number(var))
  vls <- vls %>% mutate(year = str_sub(var, 18, 21))
  vls <- vls %>% mutate(year = as.numeric(year))
  mxm <- max(vls$value)
  vls <- vls %>% mutate(max = mxm)
  vls <- vls %>% mutate(period = if_else(mxm > 60, 'Yearly', 'Season'))
  vls <- vls %>% mutate(year = factor(year, levels = 1983:2016))
  saveRDS(object = vls, file = paste0('../rds/indicators/heat_stress_daily/daily_tbl_', cr, '.rds'))
  print('Done!')
  return(vls)
}
make_map <- function(fle, yrs, nme){
  
  # Proof
  # fle <- vls_all[10]
  # yrs <- prd[1]
  # nme <- 'barley s1'
  
  print(nme)
  tbl <- readRDS(fle)
  tpe <- unique(tbl$period)
  yrs <- as.numeric(yrs[[1]])
  yrs <- as.character(yrs)
  tbl <- tbl %>% filter(year %in% yrs)
  
  y1 <- yrs[1]
  y2 <- yrs[length(yrs)]
  
  if(tpe == 'Season'){
    
    g <- ggplot(tbl)  +
      geom_tile(aes(x = x, y =  y, fill = value)) +
      facet_wrap(~ year) +
      scale_fill_gradientn(colours = RColorBrewer::brewer.pal(n = 8, name = "YlOrRd"), 
                           na.value = 'white', limits = c(0, 60)) +
      theme_bw() +
      geom_polygon(data = shp, aes(x = long, y = lat, group = group), color = 'grey', fill = 'NA') +
      coord_equal() +
      labs(title = paste0('Heat stress ', nme), fill = 'days',  x = 'Longitude', y = 'Latitude') +
      theme(legend.position = 'bottom',
            plot.title = element_text(hjust = 0.5),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            legend.key.width = unit(5, 'line')) +
      guides(shape = guide_legend(override.aes = list(size = 10)))
    ggsave(plot = g, 
           filename = paste0('../png/maps/heat_stress_daily/map_africa_', nme, '_', y1, '_', y2, '.png'),
           units = 'in',
           width = 12,
           height = 10,
           dpi = 300)
    
  } else {
    
    g <- ggplot(tbl)  +
      geom_tile(aes(x = x, y =  y, fill = value)) +
      facet_wrap(~ year) +
      scale_fill_gradientn(colours = RColorBrewer::brewer.pal(n = 8, name = "YlOrRd"), 
                           na.value = 'white', limits = c(0, 365)) +
      theme_bw() +
      geom_polygon(data = shp, aes(x = long, y = lat, group = group), color = 'grey', fill = 'NA') +
      coord_equal() +
      labs(title = paste0('Heat stress ', nme), fill = 'days',  x = 'Longitude', y = 'Latitude') +
      theme(legend.position = 'bottom',
            plot.title = element_text(hjust = 0.5),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            legend.key.width = unit(5, 'line')) +
      guides(shape = guide_legend(override.aes = list(size = 10)))
    ggsave(plot = g, 
           filename = paste0('../png/maps/heat_stress_daily/map_africa_', nme, '_', y1, '_', y2, '.png'),
           units = 'in',
           width = 12,
           height = 10,
           dpi = 300)
    
  }
  
}

# Load data ---------------------------------------------------------------
pth <- '../raster/indicators/heat_stress'
mss <- c('maiz_s2', 'maiz_s2')

fls <- list.files(pth, full.names = T, pattern = '.tif$')
fls <- grep('Africa', fls, value = TRUE)

shp <- shapefile('../data/shp/base/continents_1.shp')
shp <- shp[shp@data$CONTINENT == 'Africa',]
prd <- list(1983:1994, 1995:2006, 2007:2016)

crp <- unique(str_sub(basename(fls), start = 11, end = nchar(basename(fls)) - 16))

# Extract tables ----------------------------------------------------------
cl <- makeCluster(10)
registerDoSNOW(cl)

vls_all <- foreach(i = 1:length(crp), .verbose = TRUE) %dopar%{
  library(tidyverse)
  library(raster)
  library(rgdal)
  library(rgeos)
  library(stringr)
  extract_table(cr = crp[i])
}

# Read the results tables -------------------------------------------------
vls_all <- list.files('../rds/indicators/heat_stress_daily', full.names = TRUE, pattern = '.rds')
vls_all <- grep(paste0(mss, collapse = '|'), vls_all, value = TRUE)
yrs_all <- prd
nme_all <- gsub('daily_tbl_', '', basename(vls_all))
nme_all <- gsub('.rds', '', nme_all)
nme_all<- 'maize_s1'

# To make the maps --------------------------------------------------------
for(i in 1:length(nme_all)){
  for(j in 1:length(yrs_all)){
    # print(nme_all[i])
    print(yrs_all[j])
    make_map(fle = vls_all[i], yrs = yrs_all[j], nme = 'cotton')
  }
}








