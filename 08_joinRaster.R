

## 
require(pacman)
pacman::p_load(raster, rgdal, rgeos, stringr, velox, sf, tidyverse)
rm(list = ls())
options(scipen = 999, stringsAsFactors = FALSE)

# Functions
mosaicRaster <- function(x, y){
  # x <- 't1_'
  # y <- 981
  
  fl <- grep(x, fls, value = T) %>% 
    grep(y, ., value = T)
  tb <- map(.x = fl, .f = readRDS)
  
  ly <- list()
  for(i in 1:length(tb)){ly[[i]] <- rasterFromXYZ(tb[[i]])}
  
  print(length(ly))
  
  ms1 <- mosaic(x = ly[[1]], y = ly[[2]], fun = 'mean')
  ms2 <- mosaic(x = ly[[3]], y = ly[[4]], fun = 'mean')
  ms3 <- mosaic(x = ms1, y = ms2, fun = 'mean')
  ms4 <- mosaic(x = ms3, y = ly[[5]], fun = 'mean')
  
  print('To write')
  writeRaster(ms4, paste0('../raster/indicators/dryDays/mosaic/dryDays_', y, '_', x, '.tif'), overwrite = TRUE)
  print('Done!!!')
}

# Load data

fls <- list.files('../rds/chirps/dryDays', full.names = TRUE, pattern = '.rds')

vec <- str_split(string = fls, pattern = '_')

mnt <- sapply(1:length(vec), function(k) vec[[k]][2])
yrs <- sapply(1:length(vec), function(k) vec[[k]][5]) %>% gsub('.rds', '', .)
cnt <- sapply(1:length(vec), function(k) vec[[k]][3])

dfm <- data.frame(year = yrs, month = mnt, continent = cnt)
dfm <- dfm %>% 
  mutate(year = str_sub(string = year, start = nchar(year) - 2, end = nchar(year)),
         month = parse_number(month))

dfm %>% 
  group_by(continent) %>% 
  dplyr::summarise(count = n()) %>% 
  ungroup()

dfm %>% 
  filter(continent == 'North America') %>% 
  pull(year) %>% 
  unique()

sub <- dfm %>% 
  filter(year %in% c('981', '982', '983', '984'))


yrs_sub <- c('981', '982', '983', '984')
for(i in 1:12){
  for(j in 1:length(yrs_sub)){
    mosaicRaster(x = i, y = yrs_sub[[j]])
  }
}

