
# Load libraries ------------------------------------------------------------
require(pacman)
pacman::p_load(raster, rgdal, rgeos, stringr, terra, sf, tidyverse, glue, gtools, hablar, readxl)
rm(list = ls())

# Functions to use --------------------------------------------------------
dr_stress <- function(PREC){
  runs <- rle(PREC < 1)
  cons_days <- max(runs$lengths[runs$values==1], na.rm=TRUE)
  return(cons_days)
}
dryDays <- function(cn, yr){
  print(yr)
  fl <- grep(cn, fls, value = T) %>% 
    grep(yr, ., value = T)
  tb <- readRDS(fl) %>% 
    mutate(id = 1:nrow(.)) %>% 
    gather(var, value, -x, -y, -id) %>% 
    as_tibble() %>% 
    mutate(year = str_sub(var, 3, 6),
           month = str_sub(var, 8, 9),
           day = str_sub(var, 11, 12)) 
  cr <- tb %>% distinct(id, x, y)
  dy <- tb %>% 
    group_by(id) %>% 
    summarise(n_days = dr_stress(value)) %>% 
    ungroup()
  rs <- inner_join(dy, cr, by = 'id') %>% 
    mutate(year = yr)
  print(' Done!!! ')
  return(rs)
}

# Load data ---------------------------------------------------------------
fls <- list.files('../rds/chirps', full.names = T, pattern = '.rds')
cnt <- c('Africa', 'Asia', 'North America', 'Oceania', 'South America')
yrs <- as.character(1981:2019)

# Apply the function ------------------------------------------------------
afr <- lapply(1:length(yrs), function(k)dryDays(yr = yrs[k], cn = cnt[1]))
afr <- bind_rows(afr)
saveRDS(afr, '../rds/chirps/consDaysYearly/dryConsDaysYear_Africa.rds')

oce <- lapply(1:length(yrs), function(k)dryDays(yr = yrs[k], cn = cnt[4]))
oce <- bind_rows(oce)
saveRDS(oce, '../rds/chirps/consDaysYearly/dryConsDaysYear_Oceania.rds')
