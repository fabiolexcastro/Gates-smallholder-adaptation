

# Load the libraries ------------------------------------------------------
require(pacman)
pacman::p_load(raster, rgdal, rgeos, stringr, terra, sf, tidyverse, glue, gtools, hablar, readxl)
rm(list = ls())
source('./02_03_functions.R')

# Load data ---------------------------------------------------------------
fls <- list.files('../rds/chirps', full.names = T, pattern = '.rds')
cnt <- c('Africa', 'Asia', 'North America', 'Oceania', 'South America')
lbl <- read_excel('../tbl/seasons.xlsx') %>% 
  mutate(mnt = paste(m1, m2, m3, sep = '-')) %>% 
  dplyr::select(season, mnt)
yrs <- as.character(1981:2019)

# Processing --------------------------------------------------------------
lapply(15:length(yrs), function(k)drySeason(yr = yrs[k], cn = cnt[1]))

