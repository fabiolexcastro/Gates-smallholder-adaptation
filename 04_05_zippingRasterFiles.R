
# Load libraries ----------------------------------------------------------
require(pacman)
pacman::p_load(raster, rgdal, rgeos, stringr, sf, tidyverse, zip)

# Functions to use --------------------------------------------------------
my_zip <- function(crp){
  print(crp)
  fls <- grep(crp, fls, value = TRUE)
  zip::zipr(zipfile = paste0('../zip/', crp, '.zip'), files = fls)
  print('Done!')
}

# Load data ---------------------------------------------------------------
fls <- list.files('../raster/indicators/heat_stress', full.names = TRUE, pattern = '.tif')
fls <- grep('Africa', fls, value = TRUE)

nms <- grep('1983', fls, value = T)
nms <- str_sub(basename(nms), start = 11, end = nchar(basename(nms)) - 16)

map(.x = nms, .f = my_zip)
