
# 
require(pacman)
pacman::p_load(raster, rgdal, rgeos, stringr, velox, sf, tidyverse)
rm(list = ls())


# Read
fls <- list.files('../rds/chirps/complet')

# Proof
tbl <- readRDS(fls[1])