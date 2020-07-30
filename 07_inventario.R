
#
require(pacman)
pacman::p_load(raster, rgdal, rgeos, stringr, velox, tidyverse)

fls <- list.files('../rds/chirps', full.names = F, pattern = '.rds')
fls <- str_split(string = fls, pattern = '_')

dfm <- 
data.frame(continent = sapply(1:length(fls), function(k) fls[[k]][1]),
           year = sapply(1:length(fls), function(k) fls[[k]][3]) %>% gsub('.rds', '', .))

write.csv(dfm, '../tbl/inventario_tables.csv', row.names = FALSE)
