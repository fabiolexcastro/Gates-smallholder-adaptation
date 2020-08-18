
# Load libraries ----------------------------------------------------------
require(pacman)
pacman::p_load(raster, rgdal, doSNOW, foreach, parallel, magrittr, rgeos, gtools, stringr, sf, glue, tidyverse, terra, gtools)
g <- gc(reset = TRUE)
rm(list = ls())

# Functions to use --------------------------------------------------------
dry_days <- function(rcp, prd, yrs){
  
  rcp <- 'rcp45'
  prd <- '2020_2049'
  yrs <- years[1]
  
  ftr <- grep(rcp, fld_ftr, value = TRUE) %>% 
    grep(prd, ., value = TRUE) %>% 
    mixedsort() %>% 
    list.files(., full.names = TRUE, pattern = '.tif$') 
  ftr <- ftr[grep(paste0('_', yrs, '_'), basename(ftr), value = F)]
  
  # yyy <- str_sub(basename(ftr), start = 5, end = nchar(basename(ftr)) - 4)
  # yyy <- unique(yyy)
  
  rsl <- lapply(1:12, function(k){
    
    print(k)
    
    x <- ifelse(k < 10, paste0('0', k), k)
    x <- paste0('_', x, '_')
    x <- grep(x, ftr, value = TRUE)
    x <- stack(x)
    x <- rasterToPoints(x)
    x <- as_tibble(x)
    x$gid <- 1:nrow(x)
    x <- x %>% gather(var, value, -x, -y, -gid)
    x$year <- str_sub(x$var, start = 5, end = 8)
    x$mont <- str_sub(x$var, start = 10, end = 11)
    x$day  <- str_sub(x$var, start = 13, end = 14)  
    
    y <- x %>% 
      mutate(binary = if_else(value == 0, 1, 0)) %>% 
      group_by(gid, x, y, year, mont) %>% 
      summarise(count = sum(binary)) %>% 
      ungroup()
    print('Done!')
    return(y)
    
  })
  
  rsl <- bind_rows(rsl)
  rsl <- rsl %>% mutate(mont = as.numeric(mont))
  
  rst <- lapply(1:12, function(k){
    print(k)
    rs <- rsl %>% filter(mont == k)
    rs <- rasterFromXYZ(rs[,c(2, 3, 6)])
    print('Done!')
    return(rs)
  })
  
  Map('writeRaster', x = rst, filename = paste0('../tif/dryDays/', rcp, '_', pr, '/', 'dryDays_', yr, '_', 1:12, '.tif'), overwrite = TRUE)
  
  
  
}

# Load data ---------------------------------------------------------------
fld_ftr <- list.files('../tif/future', full.names = TRUE) %>% grep('rcp', ., value = TRUE)
rcps <- c('rcp45', 'rcp45','rcp85', 'rcp85')
prds <- c('2020_2049', '2040_2069')
years <- 2020:2049

#  Apply the function -----------------------------------------------------
for(i in 2:length(yrs)){
  print(yrs[i])
  dry_days(rcp = rcps[1], pr = prd[1], yr = yrs[i])
}

done <- list.files('./tif/chirps/dryDays/future')
done <- str_sub(done, 19, 22)
done <- unique(done)
