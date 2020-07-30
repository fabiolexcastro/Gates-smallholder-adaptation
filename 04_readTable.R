
require(tidyverse)
require(stringr)

dryDays <- function(pth){
  # pth <- fls[1]
  tbl <- readRDS(pth)
  tbl <- tbl %>% 
    mutate(id = 1:nrow(.)) %>% 
    gather(var, value, -id, -x, -y)
  tbl <- tbl %>% 
    mutate(year = str_sub(var, 3, 6),
           month = str_sub(var, 8, 9),
           day = str_sub(var, 11, 12)) %>% 
    dplyr::select(id, x, y, year, month, day, value)
  
  ##
  smm <- tbl %>% 
    retype() %>% 
    filter(value < 1) %>% 
    group_by(id, x, y, year, month) %>% 
    dplyr::summarise(count = n()) %>% 
    ungroup()
  mnt <- unique(smm$month)
  
  yrs2rst <- function(m){
    tbl <- smm %>% 
      filter(month == m) %>% 
      dplyr::select(x, y, count)
    # rst <- rasterFromXYZ(tbl) 
    saveRDS(object = tbl, file = paste0('../rds/chirps/dryDays/dry_mnt', m, '_', basename(pth)))
  }
  
  lapply(1:12, function(k) yrs2rst(m = k))
  
  print('Done!')
}


fls <- list.files('../rds/chirps', full.names = T, pattern = '.rds$')
# fle <- '../rds/chirps/South America_chirps_1982.rds'
map(.x = fls, .f = dryDays)




