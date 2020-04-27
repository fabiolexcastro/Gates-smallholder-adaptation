

# Load libraries ------------------------------------------------------------
require(pacman)
pacman::p_load(raster, rgdal, rgeos, stringr, terra, sf, tidyverse, glue, gtools, hablar, readxl)
rm(list = ls())


# Load Data -----------------------------------------------------------------
cnt <- c('Africa', 'Asia', 'North America', 'Oceania', 'South America')
fls <- '../rds/chirps/seasonDry' %>% list.files(., full.names = T)

tbl2rst <- function(cn){
  cn <- cnt[1]
  
  print(cn)
  fl <- grep(cn, fls, value = T)
  
  tb <- map(.x = fl, .f = readRDS)
  rs <- lapply(1:12, function(k){
    d <- tb[[k]]
    p1 <- is.infinite(d$n)
    p2 <- is.nan(d$n)
    p1 <- which(p1 == TRUE)
    p2 <- which(p2 == TRUE)
    d <- d[-p1,]
    d <- d[-p2,]
    a <- d %>% pull(n) %>% mean(., na.rm = T)
    return(a)
  })
  
  
}


