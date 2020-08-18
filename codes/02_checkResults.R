

library(tidyverse)
library(stringr)
library(hablar)

# Check results -----------------------------------------------------------
fldr <- list.files('../tif/future', full.names = TRUE)

date <- data.frame(period = c(rep(c('2020_2049'), length(1981:2010)), rep(c('2040_2069'), length(1981:2010))),
                   year_1 = c(1981:2010, 1981:2010),
                   year_2 = c(2020:2049, 2040:2069))


checking <- function(fld){
  
  fld <- fldr[1]
  dfm <- data_frame(path = list.files(fld, full.names = FALSE, pattern = '.tif$'))
  prt <- str_split(string = dfm$path, pattern = '_')
  
  dfm <- dfm %>% 
    mutate(year = sapply(1:length(prt), function(x) prt[[x]][2]),
           month = sapply(1:length(prt), function(x) prt[[x]][3]),
           day = gsub('.tif', '', sapply(1:length(prt), function(x) prt[[x]][4]))) %>% 
    retype()
  
  smm <- dfm %>%
      group_by(year) %>%
      summarise(count = n()) %>%
      ungroup()

  smm <- smm %>%
      filter(count < 365)
  
  return(smm)
  
}

rslt <- lapply(1:4, function(x) checking(fld = fldr[x]))
rslt

# RCP 45 2020 - 2069
rcp45_2030 <- checking(fld = fldr[1])
rcp45_2030 <- drop_na(rcp45_2030)
rcp45_2030 <- rcp45_2030 %>% filter(year == 2051)
rcp45_2030 %>% 
  group_by(month) %>% 
  summarize(count = n())

# RCP 45 2040 - 2069
rcp45_2050 <- checking(fld = fldr[2])
rcp45_2050 <- drop_na(rcp45_2050)
rcp45_2050 <- rcp45_2050 %>% filter(year == 2063)
rcp45_2050 %>% 
  group_by(month) %>% 
  summarize(count = n())

# RCP 85 2020 - 2049
rcp85_2030 <- checking(fld = fldr[3])
rcp85_2030 <- drop_na(rcp85_2030)
rcp85_2030 <- rcp85_2030 %>% filter(year == 2020)
rcp85_2030 %>% 
  group_by(month) %>% 
  summarize(count = n())

# RCP 85 2040 - 2069
rcp85_2050 <- checking(fld = fldr[4])
rcp85_2050 <- drop_na(rcp85_2050)
rcp85_2050 <- rcp85_2050 %>% filter(year == 2063)
rcp85_2050 %>% 
  group_by(month) %>% 
  summarize(count = n())
