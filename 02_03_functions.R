
dr_stress <- function(PREC, p_thresh = 1){
  runs <- rle(PREC < p_thresh)
  cons_days <- max(runs$lengths[runs$values==1], na.rm=TRUE)
  return(cons_days)
}

ssn <- function(tabla, mnt, yrs){
  # mnt <- c('01', '02', '03')
  # yrs <- c('1982', '1982', '1982')
  print(paste0('Processing ----------> ',  mnt, ' ', yrs))
  tb <- tabla %>% 
    filter(year %in% yrs & month %in% mnt[[1]] |
             year %in% yrs & month %in% mnt[[2]] |
             year %in% yrs & month %in% mnt[[3]]) %>% 
    dplyr::select(-var)
  cr <- tb %>% distinct(id, x, y)
  sm <- tb %>% 
    group_by(id) %>% 
    dplyr::summarise(n = dr_stress(value)) %>% 
    ungroup() %>% 
    inner_join(., cr, by = 'id') %>% 
    dplyr::select(id, x, y, n) 
  sm <- sm %>% 
    mutate(ssn = paste(mnt[1], mnt[2], mnt[3], sep = '-')) %>% 
    inner_join(., lbl, by = c('ssn' = 'mnt')) %>% 
    dplyr::select(id, x, y, season, n)
  print('Done')
  return(sm)
}

drySeason <- function(yr, cn){
  
  # Proof
  # yr <- yrs[1]; cn <- cnt[1]
  
  # Processing
  tb <- grep(cn, fls, value = T) %>% 
    grep(yr, ., value = T) %>% 
    readRDS() %>% 
    mutate(id = 1:nrow(.)) %>% 
    gather(var, value, -x, -y, -id) %>% 
    as_tibble() %>% 
    mutate(year = str_sub(var, 3, 6),
           month = str_sub(var, 8, 9),
           day = str_sub(var, 11, 12)) 
  
  rs <- lapply(1:10, function(k){
    r <- lbl[k,2] %>% as.character %>% str_split(., pattern = '-')
    r <- r[[1]]
    s <- ssn(tabla = tb, mnt = r, yrs = c(yr,yr,yr))
  })
  
  # Write the final table
  rs <- bind_rows(rs)
  saveRDS(rs, paste0('../rds/chirps/seasonDry/ssn_dry_', cn, '_', yr, '.rds'))
  
}