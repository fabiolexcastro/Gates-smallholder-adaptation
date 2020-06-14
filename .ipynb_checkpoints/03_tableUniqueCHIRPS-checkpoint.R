g=gc()
rm(list=ls())
require(raster);require(dplyr); require(foreach); require(parallel)
path <- '//dapadfs/data_cluster_4/observed/gridded_products/chirps/daily/32bits/'
base <- raster(paste0(path, "chirps-v2.0.1981.01.01.tif"))
fls <- list.files(path = path, pattern = ".tif",full.names = TRUE)
year <- 1985:2019
year <- as.character(year)

n_a <- shapefile("//dapadfs/workspace_cluster_13/GATES/data/shp/base/continents.shp")
proj4string(n_a) <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
n_a <- n_a[n_a@data$CONTINENT == "North America",]


library(foreach)
library(parallel)
library(doParallel)
cores<- detectCores()
cl<- makeCluster(cores-10)
registerDoParallel(cl)

system.time(prec<-foreach(i=1:length(year)) %dopar% {
  mylist <- list()
  require(raster)
  fls_year <- fls[grep(year[i], fls)]
  r <- raster:::stack(fls_year)
  df <- lapply(1:nlayers(r) , function(j){
    raster <- r[[j]]
    
    raster <- crop(raster, c(-179.1411, -40, 5.499059, 83.62746))
    n_a <- crop(n_a, c(-179.1411, -40, 5.499059, 83.62746))
    raster <- crop(raster,c(-179.1411, -40, 5.499059, 83.62746))
    raster <- mask(raster,n_a)
    df <- as.data.frame(rasterToPoints(r))
    col <- colnames(df)
    col <- col[3] 
    #   dte <- stringr::str_sub(col, start = 13, end = nchar(col))
    #   dte <- gsub('\\.', '-', dte)
    #   colnames(df) <- c("x","y",dte)
    #   cellID <- cellFromXY(base, df[,1:2])
    #   df <- cbind(cellID, df)
    #   df$x <- NULL 
    #   df$y <- NULL
    
    
  })
  
  
  system.time(df <- as.data.frame(rasterToPoints(r)))
  df[,3][which(df[,3] < 0)] <- NA 
  df <- na.omit(df)
  cellID <- cellFromXY(base, df[,1:2])
  df <- cbind(cellID,df)
  df[,2:3]<- NULL
  mylist[[i]] <- df
 })
stopCluster(cl)

tabla <- do.call(cbind, prec)
n<- seq(3,5000,by=2)
tabla1 <- tabla[,-n]
saveRDS(tabla1, "//dapadfs/workspace_cluster_13/GATES/rds/chirps/complet/prec_5.rds")
  

