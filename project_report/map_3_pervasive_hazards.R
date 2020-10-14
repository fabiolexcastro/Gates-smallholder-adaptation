#JRV / pervasive hazard map

#libraries
library(raster)
library(tidyverse)
library(rgdal)

#working directories
wd <- "~/work/hazard_layers"

#output maps
outdir <- paste(wd, "/report_maps", sep="")
if (!file.exists(outdir)) {dir.create(outdir)}

#shapefile of Africa
sh_ctry <- readOGR(paste("~/work/Africa_shp/African_continet.shp",sep=""))

#mask
msk <- raster(paste(wd,"/chirps_cv/cvr_africa.tif",sep=""))
msk[which(!is.na(msk[]))] <- 1

#spam crop mask
spam <- list.files(paste(wd,"/spam_2017_mask/",sep=""), pattern="\\_A.tif")
spam <- stack(paste(wd,"/spam_2017_mask/",spam,sep="")) %>%
            crop(., msk)
spam <- sum(spam, na.rm=TRUE)
spam[which(spam[] == 0)] <- NA
spam <- resample(spam, msk)
spam[which(!is.na(spam[]))] <- 1

#load hazard layers
#aridity
hz1 <- raster(paste(wd,"/aridity/aridity_thornthwaite_fut_mmm_2030_rcp8.5.tif",sep="")) %>%
        raster::crop(., msk) %>%
        raster::resample(., msk)
hz1_rc <- raster(hz1)
hz1_rc[which(hz1[] <= 40)] <- 0
hz1_rc[which(hz1[] > 40 & hz1[] <= 60)] <- 1
hz1_rc[which(hz1[] > 60 & hz1[] <= 80)] <- 2
hz1_rc[which(hz1[] > 80)] <- 3

#chrips cv
hz2 <- raster(paste(wd,"/chirps_cv/cv_rcp85_30s.tif",sep="")) %>%
        raster::crop(., msk) %>%
        raster::resample(., msk)
hz2_rc <- raster(hz2)
hz2_rc[which(hz2[] <= 15)] <- 0
hz2_rc[which(hz2[] > 15 & hz2[] <= 25)] <- 1
hz2_rc[which(hz2[] > 25 & hz2[] <= 35)] <- 2
hz2_rc[which(hz2[] > 35)] <- 3

#dry days
hz3 <- raster(paste(wd,"/dry_days_future/yearly/dry_days_rcp85_30s/dry_days_rcp85_30s_2020_2049_mean.tif",sep="")) %>%
        raster::crop(., msk) %>%
        raster::resample(., msk)
hz3_rc <- raster(hz3)
hz3_rc[which(hz3[] <= 15)] <- 0
hz3_rc[which(hz3[] > 15 & hz3[] <= 20)] <- 1
hz3_rc[which(hz3[] > 20 & hz3[] <= 25)] <- 2
hz3_rc[which(hz3[] > 25)] <- 3

#heat stress crops
hz4 <- raster(paste(wd,"/heat_stress_days_future/heat_stress_days_rcp85_30s/generic/heat_crop_generic_2020_2049_Africa_mean.tif",sep="")) %>%
        raster::crop(., msk) %>%
        raster::resample(., msk)
hz4_rc <- raster(hz4)
hz4_rc[which(hz4[] == 0)] <- 0
hz4_rc[which(hz4[] > 0 & hz4[] <= 5)] <- 1
hz4_rc[which(hz4[] > 5 & hz4[] <= 10)] <- 2
hz4_rc[which(hz4[] > 10)] <- 3

#heat stress cattle
hz5 <- raster(paste(wd,"/thi/thi_max_2030_rcp85.asc",sep="")) %>%
        raster::crop(., msk) %>%
        raster::resample(., msk)
hz5_rc <- raster(hz5)
hz5_rc[which(hz5[] <= 72)] <- 0
hz5_rc[which(hz5[] > 72 & hz5[] <= 78)] <- 1
hz5_rc[which(hz5[] > 78 & hz5[] <= 89)] <- 2
hz5_rc[which(hz5[] > 89)] <- 3

#heat stress humans
hz6 <- raster(paste(wd,"/hi/hi_max_2030_rcp85.asc",sep="")) %>%
        raster::crop(., msk) %>%
        raster::resample(., msk)
hz6_rc <- raster(hz5)
hz6_rc[which(hz6[] <= 26)] <- 0
hz6_rc[which(hz6[] > 26 & hz6[] <= 32)] <- 1
hz6_rc[which(hz6[] > 32 & hz6[] <= 41)] <- 2
hz6_rc[which(hz6[] > 41)] <- 3

#total hazard score
hz_score <- mean(stack(hz1_rc, hz2_rc, hz3_rc, hz4_rc, hz5_rc, hz6_rc))
hz_score <- mask(hz_score, spam)

#plot(hz_score, col=rev(heat.colors(9))); plot(sh_ctry, add=TRUE)

#number of hazards in categories 2 or 3
for (i in 1:6) {
    hz_lyr <- get(paste("hz",i,"_rc",sep=""))
    hzx <- raster(hz_lyr)
    hzx[which(hz_lyr[] < 2)] <- 0
    hzx[which(hz_lyr[] >= 2)] <- 1
    if (i == 1) {hz_count <- hzx} else {hz_count <- hz_count + hzx}
    rm(list=c("hz_lyr","hzx"))
}
hz_count <- mask(hz_count, spam)
#plot(hz_count, col=rev(heat.colors(3))); plot(sh_ctry, add=TRUE)

#check correlation between the two plotting variables
rs_vls <- rasterToPoints(hz_score) %>% 
        as_tibble() %>% 
        setNames(c('x', 'y', 'score')) %>%
        dplyr::mutate(., count=raster::extract(hz_count, .[,c('x','y')]))
plot(rs_vls$score, rs_vls$count, pch=20, cex=0.5)

#save outputs
writeRaster(hz_score, paste(outdir, "/hazard_score.tif", sep=""))
writeRaster(hz_count, paste(outdir, "/hazard_count.tif", sep=""))
