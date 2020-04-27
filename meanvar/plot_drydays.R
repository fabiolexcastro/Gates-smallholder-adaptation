#JRV - plot mean and variability indicators
#April 2020

#working directory
wd <- "~/work"
hdir <- paste(wd,"/hazard_layers",sep="")
sdir <- paste(wd,"/spam2010V1r1",sep="")

#libraries
library(raster)
library(rgdal)
library(stringr)
library(sf)
library(tidyverse)
library(RColorBrewer)
library(glue)
library(ggspatial)
#library(cowplot)
#library(ggpubr)
#library(shadowtext)

#layer name
lname <- "dry_days"
period <- "hist"
toplot <- "mean" #c("mean","cv","95p","05p")

#crop name
#list of crops: ACOF, BANA, BARL, BEAN, CASS, CHIC, CNUT, COCO, COTT, COWP, GROU, LENT, MAIZ
#               OCER, OFIB, OILP, OOIL, OPUL, ORTS, PIGE, PLNT, PMIL, POTA, RAPE, RCOF, REST
#               RICE, SESA, SMIL, SORG, SOYB, SUGB, SUGC, SUNF, SWPO, TEAS, TEMF, TOBA, TROF
#               VEGE, WHEA, YAMS
cropname <- "MAIZ"

#name of crop calendar
calname <- "Maize"

#technology, list of technologies
#*_A	all technologies together, ie complete crop
#*_H	rainfed high inputs portion of crop
#*_I	irrigated portion of crop
#*_L	rainfed low inputs portion of crop
#*_R	rainfed portion of crop (= A - I, or H + L + S)
#*_S	rainfed subsistence portion of crop
technol <- "S"

#years
yi <- 1981
yf <- 2019

#load world shapefile
sh <- shapefile(paste(wd,"/world_shp/ne_50m_admin_0_countries.shp",sep=""))
sh_st <- st_as_sf(sh)

#load hazard layer
rsh <- raster(paste(hdir,"/",lname,"_",period,"/yearly/",lname,"_",yi,"_",yf,"_",toplot,".asc",sep=""))

#load harvested area layer
rsa <- raster(paste(sdir,"/spam2010V1r1_global_H_",cropname,"_",technol,".tif",sep=""))

#resample, crop harvested area into hazard layer, and create 0,1 layer
rsa <- resample(rsa, rsh)
rsa <- crop(rsa, rsh)
rsa[which(rsa[] < 0.01)] <- NA
rsa[which(!is.na(rsa[]))] <- 1

#to points
rs_vls <- rasterToPoints(rsh) %>% 
  as_tibble() %>% 
  setNames(c('x', 'y', 'value'))

#plot details (colorbar, limits, breaks)
plt <- RColorBrewer::brewer.pal(n = 9, name = "YlOrRd")
lims <- c(0, 365)
brks <- round(seq(0,365,by=36.5),0)

#plot
g1 <- ggplot(rs_vls)  +
  geom_tile(aes(x = x, y =  y, fill = value)) +
  scale_fill_gradientn(colours = plt, 
                       na.value = 'white',
                       limits=c(0,1),
                       breaks=brks) +
  geom_polygon(data=sh, aes(x=long, y=lat, group=group), 
               fill=NA,color="black", size=0.5) +
  theme_bw() +
  coord_fixed(ratio=1,xlim = c(-150,150), ylim = extent(rsh)[3:4]) +
  labs(x = 'Longitude', y = 'Latitude', fill = 'days', caption = 'Smallholder Adaptation Atlas') +
  theme(legend.position = 'bottom',
        plot.title = element_text(size = 20, hjust = 0.5, face = 'bold'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.key.width = unit(15, 'line'),
        strip.background = element_blank(),
        strip.text = element_blank()) +
  guides(shape = guide_legend(override.aes = list(size = 10)))
ggsave(plot = g1, filename = paste(hdir,"/",lname,"_",period,"/plot_",lname,"_",period,".png",sep=""), units = 'in', width = 20, height = 9, dpi = 300)


#calculate dry days for cropping season
#load cropping season data
rsc_p <- raster(paste(wd,"/crop_calendar/",calname,".crop.calendar.fill.nc",sep=""),varname="plant")
rsc_p <- resample(rsc_p, rsh)
rsc_p <- crop(rsc_p, rsh)

rsc_h <- raster(paste(wd,"/crop_calendar/",calname,".crop.calendar.fill.nc",sep=""),varname="harvest")
rsc_h <- resample(rsc_h, rsh)
rsc_h <- crop(rsc_h, rsh)

#load monthly hazard layer data
rsh <- stack(paste(hdir,"/",lname,"_",period,"/monthly/",lname,"_",yi,"_",yf,"_",1:12,"_",toplot,".asc",sep=""))
rsh <- stack(rsh, rsc_p, rsc_h)
rsh <- readAll(rsh)
rsh <- mask(rsh, rsa)

#to test function
x <- as.numeric(raster::extract(rsh, data.frame(lon=-75,lat=3.75)))

#function to calculate total of index
calc_gsi <- function(x) {
    if (NA %in% x) {
        y <- NA
    } else {
        #gsi <- x[13]; gsf <- x[14]
        #if (gsi <= gsf) {gsl <- length(gsi:gsf)} else {gsl <- length(c(gsi:365,1:gsf))}
        mondata <- x[1:12]
        gsi <- round(x[13] / 30, 0); if (gsi == 0) {gsi <- 12}
        gsf <- round(x[14] / 30, 0); if (gsf == 13) {gsf <- 1}
        if (gsi <= gsf) {
            gsl <- length(gsi:gsf) * 30
            y <- sum(mondata[gsi:gsf],na.rm=T) / gsl
        } else if (gsi > gsf) {
            gsl <- length(c(gsi:12,1:gsf)) * 30
            y <- (sum(mondata[gsi:12],na.rm=T) + sum(mondata[1:gsf],na.rm=T)) / gsl
        }
    }
    return(y)
}

#run calculation
rsout <- calc(rsh, fun=calc_gsi, forceapply=T)
rsout[which(rsout[] > 1)] <- 1

#to points
rs_vls <- rasterToPoints(rsout) %>% 
  as_tibble() %>% 
  setNames(c('x', 'y', 'value'))

#update limits and breaks
lims <- c(0, 1)
brks <- seq(0,1,by=0.1)

#plot
g1 <- ggplot(rs_vls)  +
  geom_tile(aes(x = x, y =  y, fill = value)) +
  scale_fill_gradientn(colours = plt, 
                       na.value = 'white',
                       limits=lims,
                       breaks=brks) +
  geom_polygon(data=sh, aes(x=long, y=lat, group=group), 
               fill=NA,color="black", size=0.5) +
  theme_bw() +
  coord_fixed(ratio=1,xlim = c(-150,150), ylim = extent(rsh)[3:4]) +
  labs(x = 'Longitude', y = 'Latitude', fill = 'days', caption = 'Smallholder Adaptation Atlas') +
  theme(legend.position = 'bottom',
        plot.title = element_text(size = 20, hjust = 0.5, face = 'bold'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.key.width = unit(15, 'line'),
        strip.background = element_blank(),
        strip.text = element_blank()) +
  guides(shape = guide_legend(override.aes = list(size = 10)))
ggsave(plot = g1, filename = paste(hdir,"/",lname,"_",period,"/plot_",lname,"_",period,"_gseason_frac.png",sep=""), units = 'in', width = 20, height = 9, dpi = 300)

