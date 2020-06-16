#JRV - mockup tests
#June 2020

#working directory
wd <- "~/work"
hdir <- paste(wd,"/hazard_layers/aridity",sep="")
odir <- paste(wd,"/hazard_layers/mockup_tests",sep="")
if (!file.exists(odir)) {dir.create(odir)}

#libraries
library(raster)
library(rgdal)
library(tidyverse)
library(RColorBrewer)

#load shapefile of Africa
sh_ctry <- readOGR(paste(wd,"/Africa_shp/African_continet.shp",sep=""))
sh_ctry <- sh_ctry[which(sh_ctry$ADM0_NAME %in% "Kenya"),]
sh_xt <- extent(sh_ctry)
sh_xt@xmin <- sh_xt@xmin-1; sh_xt@ymin <- sh_xt@ymin-1; sh_xt@xmax <- sh_xt@xmax+1; sh_xt@ymax <- sh_xt@ymax+1

#given classes
m <- c(0, 40, 1,   #no significant
       40, 60, 2,  #moderate
       60, 80, 3,  #severe
       80, 101, 4) #extreme
rclmat <- matrix(m, ncol=3, byrow=TRUE)

#load and crop historical layer
rsh <- stack(c(paste(hdir,"/aridity_thornthwaite_hist.tif",sep=""),
              list.files(hdir, pattern="aridity_thornthwaite_fut_mmm_", full.names=TRUE)))
rsh <- crop(rsh, sh_xt)

#create mask from shapefile
msk <- rasterize(sh_ctry, rsh)

#mask layer to country
rsh <- mask(rsh, msk)

#calculate frequency table
ftable <- data.frame()
for (k in 1:nlayers(rsh)) {
    #k <- 1
    cat("processing",names(rsh)[k],"\n")
    rst <- rsh[[k]]; rst1 <- raster(rst)
    
    #reclassify raster
    for (i in 1:nrow(rclmat)) {
        #i <- 1
        rst1[which(rst[] >= rclmat[i,1] & rst[] < rclmat[i,2])] <- rclmat[i,3]
        val_i <- length(which(rst1[] == rclmat[i,3])) / length(which(!is.na(rst[])))
        out_row <- data.frame(layer=k, 
                              layer_name=names(rsh)[k], 
                              class=rclmat[i,3],
                              areafrac=val_i)
        ftable <- rbind(ftable, out_row)
    }
    
    #to points
    rs_vls <- rasterToPoints(rst1) %>% 
      as_tibble() %>% 
      setNames(c('x', 'y', 'value'))
    rs_vls$class <- factor(paste(rs_vls$value), 
                           labels=c("no sig. stress", "moderate", "severe", "extreme"))
    
    #plot details (colorbar, limits, breaks)
    plt <- RColorBrewer::brewer.pal(n = 9, name = "YlOrRd")
    brks <- c(1:4)
    
    #plot
    g1 <- ggplot(rs_vls)  +
      geom_tile(aes(x = x, y =  y, fill = class)) +
      scale_fill_manual(values = plt[c(1,3,6,9)]) +
      geom_polygon(data=sh_ctry, aes(x=long, y=lat, group=group), 
                   fill=NA,color="black", size=0.5) +
      theme_bw() +
      coord_fixed(ratio=1, xlim = extent(rst1)[1:2], ylim = extent(rst1)[3:4]) +
      labs(x = '', y = '', fill = '', caption = '') +
      theme(legend.position = 'top',
            plot.title = element_text(size = 20, hjust = 0.5, face = 'bold'),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            legend.key.width = unit(1.5, 'line'),
            strip.background = element_blank(),
            strip.text = element_blank())
    #print(g1)
    if (!file.exists(paste(odir,"/",names(rsh)[k],".png",sep=""))) {
        ggsave(plot = g1, 
               filename = paste(odir,"/",names(rsh)[k],".png",sep=""), 
               units = 'in', width = 10, height = 10, dpi = 300)
    }
}

#organize ftable
ftable$class_name <- factor(paste(ftable$class), levels=paste(4:1),
                           labels=rev(c("no sig. stress", "moderate", "severe", "extreme")))
ftable$scenario <- gsub("aridity_thornthwaite_", "", ftable$layer_name)
ftable$scenario <- gsub("fut_mmm_", "", ftable$scenario)
ftable$period <- unlist(lapply(strsplit(ftable$scenario, split="_", fixed=TRUE), FUN=function(x) {return(x[1])}))
ftable$scenario <- unlist(lapply(strsplit(ftable$scenario, split="_", fixed=TRUE), FUN=function(x) {return(x[2])}))
hpart1 <- hpart2 <- ftable[which(is.na(ftable$scenario)),]
hpart1$scenario <- "rcp4.5"; hpart2$scenario <- "rcp8.5"
ftable <- ftable[which(!is.na(ftable$scenario)),]
ftable <- rbind(hpart1,hpart2,ftable)
ftable$period <- factor(ftable$period, levels=c("hist", "2030", "2050"), labels=c("hist", "2030", "2050"))
ftable$scenario <- factor(ftable$scenario, levels=c("rcp4.5", "rcp8.5"), labels=c("RCP 4.5", "RCP 8.5"))

#stacked barplot
fplot <- ggplot(ftable, aes(fill=class_name, y=areafrac, x=period)) + 
            geom_bar(stat="identity") +
            scale_fill_manual(values = rev(plt[c(1,3,6,9)])) +
            labs(x = '', y = 'Cumulative percentage of area [%]', fill = '', caption = '') +
            theme(legend.position = 'bottom',
                plot.title = element_text(size = 20, hjust = 0.5, face = 'bold'),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                legend.key.width = unit(1.5, 'line'),
                strip.background = element_blank(),
                strip.text = element_text(size=10, hjust=0.5, face='bold')) +
            facet_wrap(~scenario)
print(fplot)
ggsave(plot = fplot, 
           filename = paste(odir,"/aridity_thornthwaite_freqplot.png",sep=""), 
           units = 'in', width = 6.5, height = 5, dpi = 300)

#write csv with values
ftowrite <- ftable; ftowrite$scenario <- paste(ftowrite$scenario); ftowrite$period <- paste(ftowrite$period)
ftowrite <- ftowrite[-which(ftable$scenario == "RCP 4.5" & ftable$period == "hist"),]
ftowrite$scenario[which(ftowrite$period == "hist")] <- "hist"
write.csv(ftowrite, paste(odir,"/aridity_thornthwaite.csv",sep=""), row.names=FALSE)
