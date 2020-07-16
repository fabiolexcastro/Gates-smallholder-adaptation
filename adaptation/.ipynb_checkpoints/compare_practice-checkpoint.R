#JRV - ERA analyses for adaptation atlas

#working directory
wd <- "~/work"
cimdir <- paste(wd,"/ERA_comparison",sep="")
soldir <- paste(wd,"/soilgrids_global/5km",sep="")
hzdir <- paste(wd,"/hazard_layers",sep="")

#analysis version
vr <- 3

#load library
library(analogues)
library(rgdal)
library(raster)
library(maptools); data(wrld_simpl)
library(data.table)
library(Hmisc)
library(tidyverse)
library(ggplot2)
library(metR)
library(pROC)

#create a mask
load(paste(cimdir,"/wc_prec.rda",sep=""))
msk <- wc_prec[[1]]; rm(wc_prec)
msk[which(!is.na(msk[]))] <- 1

#Practice (n) as follows
#   1. Reduced Tillage (47)
#   2. Mulch (42)
#   3. Agroforestry Pruning (38)
#   4. Green Manure (37)
#   5. Water Harvesting (35)
pr_list <- c("Mulch","Reduced Tillage","Green Manure","Water Harvesting","Agroforestry Pruning")

#read ERA data file
all_data <- read.csv(paste(cimdir,"/ERA_Test_Data.csv",sep="")) %>%
                dplyr::filter(., Product.Simple == "Maize") %>%
                dplyr::filter(., PrName %in% pr_list) %>%
                droplevels(.)

#select practice and organize practice data
#aggregate data up to site level
rnd <- 4
agg_by <- c("Site.ID","PrName")

#Identify outliers
outlier_calc <- function(x){
    return((x < quantile(x)[2] - 3 *  IQR(x)  | x > quantile(x)[4] + 3 * IQR(x)))
}

#determine outliers
all_data <- as.data.table(all_data)
outliers <- unlist(all_data[,R:=1:nrow(all_data)][,list(outliers=list(R[outlier_calc(yi)])), by=agg_by][,outliers])

#Calculate Weighted mean effect size by practice and site combination
wt_grp <- unique(c("Code",agg_by))
data_sites <- all_data[!outliers
                      ][,N.Obs.Study:=.N,by=wt_grp # Recalculate Weightings by study within obervations grouping 
                       ][,Weight.Study:=(Rep^2/(2*Rep))/N.Obs.Study # Recalculate Weightings
                        ][,list(Observations=.N,
                                Studies=length(unique(Code)),
                                Sites=length(unique(Site.ID)),
                                Lat=mean(Latitude),
                                Lon=mean(Longitude),
                                RR=round(weighted.mean(yi,Weight.Study,na.rm=T),rnd),
                                RR.sd=suppressWarnings(round(abs(wtd.var(yi,Weight.Study,na.rm=T))^0.5,rnd))),
                          by=agg_by  # this 
                         ][,pc:=round(100*exp(RR)-100,rnd-1)
                          ][,pc.sd:=round(100*exp(RR.sd)-100,rnd-1)]
data_sites$Label <- paste(data_sites$Site.ID, data_sites$PrName)

#back to data.frame
data_sites <- as.data.frame(data_sites)

#loop practices
for (prname in pr_list) {
    #prname <- pr_list[1]
    cat("...processing practice=", prname, "\n")
    
    #output directory
    pr_odir <- paste(cimdir,"/",gsub(" ", "_", prname, fixed=T),"_v",vr,sep="")
    
    ####
    pdata <- data_sites %>%
                dplyr::filter(., PrName == prname) %>%
                tidyr::drop_na(., Lat, Lon) %>%
                droplevels(.) %>%
                select(., Lon, Lat, RR, PrName)
    
    #cleanup any data points outside African continent, and add positive and negative columns
    pdata <- pdata %>%
                dplyr::mutate(., isvalid=raster::extract(msk, .[,c("Lon","Lat")])) %>%
                tidyr::drop_na(., isvalid) %>%
                dplyr::select(., -isvalid) %>%
                dplyr::mutate(., class=ifelse(.[,'RR'] >= 0.15, "positive", "negative"))

    #load CL's output
    if (prname == "Agroforestry Pruning") {
        rs_cl <- paste(cimdir, "/niche_layers/PRELIMINARY/Current Suitability/Agroforestry-Current-ChanceSuitable.tif",sep="") %>%
                    raster(.)
    } else {
        rs_cl <- paste(cimdir, "/niche_layers/PRELIMINARY/Current Suitability/", 
                       gsub(" ", "", prname, fixed=TRUE), "-Current-ChanceSuitable.tif",sep="") %>%
                    raster(.)
    }
    
    #load max pos and neg
    rs_ps <- raster(paste(pr_odir, "/max_similarity_pos.tif", sep="")) %>%
                raster::crop(., rs_cl) %>%
                raster::resample(., rs_cl, method='bilinear') %>%
                raster::mask(., rs_cl)
    rs_ng <- raster(paste(pr_odir, "/max_similarity_neg.tif", sep="")) %>%
                raster::crop(., rs_cl) %>%
                raster::resample(., rs_cl, method='bilinear') %>%
                raster::mask(., rs_cl)
    
    #correction with negative outcome similarity
    #options [a] p/(1+n); [b] p/exp(n); [c] (p-n); [d] p/n; [e] p*(1-n)
    if (!file.exists(paste(pr_odir, "/f1_similarity.tif", sep=""))) {
        rs_f1 <- rs_ps / (1 + rs_ng)
        rs_f1 <- (rs_f1 - min(rs_f1[], na.rm=T)) / (max(rs_f1[], na.rm=T) - min(rs_f1[], na.rm=T))
        rs_f1 <- writeRaster(rs_f1, paste(pr_odir, "/f1_similarity.tif", sep=""))
    } else {
        rs_f1 <- raster(paste(pr_odir, "/f1_similarity.tif", sep=""))
    }
    
    if (!file.exists(paste(pr_odir, "/f2_similarity.tif", sep=""))) {
        rs_f2 <- rs_ps / exp(rs_ng)
        rs_f2 <- (rs_f2 - min(rs_f2[], na.rm=T)) / (max(rs_f2[], na.rm=T) - min(rs_f2[], na.rm=T))
        rs_f2 <- writeRaster(rs_f2, paste(pr_odir, "/f2_similarity.tif", sep=""))
    } else {
        rs_f2 <- raster(paste(pr_odir, "/f2_similarity.tif", sep=""))
    }
    
    if (!file.exists(paste(pr_odir, "/f3_similarity.tif", sep=""))) {
        rs_f3 <- rs_ps - rs_ng
        rs_f3 <- (rs_f3 - min(rs_f3[], na.rm=T)) / (max(rs_f3[], na.rm=T) - min(rs_f3[], na.rm=T))
        rs_f3 <- writeRaster(rs_f3, paste(pr_odir, "/f3_similarity.tif", sep=""))
    } else {
        rs_f3 <- raster(paste(pr_odir, "/f3_similarity.tif", sep=""))
    }
    
    if (!file.exists(paste(pr_odir, "/f4_similarity.tif", sep=""))) {
        rs_f4 <- rs_ps / rs_ng
        rs_f4 <- (rs_f4 - min(rs_f4[], na.rm=T)) / (max(rs_f4[], na.rm=T) - min(rs_f4[], na.rm=T))
        rs_f4 <- writeRaster(rs_f4, paste(pr_odir, "/f4_similarity.tif", sep=""))
    } else {
        rs_f4 <- raster(paste(pr_odir, "/f4_similarity.tif", sep=""))
    }
    
    if (!file.exists(paste(pr_odir, "/f5_similarity.tif", sep=""))) {
        rs_f5 <- rs_ps * (1 - rs_ng)
        rs_f5 <- (rs_f5 - min(rs_f5[], na.rm=T)) / (max(rs_f5[], na.rm=T) - min(rs_f5[], na.rm=T))
        rs_f5 <- writeRaster(rs_f5, paste(pr_odir, "/f5_similarity.tif", sep=""))
    } else {
        rs_f5 <- raster(paste(pr_odir, "/f5_similarity.tif", sep=""))
    }
    
    #extract data for the locations in question
    pdata_cl <- pdata %>%
                    dplyr::mutate(., 
                                  chance_suit=raster::extract(rs_cl, .[,c("Lon","Lat")]),
                                  simp=raster::extract(rs_ps, .[,c("Lon","Lat")]),
                                  simf1=raster::extract(rs_f1, .[,c("Lon","Lat")]),
                                  simf2=raster::extract(rs_f2, .[,c("Lon","Lat")]),
                                  simf3=raster::extract(rs_f3, .[,c("Lon","Lat")]),
                                  simf4=raster::extract(rs_f4, .[,c("Lon","Lat")]),
                                  simf5=raster::extract(rs_f5, .[,c("Lon","Lat")])) %>%
                    tidyr::drop_na(., chance_suit) %>%
                    dplyr::mutate(., classnum=ifelse(class == 'positive', 1, 0))
    
    roc1 <- roc(pdata_cl$classnum, pdata_cl$chance_suit)
    png(paste(pr_odir,"/plot_auc_niches.png",sep=""),width=5,height=5,units='in',res=300,pointsize=12)
    plot(roc1); text(0.8,0.95,paste("AUC=",round(roc1$auc,3), sep=""), cex=1.5)
    dev.off()
    
    roc2 <- roc(pdata_cl$classnum, pdata_cl$simp)
    png(paste(pr_odir,"/plot_auc_similarity.png",sep=""),width=5,height=5,units='in',res=300,pointsize=12)
    plot(roc2); text(0.8,0.95,paste("AUC=",round(roc2$auc,3), sep=""), cex=1.5)
    dev.off()
    
    #extract values into matrix
    rs_vls <- rasterToPoints(rs_cl) %>% 
        as_tibble() %>% 
        setNames(c('x', 'y', 'icraf')) %>%
        dplyr::mutate(., 
                      pos=raster::extract(rs_ps, .[,c('x','y')]),
                      neg=raster::extract(rs_ng, .[,c('x','y')]),
                      f1=raster::extract(rs_f1, .[,c('x','y')]),
                      f2=raster::extract(rs_f2, .[,c('x','y')]),
                      f3=raster::extract(rs_f3, .[,c('x','y')]),
                      f4=raster::extract(rs_f4, .[,c('x','y')]),
                      f5=raster::extract(rs_f5, .[,c('x','y')])) %>%
        tidyr::drop_na(.)
    
    #plot(rs_vls[sample(1:nrow(rs_vls), size=5000), c('pos', 'neg')], pch=20, xlim=c(0,1), ylim=c(0,1))
    #abline(h=0.5, v=0.5, col='red')
    
    #correlation matrix
    cormat <- as.data.frame(cor(rs_vls[,c(3:ncol(rs_vls))], method='pearson'))
    write.csv(cormat, paste(pr_odir, "/cor_matrix.csv", sep=""))
    
    #correlation matrix for pos <= 0.5
    rs_vls_filt <- dplyr::filter(rs_vls, pos <= 0.5) %>%
                    dplyr::select(., -x, -y)
    cormat_filt <- as.data.frame(cor(rs_vls_filt, method='pearson'))
    write.csv(cormat_filt, paste(pr_odir, "/cor_matrix_midlow.csv", sep=""))
    
    #correlation matrix for pos >= 0.5
    rs_vls_filt <- dplyr::filter(rs_vls, pos >= 0.5) %>%
                    dplyr::select(., -x, -y)
    cormat_filt <- as.data.frame(cor(rs_vls_filt, method='pearson'))
    write.csv(cormat_filt, paste(pr_odir, "/cor_matrix_midupp.csv", sep=""))
    
    #correlation matrix for pos >= 0.8
    rs_vls_filt <- dplyr::filter(rs_vls, pos >= 0.8) %>%
                    dplyr::select(., -x, -y)
    cormat_filt <- as.data.frame(cor(rs_vls_filt, method='pearson'))
    write.csv(cormat_filt, paste(pr_odir, "/cor_matrix_upp.csv", sep=""))
}
