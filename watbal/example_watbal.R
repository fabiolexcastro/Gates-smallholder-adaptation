#JRV, Mar 2020

#load libraries
library(GSIF)

#working directory
wd <- "~/Documents/Gates-Ag-Adapt-Atlas/watbal"

#source watbal functions
source(paste(wd,"/watbal.R",sep=""))

#load example weather file
#srad in MJ/day; tmax, tmin in Celsius; and rain in mm
wth <- read.table(paste(wd,"/UFON9911.WTH",sep=""), skip=5, colClasses=c('character',rep('numeric',4)))
names(wth) <- c("DATE","SRAD","TMAX","TMIN","RAIN")

# get example soil
data(afsp) #load data
sel <- afsp$horizons$SOURCEID=="NG 28440_Z5" #define soil profile of interest
hor <- afsp$horizons[sel,] #select this soil profile
BLDf <- ifelse(is.na(hor$BLD), mean(hor$BLD, na.rm=TRUE), hor$BLD) # replace missing values by average BLD

#note i use h3=-33 because at field capacity the pressure is -33 kPa. Hence for ASW we will use AWCh3
#note BLDf (or BLD) needs to be multiplied by 1000
hor <- cbind(hor, AWCPTF(hor$SNDPPT, hor$SLTPPT, 
                         hor$CLYPPT, hor$ORCDRC, BLD=BLDf*1000, hor$CEC, 
                         hor$PHIHOX, h1=-10, h2=-20, h3=-33))

#now calculate the ASW in mm for each soil horizon
hor$tetaFC <- hor$WWP + hor$AWCh3 #volumetric water content at field capacity (fraction)
hor$AWSat <- hor$tetaS - hor$tetaFC

#soilcap for a 60cm rooting depth
#minval and maxval are in cm, asw values per depth are in %
scp <- soilcap_calc(x=hor$AWCh3, y=hor$LHDICM, rdepth=60, minval=45, maxval=100)
ssat <- soilcap_calc(x=hor$AWSat, y=hor$LHDICM, rdepth=60, minval=45, maxval=100)

#run watbal
wtbal <- watbal_wrapper(wth, scp, ssat)

#plots to understand water balance (no need to run this part every time)
plot(wtbal$ERATIO[1:365],ty='l',ylim=c(0,1.2)) #dimensionless (fraction)
lines(wtbal$RAIN[1:365]*0.01,col='blue') #in mm, but multiplied by 0.01 to put in the same plot
lines(wtbal$AVAIL[1:365]*0.01,col='green') #in mm, but multiplied by 0.01 to put in the same plot
lines(wtbal$LOGGING[1:365]*0.01,col='orange') #in mm, but multiplied by 0.01 to put in the same plot
lines(wtbal$RUNOFF[1:365]*0.01,col='red') #in mm, but multiplied by 0.01 to put in the same plot
abline(h=scp*0.01) #in mm, but multiplied by 0.01 to put in the same plot
abline(h=ssat*0.01) #in mm, but multiplied by 0.01 to put in the same plot

#split year and doy from field DATE (only needed because i used a DSSAT file)
wtbal$YEAR <- as.numeric(substr(wtbal$DATE,1,2))
wtbal$dumm <- as.numeric(wtbal$YEAR) %% 100 > 1950%%100
wtbal$YEAR[which(wtbal$dumm)] <- wtbal$YEAR[which(wtbal$dumm)] + 1900
wtbal$YEAR[which(!wtbal$dumm)] <- wtbal$YEAR[which(!wtbal$dumm)] + 2000
wtbal$dumm <- NULL
wtbal$DOY <- as.numeric(substr(wtbal$DATE,3,5))

#calculate number of waterlogging days for year 1999
nd_wlog1 <- nrow(wtbal[which(wtbal$YEAR == 1999 & wtbal$LOGGING > 0),])

#calculate number of days with more than 50% waterlogging (water content at 50% of saturation)
nd_wlog2 <- nrow(wtbal[which(wtbal$YEAR == 1999 & wtbal$LOGGING > (ssat*0.5)),])

#calculate number of days with more than 90% waterlogging (water content at saturation)
nd_wlog3 <- nrow(wtbal[which(wtbal$YEAR == 1999 & wtbal$LOGGING >= ssat),])

#calculate max. number of continuous logging days for year 1999
wtbal_yr <- wtbal[which(wtbal$YEAR == 1999),]
ndays <- 0; maxndays <- c()
for (i in 1:nrow(wtbal_yr)) {
  if (wtbal_yr$LOGGING[i] > 0) {ndays <- ndays+1} else {maxndays <- c(maxndays,ndays); ndays <- 0}
}
maxndays <- max(maxndays)


