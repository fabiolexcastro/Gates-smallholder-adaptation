#JRV, April 2020

library(lubridate)

#working directory
wd <- "~/Repositories/Gates-smallholder-adaptation/watbal"

#source watbal functions
source(paste(wd,"/wgen_srad.R",sep=""))

#load example weather data
wth <- read.table(paste(wd,"/UFON9911.WTH",sep=""), skip=5, colClasses=c('character',rep('numeric',4)))
names(wth) <- c("DATE","SRAD","TMAX","TMIN","RAIN")

#daily_hist: data.frame with variables date, year, prate, srad, and month (from AgMERRA)
#tdates: data.frame with daily data with the variables "date","prec","tasmin","tasmax"
#wgen_srad(daily_hist,dates,year,lon,lat)

#getting the right type of matrix for 'daily_hist'
#NOTE: this object is all daily data from AgMERRA for the site
x <- wth[,c("DATE","SRAD","RAIN")]
x$YEAR <- as.numeric(substr(x$DATE,1,2))
x$dumm <- as.numeric(x$YEAR) %% 100 > 1950%%100
x$YEAR[which(x$dumm)] <- x$YEAR[which(x$dumm)] + 1900
x$YEAR[which(!x$dumm)] <- x$YEAR[which(!x$dumm)] + 2000
x$dumm <- NULL
x$DOY <- as.numeric(substr(x$DATE,3,5))
x$month <- month(as.Date((x$DOY-1), origin = paste(x$YEAR,"-01-01",sep="")))
x$date1 <- paste(as.Date((x$DOY-1), origin = paste(x$YEAR,"-01-01",sep="")))
x <- x[,c("date1","RAIN","SRAD","month")]
names(x) <- c("date","prate","srad","month")

#this object is all daily data from CHIRPS and CHIRTS
#getting the right type of matrix for 'tdates'
y <- wth[,c("DATE","RAIN","TMIN","TMAX")]
y$YEAR <- as.numeric(substr(y$DATE,1,2))
y$dumm <- as.numeric(y$YEAR) %% 100 > 1950%%100
y$YEAR[which(y$dumm)] <- y$YEAR[which(y$dumm)] + 1900
y$YEAR[which(!y$dumm)] <- y$YEAR[which(!y$dumm)] + 2000
y$dumm <- NULL
y$DOY <- as.numeric(substr(y$DATE,3,5))
y$date1 <- paste(as.Date((y$DOY-1), origin = paste(y$YEAR,"-01-01",sep="")))
y <- y[,c("date1","RAIN","TMIN","TMAX","YEAR")]
names(y) <- c("date","prec","tasmin","tasmax","year")

#now calculate radiation, per year
#lat=27.398; lon=-81.940
srad_all <- data.frame()
for (yr in unique(y$year)) {
  #yr <- unique(y$year)[1]
  cat(yr,"\n")
  yi <- y[which(y$year == yr),]
  yi$year <- NULL
  srad <- wgen_srad(x,yi,yr,lon=-81.940,lat=27.398)
  srad_all <- rbind(srad_all, srad)
}

#now compare the estimated and the original one
par(mar=c(5,5,1,1))
plot(srad_all$srad[1:365], ty='l', ylab="Solar radiation [MJ / day]", xlab="Day of year")
lines(wth$SRAD[1:365], col='red')
legend(x=250,y=28,lty=c(1,1),cex=0.75,col=c("black","red"),legend=c("estimated","original"))


