### calculate soilcap in mm
# x: awc in fraction
# y: lower limit of horizon, in cm
# rdepth: rooting depth (in cm) for given crop. If no idea use 60cm
# minval, maxval: are the (theoretical) limits for rooting depth (taken from a paper)
soilcap_calc <- function(x, y, rdepth=60, minval, maxval) {
  if (length(x) != length(y)) {stop("length of x and y must be the same")}
  rdepth <- max(c(rdepth,minval)) #cross check
  rdepth <- min(c(rdepth,maxval)) #cross-check
  wc_df <- data.frame(depth=y,wc=x)
  if (!rdepth %in% wc_df$depth) {
    wc_df1 <- wc_df[which(wc_df$depth < rdepth),]
    wc_df2 <- wc_df[which(wc_df$depth > rdepth),]
    y1 <- wc_df1$wc[nrow(wc_df1)]; y2 <- wc_df2$wc[1]
    x1 <- wc_df1$depth[nrow(wc_df1)]; x2 <- wc_df2$depth[1]
    ya <- (rdepth-x1) / (x2-x1) * (y2-y1) + y1
    wc_df <- rbind(wc_df1,data.frame(depth=rdepth,wc=ya),wc_df2)
  }
  wc_df <- wc_df[which(wc_df$depth <= rdepth),]
  wc_df$soilthick <- wc_df$depth - c(0,wc_df$depth[1:(nrow(wc_df)-1)])
  wc_df$soilcap <- wc_df$soilthick * wc_df$wc
  soilcp <- sum(wc_df$soilcap) * 10 #in mm
  return(soilcp)
}

#potential evapotranspiration
peest <- function(srad,tmin,tmax) {
  #constants
  albedo <- 0.2
  vpd_cte <- 0.7
  
  #soil heat flux parameters
  a_eslope=611.2
  b_eslope=17.67
  c_eslope=243.5
  
  #input parameters
  tmean <- (tmin+tmax)/2
  
  #net radiation
  rn = (1-albedo) * srad
  
  #soil heat flux
  eslope=a_eslope*b_eslope*c_eslope/(tmean+c_eslope)^2*exp(b_eslope*tmean/(tmean+c_eslope))
  
  #estimate vpd
  esat_min=0.61120*exp((17.67*tmin)/(tmin+243.5))
  esat_max=0.61120*exp((17.67*tmax)/(tmax+243.5))
  vpd=vpd_cte*(esat_max-esat_min) #kPa
  
  #Priestley-Taylor
  pt_const=1.26
  pt_fact=1
  vpd_ref=1
  psycho=62
  rho_w=997
  rlat_ht=2.26E6
  
  pt_coef=pt_fact*pt_const
  pt_coef = 1 + (pt_coef-1) * vpd / vpd_ref
  
  #*10^6? To convert fluxes MJ to J
  #rlat_ht? Latent heat flux to water flux
  #100/rho_w? Kg/m^2 to cm
  et_max=(pt_coef * rn * eslope/(eslope+psycho) * 10^6 / rlat_ht * 100/rho_w)*10 #in mm
  return(et_max)
}


#the two functions below estimate the ea/ep
#based on Jones (1987)
#ea/ep: actual to potential evapotranspiration ratio
eabyep_calc <- function(soilcp=100,soilsat=100,cropfc=1,avail=50,rain,evap) {
  avail <- min(c(avail,soilcp))
  eratio <- eabyep(soilcp,avail)
  demand <- eratio*cropfc*evap
  result <- avail + rain - demand
  logging <- result - soilcp
  logging <- max(c(logging,0))
  logging <- min(c(soilsat,logging))
  runoff <- result - (logging+soilcp)
  avail <- min(c(soilcp,result))
  avail <- max(c(avail,0))
  runoff <- max(c(runoff,0))
  out <- data.frame(AVAIL=avail,DEMAND=demand,ERATIO=eratio,RAIN=rain,LOGGING=logging,RUNOFF=runoff)
  return(out)
}


#ea/ep function
eabyep <- function(soilcp,avail) {
  percwt <- min(c(100,avail/soilcp*100))
  percwt <- max(c(1,percwt))
  eratio <- min(c(percwt/(97-3.868*sqrt(soilcp)),1))
  return(eratio)
}


#wrapper to calculate the water balance modeling variables
watbal_wrapper <- function(out_all, soilcp, soilsat)  {
  out_all$ETMAX <- out_all$AVAIL <- out_all$ERATIO <- out_all$LOGGING <- out_all$RUNOFF <- out_all$DEMAND <- out_all$CUM_RAIN <- NA
  for (d in 1:nrow(out_all)) {
    out_all$ETMAX[d] <- peest(out_all$SRAD[d],out_all$TMIN[d],out_all$TMAX[d])
    
    if (d==1) {
      out_all$CUM_RAIN[d] <- out_all$RAIN[d]
      sfact <- eabyep_calc(soilcp=soilcp,soilsat=soilsat,cropfc=1,avail=0,rain=out_all$RAIN[d],evap=out_all$ETMAX[d])
      out_all$AVAIL[d] <- sfact$AVAIL
      out_all$ERATIO[d] <- sfact$ERATIO
      out_all$LOGGING[d] <- sfact$LOGGING
      out_all$RUNOFF[d] <- sfact$RUNOFF
      out_all$DEMAND[d] <- sfact$DEMAND
      
    } else {
      out_all$CUM_RAIN[d] <- out_all$CUM_RAIN[d-1] + out_all$RAIN[d]
      sfact <- eabyep_calc(soilcp=soilcp,soilsat=soilsat,cropfc=1,avail=out_all$AVAIL[d-1],rain=out_all$RAIN[d],evap=out_all$ETMAX[d])
      out_all$AVAIL[d] <- sfact$AVAIL
      out_all$ERATIO[d] <- sfact$ERATIO
      out_all$LOGGING[d] <- sfact$LOGGING
      out_all$RUNOFF[d] <- sfact$RUNOFF
      out_all$DEMAND[d] <- sfact$DEMAND
    }
  }
  return(out_all)
}
