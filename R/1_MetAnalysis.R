# -------------------------------------------
# Assessing the scales and magnitude of climate change variability over the past millennium in the MIP drivers

# 1. Stastical detection of significant change
#    - Use loess or TP regression splines
# 2. Comparison with paleoclimate reconstructions
# -------------------------------------------
rm(list=ls())


# -------------------------------------------
# Load libraries; set file paths
# -------------------------------------------
library(ncdf4)
library(ggplot2)
library(mgcv)
setwd("~/Desktop/Research/PalEON_CR/PalEON_MIP_Site/Analyses/Change-and-Stability")
path.RFunc <- "~/Desktop/Research/R_Functions/"

raw.dir <- "raw_data"
dat.out <- "Data/Met"
fig.out <- "Figures/Met"

if(!dir.exists(dat.out)) dir.create(dat.out)
if(!dir.exists(fig.out)) dir.create(fig.out)
# -------------------------------------------

# -------------------------------------------
# Define some useful variables
# -------------------------------------------
met.vars <- c("tair", "precipf", "swdown", "lwdown", "psurf", "qair", "wind")
hips <- data.frame(Site = c( "PHA",  "PHO",  "PUN",  "PBL",  "PDL",  "PMB"),
                   lat  = c( 42.54,  45.25,  46.22,  46.28,  47.17,  43.61),
                   lon  = c(-72.18, -68.73, -89.53, -94.58, -95.17, -82.83),
                   lat2 = c( 42.75,  45.25,  46.25,  46.25,  47.25,  43.75),
                   lon2 = c(-72.25, -68.75, -89.75, -94.75, -95.25, -82.75))
years <- 850:2010
# -------------------------------------------

# -------------------------------------------
# Covert raw monthly drivers to annual
# Base dataset = phase2_met_regional_v2_monthly (available on iPlant)
# -------------------------------------------
{
# years <- 850:2010
# months.years <- vector()
# for(y in years){
# 	months.years <- c(months.years, rep(y, 12))
# }
# length(months.years)
# months.years[1:25]
# 
# 
# 
# for(v in met.vars){
# 	nc.v <- nc_open(file.path(raw.dir, "phase2_met_regional_v2_monthly", paste0(v, ".nc")))
# 	dat.v <- ncvar_get(nc.v, v)
# 	lat   <- ncvar_get(nc.v, "lat" )
# 	lon   <- ncvar_get(nc.v, "lon" )
# 	nc_close(nc.v)
# 	dim(dat.v)
# 	
# 	# make a new array fo rthe year
# 	dat.new <- array(dim=c(dim(dat.v)[1], dim(dat.v)[2], length(unique(months.years))))
# 	
# 	# get the mean for each year
# 	for(y in 1:length(years)){
# 		dat.new[,,y] <- apply(dat.v[,,(y*12-11):(y*12)], c(1,2), FUN=mean)
# 	}
# 
# 	# Save it as a new .nc file
# 	nc_time_units <- paste('Year', sep='')
# 	dim.time      <- ncdim_def("time",nc_time_units,years,unlim=TRUE)
# 	
# 	# Print correct units 
# 	if (v == 'lwdown') {
# 	  nc_variable_long_name=paste('Incident (downwelling) longwave ',
# 	                              'radiation averaged over the time step of the forcing data', sep='')
# 	  nc_variable_units='W m-2'
# 	} else if (v == 'precipf') {
# 	  nc_variable_long_name=paste('The per unit area and time ',
# 	                              'precipitation representing the sum of convective rainfall, ',
# 	                              'stratiform rainfall, and snowfall', sep='')
# 	  nc_variable_units='kg m-2 s-1'
# 	} else if (v == 'psurf') {
# 	  nc_variable_long_name='Pressure at the surface'
# 	  nc_variable_units='Pa'
# 	} else if (v == 'qair') {
# 	  nc_variable_long_name=
# 	    'Specific humidity measured at the lowest level of the atmosphere'
# 	  nc_variable_units='kg kg-1'
# 	} else if (v == 'swdown') {
# 	  nc_variable_long_name=paste('Incident (downwelling) radiation in ',
# 	                              'the shortwave part of the spectrum averaged over the time ',
# 	                              'step of the forcing data', sep='')
# 	  nc_variable_units='W m-2'
# 	} else if (v == 'tair') {
# 	  nc_variable_long_name='2 meter air temperature'
# 	  nc_variable_units='K'
# 	} else if (v == 'wind') {
# 	  nc_variable_long_name='Wind speed'
# 	  nc_variable_units='m s-1'
# 	}
# 	
# 	# Make a few dimensions we can use
# 	dimY <- ncdim_def( "lat", "longitude: degrees", lat )
# 	dimX <- ncdim_def( "lon", "latitude: degrees", lon )
# 	
# 	nc_var  <- ncvar_def(v,nc_variable_units,
# 	                     list(dimX,dimY,dim.time), missval=NA, longname=nc_variable_long_name,prec="double")
# 	
# 	ofname  <- paste(file.path(dat.out,paste0(v,'.nc')))
# 	newfile <- nc_create( ofname, nc_var ) # Initialize file 
# 	
# 	ncatt_put( newfile, nc_var, 'days since 850', years)
# 	ncatt_put( newfile, 0, 'description',"PalEON formatted Phase 1 met driver")
# 	
# 	ncvar_put(newfile, nc_var, dat.new) # Write netCDF file
# 	
# 	nc_close(newfile)  
# 
# }

}
# -------------------------------------------

# -------------------------------------------
# Graph the range and HIPS site met
# -------------------------------------------
met.region <- merge(years, met.vars)
names(met.region) <- c("Year", "MetVar")
summary(met.region)

met.sites <- data.frame(merge(hips, years))
names(met.sites) <- c(names(hips), "Year")

summary(met.region)
summary(met.sites)

# Adding lat/lon to the hips data frame
for(v in met.vars){
  nc.v <- nc_open(file.path(dat.out, paste0(v, ".nc")))
  met.region[met.region$MetVar==v, paste0("mean")] <- apply(ncvar_get(nc.v, v), 3, mean, na.rm=T)
  met.region[met.region$MetVar==v, paste0( "min")] <- apply(ncvar_get(nc.v, v), 3,  min, na.rm=T)
  met.region[met.region$MetVar==v, paste0( "max")] <- apply(ncvar_get(nc.v, v), 3,  max, na.rm=T)
  
  for(s in unique(met.sites$Site)){
    x.use <- which(ncvar_get(nc.v, "lat")==mean(met.sites[met.sites$Site==s, "lat2"]))
    y.use <- which(ncvar_get(nc.v, "lon")==mean(met.sites[met.sites$Site==s, "lon2"]))
    met.sites[met.sites$Site==s,v] <- ncvar_get(nc.v, v)[y.use, x.use, ]
  }
  nc_close(nc.v)
}

summary(met.sites)
summary(met.region)

met.sites2 <- stack(met.sites[,met.vars])
names(met.sites2) <- c("value", "MetVar")
met.sites2[,c("Site", "Year")] <- met.sites[,c("Site", "Year")]
summary(met.sites2)


# pdf(file.path(fig.out, "MetVars_Regional_vs_Site_Smoother.pdf"))
# for(v in met.vars){
# print(
#   ggplot(data=met.region[met.region$MetVar==v,]) + 
# #   facet_wrap(~Site, scales="fixed") +
#   #   geom_ribbon(aes(x=Year, ymin=min, ymax=max), alpha=0.5) +
#     geom_line(data=met.sites2[met.sites2$MetVar==v,], aes(x=Year, y=value, color=Site), alpha=0.5) +
#   geom_line(aes(x=Year, y=mean), alpha=0.5) +
#   stat_smooth(aes(x=Year, y=mean), color="black", fill="black", size=2) +
#   stat_smooth(data=met.sites2[met.sites2$MetVar==v,], aes(x=Year, y=value, color=Site, fill=Site), size=2) +
#   ggtitle(v)
#   )
# }
# dev.off()

# Met Splines (temporal trends) for sites & region
summary(met.region)
summary(met.sites)

# reshaping met.region a bit to make it easy to merge with the sites
# Note: Here we're doing the regression on the mean temp, not the mean shape of the deviations
library(reshape2)
met.region$Year <- as.factor(met.region$Year)
met.region2 <- recast(met.region[,c("Year", "MetVar", "mean")], Year ~ MetVar)
met.region2$Site <- as.factor("Region")

# Put year back as numeric in both data frames!
met.region$Year <- as.numeric(paste(met.region$Year))
met.region2$Year <- as.numeric(paste(met.region2$Year))
summary(met.region)
summary(met.region2)

met.sites <- merge(met.sites, met.region2, all.x=T, all.y=T)
summary(met.sites)

# Fit some simple gams with temporal thin-plate regression splines on year so we can compare sites and the regional mean
gam.tair    <- gam(   tair ~ s(Year,by=Site, k=46) + Site, data=met.sites)
gam.precipf <- gam(precipf ~ s(Year,by=Site, k=46) + Site, data=met.sites)
gam.swdown  <- gam( swdown ~ s(Year,by=Site, k=46) + Site, data=met.sites)
gam.lwdown  <- gam( lwdown ~ s(Year,by=Site, k=46) + Site, data=met.sites)
gam.psurf   <- gam(  psurf ~ s(Year,by=Site, k=46) + Site, data=met.sites)
gam.qair    <- gam(   qair ~ s(Year,by=Site, k=46) + Site, data=met.sites)
gam.wind    <- gam(   wind ~ s(Year,by=Site, k=46) + Site, data=met.sites)

gam.tair1    <- gam(   tair ~ s(Year,by=Site, k=40) + Site, data=met.sites[met.sites$Year<1850,])
gam.precipf1 <- gam(precipf ~ s(Year,by=Site, k=40) + Site, data=met.sites[met.sites$Year<1850,])
gam.swdown1  <- gam( swdown ~ s(Year,by=Site, k=40) + Site, data=met.sites[met.sites$Year<1850,])
gam.lwdown1  <- gam( lwdown ~ s(Year,by=Site, k=40) + Site, data=met.sites[met.sites$Year<1850,])
gam.psurf1   <- gam(  psurf ~ s(Year,by=Site, k=40) + Site, data=met.sites[met.sites$Year<1850,])
gam.qair1    <- gam(   qair ~ s(Year,by=Site, k=40) + Site, data=met.sites[met.sites$Year<1850,])
gam.wind1    <- gam(   wind ~ s(Year,by=Site, k=40) + Site, data=met.sites[met.sites$Year<1850,])

gam.tair2    <- gam(   tair ~ s(Year,by=Site, k=6) + Site, data=met.sites[met.sites$Year>1900,])
gam.precipf2 <- gam(precipf ~ s(Year,by=Site, k=6) + Site, data=met.sites[met.sites$Year>1900,])
gam.swdown2  <- gam( swdown ~ s(Year,by=Site, k=6) + Site, data=met.sites[met.sites$Year>1900,])
gam.lwdown2  <- gam( lwdown ~ s(Year,by=Site, k=6) + Site, data=met.sites[met.sites$Year>1900,])
gam.psurf2   <- gam(  psurf ~ s(Year,by=Site, k=6) + Site, data=met.sites[met.sites$Year>1900,])
gam.qair2    <- gam(   qair ~ s(Year,by=Site, k=6) + Site, data=met.sites[met.sites$Year>1900,])
gam.wind2    <- gam(   wind ~ s(Year,by=Site, k=6) + Site, data=met.sites[met.sites$Year>1900,])

source(file.path(path.RFunc, "Calculate_GAMM_Posteriors.R"))
summary(met.sites)
tair.predict    <- post.distns(model.gam=gam.tair   , model.name="tair"   , n=100, newdata=met.sites, vars="Year", terms=F)
precipf.predict <- post.distns(model.gam=gam.precipf, model.name="precipf", n=100, newdata=met.sites, vars="Year", terms=F)
swdown.predict  <- post.distns(model.gam=gam.swdown , model.name="swdown" , n=100, newdata=met.sites, vars="Year", terms=F)
lwdown.predict  <- post.distns(model.gam=gam.lwdown , model.name="lwdown" , n=100, newdata=met.sites, vars="Year", terms=F)
psurf.predict   <- post.distns(model.gam=gam.psurf  , model.name="psurf"  , n=100, newdata=met.sites, vars="Year", terms=F)
qair.predict    <- post.distns(model.gam=gam.qair   , model.name="qair"   , n=100, newdata=met.sites, vars="Year", terms=F)
wind.predict    <- post.distns(model.gam=gam.wind   , model.name="wind"   , n=100, newdata=met.sites, vars="Year", terms=F)

tair.predict1    <- post.distns(model.gam=gam.tair1   , model.name="tair"   , n=100, newdata=met.sites, vars="Year", terms=F)
precipf.predict1 <- post.distns(model.gam=gam.precipf1, model.name="precipf", n=100, newdata=met.sites, vars="Year", terms=F)
swdown.predict1  <- post.distns(model.gam=gam.swdown1 , model.name="swdown" , n=100, newdata=met.sites, vars="Year", terms=F)
lwdown.predict1  <- post.distns(model.gam=gam.lwdown1 , model.name="lwdown" , n=100, newdata=met.sites, vars="Year", terms=F)
psurf.predict1   <- post.distns(model.gam=gam.psurf1  , model.name="psurf"  , n=100, newdata=met.sites, vars="Year", terms=F)
qair.predict1    <- post.distns(model.gam=gam.qair1   , model.name="qair"   , n=100, newdata=met.sites, vars="Year", terms=F)
wind.predict1    <- post.distns(model.gam=gam.wind1   , model.name="wind"   , n=100, newdata=met.sites, vars="Year", terms=F)

tair.predict2    <- post.distns(model.gam=gam.tair2   , model.name="tair"   , n=100, newdata=met.sites, vars="Year", terms=F)
precipf.predict2 <- post.distns(model.gam=gam.precipf2, model.name="precipf", n=100, newdata=met.sites, vars="Year", terms=F)
swdown.predict2  <- post.distns(model.gam=gam.swdown2 , model.name="swdown" , n=100, newdata=met.sites, vars="Year", terms=F)
lwdown.predict2  <- post.distns(model.gam=gam.lwdown2 , model.name="lwdown" , n=100, newdata=met.sites, vars="Year", terms=F)
psurf.predict2   <- post.distns(model.gam=gam.psurf2  , model.name="psurf"  , n=100, newdata=met.sites, vars="Year", terms=F)
qair.predict2    <- post.distns(model.gam=gam.qair2   , model.name="qair"   , n=100, newdata=met.sites, vars="Year", terms=F)
wind.predict2    <- post.distns(model.gam=gam.wind2   , model.name="wind"   , n=100, newdata=met.sites, vars="Year", terms=F)

gams.out <- rbind(tair.predict, precipf.predict, swdown.predict, lwdown.predict, psurf.predict, qair.predict, wind.predict)
# gams.out$Year <- gams.out$x

gams.out1 <- rbind(tair.predict1, precipf.predict1, swdown.predict1, lwdown.predict1, psurf.predict1, qair.predict1, wind.predict1)
# gams.out1$Year <- gams.out1$x

gams.out2 <- rbind(tair.predict2, precipf.predict2, swdown.predict2, lwdown.predict2, psurf.predict2, qair.predict2, wind.predict2)
# gams.out2$Year <- gams.out2$x

summary(gams.out)

# # Relativizing things to be able to compare the relative shifts at the regional scale
# for(v in unique(gams.out$Model)){
#   var.mean <- mean(met.sites[met.sites$Site=="Region",v], na.rm=T)
#   gams.out[gams.out$Model==v, "mean.rel"] <- gams.out[gams.out$Model==v, "mean"]/var.mean
#   gams.out[gams.out$Model==v, "lwr.rel"] <- gams.out[gams.out$Model==v, "lwr"]/var.mean
#   gams.out[gams.out$Model==v, "upr.rel"] <- gams.out[gams.out$Model==v, "upr"]/var.mean
# }
# summary(gams.out)




# Calculating the derivative 
source(file.path(path.RFunc, "Calculate_GAMM_Derivs.R"))
deriv.tair    <- calc.derivs(gam.tair   , newdata=met.sites, vars="Year")
deriv.precipf <- calc.derivs(gam.precipf, newdata=met.sites, vars="Year")
deriv.swdown  <- calc.derivs(gam.swdown , newdata=met.sites, vars="Year")
deriv.lwdown  <- calc.derivs(gam.lwdown , newdata=met.sites, vars="Year")
deriv.qair    <- calc.derivs(gam.qair   , newdata=met.sites, vars="Year")
deriv.psurf   <- calc.derivs(gam.psurf  , newdata=met.sites, vars="Year")
deriv.wind    <- calc.derivs(gam.wind   , newdata=met.sites, vars="Year")

names(deriv.tair   )[1] <- "x"
names(deriv.precipf)[1] <- "x"
names(deriv.swdown )[1] <- "x"
names(deriv.lwdown )[1] <- "x"
names(deriv.qair   )[1] <- "x"
names(deriv.psurf  )[1] <- "x"
names(deriv.wind   )[1] <- "x"

deriv.tair   $var <- "tair"
deriv.precipf$var <- "precipf"
deriv.swdown $var <- "swdown"
deriv.lwdown $var <- "lwdown"
deriv.qair   $var <- "qair"
deriv.psurf  $var <- "psurf"
deriv.wind   $var <- "wind"

# Looking at the 850-1850 window only
deriv.tair1    <- calc.derivs(gam.tair1   , newdata=met.sites[met.sites$Year<1850,], vars="Year")
deriv.precipf1 <- calc.derivs(gam.precipf1, newdata=met.sites[met.sites$Year<1850,], vars="Year")
deriv.swdown1  <- calc.derivs(gam.swdown1 , newdata=met.sites[met.sites$Year<1850,], vars="Year")
deriv.lwdown1  <- calc.derivs(gam.lwdown1 , newdata=met.sites[met.sites$Year<1850,], vars="Year")
deriv.qair1    <- calc.derivs(gam.qair1   , newdata=met.sites[met.sites$Year<1850,], vars="Year")
deriv.psurf1   <- calc.derivs(gam.psurf1  , newdata=met.sites[met.sites$Year<1850,], vars="Year")
deriv.wind1    <- calc.derivs(gam.wind1   , newdata=met.sites[met.sites$Year<1850,], vars="Year")

names(deriv.tair1   )[1] <- "x"
names(deriv.precipf1)[1] <- "x"
names(deriv.swdown1 )[1] <- "x"
names(deriv.lwdown1 )[1] <- "x"
names(deriv.qair1   )[1] <- "x"
names(deriv.psurf1  )[1] <- "x"
names(deriv.wind1   )[1] <- "x"

deriv.tair1   $var <- "tair"
deriv.precipf1$var <- "precipf"
deriv.swdown1 $var <- "swdown"
deriv.lwdown1 $var <- "lwdown"
deriv.qair1   $var <- "qair"
deriv.psurf1  $var <- "psurf"
deriv.wind1   $var <- "wind"

# Looking at the post-1900 window only
deriv.tair2    <- calc.derivs(gam.tair2   , newdata=met.sites[met.sites$Year>1900,], vars="Year")
deriv.precipf2 <- calc.derivs(gam.precipf2, newdata=met.sites[met.sites$Year>1900,], vars="Year")
deriv.swdown2  <- calc.derivs(gam.swdown2 , newdata=met.sites[met.sites$Year>1900,], vars="Year")
deriv.lwdown2  <- calc.derivs(gam.lwdown2 , newdata=met.sites[met.sites$Year>1900,], vars="Year")
deriv.qair2    <- calc.derivs(gam.qair2   , newdata=met.sites[met.sites$Year>1900,], vars="Year")
deriv.psurf2   <- calc.derivs(gam.psurf2  , newdata=met.sites[met.sites$Year>1900,], vars="Year")
deriv.wind2    <- calc.derivs(gam.wind2   , newdata=met.sites[met.sites$Year>1900,], vars="Year")

names(deriv.tair2   )[1] <- "x"
names(deriv.precipf2)[1] <- "x"
names(deriv.swdown2 )[1] <- "x"
names(deriv.lwdown2 )[1] <- "x"
names(deriv.qair2   )[1] <- "x"
names(deriv.psurf2  )[1] <- "x"
names(deriv.wind2   )[1] <- "x"

deriv.tair2   $var <- "tair"
deriv.precipf2$var <- "precipf"
deriv.swdown2 $var <- "swdown"
deriv.lwdown2 $var <- "lwdown"
deriv.qair2   $var <- "qair"
deriv.psurf2  $var <- "psurf"
deriv.wind2   $var <- "wind"


derivs <- rbind(deriv.tair, deriv.precipf, deriv.swdown, deriv.lwdown, deriv.qair, deriv.psurf, deriv.wind)
derivs$var <- as.factor(derivs$var)
derivs$mean.sig <- ifelse(derivs$sig=="*", derivs$mean, NA)
summary(derivs)

derivs1 <- rbind(deriv.tair1, deriv.precipf1, deriv.swdown1, deriv.lwdown1, deriv.qair1, deriv.psurf1, deriv.wind1)
derivs1$var <- as.factor(derivs1$var)
derivs1$mean.sig <- ifelse(derivs1$sig=="*", derivs1$mean, NA)
summary(derivs1)

derivs2 <- rbind(deriv.tair2, deriv.precipf2, deriv.swdown2, deriv.lwdown2, deriv.qair2, deriv.psurf2, deriv.wind2)
derivs2$var <- as.factor(derivs2$var)
derivs2$mean.sig <- ifelse(derivs2$sig=="*", derivs2$mean, NA)
summary(derivs2)

# Merging the derivs in with the raw gamm output
gams.out$var <- gams.out$Model
names(derivs)[which(names(derivs) %in% c("mean", "lwr", "upr"))] <- paste0("deriv.", c("mean", "lwr", "upr"))
summary(derivs)
summary(gams.out)

gams.out <- merge(gams.out, derivs, all.x=T, all.y=T)
gams.out$mean.sig <- ifelse(gams.out$sig=="*", gams.out$mean, NA)
summary(gams.out)


gams.out1$var <- gams.out1$Model
names(derivs1)[which(names(derivs1) %in% c("mean", "lwr", "upr"))] <- paste0("deriv.", c("mean", "lwr", "upr"))
summary(derivs1)
summary(gams.out1)

gams.out1 <- merge(gams.out1, derivs1, all.x=T, all.y=T)
gams.out1$mean.sig <- ifelse(gams.out1$sig=="*", gams.out1$mean, NA)
summary(gams.out1)

gams.out2$var <- gams.out2$Model
names(derivs2)[which(names(derivs2) %in% c("mean", "lwr", "upr"))] <- paste0("deriv.", c("mean", "lwr", "upr"))
summary(derivs2)
summary(gams.out2)

gams.out2 <- merge(gams.out2, derivs2, all.x=T, all.y=T)
gams.out2$mean.sig <- ifelse(gams.out2$sig=="*", gams.out2$mean, NA)
summary(gams.out2)

# -------------------------------------------
met.sites2$Site <- droplevels(met.sites2$Site)
met.sites2$Site <- factor(met.sites2$Site, levels=c("PDL", "PBL", "PUN", "PMB", "PHA", "PHO"))
gams.out$Site <- factor(gams.out$Site, levels=c("PDL", "PBL", "PUN", "PMB", "PHA", "PHO", "Region"))
colors.sites <- c("firebrick2", "darkorange2", "goldenrod3", "forestgreen", "dodgerblue4", "darkorchid4", "black")
# colors.sites <- c("darkorchid4", "dodgerblue4", "forestgreen", "goldenrod3", "darkorange2", "firebrick2", "black")
colors.met <- c("red3", "blue3", "goldenrod2", "darkorange2", "darkolivegreen3", "cadetblue1", "gray50")

names(gams.out)[which(names(gams.out)=="var")] <- "MetVar"
summary(gams.out)
summary(met.sites2)

met.sites3 <- merge(met.sites2, gams.out, all.x=T, all.y=T)
dim(met.sites3)

pdf(file.path(fig.out, "MetVars_GAMM_Continuous.pdf"))
for(v in unique(met.sites2$MetVar)){
  print(
    ggplot(data=met.sites3[met.sites3$MetVar==v,]) + 
      geom_line(aes(x=Year, y=value, color=Site), size=0.5, alpha=0.3) +
      geom_ribbon(aes(x=Year, ymin=lwr, ymax=upr, fill=Site), alpha=0.5) +
      geom_line(aes(x=Year, y=mean, color=Site), size=1, alpha=0.3) +
      geom_line(aes(x=Year, y=mean.sig, color=Site), size=1.5, alpha=1) +
      geom_vline(xintercept=1850, linetype="dashed") +
      scale_x_continuous(expand=c(0,0), name="Year") +
      scale_y_continuous(expand=c(0,0), name="Temporal Trend") +
      scale_color_manual(values=colors.sites) +
      scale_fill_manual(values=colors.sites) + 
      ggtitle(v) +
      theme_bw()
  )
}
dev.off()

names(gams.out1)[which(names(gams.out1)=="var")] <- "MetVar"
names(gams.out2)[which(names(gams.out2)=="var")] <- "MetVar"
names(gams.out1)[which(names(gams.out1)=="x")] <- "value"
names(gams.out2)[which(names(gams.out2)=="x")] <- "value"

gams.out1b <- gams.out1[gams.out1$Year<1850,!names(gams.out1) %in% c("Model", "value")]
gams.out2b <- gams.out2[gams.out2$Year>1900,!names(gams.out2) %in% c("Model", "value")]

gams.out.b <- rbind(gams.out1b, gams.out2b)

met.sites4 <- merge(met.sites2, gams.out.b, all.x=T, all.y=T)
met.sites4$mean.sig <- ifelse(met.sites4$sig=="*", met.sites4$mean, NA)
summary(met.sites4)
dim(met.sites4)

pdf(file.path(fig.out, "MetVars_GAMM_Separate.pdf"))
for(v in unique(met.sites2$MetVar)){
  print(
    ggplot(data=met.sites4[met.sites4$MetVar==v,]) + 
      geom_line(aes(x=Year, y=value, color=Site), size=0.5, alpha=0.3) +
      geom_ribbon(data=met.sites4[met.sites4$MetVar==v & met.sites4$Year<1850,], aes(x=Year, ymin=lwr, ymax=upr, fill=Site), alpha=0.5) +
      geom_ribbon(data=met.sites4[met.sites4$MetVar==v & met.sites4$Year>1900,], aes(x=Year, ymin=lwr, ymax=upr, fill=Site), alpha=0.5) +
      geom_line(aes(x=Year, y=mean, color=Site), size=1, alpha=0.3) +
      geom_line(aes(x=Year, y=mean.sig, color=Site), size=1.5, alpha=1) +
      geom_vline(xintercept=1850, linetype="dashed") +
      scale_x_continuous(expand=c(0,0), name="Year") +
      scale_y_continuous(expand=c(0,0), name="Temporal Trend") +
      scale_color_manual(values=colors.sites) +
      scale_fill_manual(values=colors.sites) + 
      ggtitle(v) +
      theme_bw()
  )
}
dev.off()
