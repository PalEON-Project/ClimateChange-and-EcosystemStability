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


pdf(file.path(fig.out, "MetVars_Regional_vs_Site_Smoother.pdf"))
for(v in met.vars){
print(
  ggplot(data=met.region[met.region$MetVar==v,]) + 
#   facet_wrap(~Site, scales="fixed") +
  #   geom_ribbon(aes(x=Year, ymin=min, ymax=max), alpha=0.5) +
    geom_line(data=met.sites2[met.sites2$MetVar==v,], aes(x=Year, y=value, color=Site), alpha=0.5) +
  geom_line(aes(x=Year, y=mean), alpha=0.5) +
  stat_smooth(aes(x=Year, y=mean), color="black", fill="black", size=2) +
  stat_smooth(data=met.sites2[met.sites2$MetVar==v,], aes(x=Year, y=value, color=Site, fill=Site), size=2) +
  ggtitle(v)
  )
}
dev.off()

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
gam.tair    <- gam(   tair ~ s(Year,by=Site) + Site, data=met.sites)
gam.precipf <- gam(precipf ~ s(Year,by=Site) + Site, data=met.sites)
gam.swdown  <- gam( swdown ~ s(Year,by=Site) + Site, data=met.sites)
gam.lwdown  <- gam( lwdown ~ s(Year,by=Site) + Site, data=met.sites)
gam.psurf   <- gam(  psurf ~ s(Year,by=Site) + Site, data=met.sites)
gam.qair    <- gam(   qair ~ s(Year,by=Site) + Site, data=met.sites)
gam.wind    <- gam(   wind ~ s(Year,by=Site) + Site, data=met.sites)


source("~/Desktop/Research/R_Functions/Calculate_GAMM_Posteriors.R")
summary(met.sites)
tair.predict    <- post.distns(model.gam=gam.tair   , model.name="tair"   , n=100, newdata=met.sites, vars="Year", terms=T)
precipf.predict <- post.distns(model.gam=gam.precipf, model.name="precipf", n=100, newdata=met.sites, vars="Year", terms=T)
swdown.predict  <- post.distns(model.gam=gam.swdown , model.name="swdown" , n=100, newdata=met.sites, vars="Year", terms=T)
lwdown.predict  <- post.distns(model.gam=gam.lwdown , model.name="lwdown" , n=100, newdata=met.sites, vars="Year", terms=T)
psurf.predict   <- post.distns(model.gam=gam.psurf  , model.name="psurf"  , n=100, newdata=met.sites, vars="Year", terms=T)
qair.predict    <- post.distns(model.gam=gam.qair   , model.name="qair"   , n=100, newdata=met.sites, vars="Year", terms=T)
wind.predict    <- post.distns(model.gam=gam.wind   , model.name="wind"   , n=100, newdata=met.sites, vars="Year", terms=T)

gams.out <- rbind(tair.predict, precipf.predict, swdown.predict, lwdown.predict, psurf.predict, qair.predict, wind.predict)
gams.out$Year <- gams.out$x
summary(gams.out)

# # Relativizing things to be able to compare the relative shifts at the regional scale
# for(v in unique(gams.out$Model)){
#   var.mean <- mean(met.sites[met.sites$Site=="Region",v], na.rm=T)
#   gams.out[gams.out$Model==v, "mean.rel"] <- gams.out[gams.out$Model==v, "mean"]/var.mean
#   gams.out[gams.out$Model==v, "lwr.rel"] <- gams.out[gams.out$Model==v, "lwr"]/var.mean
#   gams.out[gams.out$Model==v, "upr.rel"] <- gams.out[gams.out$Model==v, "upr"]/var.mean
# }
# summary(gams.out)

colors.sites <- c("darkorchid4", "dodgerblue4", "forestgreen", "goldenrod3", "darkorange2", "firebrick2", "black")
colors.met <- c("red3", "blue3", "goldenrod2", "darkorange2", "darkolivegreen3", "cadetblue1", "gray50")

pdf(file.path(fig.out, "MetVars_TemporalTrends.pdf"), height=11, width=8.5)
{
print(
ggplot(data=gams.out) + 
  facet_grid(Model~., scales="free_y") +
  geom_ribbon(aes(x=Year, ymin=lwr, ymax=upr, fill=Site), alpha=0.5) +
  geom_line(aes(x=Year, y=mean, color=Site), size=1) +
  scale_x_continuous(expand=c(0,0), name="Year") +
  scale_y_continuous(expand=c(0,0), name="Temporal Trend") +
  scale_color_manual(values=colors.sites) +
  scale_fill_manual(values=colors.sites) + 
  ggtitle("Temporal Trends by Site") +
  theme_bw()
)
print(
ggplot(data=gams.out[gams.out$Site=="Region",]) +
  facet_grid(Model~., scales="free_y") +
  geom_ribbon(aes(x=Year, ymin=lwr, ymax=upr, fill=Model), alpha=0.5) +
  geom_line(aes(x=Year, y=mean, color=Model), size=2) +
  scale_x_continuous(expand=c(0,0), name="Year") +
  scale_y_continuous(expand=c(0,0), name="Temporal Trend") +
  scale_fill_manual(values=colors.met) +
  scale_color_manual(values=colors.met)
)
}
dev.off()

png(file.path(fig.out, "MetVars_TemporalTrends.png"), height=11, width=8.5, units="in", res=180)
{
  ggplot(data=gams.out) + 
      facet_grid(Model~., scales="free_y") +
      geom_ribbon(aes(x=Year, ymin=lwr, ymax=upr, fill=Site), alpha=0.5) +
      geom_line(aes(x=Year, y=mean, color=Site), size=1) +
      scale_x_continuous(expand=c(0,0), name="Year") +
      scale_y_continuous(expand=c(0,0), name="Temporal Trend") +
      scale_color_manual(values=colors.sites) +
      scale_fill_manual(values=colors.sites) + 
      ggtitle("Temporal Trends by Site") +
      theme_bw()
}
dev.off()




summary(met.gam.region)
par(mfrow=c(3,3)); plot(met.gam.region)
# -------------------------------------------
