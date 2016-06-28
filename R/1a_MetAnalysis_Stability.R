# -------------------------------------------
# Assessing the scales and magnitude of climate change variability over the past millennium in the MIP drivers
# Author: Christy Rollinson, crollinson@gmail.com

# 1. Stastical detection of significant change
#    - Use loess or TP regression splines
# 2. Comparison with paleoclimate reconstructions
# -------------------------------------------
rm(list=ls())


# -------------------------------------------
# Load libraries; set file paths
# -------------------------------------------
library(ncdf4)
library(ggplot2); library(gridExtra)
library(mgcv)
library(plyr); library(parallel)
setwd("~/Dropbox/PalEON_CR/PalEON_MIP_Site/Analyses/Change-and-Stability")
path.gamm.func <- "~/Desktop/R_Functions/"  # Path to github repository of my GAMM helper functions: https://github.com/crollinson/R_Functions.git
mip.utils <- "~/Dropbox/Research/PalEON_CR/MIP_Utils/" # Path to PalEON MIP Utility repository: https://github.com/PalEON-Project/MIP_Utils.git

raw.dir <- "raw_data"
out.dir <- "Data/Met"
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
# 1. Calculating and Mapping met stability
# -------------------------------------------
{
source("R/0_TimeAnalysis.R")
source(file.path(path.gamm.func, "Calculate_GAMM_Derivs.R"))

met.region <- merge(years, met.vars)
names(met.region) <- c("Year", "MetVar")
summary(met.region)

met.sites <- data.frame(merge(hips, years))
names(met.sites) <- c(names(hips), "Year")

summary(met.region)
summary(met.sites)


calc.stability <- function(x){
  dat.tmp <- data.frame(Y=x, Year=1:length(x))
  k.use=round(length(years)/25, 0)
  mod.gam <- gam(Y ~ s(Year, k=k.use), data=dat.tmp)
  mod.deriv <- calc.derivs(mod.gam, newdata=dat.tmp, vars="Year")
  return(mod.deriv$mean)
}

# Adding lat/lon to the hips data frame
for(v in met.vars){
  nc.v <- nc_open(file.path(dat.out, paste0(v, ".nc")))

  # Running The Derivative calculation on each metvar
  years <- 850:2010

  # Extract the regional data
  lon <- ncvar_get(nc.v, "lon")
  lat <- ncvar_get(nc.v, "lat")
  df.met <- ncvar_get(nc.v, v)
  dimnames(df.met) <- list(lon=lon, lat=lat, Year=years)
  
  # Extract the sites; will calc stability separately so we can assess stat. sig.
  for(s in unique(met.sites$Site)){
    y.use <- which(lat==mean(met.sites[met.sites$Site==s, "lat2"]))
    x.use <- which(lon==mean(met.sites[met.sites$Site==s, "lon2"]))
    met.sites[met.sites$Site==s,v] <- ncvar_get(nc.v, v)[x.use, y.use, ]
  }
  
  nc_close(nc.v)

  # ---------------------------
  # Calculating the rate of change for the region 
  #
  # I don't know a better way to get parallel working, so we'll very quickly make 
  # things a list and then us mclapply
  # ---------------------------
  dat.list1 <- list()
  dat.list2 <- list()
  for(x in 1:length(lon)){
    for(y in 1:length(lat)){
      if(is.na(max(df.met[x,y,]))) next
      dat.list1[[paste0("lat", lat[y], "lon", lon[x])]] <- df.met[x,y,which(years<1850)]
      dat.list2[[paste0("lat", lat[y], "lon", lon[x])]] <- df.met[x,y,which(years>1900)]
    }
  }

  deriv.out1 <- mclapply(dat.list1, calc.stability, mc.cores=14)
  deriv.out2 <- mclapply(dat.list2, calc.stability, mc.cores=14)
  
  # Plugging things back into the appropriate array
  df.deriv1 <- array(dim=c(dim(df.met)[1:2],length(which(years<1850))))
  df.deriv2 <- array(dim=c(dim(df.met)[1:2],length(which(years>1900))))
  dimnames(df.deriv1)[[1]] <- dimnames(df.deriv2)[[1]] <- lon
  dimnames(df.deriv1)[[2]] <- dimnames(df.deriv2)[[2]] <- lat
  dimnames(df.deriv1)[[3]] <- years[which(years<1850)]
  dimnames(df.deriv2)[[3]] <- years[which(years>1900)]
  names(dimnames(df.deriv1)) <- names(dimnames(df.deriv2)) <- c("lon", "lat", "Year")

  for(x in 1:length(lon)){
    for(y in 1:length(lat)){
      if(is.na(max(df.met[x,y,]))) next
      df.deriv1[x,y,] <- deriv.out1[[paste0("lat", lat[y], "lon", lon[x])]]
      df.deriv2[x,y,] <- deriv.out2[[paste0("lat", lat[y], "lon", lon[x])]]
    }
  }

  # Finding and graphing the mean derivative
  deriv.mean1 <- apply(abs(df.deriv1), c(1,2), FUN=mean)
  deriv.mean2 <- apply(abs(df.deriv2), c(1,2), FUN=mean)
  abs.mean1   <- apply(df.met[,,which(years<1850)], c(1,2), FUN=mean)
  abs.mean2   <- apply(df.met[,,which(years>1900)], c(1,2), FUN=mean)
  summary(deriv.mean1)
  
  deriv.stack <- stack(data.frame(deriv.mean1))    
  names(deriv.stack) <- c("deriv.pre1850", "lat")
  deriv.stack$lon <- as.numeric(paste(dimnames(df.deriv1)[[1]]))
  deriv.stack$lat <- as.numeric(substr(deriv.stack$lat, 2, nchar(paste(deriv.stack$lat))))
  deriv.stack$deriv.post1900 <- stack(data.frame(deriv.mean2))[,1]    
  deriv.stack$mean.pre1850   <- stack(data.frame(abs.mean1))[,1]    
  deriv.stack$mean.post1900  <- stack(data.frame(abs.mean2))[,1]    
  summary(deriv.stack)
  # ---------------------------

  # ---------------------------
  # Graphing the met var & change
  # ---------------------------
  # ------------
  # Means
  # ------------
  mean1 <- ggplot(data=deriv.stack) +
    geom_raster(aes(x=lon, y=lat, fill=mean.pre1850)) +
    geom_point(data=hips, aes(x=lon, y=lat), color="red") +
    scale_x_continuous(expand=c(0,0)) +
    scale_y_continuous(expand=c(0,0)) +
    scale_fill_gradient(name=paste0(v), limits=range(deriv.stack[,c("mean.pre1850", "mean.post1900")], na.rm=T)) +
    ggtitle(paste0(v, " pre-1850 mean")) +
    coord_equal(ratio=1)
  mean2 <- ggplot(data=deriv.stack) +
    geom_raster(aes(x=lon, y=lat, fill=mean.post1900)) +
    geom_point(data=hips, aes(x=lon, y=lat), color="red") +
    scale_x_continuous(expand=c(0,0)) +
    scale_y_continuous(expand=c(0,0)) +
    scale_fill_gradient(name=paste0(v), limits=range(deriv.stack[,c("mean.pre1850", "mean.post1900")], na.rm=T)) +
    ggtitle(paste0(v, " post-1900 mean")) +
    coord_equal(ratio=1)

  png(file.path(fig.out, paste0(v, "_means.png")), height=8, width=11, units="in", res=200)
  grid.arrange(mean1, mean2, ncol=1)
  dev.off()
  # ------------
  
  # ------------
  # Rates of change
  # ------------
  # ggplot(data=deriv.stack) +
  #   geom_raster(aes(x=lon, y=lat, fill=deriv.pre1850)) +
  #   geom_point(data=hips, aes(x=lon, y=lat), color="red") +
  #   scale_x_continuous(expand=c(0,0)) +
  #   scale_y_continuous(expand=c(0,0)) +
  #   scale_fill_gradient(name=paste0(v, " deriv")) +
  #   ggtitle(paste0(v, "Mean Rate of Change (absolute), pre-1850 (free scale)")) +
  #   coord_equal(ratio=1)
  deriv1 <- ggplot(data=deriv.stack) +
    geom_raster(aes(x=lon, y=lat, fill=deriv.pre1850)) +
    geom_point(data=hips, aes(x=lon, y=lat), color="red") +
    scale_x_continuous(expand=c(0,0)) +
    scale_y_continuous(expand=c(0,0)) +
    scale_fill_gradient(name=paste0(v, " deriv"), limits=range(deriv.stack[,c("deriv.pre1850", "deriv.post1900")], na.rm=T)) +
    ggtitle(paste0(v, " -- Mean Rate of Change (absolute), pre-1850")) +
    coord_equal(ratio=1)
  deriv2 <- ggplot(data=deriv.stack) +
    geom_raster(aes(x=lon, y=lat, fill=deriv.post1900)) +
    geom_point(data=hips, aes(x=lon, y=lat), color="red") +
    scale_x_continuous(expand=c(0,0)) +
    scale_y_continuous(expand=c(0,0)) +
    scale_fill_gradient(name=paste0(v, " deriv"), limits=range(deriv.stack[,c("deriv.pre1850", "deriv.post1900")], na.rm=T)) +
    ggtitle(paste0(v, " -- Mean Rate of Change (absolute), Post-1900")) +
    coord_equal(ratio=1)

  png(file.path(fig.out, paste0(v, "_derivs.png")), height=8, width=11, units="in", res=200)
  grid.arrange(deriv1, deriv2, ncol=1)
  dev.off()
  # ------------
  # ---------------------------
    
  # ---------------------------
  # Saving Output
  # ---------------------------
  save(df.deriv1, df.deriv2, file=file.path(out.dir, paste0(v, "_derivs.Rdata")))
  # ---------------------------
  
}  
# Save HIPS met + derivs
write.csv(met.sites, file=file.path(out.dir, "Met_Sites_Annual.csv"), row.names=F)
}
# -------------------------------------------


# -------------------------------------------
# Analyzing Met stability and comparing with Paleo drivers
# -------------------------------------------
source("R/0_TimeAnalysis.R")

# Note: I've saved the 
met.sites <- read.csv(file.path(out.dir, "Met_Sites_Annual.csv"))
summary(met.sites)

# Package the data -- Variables in separate layers to parallelize
dat.list <- list()
dat.list2 <- list()
for(v in met.vars){
    dat.list[[v]] <- met.sites[met.sites$Year<1850, c("Site", "Year", v)]
    dat.list2[[v]] <- met.sites[met.sites$Year>1900, c("Site", "Year", v)]

    # Making a dummy variable of "v" to make the parallizing work 
    dat.list[[v]][,"v"] <- dat.list[[v]][,v]
    dat.list2[[v]][,"v"] <- dat.list2[[v]][,v]
}
  
# Run the stats
cores.use <- min(12, length(dat.list))
dat.out  <- mclapply(dat.list, analyze.time, mc.cores=cores.use, Y="v", fac.fit="Site", k.freq=25, path.gamm.func=path.gamm.func)
dat.out2 <- mclapply(dat.list2, analyze.time, mc.cores=cores.use, Y="v", fac.fit="Site", k.freq=25, path.gamm.func=path.gamm.func)
 
# Format the output
for(m in names(dat.out)){
  dat.out[[m]]$out$Model <- as.factor(m)
  dat.out[[m]]$out$var   <- as.factor(m)
  dat.out2[[m]]$out$Model <- as.factor(m)
  dat.out2[[m]]$out$var   <- as.factor(m)
  if(m == names(dat.out)[1]){
    dat.out3 <- rbind(dat.out[[m]]$out, dat.out2[[m]]$out)
  } else{
    dat.out3 <- rbind(dat.out3, rbind(dat.out[[m]]$out, dat.out2[[m]]$out))
  }
} # End model formatting

# dat.out3 <- dat.out3[!is.na(dat.out3$Y),]
summary(dat.out3)

ggplot(dat.out3[,]) +
  facet_wrap(~Model, scales="free_y") +
  geom_line(aes(x=Year, y=Y, color=Site), size=0.2)
# geom_point(aes(x=Year, y=mean, color=Site), size=0.2)

# Adding the var to met sites to make the merge happen correctly
met.sites2 <- stack(met.sites[,met.vars])
names(met.sites2) <- c("Y", "var")
met.sites2$Year <- met.sites$Year
met.sites2$Site <- met.sites$Site
met.sites2$Model <- met.sites2$var
summary(met.sites2)

dat.out4 <- merge(met.sites2, dat.out3, all.x=T, all.y=T)
summary(dat.out4)

## Some additional formatting
dat.out4$mean.sig <- ifelse(dat.out4$sig=="*", dat.out4$mean, NA)
summary(dat.out4)

# Save the output
write.csv(dat.out4, file.path(out.dir, paste0("StabilityCalcs_SiteMet_850-1850_1900-2010.csv")), row.names=F)
save(dat.out, dat.out2, file=file.path(out.dir, paste0("StabilityCalcs_SiteMet_850-1850_1900-2010.RData")))
  
png(file.path(fig.out, paste0("Stability_SiteMet_850-1850_1900-2010.png")), height=11, width=8.5, units="in", res=180)
print(
  ggplot(data=dat.out4[,]) + 
    facet_grid(var~., scales="free_y") +
    geom_line(aes(x=Year, y=Y, color=Site), size=0.5, alpha=0.3) +
    geom_ribbon(data=dat.out4[dat.out4$Year<1850,], aes(x=Year, ymin=lwr, ymax=upr, fill=Site), alpha=0.3) +
    geom_ribbon(data=dat.out4[dat.out4$Year>1900,], aes(x=Year, ymin=lwr, ymax=upr, fill=Site), alpha=0.3) +
    geom_line(data=dat.out4[dat.out4$Year<1850,], aes(x=Year, y=mean, color=Site), size=1, alpha=0.2) +
    geom_line(data=dat.out4[dat.out4$Year>1900,], aes(x=Year, y=mean, color=Site), size=2, alpha=1) +
    geom_line(aes(x=Year, y=mean.sig, color=Site), size=2, alpha=1) +
    geom_vline(xintercept=1850, linetype="dashed") +
    geom_vline(xintercept=1900, linetype="dashed") +
    scale_x_continuous(expand=c(0,0), name="Year") +
    scale_y_continuous(expand=c(0,0), name=v) +
    # scale_color_manual(values=col.model) +
    # scale_fill_manual(values=col.model) + 
    theme_bw()
)
dev.off()
  

# -------------------------------------------
