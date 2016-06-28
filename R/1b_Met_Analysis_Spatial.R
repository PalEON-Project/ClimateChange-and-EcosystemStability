# -------------------------------------------
# 1. Spatial descriptions of meteorology and change from initial state
# Author: Christy Rollinson, crollinson@gmail.com
# -------------------------------------------
rm(list=ls())


# -------------------------------------------
# Load libraries; set file paths
# -------------------------------------------
library(ncdf4)
library(ggplot2); library(gridExtra); library(grid); library(scales)
library(mgcv)
library(plyr); library(parallel)
setwd("~/Dropbox/PalEON_CR/PalEON_MIP_Site/Analyses/Change-and-Stability")
path.gamm.func <- "~/Desktop/R_Functions/"  # Path to github repository of my GAMM helper functions: https://github.com/crollinson/R_Functions.git
mip.utils <- "~/Dropbox/PalEON_CR/MIP_Utils/" # Path to PalEON MIP Utility repository: https://github.com/PalEON-Project/MIP_Utils.git

raw.dir <- "raw_data"
dat.out <- "Data/Met"
fig.out <- "Figures/Met"

if(!dir.exists(dat.out)) dir.create(dat.out)
if(!dir.exists(fig.out)) dir.create(fig.out)
# -------------------------------------------


# -------------------------------------------
# Define some useful variables & Functions
# -------------------------------------------
met.vars <- c("tair", "precipf", "swdown", "lwdown", "psurf", "qair", "wind")
hips <- data.frame(Site = c( "PHA",  "PHO",  "PUN",  "PBL",  "PDL",  "PMB"),
                   lat  = c( 42.54,  45.25,  46.22,  46.28,  47.17,  43.61),
                   lon  = c(-72.18, -68.73, -89.53, -94.58, -95.17, -82.83),
                   lat2 = c( 42.75,  45.25,  46.25,  46.25,  47.25,  43.75),
                   lon2 = c(-72.25, -68.75, -89.75, -94.75, -95.25, -82.75))
years <- 850:2010

calc.stability <- function(x){
  dat.tmp <- data.frame(Y=x, Year=1:length(x))
  k.use=round(length(years)/25, 0)
  mod.gam <- gam(Y ~ s(Year, k=k.use), data=dat.tmp)
  mod.deriv <- calc.derivs(mod.gam, newdata=dat.tmp, vars="Year")
  return(mod.deriv$mean)
}

# -------------------------------------------


# -------------------------------------------
# Extract Met means: 1830-1849 (GCM only); 1991-2010 (CRUNCEP)
# -------------------------------------------
source("R/0_TimeAnalysis.R")
source(file.path(path.gamm.func, "Calculate_GAMM_Derivs.R"))

met.sites <- data.frame(merge(hips, years))
names(met.sites) <- c(names(hips), "Year")

summary(met.sites)

for(v in met.vars){
  print(paste0(" **** ", v, " **** "))
  nc.v <- nc_open(file.path(dat.out, paste0(v, ".nc")))
  
  # Running The Derivative calculation on each metvar
  years <- 850:2010
  
  # Extract the regional data
  lon <- ncvar_get(nc.v, "lon")
  lat <- ncvar_get(nc.v, "lat")
  df.met <- ncvar_get(nc.v, v)
  dimnames(df.met) <- list(lon=lon, lat=lat, Year=years)
  
  # Extract the entire time series for sites
  for(s in unique(met.sites$Site)){
    y.use <- which(lat==mean(met.sites[met.sites$Site==s, "lat2"]))
    x.use <- which(lon==mean(met.sites[met.sites$Site==s, "lon2"]))
    met.sites[met.sites$Site==s,v] <- ncvar_get(nc.v, v)[x.use, y.use, ]
  }
  
  # Extract the mean for the spinup period (first 20 years) for the entire region
  met0 <- apply(df.met[,,which(years>=850 & years<=869)], c(1,2), mean, na.rm=F)
  met.x0 <- stack(data.frame(met0))
  names(met.x0) <- c(v, "lat")
  met.x0$lat <- as.numeric(substr(met.x0$lat, 2, 6))
  met.x0$lon <- as.numeric(dimnames(df.met)[[1]])
  met.x0$Time <- as.factor("Spinup")
  
  # Extract the mean for 1830-1849 (20 years before 1850) for the entire region
  met1 <- apply(df.met[,,which(years>=1830 & years<=1850)], c(1,2), mean, na.rm=F)
  met.x1 <- stack(data.frame(met1))
  names(met.x1) <- c(v, "lat")
  met.x1$lat <- as.numeric(substr(met.x1$lat, 2, 6))
  met.x1$lon <- as.numeric(dimnames(df.met)[[1]])
  met.x1$Time <- as.factor("Settlement")

  # Extract the mean for 1991-2010 (most recent 20 years) for the entire region
  met2 <- apply(df.met[,,which(years>=1991 & years<=2010)], c(1,2), mean, na.rm=F)
  met.x2 <- stack(data.frame(met2))
  names(met.x2) <- c(v, "lat")
  met.x2$lat <- as.numeric(substr(met.x2$lat, 2, 6))
  met.x2$lon <- as.numeric(dimnames(df.met)[[1]])
  met.x2$Time <- as.factor("Modern")
  
  # Extract the difference between the modern and pre-set avgs (most recent 20 years) for the entire region
  met.diff <- met2 - met1
  met.x3 <- stack(data.frame(met.diff))
  names(met.x3) <- c(v, "lat")
  met.x3$lat <- as.numeric(substr(met.x3$lat, 2, 6))
  met.x3$lon <- as.numeric(dimnames(df.met)[[1]])
  met.x3$Time <- as.factor("diff.Settlement")
  
  # Extract the difference between the modern and pre-set avgs (most recent 20 years) for the entire region
  met.diff2 <- met2 - met0
  met.x4 <- stack(data.frame(met.diff2))
  names(met.x4) <- c(v, "lat")
  met.x4$lat <- as.numeric(substr(met.x4$lat, 2, 6))
  met.x4$lon <- as.numeric(dimnames(df.met)[[1]])
  met.x4$Time <- as.factor("diff.Spinup")
  
  
  # Bind the data together
  if(v==met.vars[1]){
    met.region <- rbind(met.x0, met.x1, met.x2, met.x3, met.x4)
  } else {
    met.region[,v] <- rbind(met.x0, met.x1, met.x2, met.x3, met.x4)[,v]
  }
  
  nc_close(nc.v)
}

summary(met.sites)
summary(met.region)

# Convert some units
# Precip: kg/m2/s (== mm/s) to mm/yr
sec2yr <- 60*60*24*365
met.sites $precipf <- met.sites $precipf*sec2yr
met.region$precipf <- met.region$precipf*sec2yr
met.sites $tair    <- met.sites $tair - 273.15
met.region[!substr(met.region$Time, 1,4)=="diff", "tair"]    <- met.region[!substr(met.region$Time, 1,4)=="diff", "tair"] - 273.15
met.region <- met.region[complete.cases(met.region),]


summary(met.sites[,])
summary(met.region[met.region$Time=="Modern",])
summary(met.region[met.region$Time=="diff.Spinup",])
summary(met.region[met.region$Time=="diff.Settlement",])

write.csv(met.region, file.path(dat.out, "Met_Region_Annual.csv"), row.names=F)
# -------------------------------------------
