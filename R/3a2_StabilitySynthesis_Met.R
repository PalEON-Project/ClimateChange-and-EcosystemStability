# Looking at Spatio-Temporal patterns of Climate Change in the regional met drivers
# This is partially for the stability manuscript & partially for site selection of the regional drivers

# Figures to make:
# 1. Spatial Map: Mean change in Tair & Precip by century
# 2. Climate Space: Absolute means, derivs, w/ sites
# 3. 

# -------------------------------------------
# Load libraries; set file paths
# -------------------------------------------
library(ncdf4)
library(ggplot2); library(gridExtra); library(scales)
library(mgcv)
library(plyr); library(parallel)
library(stringr)
setwd("~/Dropbox/PalEON_CR/PalEON_MIP_Site/Analyses/Change-and-Stability")
# path.gamm.func <- "~/Desktop/R_Functions/"  # Path to github repository of my GAMM helper functions: https://github.com/crollinson/R_Functions.git
path.gamm.func <- "~/Desktop/Research/R_Functions/"
# mip.utils <- "~/Dropbox/Research/PalEON_CR/MIP_Utils/" # Path to PalEON MIP Utility repository: https://github.com/PalEON-Project/MIP_Utils.git
mip.utils <- "~/Desktop/Research/PalEON_CR/MIP_Utils/" # Path to PalEON MIP Utility repository: https://github.com/PalEON-Project/MIP_Utils.git

raw.dir <- "raw_data"
out.dir <- "Data/Met"
fig.out <- "Figures/Met"
dir.co2 <- "~/Dropbox/PalEON_CR/env_regional/env_paleon/co2/paleon_annual_co2.nc"


source("R/0_TimeAnalysis.R")
source(file.path(path.gamm.func, "Calculate_GAMM_Derivs.R"))
met.vars <- c("tair", "precipf", "swdown", "lwdown", "psurf", "qair", "wind", "co2")
hips <- data.frame(Site = c( "PHA",  "PHO",  "PUN",  "PBL",  "PDL",  "PMB"),
                   lat  = c( 42.54,  45.25,  46.22,  46.28,  47.17,  43.61),
                   lon  = c(-72.18, -68.73, -89.53, -94.58, -95.17, -82.83),
                   lat2 = c( 42.75,  45.25,  46.25,  46.25,  47.25,  43.75),
                   lon2 = c(-72.25, -68.75, -89.75, -94.75, -95.25, -82.75))
years <- 850:2010

if(!dir.exists(out.dir)) dir.create(out.dir)
if(!dir.exists(fig.out)) dir.create(fig.out)
# -------------------------------------------


# -------------------------------------------
# Loading and organizing the data
# -------------------------------------------

# ------------
# Loading the raw met drivers
# ------------
met.ann <- list()

# 1. Temperature
nc.v <- nc_open(file.path(out.dir, paste0("tair", ".nc")))
# Extract the regional data
lon <- ncvar_get(nc.v, "lon")
lat <- ncvar_get(nc.v, "lat")
met.ann[["tair"]] <- ncvar_get(nc.v, "tair")
dimnames(met.ann$tair) <- list(lon=lon, lat=lat, Year=years)
nc_close(nc.v)

# 2. Precipitation
# Temperature first
nc.v <- nc_open(file.path(out.dir, paste0("precipf", ".nc")))
# Extract the regional data
lon <- ncvar_get(nc.v, "lon")
lat <- ncvar_get(nc.v, "lat")
met.ann[["precipf"]] <- ncvar_get(nc.v, "precipf")
dimnames(met.ann$precipf) <- list(lon=lon, lat=lat, Year=years)
nc_close(nc.v)

summary(met.ann)
# ------------

# ------------
# Now, load the derivatives
# ------------
met.derivs <- list()

load(file.path(out.dir, paste0("tair", "_derivs.Rdata")))
met.derivs[["tair"]] <- list(pre1850=df.deriv1, post1900=df.deriv2)

load(file.path(out.dir, paste0("precipf", "_derivs.Rdata")))
met.derivs[["precipf"]] <- list(pre1850=df.deriv1, post1900=df.deriv2)
# ------------
# -------------------------------------------


# -------------------------------------------
# Doing some maps by century for Andria & Jack et al.
# 900-1800 using 50 years as midpoints (a la STEPPS) + 1900-2000
# 4-panel figures:
#  Rows: Temp, Precip
#  Columns: Means Values, Mean Change
# -------------------------------------------
# Just getting the range of centennial temp & precip means for color mapping
range.t.abs <- NA
range.p.abs <- NA
range.temp <- NA
range.precip <- NA
for(i in seq(from=900, to=1800, by=100)){
  temp.t <- range(apply(met.ann$tair[,,which(years>=i & years<=i+100) ], c(1,2), FUN=mean), na.rm=T)
  temp.p <- range(apply(met.ann$precipf[,,which(years>=i & years<=i+100) ], c(1,2), FUN=mean), na.rm=T)
  
  range.t.abs <- range(range.t.abs, temp.t, na.rm=T)
  range.p.abs <- range(range.p.abs, temp.p, na.rm=T)
  
  if(i == 1800) next
  temp.t <- range(apply(met.derivs$tair   $pre1850[,,which(years>=i & years<=i+100) ], c(1,2), FUN=mean), na.rm=T)
  temp.p <- range(apply(met.derivs$precipf$pre1850[,,which(years>=i & years<=i+100) ], c(1,2), FUN=mean), na.rm=T)
  
  range.temp <- range(range.temp, temp.t, na.rm=T)
  range.precip <- range(range.precip, temp.p, na.rm=T)
}
range.t.abs
range.p.abs

for(i in seq(from=900, to=1999, by=100)){
  if(i == 1800) next
  if(i<1900){
    tair.range <- range.temp
    precip.range <- range.precip
  } else {
    tair.range   <- range(apply(met.derivs$tair   $post1900, c(1,2), FUN=mean), na.rm=T)
    precip.range <- range(apply(met.derivs$precipf$post1900, c(1,2), FUN=mean), na.rm=T)
  }
  
  # Putting things into arrays --> finding the means for each time period of interest
  tair.mean    <- apply(met.ann$tair   [,,which(years>=i & years<=i+100) ], c(1,2), FUN=mean)
  precipf.mean <- apply(met.ann$precipf[,,which(years>=i & years<=i+100) ], c(1,2), FUN=mean)
  
  if(i < 1800){
    tair.deriv    <- apply(met.derivs$tair   $pre1850[,,which(years>=i & years<=i+100) ], c(1,2), FUN=mean)
    precipf.deriv <- apply(met.derivs$precipf$pre1850[,,which(years>=i & years<=i+100) ], c(1,2), FUN=mean)
  } else if(i==1900) {
    tair.deriv    <- apply(met.derivs$tair   $post1900[,,1:100], c(1,2), FUN=mean)
    precipf.deriv <- apply(met.derivs$precipf$post1900[,,1:100], c(1,2), FUN=mean)
  } else {
    tair.deriv    <- array(NA, dim=dim(tair.mean))
    precipf.deriv <- array(NA, dim=dim(tair.mean))
  }
  
  
  # Finding and graphing the mean derivative
  dat.out <- stack(data.frame(tair.mean))
  names(dat.out) <- c("tair", "lat")
  dat.out$lon <- as.numeric(paste(dimnames(tair.mean)[[1]]))
  dat.out$lat <- as.numeric(substr(dat.out$lat, 2, nchar(paste(dat.out$lat))))
  dat.out$precipf <- stack(data.frame(precipf.mean))[,1]
  dat.out$tair.deriv    <- stack(data.frame(tair.deriv))[,1]
  dat.out$precipf.deriv <- stack(data.frame(precipf.deriv))[,1]
  summary(dat.out)
  
  # Creating the graphs
  plot.tair <- ggplot(data=dat.out) +
    geom_raster(aes(x=lon, y=lat, fill=tair)) +
    geom_point(data=hips, aes(x=lon, y=lat), color="gray80", size=2) +
    scale_x_continuous(expand=c(0,0)) +
    scale_y_continuous(expand=c(0,0)) +
    scale_fill_gradient(name=paste0("K"), high="red", limits=range.t.abs) +
    ggtitle(paste0("Temp: ", i, " - ", i+100, " A.D.")) +
    theme(panel.border=element_rect(color="black", fill=NA)) +
    coord_equal(ratio=1)
  plot.precipf <- ggplot(data=dat.out) +
    geom_raster(aes(x=lon, y=lat, fill=precipf)) +
    geom_point(data=hips, aes(x=lon, y=lat), color="gray80", size=2) +
    scale_x_continuous(expand=c(0,0)) +
    scale_y_continuous(expand=c(0,0)) +
    scale_fill_gradient(name=paste0("kg s-1"), limits=range.p.abs) +
    ggtitle(paste0("Precip: ", i, " - ", i+100, " A.D.")) +
    theme(panel.border=element_rect(color="black", fill=NA)) +
    coord_equal(ratio=1)

  plot.tair.deriv <- ggplot(data=dat.out) +
    geom_raster(aes(x=lon, y=lat, fill=tair.deriv)) +
    geom_point(data=hips, aes(x=lon, y=lat), color="black", size=2) +
    scale_x_continuous(expand=c(0,0)) +
    scale_y_continuous(expand=c(0,0)) +
    # scale_fill_gradient(name=paste0("K yr-1"), high="red") +
    scale_fill_gradient2(name=paste0("K yr-1"), limits=tair.range, low=muted("blue"), mid="white", high=muted("red"), midpoint=0) +
    ggtitle(paste0("Change in Temp: ", i, " - ", i+100, " A.D.")) +
    theme(panel.border=element_rect(color="black", fill=NA)) +
    coord_equal(ratio=1)
  plot.precipf.deriv <- ggplot(data=dat.out) +
    geom_raster(aes(x=lon, y=lat, fill=precipf.deriv)) +
    geom_point(data=hips, aes(x=lon, y=lat), color="black", size=2) +
    scale_x_continuous(expand=c(0,0)) +
    scale_y_continuous(expand=c(0,0)) +
    # scale_fill_gradient(name=paste0("kg s-1 yr-1")) +
    scale_fill_gradient2(name=paste0("kg s-1 yr-1"), limits=precip.range, low=muted("blue"), mid="white", high=muted("red"), midpoint=0) +
    ggtitle(paste0("Change in Precip: ", i, " - ", i+100, " A.D.")) +
    theme(panel.border=element_rect(color="black", fill=NA)) +
    coord_equal(ratio=1)
  
  png(file.path(fig.out, paste0("ClimateChange_means_", str_pad(i, 4, pad="0"), "-", i+100, ".png")), height=5, width=11, unit="in", res=200)
  grid.arrange(plot.tair, plot.tair.deriv, plot.precipf, plot.precipf.deriv, ncol=2)
  dev.off()
  
}
# -------------------------------------------


# -------------------------------------------
# Mapping Climate Space (T vs. P)
# Since spatial patterns are fairly stable through time, just going to use the 1900-2010 mean
# -------------------------------------------
tair.mean    <- apply(met.ann$tair   [,,which(years>1900)], c(1,2), FUN=mean)
precipf.mean <- apply(met.ann$precipf[,,which(years>1900) ], c(1,2), FUN=mean)
tair.deriv.max    <- apply(met.derivs$tair   $pre1850[,,], c(1,2), FUN=function(x){x[which(abs(x)==max(abs(x)))[1]]})
precipf.deriv.max <- apply(met.derivs$precipf$pre1850[,,], c(1,2), FUN=function(x){x[which(abs(x) == max(abs(x)))[1]]})

# Finding and graphing the mean derivative
dat.out <- stack(data.frame(tair.mean))
names(dat.out) <- c("tair.mean", "lat")
dat.out$lon <- as.numeric(paste(dimnames(tair.mean)[[1]]))
dat.out$lat <- as.numeric(substr(dat.out$lat, 2, nchar(paste(dat.out$lat))))
dat.out$precipf.mean <- stack(data.frame(precipf.mean))[,1]
dat.out$tair.deriv.max    <- stack(data.frame(tair.deriv.max))[,1]
dat.out$precipf.deriv.max <- stack(data.frame(precipf.deriv.max))[,1]
summary(dat.out)

# Loading the ED coarse grid to see how that gets our environmental space
ed.status <- read.csv("~/Dropbox/PalEON_CR/ED_PalEON/MIP2_Region/0_setup/Paleon_MIP_Phase2_ED_Order_Status.csv")
summary(ed.status)


#  ---------------
# Extracting Temp & Precip for HIPS sites
#  ---------------
for(i in 1:nrow(hips)){
  row.dat <-  which(dat.out$lat-0.25<=hips[i,"lat"] & dat.out$lat+0.25>=hips[i,"lat"] & dat.out$lon-0.25<=hips[i,"lon"] & dat.out$lon+0.25>=hips[i,"lon"])
  
  hips[i,"tair.mean"] <- dat.out[row.dat,"tair.mean"]
  hips[i,"precipf.mean"] <- dat.out[row.dat,"precipf.mean"]
  hips[i,"tair.deriv.max"] <- dat.out[row.dat,"tair.deriv.max"]
  hips[i,"precipf.deriv.max"] <- dat.out[row.dat,"precipf.deriv.max"]
}

# Subset everything coarser than 5 degrees for ED Coarse Grid
summary(ed.status$notes)
grid.5degree <- ed.status[ed.status$notes %in% c("10-degree", "5-degree"),]
nrow(grid.5degree)
for(i in 1:nrow(grid.5degree)){
  row.dat <-  which(dat.out$lat-0.25<=grid.5degree[i,"lat"] & dat.out$lat+0.25>=grid.5degree[i,"lat"] & dat.out$lon-0.25<=grid.5degree[i,"lon"] & dat.out$lon+0.25>=grid.5degree[i,"lon"])
  
  grid.5degree[i,"tair.mean"] <- dat.out[row.dat,"tair.mean"]
  grid.5degree[i,"precipf.mean"] <- dat.out[row.dat,"precipf.mean"]
  grid.5degree[i,"tair.deriv.max"] <- dat.out[row.dat,"tair.deriv.max"]
  grid.5degree[i,"precipf.deriv.max"] <- dat.out[row.dat,"precipf.deriv.max"]
}


sec2yr <- 60*60*24*365
plot.means <- ggplot(data=dat.out) +
  stat_bin2d(aes(x=tair.mean, y=precipf.mean*sec2yr)) +
  geom_point(data=grid.5degree, aes(x=tair.mean, y=precipf.mean*sec2yr), col="gray50", size=5) +
  geom_text(data=hips, aes(x=tair.mean, y=precipf.mean*sec2yr, label=Site), fontface="bold", col="indianred", size=8) +
  labs(x="Tmean (K, 1900-2000)", y="Precip (mm/yr, 1900-2000") +
  ggtitle("Modern Climate Space") +
  geom_hline(yintercept=mean(dat.out$precipf.mean*sec2yr, na.rm=T), linetype="dashed", color="red") +
  geom_vline(xintercept=mean(dat.out$tair.mean, na.rm=T), linetype="dashed", color="red") +
  # scale_fill_gradient2(low="black", mid="gray50", high="gray80") +
  theme_bw()

plot.derivs <- ggplot(data=dat.out) +
  stat_bin2d(aes(x=tair.deriv.max, y=precipf.deriv.max*sec2yr)) +
  geom_point(data=grid.5degree, aes(x=tair.deriv.max, y=precipf.deriv.max*sec2yr), col="gray50", size=5) +
  geom_text(data=hips, aes(x=tair.deriv.max, y=precipf.deriv.max*sec2yr, label=Site), fontface="bold", col="indianred", size=8) +
  labs(x="Tmean (K yr-1, 1900-2000)", y="Precip (mm/yr yr-1, 1900-2000") +
  ggtitle("Pre-1850 Maximum Change") +
  geom_hline(yintercept=0, linetype="dashed", color="red") +
  geom_vline(xintercept=0, linetype="dashed", color="red") +
  # scale_fill_gradient2(low="black", mid="gray50", high="gray80") +
  theme_bw()

plot.tair <-  ggplot(data=dat.out) +
  stat_bin2d(aes(x=tair.mean, y=tair.deriv.max)) +
  geom_point(data=grid.5degree, aes(x=tair.mean, y=tair.deriv.max), col="gray50", size=5) +
  geom_text(data=hips, aes(x=tair.mean, y=tair.deriv.max, label=Site), fontface="bold", col="indianred", size=8) +
  labs(x="Tmean (K yr-1, 1900-2000)", y="Tmean Change (K yr-1, 1900-2000)") +
  ggtitle("Tair: Climate versus Pre-1850 Maximum Change") +
  geom_hline(yintercept=0, linetype="dashed", color="red") +
  geom_vline(xintercept=mean(dat.out$tair.mean, na.rm=T), linetype="dashed", color="red") +
  # scale_fill_gradient2(low="black", mid="gray50", high="gray80") +
  theme_bw()

plot.precipf <-  ggplot(data=dat.out) +
  stat_bin2d(aes(x=precipf.mean*sec2yr, y=precipf.deriv.max*sec2yr)) +
  geom_point(data=grid.5degree, aes(x=precipf.mean*sec2yr, y=precipf.deriv.max*sec2yr), col="gray50", size=5) +
  geom_text(data=hips, aes(x=precipf.mean*sec2yr, y=precipf.deriv.max*sec2yr, label=Site), fontface="bold", col="indianred", size=8) +
  labs(x="Precip (mm/yr, 1900-2000", y="Precip Change (mm/yr yr-1, 1900-2000") +
  ggtitle("precipf: Climate versus Pre-1850 Maximum Change") +
  geom_hline(yintercept=0, linetype="dashed", color="red") +
  geom_vline(xintercept=mean(dat.out$precipf.mean*sec2yr, na.rm=T), linetype="dashed", color="red") +
  # scale_fill_gradient2(low="black", mid="gray50", high="gray80") +
  theme_bw()

png(file.path(fig.out, paste0("ClimateSpace_HIPS_5degreeGrid.png")), height=8, width=11, units="in", res=200)
grid.arrange(plot.means, plot.derivs, ncol=1)
dev.off()
#  ---------------

# -------------------------------------------
