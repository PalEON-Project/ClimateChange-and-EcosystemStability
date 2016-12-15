# -------------------------------------------
# Assessing the scales and magnitude of change in ecosystem state
# Author: Christy Rollinson, crollinson@gmail.com

# 1. Stastical detection of significant change
#    - Use TP regression splines
# -------------------------------------------
rm(list=ls())


# -------------------------------------------
# Load libraries; set file paths
# -------------------------------------------
library(car)
library(ggplot2)
library(mgcv)

library(gridExtra); library(gtable)
library(cowplot)

setwd("~/Desktop/Research/Conferences etc/AGU 2016/AGU2016_code")
path.gamm.func <- "~/Desktop/R_Functions/"  # Path to github repository of my GAMM helper functions: https://github.com/crollinson/R_Functions.git
inputs    <- "~/Desktop/Research/PalEON_CR/PalEON_MIP_Site/phase1a_output_variables/" # Path to my cleaned model output
mip.utils <- "~/Desktop/Research/PalEON_CR/MIP_Utils/" # Path to PalEON MIP Utility repository: https://github.com/PalEON-Project/MIP_Utils.git
path.raw <- "~/Desktop/Research/PalEON_CR/PalEON_MIP_Site/phase1a_model_output/" # Path to raw model output

analy.out <- "~/Desktop/Research/PalEON_CR/PalEON_MIP_Site/Analyses/Change-and-Stability/Data/EcosysChange/"

out.dir <- "data/" # Path to where the analysis output should go
fig.dir <- "figures/" # Path to where figures should go

if(!dir.exists(out.dir)) dir.create(out.dir)
if(!dir.exists(fig.dir)) dir.create(fig.dir)
# -------------------------------------------


# -------------------------------------------
# Reading in raw ecosystem model output & Adding the PFT info
# -------------------------------------------
{
  ecosys <- read.csv("~/Desktop/Research/PalEON_CR/PalEON_MIP_Site/Analyses/Change-and-Stability/Data/PalEON_MIP_Yearly_withFcomp.csv")
  ecosys$Model.Order <- recode(ecosys$Model, "'clm.bgc'='01'; 'clm.cn'='02'; 'ed2'='03'; 'ed2.lu'='04';  'jules.stat'='05'; 'jules.triffid'='06'; 'linkages'='07'; 'lpj.guess'='08'; 'lpj.wsl'='09'; 'sibcasa'='10'")
  levels(ecosys$Model.Order) <- c("CLM-BGC", "CLM-CN", "ED2", "ED2-LU", "JULES-STATIC", "JULES-TRIFFID", "LINKAGES", "LPJ-GUESS", "LPJ-WSL", "SiBCASA")
  summary(ecosys)
  
  # Colors used for graphing
  model.colors <- read.csv("~/Desktop/Research/PalEON_CR/PalEON_MIP_Site/Model.Colors.csv")
  model.colors $Model.Order <- recode(model.colors$Model, "'CLM4.5-BGC'='01'; 'CLM4.5-CN'='02'; 'ED2'='03'; 'ED2-LU'='04';  'JULES-STATIC'='05'; 'JULES-TRIFFID'='06'; 'LINKAGES'='07'; 'LPJ-GUESS'='08'; 'LPJ-WSL'='09'; 'SiBCASA'='10'")
  levels(model.colors$Model.Order)[1:10] <- c("CLM-BGC", "CLM-CN", "ED2", "ED2-LU", "JULES-STATIC", "JULES-TRIFFID", "LINKAGES", "LPJ-GUESS", "LPJ-WSL", "SiBCASA")
  model.colors <- model.colors[order(model.colors$Model.Order),]
  model.colors
  
  
  model.names <- data.frame(Model=unique(ecosys$Model), Model.Order=unique(ecosys$Model.Order))
  vars.agu <- c("GPP","AGB", "Fcomp") # Add dominant PFT
}
# -------------------------------------------
 
# -------------------------------------------
# Plotting PBL AGB for Ann (& me)
# -------------------------------------------
{
  summary(ecosys)
  models.agb <- unique(ecosys[ecosys$AGB>0 & !is.na(ecosys$AGB),"Model.Order"])
  colors.agb <- model.colors[model.colors$Model.Order %in% models.agb, "color"]
  
  png(file.path(fig.dir, "AGB_PBL.png"), height=9, width=14, units="in", res=180)
  print(
  ggplot(data=ecosys[ecosys$Site=="PBL" & ecosys$Model.Order %in% models.agb,]) +
    geom_line(aes(x=Year, y=AGB, color=Model.Order), size=2) +
    scale_x_continuous(expand=c(0,0)) +
    scale_y_continuous(name=expression(bold(paste("AGB (kg C m"^"-2", ")")))) +
    scale_color_manual(values=paste(colors.agb)) +
    guides(col=guide_legend(nrow=2, title="Model"), fill=F) +
    theme(legend.position="top") +
    # 	theme(plot.title=element_text(face="bold", size=rel(3))) + 
    theme(legend.text=element_text(size=24), 
          legend.title=element_text(size=28, face="bold"),
          legend.key=element_blank(),
          legend.key.size=unit(2, "lines")) + 
    theme(axis.line.y=element_line(color="black", size=0.5), 
          axis.line.x=element_line(color="black", size=0.5), 
          panel.grid.major=element_blank(), 
          panel.grid.minor=element_blank(), 
          panel.border=element_blank(), 
          panel.background=element_blank(), 
          axis.text.x=element_text(angle=0, color="black", size=24), 
          axis.text.y=element_text(color="black", size=24), 
          axis.title.x=element_text(face="bold", vjust=-1, size=28),  
          axis.title.y=element_text(face="bold", vjust=1, size=28))
  )
  dev.off()
}
# -------------------------------------------


# -------------------------------------------
# Map showing PBL on Prairie-Forest Boundary
# -------------------------------------------
# setwd("~/Desktop/ITRDB Map/")
# rm(list=ls())
library(ggplot2)
library(raster)
library(maps)

nbcd <- raster("~/Desktop/Research/PalEON_CR/PalEON_MIP_Site/Analyses/Change-and-Stability/raw_data/Benchmarks/NBCD_countrywide_biomass_240m_raster/NBCD_countrywide_biomass_mosaic.tif")
nbcd

# Aggregating nbcd just to make things a little faster for testing
nbcd2 <- aggregate(nbcd, fac=8) # This should make it approx 1 km res.
nbcd2

hips <- data.frame(Site = c( "PHA",  "PHO",  "PUN",  "PBL",  "PDL",  "PMB"),
                   Name = c("Harvard", "Howland", "UNDERC", "Billy's", "Demming", "Minden"),
                   lat  = c( 42.54,  45.25,  46.22,  46.28,  47.17,  43.61),
                   lon  = c(-72.18, -68.73, -89.53, -94.58, -95.17, -82.83)
                   )
                   
hips.sp <- SpatialPointsDataFrame(hips[,c("lon", "lat")], data=hips,  proj4string=CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
hips.aea <- spTransform(hips.sp, projection(nbcd))
hips <- data.frame(hips.aea)
names(hips)[names(hips) %in% c("lon.1", "lat.1")] <- c("x", "y")
hips

usa <- map_data("usa")
usa2 <- SpatialPointsDataFrame(usa[,c("long", "lat")], data=usa, proj4string=CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
usa2 <- spTransform(usa2, projection(nbcd2))
usa3 <- data.frame(usa2)
names(usa3)[names(usa3) %in% c("long.1", "lat.1")] <- c("x", "y")
summary(usa3)

states <- map_data("state")
states2 <- SpatialPointsDataFrame(states[,c("long", "lat")], data=states, proj4string=CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
states2 <- spTransform(states2, projection(nbcd2))
states3 <- data.frame(states2)
names(states3)[names(states3) %in% c("long.1", "lat.1")] <- c("x", "y")
summary(states3)


nbcd.df <- data.frame(coordinates(nbcd2))
nbcd.df$Biomass <- as.data.frame(nbcd2)[,1]
summary(nbcd.df)
dim(nbcd.df)

png("PalEON_HIPS_NBCD_PBL.png", height=9, width=16, units="in", res=180)
ggplot(data=nbcd.df[nbcd.df$Biomass>0,]) +
  coord_equal() +
  geom_polygon(data=usa3, aes(x=x, y=y, group=group), fill="white") +
  geom_raster(aes(x=x, y=y, fill=Biomass)) +
  geom_path(data=states3, aes(x=x, y=y, group=group), color="gray20") +
  geom_point(data=hips[,], aes(x=x, y=y), size=8, color="black") +
  geom_point(data=hips[hips$Site=="PBL",], aes(x=x, y=y), size=16, color="blue") +
  scale_fill_gradientn(colors=c("darkolivegreen2", "darkolivegreen4", "darkolivegreen", "darkgreen", "red")) +
  guides(fill=F) +
  # theme_bw() +
  # coord_cartesian(xlim=c(0,extent(nbcd)[2]),ylim=c(0,extent(nbcd)[4])) +
  coord_equal(ratio=1, xlim=c(-71000,2297000),ylim=c(58000,1498000), expand=F) +
  theme(axis.title=element_blank(),
        axis.text=element_blank(),
        axis.ticks=element_blank(),
        axis.line=element_blank(),
        panel.grid.major=element_blank(),
        panel.border=element_blank(),
        panel.grid.minor=element_blank(),
        panel.background=element_rect(fill="white"),
        plot.background=element_rect(fill="white"))
dev.off()


# -------------------------------------------




# -------------------------------------------
# Plot Met drivers for PBL: Tair, precipf, co2
# -------------------------------------------
{
  # Read met stability calculations
  met <- read.csv("~/Dropbox/PalEON_CR/PalEON_MIP_Site/Analyses/Change-and-Stability/Data/Met/StabilityCalcs_SiteMet_850-1850_1900-2010.csv")
  met$var <- factor(met$var, levels=c("tair", "precipf", "swdown", "lwdown", "qair", "psurf", "wind"))
  met$Site <- factor(met$Site, levels=c("PDL", "PBL", "PUN", "PMB", "PHA", "PHO"))
  # met[met$var=="tair",c("Y", "mean", "lwr", "upr", "mean.sig")] <- met[met$var=="tair",c("Y", "mean", "lwr", "upr", "mean.sig")]-273.15 # Temp to celcius
  # met[met$var=="precipf",c("Y", "mean", "lwr", "upr", "deriv.mean", "deriv.lwr", "deriv.upr", "mean.sig")] <- met[met$var=="precipf",c("Y", "mean", "lwr", "upr", "deriv.mean", "deriv.lwr", "deriv.upr", "mean.sig")]*60*60*24*365 # Precip to 
  summary(met)
  
  met.pbl <- met[met$Site=="PBL" & met$Model %in% c("tair", "precipf", "co2"),]
  met.pbl[met.pbl$var=="tair",c("Y", "mean", "lwr", "upr", "mean.sig")] <- met.pbl[met.pbl$var=="tair",c("Y", "mean", "lwr", "upr", "mean.sig")]-273.15 # Temp to celcius
  met.pbl[met.pbl$var=="precipf",c("Y", "mean", "lwr", "upr", "deriv.mean", "deriv.lwr", "deriv.upr", "mean.sig")] <- met.pbl[met.pbl$var=="precipf",c("Y", "mean", "lwr", "upr", "deriv.mean", "deriv.lwr", "deriv.upr", "mean.sig")]*60*60*24*365 # Precip to
  
  summary(met.pbl)
  
  met.pbl$Model <- factor(met.pbl$Model, levels=c("tair", "precipf", "co2"))
  levels(met.pbl$Model) <- c("Temp", "Precip", "CO2")
  # met.pbl[met.pbl$Model=="tair", ]
  
  
  png("figures/Met_Raw_PBL.png", height=9, width=15, units="in", res=180)
  print(
  ggplot(data=met.pbl) +
    facet_grid(Model~., scales="free_y") +
    geom_line(aes(x=Year, y=Y, color=Model), size=0.5) +
    geom_line(data=met.pbl[met.pbl$Model=="CO2",], aes(x=Year, y=Y, color=Model), size=5) +
    # geom_ribbon(aes(x=Year, ymin=lwr, ymax=upr), alpha=0.5) +
    # geom_line(aes(x=Year, y=mean)) +
    # geom_point(data=met.pbl[met.pbl$Model!="co2",], aes(x=Year, y=mean.sig), size=1) +
    geom_vline(xintercept=c(1850, 1900), linetype="dashed") +
    scale_color_manual(values=c("red", "blue", "green4")) +
    scale_x_continuous(expand=c(0,0)) +
    # scale_y_continuous(expand=c(0,0)) +
    theme(legend.position="top",
          legend.key=element_blank(),
          legend.key.size=unit(2, "lines"),
          legend.text=element_text(size=20),
          legend.title=element_text(size=20, face="bold")) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line=element_line(size=0.5),
          axis.ticks.length=unit(-0.5, "lines"),
          axis.text.x=element_text(color="black", size=20, margin=unit(c(1.5,0,0,0), "lines")),
          axis.text.y=element_text(color="black", size=20, margin=unit(c(0,1.5,0,0), "lines")),
          axis.title.x=element_text(color="black", size=24, face="bold"),
          axis.title.y=element_blank(),
          strip.text=element_text(size=20, face="bold"))
  )
  dev.off()
  
  met.pbl2 <- met.pbl
  for(v in met.pbl2$Model){
    met.trim <- quantile(met.pbl2[met.pbl2$Model==v, "Y"], c(0.025, 0.975), na.rm=T)
    met.pbl2[met.pbl2$Model==v, "Y.trim"] <- ifelse(met.pbl2[met.pbl2$Model==v, "Y"] >= met.trim[1] & met.pbl2[met.pbl2$Model==v, "Y"] <= met.trim[2], met.pbl2[met.pbl2$Model==v, "Y"], NA)
  }
  
  png("figures/Met_Change_PBL.png", height=9, width=15, units="in", res=180)
  print(
  ggplot(data=met.pbl2) +
    facet_grid(Model~., scales="free_y") +
    geom_line(aes(x=Year, y=Y.trim, color=Model), size=0.5, alpha=0.3) +
    geom_line(data=met.pbl2[met.pbl2$Model=="CO2",], aes(x=Year, y=Y.trim, color=Model), size=5, alpha=0.3) +
    geom_ribbon(aes(x=Year, ymin=lwr, ymax=upr), alpha=0.4) +
    geom_line(aes(x=Year, y=mean)) +
    geom_point(aes(x=Year, y=mean.sig, color=Model), size=1.5) +
    geom_vline(xintercept=c(1850, 1900), linetype="dashed") +
    scale_color_manual(values=c("red", "blue", "green4")) +
    scale_x_continuous(expand=c(0,0)) +
    # scale_y_continuous(expand=c(0,0)) +
    theme(legend.position="top",
          legend.key=element_blank(),
          legend.key.size=unit(2, "lines"),
          legend.text=element_text(size=20),
          legend.title=element_text(size=20, face="bold")) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line=element_line(size=0.5),
          axis.ticks.length=unit(-0.5, "lines"),
          axis.text.x=element_text(color="black", size=20, margin=unit(c(1.5,0,0,0), "lines")),
          axis.text.y=element_text(color="black", size=20, margin=unit(c(0,1.5,0,0), "lines")),
          axis.title.x=element_text(color="black", size=24, face="bold"),
          axis.title.y=element_blank(),
          strip.text=element_text(size=20, face="bold"))
  )
  dev.off()
}
# -------------------------------------------




# -------------------------------------------
# Set up stuff to plot GPP, AGB, & Fcomp
# -------------------------------------------
{
# Colors used for graphing
model.colors <- read.csv("~/Desktop/Research/PalEON_CR/PalEON_MIP_Site/Model.Colors.csv")
model.colors $Model.Order <- recode(model.colors$Model, "'CLM4.5-BGC'='01'; 'CLM4.5-CN'='02'; 'ED2'='03'; 'ED2-LU'='04';  'JULES-STATIC'='05'; 'JULES-TRIFFID'='06'; 'LINKAGES'='07'; 'LPJ-GUESS'='08'; 'LPJ-WSL'='09'; 'SiBCASA'='10'")
levels(model.colors$Model.Order)[1:10] <- c("CLM-BGC", "CLM-CN", "ED2", "ED2-LU", "JULES-STATIC", "JULES-TRIFFID", "LINKAGES", "LPJ-GUESS", "LPJ-WSL", "SiBCASA")
model.colors <- model.colors[order(model.colors$Model.Order),]
model.colors


gpp <- read.csv(file.path(analy.out, "StabilityCalcs_GPP_850-1850_1900-2010.csv"))
agb <- read.csv(file.path(analy.out, "StabilityCalcs_AGB_850-1850_1900-2010.csv"))
fcomp <- read.csv(file.path(analy.out, "StabilityCalcs_Fcomp_850-1850_1900-2010.csv"))

# Do some Unit conversions!
# GPP kgC/m2/s to kgC/m2/yr
sec2yr <- 60*60*24*365.25 # 60 s/min * 60 min/hr * 24 hr/day *365.25 days/yr
gpp[,c("Y", "mean", "lwr", "upr", "deriv.mean", "deriv.lwr", "deriv.upr", "mean.sig")] <- gpp[,c("Y", "mean", "lwr", "upr", "deriv.mean", "deriv.lwr", "deriv.upr", "mean.sig")]*sec2yr

# AGB & GPP kgC/m2 to MgC/HA
kgC2MgC <- 1e-3*(100*100)
gpp[,c("Y", "mean", "lwr", "upr", "deriv.mean", "deriv.lwr", "deriv.upr", "mean.sig")] <- gpp[,c("Y", "mean", "lwr", "upr", "deriv.mean", "deriv.lwr", "deriv.upr", "mean.sig")]*kgC2MgC
agb[,c("Y", "mean", "lwr", "upr", "deriv.mean", "deriv.lwr", "deriv.upr", "mean.sig")] <- agb[,c("Y", "mean", "lwr", "upr", "deriv.mean", "deriv.lwr", "deriv.upr", "mean.sig")]*kgC2MgC

# Make one mega-dataframe
ecosys.change <- rbind(gpp, agb, fcomp)

ecosys.change$var  <- factor(ecosys.change$var , levels=c("GPP", "AGB", "Fcomp"))
ecosys.change$Site <- factor(ecosys.change$Site, levels=c("PDL", "PBL", "PUN", "PMB", "PHA", "PHO"))
summary(ecosys.change)


# Aggregating across sites by year to try and find cohesive patterns to show
# standardizing everything to the climate of the 20-year spinup since that's what should be stable
ref.window = 850:869 
summary(ecosys.change)

# Making placeholders for our anomaly variables
ecosys.change[, c("Y.anom", "mean.anom", "lwr.anom", "upr.anom", "mean.sig.anom")] <- NA
for(v in unique(ecosys.change$var)){
  print(paste0(" **** ", v, " **** "))
  for(m in unique(ecosys.change$Model)){
    for(s in unique(ecosys.change$Site)){
      # print(paste0(" ---- ", s, " ---- "))
      
      # Find the reference (spinup) state
      ref <- mean(ecosys.change[ecosys.change$var==v & ecosys.change$Site==s & ecosys.change$Model==m & ecosys.change$Year %in% ref.window, "Y"], na.rm=T)
      ecosys.change[ecosys.change$var==v & ecosys.change$Model==m  & ecosys.change$Site==s, c("Y.anom", "mean.anom", "lwr.anom", "upr.anom", "mean.sig.anom")] <- ecosys.change[ecosys.change$var==v & ecosys.change$Model==m  & ecosys.change$Site==s, c("Y", "mean", "lwr", "upr", "mean.sig")] - ref
    }
  }
}
summary(ecosys.change)

# ------------------
# Aggregate ecosys.change up to the regional Model & site levels to see the best way to synthesis
# ------------------
# Showing Model consensus by site rather than model -- might give better indication of where/when models do/do not agree
ecosys.change.site <- aggregate(ecosys.change[,c("Y", "Y.anom", "mean.anom", "lwr.anom", "upr.anom", "deriv.mean", "deriv.lwr", "deriv.upr")],
                                by=ecosys.change[,c("Year", "var", "Site")],
                                FUN=mean, na.rm=T)
ecosys.change.site$Y.anom.min <- aggregate(ecosys.change[,c("Y.anom")],
                                           by=ecosys.change[,c("Year", "var", "Site")],
                                           FUN=min, na.rm=T)[,"x"]
ecosys.change.site$Y.anom.max <- aggregate(ecosys.change[,c("Y.anom")],
                                           by=ecosys.change[,c("Year", "var", "Site")],
                                           FUN=max, na.rm=T)[,"x"]
ecosys.change.site$n.sig <- aggregate(ecosys.change[,"sig"],
                                      by=ecosys.change[,c("Year", "var", "Site")],
                                      FUN=function(x){length(which(x == "*"))})[,"x"]

# Doing some temporal smoothing on ecosys.change.site for better graphs
library(zoo)
for(v in unique(ecosys.change.site$var)){
  for(s in unique(ecosys.change.site$Site)){
    ecosys.change.site[ecosys.change.site$var==v & ecosys.change.site$Site==s, "Y.anom.10"] <- rollapply(ecosys.change.site[ecosys.change.site$var==v & ecosys.change.site$Site==s, "Y.anom"], width=10, FUN=mean, fill=NA)
    ecosys.change.site[ecosys.change.site$var==v & ecosys.change.site$Site==s, "Y.anom.min.10"] <- rollapply(ecosys.change.site[ecosys.change.site$var==v & ecosys.change.site$Site==s, "Y.anom.min"], width=10, FUN=mean, fill=NA)
    ecosys.change.site[ecosys.change.site$var==v & ecosys.change.site$Site==s, "Y.anom.max.10"] <- rollapply(ecosys.change.site[ecosys.change.site$var==v & ecosys.change.site$Site==s, "Y.anom.max"], width=10, FUN=mean, fill=NA)
  }
}


head(ecosys.change.site)
summary(ecosys.change.site)
# ------------------


model.names <- data.frame(Model=unique(ecosys.change$Model), Model.Order=unique(ecosys.change$Model.Order))
model.colors2 <- model.colors[model.colors$Model.Order %in% model.names$Model.Order,]
model.colors2$Model.Order <- factor(model.colors2$Model.Order, levels=model.names$Model.Order[order(model.names$Model.Order)])  
model.colors2 <- model.colors2[order(model.colors2$Model.Order),]


colors.use <- paste(model.colors2$color)
}
# -------------------------------------------



# -------------------------------------------
# Plot GPP
# -------------------------------------------
{
  
  model.names <- data.frame(Model=unique(ecosys.change[ecosys.change$var=="GPP" & ecosys.change$Y>0 & !is.na(ecosys.change$Y),"Model"]), Model.Order=unique(ecosys.change[ecosys.change$var=="GPP" & ecosys.change$Y>0 & !is.na(ecosys.change$Y),"Model.Order"]))
  model.colors2 <- model.colors[model.colors$Model.Order %in% model.names$Model.Order,]
  model.colors2$Model.Order <- factor(model.colors2$Model.Order, levels=model.names$Model.Order[order(model.names$Model.Order)])  
  model.colors2 <- model.colors2[order(model.colors2$Model.Order),]
  colors.use <- paste(model.colors2$color)
  
  plot.raw <- ggplot(data=ecosys.change[ecosys.change$Site=="PBL" & ecosys.change$var=="GPP",]) +
    geom_line(aes(x=Year, y=Y, color=Model.Order), size=0.5) +
    geom_vline(xintercept=c(1850, 1900), linetype="dashed") +
    scale_color_manual(values=colors.use, name="Model") +
    scale_x_continuous(expand=c(0,0)) +
    scale_y_continuous(expand=c(0,0), name=expression(bold(paste("Mg C ha"^"-1"," yr"^"-1")))) +
    guides(color=guide_legend(nrow=2)) +
    # guide_legend(title="Model") +
    theme(legend.position="top",
          legend.background=element_rect(fill="white"),
          legend.key=element_blank(),
          legend.key.size=unit(1.8, "lines"),
          legend.text=element_text(size=16),
          legend.title=element_text(size=16, face="bold")) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line=element_line(size=0.5),
          axis.ticks.length=unit(-0.5, "lines"),
          # axis.text.x=element_text(color="black", size=20, margin=unit(c(1.5,0,0,0), "lines")),
          axis.text.x=element_blank(),
          axis.text.y=element_text(color="black", size=20, margin=unit(c(0,1.5,0,0), "lines")),
          # axis.title.x=element_text(color="black", size=24, face="bold"),
          axis.title.x=element_blank(),
          axis.title.y=element_text(color="black", size=24, face="bold"),
          strip.text=element_text(size=20, face="bold"))
  
  plot.dev <- ggplot(data=ecosys.change.site[ecosys.change.site$var=="GPP" & ecosys.change.site$Site == "PBL",]) +
    # facet_grid(var~Site, scales="fixed") +
    # geom_line(aes(x=Year, y=Y.anom), size=0.5, alpha=0.3) +
    geom_ribbon(aes(x=Year, ymin=Y.anom.min.10, ymax=Y.anom.max.10), alpha=0.25) +
    geom_ribbon(data=ecosys.change.site[ecosys.change.site$var=="GPP" & ecosys.change.site$Site == "PBL" & ecosys.change.site$Year<1850,], 
                aes(x=Year, ymin=lwr.anom, ymax=upr.anom), alpha=0.3) +
    geom_ribbon(data=ecosys.change.site[ecosys.change.site$var=="GPP" & ecosys.change.site$Site == "PBL" & ecosys.change.site$Year>1900,], 
                aes(x=Year, ymin=lwr.anom, ymax=upr.anom), alpha=0.3) +
    geom_line(data=ecosys.change.site[ecosys.change.site$var=="GPP" & ecosys.change.site$Site == "PBL" & ecosys.change.site$Year<1850,], 
              aes(x=Year, y=mean.anom), size=1, alpha=0.2) +
    geom_line(data=ecosys.change.site[ecosys.change.site$var=="GPP" & ecosys.change.site$Site == "PBL" & ecosys.change.site$Year>1900,], 
              aes(x=Year, y=mean.anom), size=1, alpha=0.2) +
    geom_point(aes(x=Year, y=mean.anom, size=as.factor(n.sig), color=abs(deriv.mean)), alpha=1) +
    geom_vline(xintercept=1850, linetype="dashed") +
    geom_vline(xintercept=1900, linetype="dashed") +
    geom_hline(yintercept=0, linetype="dashed") +
    scale_x_continuous(expand=c(0,0), name="Year") +
    scale_y_continuous(expand=c(0,0), name=paste0("Dev. from Spinup")) +
    scale_size_manual(values=seq(0, 5, length.out=11), breaks=seq(0,10, by=1), name="# Models") +
    scale_color_gradient(low="gray25", high="red", guide="colorbar", name=paste0("Rate of Change")) +
    guides(size=guide_legend(nrow=2)) +
    # theme(legend.key.height=unit(0.8, "lines")) +
    theme(legend.position=c(0.35, 0.85),
          legend.key=element_blank(),
          legend.key.height=unit(1.5, "lines"),
          legend.key.width=unit(2.5, "lines"),
          legend.box="horizontal",
          legend.direction="horizontal",
          # legend.key.height=unit()
          legend.text=element_text(size=14),
          legend.title=element_text(size=16, face="bold")) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line=element_line(size=0.5),
          axis.ticks.length=unit(-0.5, "lines"),
          # axis.text.x=element_blank(),
          axis.text.x=element_text(color="black", size=20, margin=unit(c(1.5,0,0,0), "lines")),
          axis.text.y=element_text(color="black", size=20, margin=unit(c(0,1.5,0,0), "lines")),
          # axis.title.x=element_blank(),
          axis.title.x=element_text(color="black", size=24, face="bold"),
          axis.title.y=element_text(color="black", size=24, face="bold"),
          strip.text=element_text(size=20, face="bold")) +
    theme(plot.margin=unit(c(0,1,1,1), "lines"))
  
  
  png(file.path(fig.dir, "GPP_Raw_Stability.png"), height=9, width=16, units="in", res=180)
  print(
    plot_grid(plot.raw + theme(plot.margin=unit(c(1,1,0.5,1), "lines")), 
              plot.dev +theme(plot.margin=unit(c(0.5,1,1,1.8), "lines")), 
              ncol=1, nrow=2)
  )
  dev.off()
}
# -------------------------------------------



# -------------------------------------------
# Plot AGB
# -------------------------------------------
{
  model.names <- data.frame(Model=unique(ecosys.change[ecosys.change$var=="AGB" & ecosys.change$Y>0 & !is.na(ecosys.change$Y),"Model"]), Model.Order=unique(ecosys.change[ecosys.change$var=="AGB" & ecosys.change$Y>0 & !is.na(ecosys.change$Y),"Model.Order"]))
  model.colors2 <- model.colors[model.colors$Model.Order %in% model.names$Model.Order,]
  model.colors2$Model.Order <- factor(model.colors2$Model.Order, levels=model.names$Model.Order[order(model.names$Model.Order)])  
  model.colors2 <- model.colors2[order(model.colors2$Model.Order),]
  colors.use <- paste(model.colors2$color)
  
  plot.raw <- ggplot(data=ecosys.change[ecosys.change$Site=="PBL" & ecosys.change$var=="AGB" & ecosys.change$Y>0 & !is.na(ecosys.change$Y),]) +
    geom_line(aes(x=Year, y=Y, color=Model.Order), size=1.5) +
    geom_vline(xintercept=c(1850, 1900), linetype="dashed") +
    scale_color_manual(values=colors.use, name="Model") +
    scale_x_continuous(expand=c(0,0)) +
    scale_y_continuous(expand=c(0,0), name=expression(bold(paste("Mg C ha"^"-1")))) +
    guides(color=guide_legend(nrow=2)) +
    # guide_legend(title="Model") +
    theme(legend.position="top",
          legend.background=element_rect(fill="white"),
          legend.key=element_blank(),
          legend.key.size=unit(1.8, "lines"),
          legend.text=element_text(size=16),
          legend.title=element_text(size=16, face="bold")) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line=element_line(size=0.5),
          axis.ticks.length=unit(-0.5, "lines"),
          # axis.text.x=element_text(color="black", size=20, margin=unit(c(1.5,0,0,0), "lines")),
          axis.text.x=element_blank(),
          axis.text.y=element_text(color="black", size=20, margin=unit(c(0,1.5,0,0), "lines")),
          # axis.title.x=element_text(color="black", size=24, face="bold"),
          axis.title.x=element_blank(),
          axis.title.y=element_text(color="black", size=24, face="bold"),
          strip.text=element_text(size=20, face="bold"))
  
  plot.dev <- ggplot(data=ecosys.change.site[ecosys.change.site$var=="AGB" & ecosys.change.site$Site == "PBL",]) +
    # facet_grid(var~Site, scales="fixed") +
    # geom_line(aes(x=Year, y=Y.anom), size=0.5, alpha=0.3) +
    geom_ribbon(aes(x=Year, ymin=Y.anom.min.10, ymax=Y.anom.max.10), alpha=0.25) +
    geom_ribbon(data=ecosys.change.site[ecosys.change.site$var=="AGB" & ecosys.change.site$Site == "PBL" & ecosys.change.site$Year<1850,], 
                aes(x=Year, ymin=lwr.anom, ymax=upr.anom), alpha=0.3) +
    geom_ribbon(data=ecosys.change.site[ecosys.change.site$var=="AGB" & ecosys.change.site$Site == "PBL" & ecosys.change.site$Year>1900,], 
                aes(x=Year, ymin=lwr.anom, ymax=upr.anom), alpha=0.3) +
    geom_line(data=ecosys.change.site[ecosys.change.site$var=="AGB" & ecosys.change.site$Site == "PBL" & ecosys.change.site$Year<1850,], 
              aes(x=Year, y=mean.anom), size=1, alpha=0.2) +
    geom_line(data=ecosys.change.site[ecosys.change.site$var=="AGB" & ecosys.change.site$Site == "PBL" & ecosys.change.site$Year>1900,], 
              aes(x=Year, y=mean.anom), size=1, alpha=0.2) +
    geom_point(aes(x=Year, y=mean.anom, size=as.factor(n.sig), color=abs(deriv.mean)), alpha=1) +
    geom_vline(xintercept=1850, linetype="dashed") +
    geom_vline(xintercept=1900, linetype="dashed") +
    geom_hline(yintercept=0, linetype="dashed") +
    scale_x_continuous(expand=c(0,0), name="Year") +
    scale_y_continuous(expand=c(0,0), name=paste0("Dev. from Spinup")) +
    scale_size_manual(values=seq(0, 5, length.out=11), breaks=seq(0,10, by=1), name="# Models") +
    scale_color_gradient(low="gray25", high="red", guide="colorbar", name=paste0("Rate of Change")) +
    guides(size=guide_legend(nrow=2)) +
    # theme(legend.key.height=unit(0.8, "lines")) +
    theme(legend.position=c(0.5, 0.1),
          legend.key=element_blank(),
          legend.key.height=unit(1.5, "lines"),
          legend.key.width=unit(2.5, "lines"),
          legend.box="horizontal",
          legend.direction="horizontal",
          # legend.key.height=unit()
          legend.text=element_text(size=14),
          legend.title=element_text(size=16, face="bold")) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line=element_line(size=0.5),
          axis.ticks.length=unit(-0.5, "lines"),
          # axis.text.x=element_blank(),
          axis.text.x=element_text(color="black", size=20, margin=unit(c(1.5,0,0,0), "lines")),
          axis.text.y=element_text(color="black", size=20, margin=unit(c(0,1.5,0,0), "lines")),
          # axis.title.x=element_blank(),
          axis.title.x=element_text(color="black", size=24, face="bold"),
          axis.title.y=element_text(color="black", size=24, face="bold"),
          strip.text=element_text(size=20, face="bold")) +
    theme(plot.margin=unit(c(0,1,1,1), "lines"))
  
  
  png(file.path(fig.dir, "AGB_Raw_Stability.png"), height=9, width=16, units="in", res=180)
  print(
    plot_grid(plot.raw + theme(plot.margin=unit(c(1,1,0.5,1), "lines")), 
              plot.dev +theme(plot.margin=unit(c(0.5,1,1,2.1), "lines")), 
              ncol=1, nrow=2)
  )
  dev.off()
}
# -------------------------------------------


# -------------------------------------------
# Plot Fcomp
# -------------------------------------------
{
  model.names <- data.frame(Model=unique(ecosys.change[ecosys.change$var=="Fcomp" & ecosys.change$Y>0 & !is.na(ecosys.change$Y),"Model"]), Model.Order=unique(ecosys.change[ecosys.change$var=="Fcomp" & ecosys.change$Y>0 & !is.na(ecosys.change$Y),"Model.Order"]))
  model.colors2 <- model.colors[model.colors$Model.Order %in% model.names$Model.Order,]
  model.colors2$Model.Order <- factor(model.colors2$Model.Order, levels=model.names$Model.Order[order(model.names$Model.Order)])  
  model.colors2 <- model.colors2[order(model.colors2$Model.Order),]
  colors.use <- paste(model.colors2$color)
  
  
  plot.raw <- ggplot(data=ecosys.change[ecosys.change$Site=="PBL" & ecosys.change$var=="Fcomp",]) +
    geom_line(aes(x=Year, y=Y*100, color=Model.Order), size=1.25) +
    geom_vline(xintercept=c(1850, 1900), linetype="dashed") +
    scale_color_manual(values=colors.use, name="Model") +
    scale_x_continuous(expand=c(0,0)) +
    scale_y_continuous(expand=c(0,0), name=expression(bold(paste("%")))) +
    guides(color=guide_legend(nrow=2)) +
    # guide_legend(title="Model") +
    theme(legend.position="top",
          legend.background=element_rect(fill="white"),
          legend.key=element_blank(),
          legend.key.size=unit(1.8, "lines"),
          legend.text=element_text(size=16),
          legend.title=element_text(size=16, face="bold")) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line=element_line(size=0.5),
          axis.ticks.length=unit(-0.5, "lines"),
          # axis.text.x=element_text(color="black", size=20, margin=unit(c(1.5,0,0,0), "lines")),
          axis.text.x=element_blank(),
          axis.text.y=element_text(color="black", size=20, margin=unit(c(0,1.5,0,0), "lines")),
          # axis.title.x=element_text(color="black", size=24, face="bold"),
          axis.title.x=element_blank(),
          axis.title.y=element_text(color="black", size=24, face="bold"),
          strip.text=element_text(size=20, face="bold"))
  
  plot.dev <- ggplot(data=ecosys.change.site[ecosys.change.site$var=="Fcomp" & ecosys.change.site$Site == "PBL",]) +
    # facet_grid(var~Site, scales="fixed") +
    # geom_line(aes(x=Year, y=Y.anom), size=0.5, alpha=0.3) +
    geom_ribbon(aes(x=Year, ymin=Y.anom.min.10*100, ymax=Y.anom.max.10*100), alpha=0.25) +
    geom_ribbon(data=ecosys.change.site[ecosys.change.site$var=="Fcomp" & ecosys.change.site$Site == "PBL" & ecosys.change.site$Year<1850,], 
                aes(x=Year, ymin=lwr.anom*100, ymax=upr.anom*100), alpha=0.3) +
    geom_ribbon(data=ecosys.change.site[ecosys.change.site$var=="Fcomp" & ecosys.change.site$Site == "PBL" & ecosys.change.site$Year>1900,], 
                aes(x=Year, ymin=lwr.anom*100, ymax=upr.anom*100), alpha=0.3) +
    geom_line(data=ecosys.change.site[ecosys.change.site$var=="Fcomp" & ecosys.change.site$Site == "PBL" & ecosys.change.site$Year<1850,], 
              aes(x=Year, y=mean.anom*100), size=1, alpha=0.2) +
    geom_line(data=ecosys.change.site[ecosys.change.site$var=="Fcomp" & ecosys.change.site$Site == "PBL" & ecosys.change.site$Year>1900,], 
              aes(x=Year, y=mean.anom*100), size=1, alpha=0.2) +
    geom_point(aes(x=Year, y=mean.anom*100, size=as.factor(n.sig), color=abs(deriv.mean)*100), alpha=1) +
    geom_vline(xintercept=1850, linetype="dashed") +
    geom_vline(xintercept=1900, linetype="dashed") +
    geom_hline(yintercept=0, linetype="dashed") +
    scale_x_continuous(expand=c(0,0), name="Year") +
    scale_y_continuous(expand=c(0,0), name=paste0("Dev. from Spinup")) +
    scale_size_manual(values=seq(0, 5, length.out=11), breaks=seq(0,10, by=1), name="# Models") +
    scale_color_gradient(low="gray25", high="red", guide="colorbar", name=paste0("Rate of Change")) +
    guides(size=guide_legend(nrow=2)) +
    # theme(legend.key.height=unit(0.8, "lines")) +
    theme(legend.position=c(0.3, 0.9),
          legend.key=element_blank(),
          legend.key.height=unit(1.5, "lines"),
          legend.key.width=unit(2.5, "lines"),
          legend.box="horizontal",
          legend.direction="horizontal",
          # legend.key.height=unit()
          legend.text=element_text(size=14),
          legend.title=element_text(size=16, face="bold")) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line=element_line(size=0.5),
          axis.ticks.length=unit(-0.5, "lines"),
          # axis.text.x=element_blank(),
          axis.text.x=element_text(color="black", size=20, margin=unit(c(1.5,0,0,0), "lines")),
          axis.text.y=element_text(color="black", size=20, margin=unit(c(0,1.5,0,0), "lines")),
          # axis.title.x=element_blank(),
          axis.title.x=element_text(color="black", size=24, face="bold"),
          axis.title.y=element_text(color="black", size=24, face="bold"),
          strip.text=element_text(size=20, face="bold")) +
    theme(plot.margin=unit(c(0,1,1,1), "lines"))
  
  
  png(file.path(fig.dir, "Fcomp_Raw_Stability.png"), height=9, width=16, units="in", res=180)
  print(
    plot_grid(plot.raw + theme(plot.margin=unit(c(1,1,0.5,1), "lines")), 
              plot.dev +theme(plot.margin=unit(c(0.5,1,1,1), "lines")), 
              ncol=1, nrow=2)
  )
  dev.off()
}
# -------------------------------------------


# -------------------------------------------
# Plot Tair & Stabilities together
# -------------------------------------------
tair.dev  <- ggplot(data=met.pbl2[met.pbl2$Model=="Temp",]) +
  facet_grid(Model~., scales="free_y") +
  geom_line(aes(x=Year, y=Y.trim-273.15, color=Model), size=0.5, alpha=0.3, color="red") +
  geom_ribbon(aes(x=Year, ymin=lwr-273.15, ymax=upr-273.15), alpha=0.4) +
  geom_line(aes(x=Year, y=mean-273.15)) +
  geom_point(aes(x=Year, y=mean.sig-273.15, color=Model), size=1.5, color="red") +
  geom_vline(xintercept=c(1850, 1900), linetype="dashed") +
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0), name="deg. C") +
  theme(legend.position="top",
        legend.key=element_blank(),
        legend.key.size=unit(2, "lines"),
        legend.text=element_text(size=20),
        legend.title=element_text(size=20, face="bold")) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line=element_line(size=0.5),
        axis.ticks.length=unit(-0.5, "lines"),
        # axis.text.x=element_text(color="black", size=20, margin=unit(c(1.5,0,0,0), "lines")),
        axis.text.x=element_blank(),
        axis.text.y=element_text(color="black", size=20, margin=unit(c(0,1.5,0,0), "lines")),
        # axis.title.x=element_text(color="black", size=24, face="bold"),
        axis.title.x=element_blank(),
        axis.title.y=element_text(color="black", size=24, face="bold"),
        strip.text=element_text(size=20, face="bold"))
gpp.dev   <- ggplot(data=ecosys.change.site[ecosys.change.site$var=="GPP" & ecosys.change.site$Site == "PBL",]) +
  facet_grid(var~., scales="fixed") +
  # geom_line(aes(x=Year, y=Y.anom), size=0.5, alpha=0.3) +
  geom_ribbon(aes(x=Year, ymin=Y.anom.min.10, ymax=Y.anom.max.10), alpha=0.25) +
  geom_ribbon(data=ecosys.change.site[ecosys.change.site$var=="GPP" & ecosys.change.site$Site == "PBL" & ecosys.change.site$Year<1850,], 
              aes(x=Year, ymin=lwr.anom, ymax=upr.anom), alpha=0.3) +
  geom_ribbon(data=ecosys.change.site[ecosys.change.site$var=="GPP" & ecosys.change.site$Site == "PBL" & ecosys.change.site$Year>1900,], 
              aes(x=Year, ymin=lwr.anom, ymax=upr.anom), alpha=0.3) +
  geom_line(data=ecosys.change.site[ecosys.change.site$var=="GPP" & ecosys.change.site$Site == "PBL" & ecosys.change.site$Year<1850,], 
            aes(x=Year, y=mean.anom), size=1, alpha=0.2) +
  geom_line(data=ecosys.change.site[ecosys.change.site$var=="GPP" & ecosys.change.site$Site == "PBL" & ecosys.change.site$Year>1900,], 
            aes(x=Year, y=mean.anom), size=1, alpha=0.2) +
  geom_point(aes(x=Year, y=mean.anom, size=as.factor(n.sig), color=abs(deriv.mean)), alpha=1) +
  geom_vline(xintercept=1850, linetype="dashed") +
  geom_vline(xintercept=1900, linetype="dashed") +
  geom_hline(yintercept=0, linetype="dashed") +
  scale_x_continuous(expand=c(0,0), name="Year") +
  scale_y_continuous(expand=c(0,0), name=paste0("dev. GPP")) +
  scale_size_manual(values=seq(0, 5, length.out=11), breaks=seq(0,10, by=1), name="# Models") +
  scale_color_gradient(low="gray25", high="red", guide="colorbar", name=paste0("Rate of Change")) +
  guides(size=F, colorbar=F, fill=F, color=F) +
  # theme(legend.key.height=unit(0.8, "lines")) +
  theme(legend.position=c(0.35, 0.85),
        legend.key=element_blank(),
        legend.key.height=unit(1.5, "lines"),
        legend.key.width=unit(2.5, "lines"),
        legend.box="horizontal",
        legend.direction="horizontal",
        # legend.key.height=unit()
        legend.text=element_text(size=14),
        legend.title=element_text(size=16, face="bold")) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line=element_line(size=0.5),
        axis.ticks.length=unit(-0.5, "lines"),
        # axis.text.x=element_text(color="black", size=20, margin=unit(c(1.5,0,0,0), "lines")),
        axis.text.x=element_blank(),
        axis.text.y=element_text(color="black", size=20, margin=unit(c(0,1.5,0,0), "lines")),
        # axis.title.x=element_text(color="black", size=24, face="bold"),
        axis.title.x=element_blank(),
        axis.title.y=element_text(color="black", size=24, face="bold"),
        strip.text=element_text(size=20, face="bold")) +
  theme(plot.margin=unit(c(0,1,1,1), "lines")) 
agb.dev   <- ggplot(data=ecosys.change.site[ecosys.change.site$var=="AGB" & ecosys.change.site$Site == "PBL",]) +
  facet_grid(var~., scales="fixed") +
  # geom_line(aes(x=Year, y=Y.anom), size=0.5, alpha=0.3) +
  geom_ribbon(aes(x=Year, ymin=Y.anom.min.10, ymax=Y.anom.max.10), alpha=0.25) +
  geom_ribbon(data=ecosys.change.site[ecosys.change.site$var=="AGB" & ecosys.change.site$Site == "PBL" & ecosys.change.site$Year<1850,], 
              aes(x=Year, ymin=lwr.anom, ymax=upr.anom), alpha=0.3) +
  geom_ribbon(data=ecosys.change.site[ecosys.change.site$var=="AGB" & ecosys.change.site$Site == "PBL" & ecosys.change.site$Year>1900,], 
              aes(x=Year, ymin=lwr.anom, ymax=upr.anom), alpha=0.3) +
  geom_line(data=ecosys.change.site[ecosys.change.site$var=="AGB" & ecosys.change.site$Site == "PBL" & ecosys.change.site$Year<1850,], 
            aes(x=Year, y=mean.anom), size=1, alpha=0.2) +
  geom_line(data=ecosys.change.site[ecosys.change.site$var=="AGB" & ecosys.change.site$Site == "PBL" & ecosys.change.site$Year>1900,], 
            aes(x=Year, y=mean.anom), size=1, alpha=0.2) +
  geom_point(aes(x=Year, y=mean.anom, size=as.factor(n.sig), color=abs(deriv.mean)), alpha=1) +
  geom_vline(xintercept=1850, linetype="dashed") +
  geom_vline(xintercept=1900, linetype="dashed") +
  geom_hline(yintercept=0, linetype="dashed") +
  scale_x_continuous(expand=c(0,0), name="Year") +
  scale_y_continuous(expand=c(0,0), name=paste0("dev. AGB")) +
  scale_size_manual(values=seq(0, 5, length.out=11), breaks=seq(0,10, by=1), name="# Models") +
  scale_color_gradient(low="gray25", high="red", guide="colorbar", name=paste0("Rate of Change")) +
  guides(size=F, colorbar=F, fill=F, color=F) +
  # theme(legend.key.height=unit(0.8, "lines")) +
  theme(legend.position=c(0.5, 0.1),
        legend.key=element_blank(),
        legend.key.height=unit(1.5, "lines"),
        legend.key.width=unit(2.5, "lines"),
        legend.box="horizontal",
        legend.direction="horizontal",
        # legend.key.height=unit()
        legend.text=element_text(size=14),
        legend.title=element_text(size=16, face="bold")) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line=element_line(size=0.5),
          axis.ticks.length=unit(-0.5, "lines"),
          # axis.text.x=element_text(color="black", size=20, margin=unit(c(1.5,0,0,0), "lines")),
          axis.text.x=element_blank(),
          axis.text.y=element_text(color="black", size=20, margin=unit(c(0,1.5,0,0), "lines")),
          # axis.title.x=element_text(color="black", size=24, face="bold"),
          axis.title.x=element_blank(),
          axis.title.y=element_text(color="black", size=24, face="bold"),
          strip.text=element_text(size=20, face="bold")) +
    theme(plot.margin=unit(c(0,1,1,1), "lines")) 
fcomp.dev <- ggplot(data=ecosys.change.site[ecosys.change.site$var=="Fcomp" & ecosys.change.site$Site == "PBL",]) +
  facet_grid(var~., scales="fixed") +
  # geom_line(aes(x=Year, y=Y.anom), size=0.5, alpha=0.3) +
  geom_ribbon(aes(x=Year, ymin=Y.anom.min.10*100, ymax=Y.anom.max.10*100), alpha=0.25) +
  geom_ribbon(data=ecosys.change.site[ecosys.change.site$var=="Fcomp" & ecosys.change.site$Site == "PBL" & ecosys.change.site$Year<1850,], 
              aes(x=Year, ymin=lwr.anom*100, ymax=upr.anom*100), alpha=0.3) +
  geom_ribbon(data=ecosys.change.site[ecosys.change.site$var=="Fcomp" & ecosys.change.site$Site == "PBL" & ecosys.change.site$Year>1900,], 
              aes(x=Year, ymin=lwr.anom*100, ymax=upr.anom*100), alpha=0.3) +
  geom_line(data=ecosys.change.site[ecosys.change.site$var=="Fcomp" & ecosys.change.site$Site == "PBL" & ecosys.change.site$Year<1850,], 
            aes(x=Year, y=mean.anom*100), size=1, alpha=0.2) +
  geom_line(data=ecosys.change.site[ecosys.change.site$var=="Fcomp" & ecosys.change.site$Site == "PBL" & ecosys.change.site$Year>1900,], 
            aes(x=Year, y=mean.anom*100), size=1, alpha=0.2) +
  geom_point(aes(x=Year, y=mean.anom*100, size=as.factor(n.sig), color=abs(deriv.mean)*100), alpha=1) +
  geom_vline(xintercept=1850, linetype="dashed") +
  geom_vline(xintercept=1900, linetype="dashed") +
  geom_hline(yintercept=0, linetype="dashed") +
  scale_x_continuous(expand=c(0,0), name="Year") +
  scale_y_continuous(expand=c(0,0), name=paste0("dev. Fcomp")) +
  scale_size_manual(values=seq(0, 5, length.out=11), breaks=seq(0,10, by=1), name="# Models") +
  scale_color_gradient(low="gray25", high="red", guide="colorbar", name=paste0("Rate of Change")) +
  guides(size=F, colorbar=F, fill=F, color=F) +
  # theme(legend.key.height=unit(0.8, "lines")) +
  theme(legend.position=c(0.3, 0.9),
        legend.key=element_blank(),
        legend.key.height=unit(1.5, "lines"),
        legend.key.width=unit(2.5, "lines"),
        legend.box="horizontal",
        legend.direction="horizontal",
        # legend.key.height=unit()
        legend.text=element_text(size=14),
        legend.title=element_text(size=16, face="bold")) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line=element_line(size=0.5),
        axis.ticks.length=unit(-0.5, "lines"),
        # axis.text.x=element_blank(),
        axis.text.x=element_text(color="black", size=20, margin=unit(c(1.5,0,0,0), "lines")),
        axis.text.y=element_text(color="black", size=20, margin=unit(c(0,1.5,0,0), "lines")),
        # axis.title.x=element_blank(),
        axis.title.x=element_text(color="black", size=24, face="bold"),
        axis.title.y=element_text(color="black", size=24, face="bold"),
        strip.text=element_text(size=20, face="bold")) +
  theme(plot.margin=unit(c(0,1,1,1), "lines"))


png(file.path(fig.dir, "Stability_Comparison.png"), height=9, width=16, units="in", res=180)
print(
  plot_grid(tair.dev + theme(plot.margin=unit(c(1,1,0.5,1.7), "lines")),
            gpp.dev + theme(plot.margin=unit(c(0.5,1,0.5,1.1), "lines")),
            agb.dev + theme(plot.margin=unit(c(0.5,1,0.5,0.3), "lines")), 
            fcomp.dev +theme(plot.margin=unit(c(0.5,1,1,1), "lines")), 
            ncol=1, nrow=4, rel_heights=c(1,1,1,1.5))
)
dev.off()

# -------------------------------------------
