# -------------------------------------------
# Extracting & formatting benchmarks from various sources 
#
# Author: Christy Rollinson, crollinson@gmail.com
#
# ------------------
# Objectives:
# ------------------
# A. Get benchmarks into a common, usuable format specific to
#    individual sites & models PFTs when needed (composition benchmarks)
# ------------------
#
#
# ------------------
# Workflow
# ------------------
# 1. Set up file structure, libraries, etc.
#
# 2. Extract & Format Benchmarks
#    2.01. Composition -- STEPPS 1 (contibutor/contact: Andria)
#    2.02. Composition -- SetVeg (contributor/contact: Simon? Jody?)
#    2.03. Biomass -- SetVeg  (contributor/contact: Jody? Kelly?)
#    2.04. Biomass -- NBCD (contributor/contact: Christy)
#    2.05. Carbon Fluxes - Flux Towers (contributor/contact: Dave)
#    2.06. MODIS - LAI, GPP (contributor/contact: Bethany)
#    2.07. Tree RIngs -- AGB, AGBincrement (contributor/contact: Ross, Alex)
# ------------------
# -------------------------------------------
rm(list=ls())

# -------------------------------------------
# 1. Set up file structure, libraries, etc.
# -------------------------------------------
library(car)
library(ggplot2)
# library(mgcv)
# setwd("~/Desktop/Research/PalEON_CR/PalEON_MIP_Site/Analyses/Change-and-Stability") # Path to this project github repository: https://github.com/PalEON-Project/MIP-Change-and-Stability.git
setwd("~/Dropbox/PalEON_CR/PalEON_MIP_Site/Analyses/Change-and-Stability") # Path to this project github repository: https://github.com/PalEON-Project/MIP-Change-and-Stability.git
# path.gamm.func <- "~/Desktop/R_Functions/"  # Path to github repository of my GAMM helper functions: https://github.com/crollinson/R_Functions.git
inputs    <- "Data/" # Path to my cleaned model output

mip.utils <- "~/Desktop/Research/PalEON_CR/MIP_Utils/" # Path to PalEON MIP Utility repository: https://github.com/PalEON-Project/MIP_Utils.git

out.dir <- "Data/Benchmarking" # Path to where the analysis output should go
fig.dir <- "Figures/Benchmarking" # Path to where figures should go
raw.dir <- "raw_data/Benchmarks/"

if(!dir.exists(out.dir)) dir.create(out.dir)
if(!dir.exists(fig.dir)) dir.create(fig.dir)

hips <- data.frame(Site = c( "PHA",  "PHO",  "PUN",  "PBL",  "PDL",  "PMB"),
                   lat  = c( 42.54,  45.25,  46.22,  46.28,  47.17,  43.61),
                   lon  = c(-72.18, -68.73, -89.53, -94.58, -95.17, -82.83),
                   latmin = c( 42.50,  45.00,  46.00,  46.00,  47.00,  43.50),
                   latmax = c( 43.00,  45.50,  46.50,  46.50,  47.50,  44.00),
                   lonmin = c(-72.50, -69.00, -90.00, -95.00, -95.50, -83.00),
                   lonmax = c(-72.00, -68.50, -89.50, -94.50, -95.00, -82.50))

# -------------------------------------------

# -------------------------------------------
# 2. Extract & Format Benchmarks
# -------------------------------------------
{
  # ------------------
  #    2.01. Composition -- STEPPS 1 & 2
  #     -- Finding the dominant PFT by model specifications
  # ------------------
  {
    library(car); library(reshape)
    # Read in the model-STEPPS crosswalk
    pft.xwalk <- read.csv("raw_data/STEPPS_conversion.csv")
    pft.xwalk$taxon <- toupper(pft.xwalk$STEPPS)
    pft.xwalk$taxon <- recode(pft.xwalk$taxon, "'OTHER CONIFER'='OTHER.CONIFER'; 'OTHER HARDWOOD'='OTHER.HARDWOOD'")
    
    # Reading in STEPPS 1
    
    # --------
    # Loading STEPPS 2
    # --------
    stepps2 <- readRDS(file.path(raw.dir, "STEPPS2", "r_hips_v1.0.RDS"))
    stepps2$taxon <- as.factor(stepps2$taxon)
    summary(stepps2)
    
    unique(stepps2$taxon)
    
    # Translating the STEPPS data into model PFTS & storing in a list
    stepps2 <- merge(stepps2, pft.xwalk)
    
    # Changing stepps time to a calendar year & subsetting just that in the modeling temporal extent
    stepps2$year <- 2000-stepps2$time*100
    stepps2 <- stepps2[stepps2$year>=850,]
    
    stepps2.pft <- list()
    for(m in names(pft.xwalk)[2:(ncol(pft.xwalk)-1)]){
      print(m)
      if(m=="SiBCASA") next
      stepps2$pft.tmp <- stepps2[,m]
      
      # Re-aggregate to the model-specific PFTs
      fcomp1 <- aggregate(stepps2[,c("prop")], by=stepps2[,c("site", "year", "iter", m)], FUN=sum)
      
      stepps2.pft[[m]] <- merge(data.frame(year=unique(stepps2$year)), data.frame(site=unique(stepps2$site)))
      # Go through each site & each iteration to find which is the dominant PFT on average
      for(s in unique(fcomp1$site)){
        # Finding out which PFT has the highest fcomp across time & iterations
        pft.avg <- aggregate(fcomp1[fcomp1$site==s, "x"], by=list(fcomp1[fcomp1$site==s, m]), FUN=mean)
        names(pft.avg) <- c("pft", "fcomp.avg")
        pft.dom <- pft.avg$pft[which(pft.avg$fcomp.avg == max(pft.avg$fcomp.avg) )]
        
        # Record the dominant pft
        stepps2.pft[[m]][stepps2.pft[[m]]$site==s,"pft"] <- pft.dom
        
        # creating an array to store the output 
        # dat.temp <- array(dim=c(length(unique(stepps2$time)), length(unique(stepps2$iter))))
        dat.temp <- fcomp1[fcomp1[,m]==pft.dom & fcomp1$site==s,c("x", "iter", "year")]
        dat.temp$iter <- as.ordered(dat.temp$iter)
        dat.temp$year <- as.ordered(dat.temp$year)
        
        dat.temp <- recast(dat.temp, year~iter, fun.aggregate=mean)
        
        stepps2.pft[[m]][stepps2.pft[[m]]$site==s,"mean"] <- apply(dat.temp, 1, mean)
        stepps2.pft[[m]][stepps2.pft[[m]]$site==s,"lwr"] <- apply(dat.temp, 1, quantile, 0.025)
        stepps2.pft[[m]][stepps2.pft[[m]]$site==s,"upr"] <- apply(dat.temp, 1, quantile, 0.975)
        
      }
    }
    
    pdf(file.path(fig.dir, "STEPPS2_ModelPFTs.pdf"))
    for(m in names(stepps2.pft)){
      print(
        ggplot(data=stepps2.pft[[m]]) +
          geom_ribbon(aes(x=year, ymin=lwr, ymax=upr, fill=site), alpha=0.4) +
          geom_line(aes(x=year, y=mean, color=site, linetype=pft), size=1.5) +
          ggtitle(m) +
          theme_bw()
      )
    }
    dev.off()
    
    save(stepps2.pft, file=file.path(out.dir, "STEPPS2_ModelPFTs.RData"))
    # --------
  }
  # ------------------
  
  # ------------------
  #    2.02. Composition -- SetVeg
  # ------------------
  {
    library(raster); library(rgdal); library(ggplot2)
    library(ncdf4)
    library(car); library(reshape)
    # Read in the model-STEPPS crosswalk
    pft.xwalk <- read.csv("raw_data/SetVeg_conversion.csv")
    pft.xwalk$taxon <- toupper(pft.xwalk$SetVeg)
    pft.xwalk$taxon <- recode(pft.xwalk$taxon, "'OTHER CONIFER'='OTHER.CONIFER'; 'OTHER HARDWOOD'='OTHER.HARDWOOD'")
    
    # --------
    # Loading SetVeg Data & Transforming to Lat/Lon
    # 1. Load data
    # 2. Transform to lat/lon
    # 3. For each iteration take mean for spatial area for each site(properly ascribe uncertainty)
    # 4. Extract site-level for each taxon to match STEPPS2 formatting
    # --------
    proj.albers <- "+proj=aea +lat_1=42.122774 +lat_2=49.01518 +lat_0=45.568977 +lon_0=-83.248627 +x_0=1000000 +y_0=1000000 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"
    
    sv.nc <- nc_open(file.path(raw.dir, "SetVeg_Stats", "composition_v0.4.nc"))
    x <- ncvar_get(sv.nc, "x")
    y <- ncvar_get(sv.nc, "y")
    taxa <- names(sv.nc$var)
  
    # ---------------
    # Creating a lat/lon conversion for the setveg data
    # ---------------
    sv.tmp <-  data.frame(x=rep(x, length(y)), y=rep(y, each=length(x)))
    sv.tmp$DUMMY <- stack(data.frame(ncvar_get(sv.nc, taxa[1])[,,1]))[,1]
    summary(sv.tmp)
    
    setveg.spat <- sv.tmp
    coordinates(setveg.spat) <- c("x", "y")
    projection(setveg.spat) <- proj.albers
    
    setveg.spat <- spTransform(setveg.spat, CRS('+proj=longlat'))
    setveg.longlat <- data.frame(setveg.spat)
    summary(setveg.longlat)
    
    # Add Lon/Lat to the sv.tmp data frame so we have a cross walk
    sv.tmp[,c("lon", "lat")] <- setveg.longlat[,c("x", "y")]
    summary(sv.tmp)
    
    # Adding Albers references to the hips dataframe
    for(i in 1:nrow(hips)){
      hips[i,"x"] <- mean(sv.tmp[which(sv.tmp$lon>=hips[i,"lonmin"] & sv.tmp$lat>=hips[i,"latmin"] & sv.tmp$lon<=hips[i,"lonmax"] & sv.tmp$lat<=hips[i,"latmax"]), "x"])
      hips[i,"y"] <- mean(sv.tmp[which(sv.tmp$lon>=hips[i,"lonmin"] & sv.tmp$lat>=hips[i,"latmin"] & sv.tmp$lon<=hips[i,"lonmax"] & sv.tmp$lat<=hips[i,"latmax"]), "y"])
      hips[i,"xmin"] <- min(sv.tmp[which(sv.tmp$lon>=hips[i,"lonmin"] & sv.tmp$lat>=hips[i,"latmin"] & sv.tmp$lon<=hips[i,"lonmax"] & sv.tmp$lat<=hips[i,"latmax"]), "x"]) # SW corner x
      hips[i,"ymin"] <- min(sv.tmp[which(sv.tmp$lon>=hips[i,"lonmin"] & sv.tmp$lat>=hips[i,"latmin"] & sv.tmp$lon<=hips[i,"lonmax"] & sv.tmp$lat<=hips[i,"latmax"]), "y"]) # SW corner y
      hips[i,"xmax"] <- max(sv.tmp[which(sv.tmp$lon>=hips[i,"lonmin"] & sv.tmp$lat>=hips[i,"latmin"] & sv.tmp$lon<=hips[i,"lonmax"] & sv.tmp$lat<=hips[i,"latmax"]), "x"]) # NE corner x
      hips[i,"ymax"] <- max(sv.tmp[which(sv.tmp$lon>=hips[i,"lonmin"] & sv.tmp$lat>=hips[i,"latmin"] & sv.tmp$lon<=hips[i,"lonmax"] & sv.tmp$lat<=hips[i,"latmax"]), "y"]) # NE corner x
    }
    
    # Figures to make sure the spatial transformation and my indexing works
    png(file.path(fig.dir, "SetVeg_Composition_Albers.png"), height=6, width=10, units="in", res=220)
    print(
    ggplot(data=sv.tmp) +
      coord_equal() +
      geom_tile(aes(x=x, y=y, fill=DUMMY)) +
      geom_rect(data=hips, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), color="red", fill=NA) +
      geom_point(data=hips, aes(x=x, y=y), color="red")
    )
    dev.off()
  
    png(file.path(fig.dir, "SetVeg_Composition_LonLat.png"), height=6, width=10, units="in", res=220)
    print(
    ggplot(data=sv.tmp) +
      coord_equal() +
      geom_point(aes(x=lon, y=lat, color=DUMMY)) +
      geom_rect(data=hips, aes(xmin=lonmin, xmax=lonmax, ymin=latmin, ymax=latmax), color="red", fill=NA) +
      geom_point(data=hips, aes(x=lon, y=lat), color="red")
    )
    dev.off()
  
    # ---------------
    
    setveg=NULL # set the set veg to something empty to make life easier
    names(sv.nc$var)
    for(v in names(sv.nc$var)){
      print(v)
      dat.tmp <- ncvar_get(sv.nc, v)
      iters <- dim(dat.tmp)[3]
      for(s in 1:nrow(hips)){
        # print(paste0(" -- ", hips[s,"Site"]))
        xind <- which(x>=hips[s,"xmin"] & x<=hips[s,"xmax"] ) # which rows we want
        yind <- which(y>=hips[s,"ymin"] & y<=hips[s,"ymax"] ) # which columns we want
        pft.tmp <- data.frame(Site=hips[s,"Site"], iter=1:iters, taxon=as.factor(toupper(v)), prop=apply(dat.tmp[xind,yind,], 3, mean, na.rm=T))
        
        # This is the clunky way of adding in the data, but it works
        if(is.null(setveg)){
          setveg <- pft.tmp
        } else {
          setveg <- rbind(setveg, pft.tmp)
        }
      } # end site loop
    } # end variable loop
    
    # If something has NAs (because outside of sample range), fill it with 0
    setveg[is.na(setveg$prop), "prop"] <- 0
    
    summary(setveg)
    # ---------------
    # Translating the SetVeg data into model PFTS & storing in a list
    # Using the mean 1800-1850 for the models
    # ---------------
    setveg <- merge(setveg, pft.xwalk, all.x=T)
    summary(setveg)
      
    setveg.pft <- NULL
    for(m in names(pft.xwalk)[2:(ncol(pft.xwalk)-1)]){
      print(m)
      if(m=="SiBCASA") next
      setveg$pft.tmp <- setveg[,m]
      
      # Re-aggregate to the model-specific PFTs
      fcomp1 <- aggregate(setveg[,c("prop")], by=setveg[,c("Site", "iter", m)], FUN=sum)
      
      tmp <- data.frame(Model=m, Site=unique(setveg$Site))
      # Go through each site & each iteration to find which is the dominant PFT on average
      for(s in unique(fcomp1$Site)){
        # Finding out which PFT has the highest fcomp across time & iterations
        pft.avg <- aggregate(fcomp1[fcomp1$Site==s, "x"], by=list(fcomp1[fcomp1$Site==s, m]), FUN=mean)
        names(pft.avg) <- c("pft", "fcomp.avg")
        pft.dom <- pft.avg$pft[which(pft.avg$fcomp.avg == max(pft.avg$fcomp.avg) )]
        
        # Record the dominant pft
        tmp[tmp$Site==s,"pft"] <- pft.dom
        
        # creating an array to store the output 
        # dat.temp <- array(dim=c(length(unique(setveg$time)), length(unique(setveg$iter))))
        dat.temp <- fcomp1[fcomp1[,m]==pft.dom & fcomp1$Site==s,c("x", "iter")]
              
        tmp[tmp$Site==s,"mean"] <- mean(dat.temp$x)
        tmp[tmp$Site==s,"lwr"] <- quantile(dat.temp$x, 0.025)
        tmp[tmp$Site==s,"upr"] <- quantile(dat.temp$x, 0.975)    
      } # end site loop
      if(is.null(setveg.pft)){
        setveg.pft <- tmp
      } else {
        setveg.pft <- rbind(setveg.pft, tmp)
      }
    } # End Model
    
  
    png(file.path(fig.dir, "SetVeg_ModelPFTs.png"))
    print(
    ggplot(data=setveg.pft) +
      facet_wrap(~Model) +
      geom_pointrange(aes(x=Site, y=mean, ymin=lwr, ymax=upr, color=pft), size=1.5) +
      theme_bw()
    )
    dev.off()
    
    write.csv(setveg.pft, file=file.path(out.dir, "SetVeg_ModelPFTs.csv"), row.names=F)
    # ---------------
    
  }
  # ------------------
  
  # ------------------
  #    2.03. Biomass -- SetVeg
  # ------------------
  {
    
    library(raster); library(rgdal); library(ggplot2)
    library(ncdf4)
    library(car); library(reshape)
    # Read in the model-STEPPS crosswalk
    pft.xwalk <- read.csv("raw_data/SetVeg_conversion.csv")
    pft.xwalk$taxon <- toupper(pft.xwalk$SetVeg)
    pft.xwalk$taxon <- recode(pft.xwalk$taxon, "'OTHER CONIFER'='OTHER.CONIFER'; 'OTHER HARDWOOD'='OTHER.HARDWOOD'")
    
    # --------
    # Loading SetVeg Data & Transforming to Lat/Lon
    # 1. Load data
    # 2. Transform to lat/lon
    # 3. For each iteration take mean for spatial area for each site(properly ascribe uncertainty)
    # 4. Extract site-level for each taxon to match STEPPS2 formatting
    # --------
    proj.albers <- "+proj=aea +lat_1=42.122774 +lat_2=49.01518 +lat_0=45.568977 +lon_0=-83.248627 +x_0=1000000 +y_0=1000000 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"
    
    # # Biomass
    sv.biom.mean <- read.csv(file.path(raw.dir, "SetVeg_Stats", "biomass_prediction_v0.9-7_bam.csv"))
    sv.biom.mean$TDT <- rowSums(sv.biom.mean[,c("TBDT", "TNDT")])
    # sv.biom <- sv.biom[complete.cases(sv.biom),]
    summary(sv.biom.mean); dim(sv.biom.mean)
    
    sv.biom.sd <- read.csv(file.path(raw.dir, "SetVeg_Stats", "biomass_uncertainty_v0.9-7_bam.csv"))
    summary(sv.biom.sd)
    
    sv.biom <- sv.biom.mean[,c("x", "y", "cell", "Total")]
    sv.biom$Total.SD <- sv.biom.sd$Total
    
    summary(sv.biom)
    
    # ---------------
    # Creating a lat/lon conversion for the setveg data
    # ---------------
    setveg.spat <- sv.biom
    coordinates(setveg.spat) <- c("x", "y")
    projection(setveg.spat) <- proj.albers
    
    setveg.spat <- spTransform(setveg.spat, CRS('+proj=longlat'))
    setveg.longlat <- data.frame(setveg.spat)
    summary(setveg.longlat)
    
    # Add Lon/Lat to the sv.tmp data frame so we have a cross walk
    sv.biom[,c("lon", "lat")] <- setveg.longlat[,c("x", "y")]
    summary(sv.biom)
    
    # Adding Albers references to the hips dataframe
    for(i in 1:nrow(hips)){
      hips[i,"x"] <- mean(sv.tmp[which(sv.tmp$lon>=hips[i,"lonmin"] & sv.tmp$lat>=hips[i,"latmin"] & sv.tmp$lon<=hips[i,"lonmax"] & sv.tmp$lat<=hips[i,"latmax"]), "x"])
      hips[i,"y"] <- mean(sv.tmp[which(sv.tmp$lon>=hips[i,"lonmin"] & sv.tmp$lat>=hips[i,"latmin"] & sv.tmp$lon<=hips[i,"lonmax"] & sv.tmp$lat<=hips[i,"latmax"]), "y"])
      hips[i,"xmin"] <- min(sv.tmp[which(sv.tmp$lon>=hips[i,"lonmin"] & sv.tmp$lat>=hips[i,"latmin"] & sv.tmp$lon<=hips[i,"lonmax"] & sv.tmp$lat<=hips[i,"latmax"]), "x"]) # SW corner x
      hips[i,"ymin"] <- min(sv.tmp[which(sv.tmp$lon>=hips[i,"lonmin"] & sv.tmp$lat>=hips[i,"latmin"] & sv.tmp$lon<=hips[i,"lonmax"] & sv.tmp$lat<=hips[i,"latmax"]), "y"]) # SW corner y
      hips[i,"xmax"] <- max(sv.tmp[which(sv.tmp$lon>=hips[i,"lonmin"] & sv.tmp$lat>=hips[i,"latmin"] & sv.tmp$lon<=hips[i,"lonmax"] & sv.tmp$lat<=hips[i,"latmax"]), "x"]) # NE corner x
      hips[i,"ymax"] <- max(sv.tmp[which(sv.tmp$lon>=hips[i,"lonmin"] & sv.tmp$lat>=hips[i,"latmin"] & sv.tmp$lon<=hips[i,"lonmax"] & sv.tmp$lat<=hips[i,"latmax"]), "y"]) # NE corner x
    }
    
    # Figures to make sure the spatial transformation and my indexing works
    png(file.path(fig.dir, "SetVeg_Biomass_Albers.png"), height=6, width=10, units="in", res=220)
    print(
      ggplot(data=sv.biom) +
        coord_equal() +
        geom_tile(aes(x=x, y=y, fill=Total)) +
        geom_rect(data=hips, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), color="red", fill=NA) +
        geom_point(data=hips, aes(x=x, y=y), color="red") +
        scale_x_continuous(limits=range(sv.biom$x))
    )
    dev.off()

    png(file.path(fig.dir, "SetVeg_Biomass_LonLat.png"), height=6, width=10, units="in", res=220)
    print(
      ggplot(data=sv.biom) +
        coord_equal() +
        geom_point(aes(x=lon, y=lat, color=Total)) +
        geom_rect(data=hips, aes(xmin=lonmin, xmax=lonmax, ymin=latmin, ymax=latmax), color="red", fill=NA) +
        geom_point(data=hips, aes(x=lon, y=lat), color="red")+
        scale_x_continuous(limits=range(sv.biom$lon))
    )
    dev.off()
    
    ggplot(data=sv.biom) +
      coord_equal() +
      geom_tile(aes(x=x, y=y, fill=Total)) +
      geom_rect(data=hips, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), color="red", fill=NA) +
      geom_point(data=hips, aes(x=x, y=y), color="red") +
      scale_x_continuous(limits=range(sv.biom$x))
    ggplot(data=sv.biom) +
      coord_equal() +
      geom_tile(aes(x=x, y=y, fill=Total.SD)) +
      geom_rect(data=hips, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), color="red", fill=NA) +
      geom_point(data=hips, aes(x=x, y=y), color="red") +
      scale_x_continuous(limits=range(sv.biom$x))
    # ---------------
    
    # ---------------
    # Extract the mean and the SD for the region modeled y the plot
    # ---------------
    summary(sv.biom)
    
    sv.biom.out <- hips[,c("Site", "lon", "lat")]
    for(i in 1:nrow(hips)){
      sv.biom.out[i,"Biomass"] <- mean(sv.biom[which(sv.biom$lon>=hips[i,"lonmin"] & sv.biom$lon>=hips[i,"lonmax"] & sv.biom$lat>=hips[i,"latmin"] & sv.biom$lat>=hips[i,"latmax"]), "Total"])
      sv.biom.out[i,"Biomass.SD"] <- mean(sv.biom[which(sv.biom$lon>=hips[i,"lonmin"] & sv.biom$lon>=hips[i,"lonmax"] & sv.biom$lat>=hips[i,"latmin"] & sv.biom$lat>=hips[i,"latmax"]), "Total.SD"])
    }
    
    # Approximating a 95% CI from the Std Dev (95% CI = mean +/- 2 SD)
    sv.biom.out$lwr <- sv.biom.out$Biomass - 2*sv.biom.out$Biomass.SD
    sv.biom.out$upr <- sv.biom.out$Biomass + 2*sv.biom.out$Biomass.SD
    
    sv.biom.out
    
    write.csv(sv.biom.out, file=file.path(out.dir, "SetVeg_Biomass.csv"), row.names=F)
    
    # ---------------
  }
  # ------------------
    
    
  # ------------------
  #    2.04. Biomass -- NBCD
  # ------------------
  {
    library(raster); library(rgdal); library(sp)
    
    # Crop the NBCD raster to just paleon to make it easier to deal with
    # nbcd <- raster(file.path(raw.dir, "NBCD_countrywide_biomass_240m_raster", "NBCD_countrywide_biomass_mosaic.tif"))
    # paleon <- raster("~/Desktop/Resarch/PalEON_CR/env_regional/env_paleon/domain_mask/paleon_domain.nc")
  
    # Transform the paleon domain to make life easier
    # paleon2 <- projectRaster(paleon, crs=projection(nbcd), filename=file.path(raw.dir, "NBCD_countrywide_biomass_240m_raster", "paleon_reproject"), overwrite=T)
  
    # Crop to the paleon domain
    # nbcd2 <- crop(nbcd, paleon2, filename=file.path(raw.dir, "NBCD_countrywide_biomass_240m_raster", "NBCD_paleon_domain"), overwrite=T)
    # plot(nbcd2)
  
    # Transforming the cropped grid to latlon 
    # nbcd3 <- projectRaster(nbcd2, crs=projection(paleon), filename=file.path(raw.dir, "NBCD_countrywide_biomass_240m_raster", "NBCD_paleon_latlon"), overwrite=T)
    # nbcd3
    
    nbcd3 <- raster(file.path(raw.dir, "NBCD_countrywide_biomass_240m_raster", "NBCD_paleon_latlon"))
    nbcd3
    # plot(nbcd3)
    
    hips.bm <- data.frame(Site=hips$Site)
    # Extracting the NBCD values just in the 0.5 x 0.5-degree cell for each site
    # Convert to KgC/m2 (from Mg Biomass per cell)
    for(i in 1:nrow(hips)){
      print(hips$Site[i])
      site.ext <- extent(c(hips$lonmin[i], hips$lonmax[i], hips$latmin[i], hips$latmax[i]))
      bm.site <- crop(nbcd3, extent(c(hips$lonmin[i], hips$lonmax[i], hips$latmin[i], hips$latmax[i])), filename=file.path(raw.dir, "NBCD_countrywide_biomass_240m_raster", paste0("NBCD_", hips$Site[i])), overwrite=T)
  
      hips.bm[i,"nbcd.mean"] <- mean(getValues(bm.site))/(240*240)*1e3*0.5     
      hips.bm[i,"nbcd.sd"  ] <- sd(getValues(bm.site))/(240*240)*1e3*0.5
      hips.bm[i,"nbcd.min" ] <- min(getValues(bm.site))/(240*240)*1e3*0.5
      hips.bm[i,"nbcd.max" ] <- max(getValues(bm.site))/(240*240)*1e3*0.5
             
    } # extent=xmin, xmax, ymin, ymax
    hips.bm
    
    write.csv(hips.bm, file.path(out.dir, "NBCD_sites_summary.csv"), row.names=F)
  
  }
  # ------------------
  
  # ------------------
  #    2.05. Carbon Fluxes - Flux Towers
  # ------------------
  { 
    # Some constants from Dave's L2 script 102_calculate...
    #molecular weights to convert from moles to g
    MolwtC = 12.0107
    MolwtO2 =16*2 
    MolwtCO2 = 44.01
    
    # Listing which variables I care about in the L2 data
    vars.L2 <- c("NEE", "RE", "GPP") # The flux variables I want to use
    cols.L2 <- c("SITE", "SiteName", "YEAR", "GAP", "DTIME", "DOY", "HRMIN")
    
    # Set the file path
    path.flux <- file.path(raw.dir, "fluxdata")
    
    # get list of flux data sites
    flux.sites <- dir(path.flux, "L2_data")
    flux.sites
    
    # Load annual GPP -- loads data frame AnnualGPP_MIPBenchmark
    load(file.path(path.flux, "AnnualGPP_MIPBenchmark.FLUXNET.22Jun2016.RData"))
    flux.gpp <- AnnualGPP_MIPBenchmark # change the name to something shorter
    flux.gpp$SiteName <- as.factor(flux.gpp$SiteName)
    summary(flux.gpp)
    
    # Getting the site-level data
    load(file.path(path.flux, "HowlandHo1.ameriflux.allsites.L2_data.17Jun2016.RData"))
    load(file.path(path.flux, "HarvardHa1.ameriflux.allsites.L2_data.05Mar2016.RData"))
    load(file.path(path.flux, "SylvaniaSyv.ameriflux.allsites.L2_data.17Jun2016.RData"))
    load(file.path(path.flux, "WillowCreekWCr.ameriflux.allsites.L2_data.17Jun2016.RData"))
    
    # Pulling out just the columns we care about to make life easier
    HowlandHo1$SITE <- as.factor("PHO")
    HowlandHo1$SiteName <- as.factor("Howland")
    HarvardHa1$SITE <- as.factor("PHA")
    HarvardHa1$SiteName <- as.factor("Harvard")
    SylvaniaSyv$SITE <- as.factor("PUN")
    SylvaniaSyv$SiteName <- as.factor("Sylvania")
    WillowCreekWCr$SITE <- as.factor("PUN")
    WillowCreekWCr$SiteName <- as.factor("WillowCreek")
    
    # Making everything into 1 dataframe
    flux.l2 <- rbind(HowlandHo1[,c(cols.L2, vars.L2)], HarvardHa1[,c(cols.L2, vars.L2)], SylvaniaSyv[,c(cols.L2, vars.L2)], WillowCreekWCr[,c(cols.L2, vars.L2)])
    summary(flux.l2)
    
    # Change our Missing (-9999) or not reported (-6999) values to NA
    flux.l2[flux.l2==-9999 | flux.l2==-6999] <- NA
    
    # Going ahead and converting fluxes from umol/m2/s to kgC/m2/s (MIP units)
    # Mol Weight * 1e-6 = conversion from umol to umol to mol
    flux.l2[,vars.L2] <- flux.l2[,vars.L2]*1e-6*MolwtC*1e-3
    
    # Aggregating to day before going to the year
    # ** Using na.rm=FALSE so that we can really see how big a deal the missing data is
    flux.l2.day <- aggregate(flux.l2[,vars.L2], by=flux.l2[,c("SITE", "SiteName", "YEAR", "DOY")], FUN=mean, na.rm=F)
    
    # Aggregating to the Year
    flux.l2.yr <- aggregate(flux.l2.day[,vars.L2], by=flux.l2.day[,c("SITE", "SiteName", "YEAR")], FUN=mean, na.rm=F)
    
    # Changing yearly fluxes to kgC/m2/Yr
    flux.l2.yr[,vars.L2] <- flux.l2.yr[,vars.L2]*365*24*60*60
    summary(flux.l2.yr)
    
    # QA/QC plot
    flux.stack <- stack(flux.l2.yr[,vars.L2])
    names(flux.stack) <- c("values", "Flux")
    flux.stack[,c("SITE", "SiteName", "YEAR")] <- flux.l2.yr[,c("SITE", "SiteName", "YEAR")]
    
    png(file.path(fig.dir, "FluxData_Year.png"), height=8, width=8, units="in", res=220)
    ggplot(data=flux.stack) +
      facet_grid(Flux~., scales="free_y") +
      geom_line(aes(x=YEAR, y=values, color=SiteName), size=2) +
      scale_x_continuous(expand=c(0,0)) +
      scale_y_continuous(expand=c(0,0), name="kgC/m2/yr") +
      theme_bw() +
      theme(legend.position="top")
    dev.off()
    
    # Write the data out
    write.csv(flux.l2.yr, file.path(out.dir, "FluxData_Year.csv"), row.names=F)
    
  }
  # ------------------
  
  # ------------------
  #    2.06. MODIS -- LAI, GPP
  # ------------------
  {
    modis.lai <- read.csv(file.path(raw.dir, "MODIS", "Benchmarking_LAI_MODIS_SummaryStats.csv"))
    modis.gpp <- read.csv(file.path(raw.dir, "MODIS", "Benchmarking_GPP_MODIS_SummaryStats.csv"))
    modis.lai$Var = as.factor("LAI")
    modis.gpp$Var = as.factor("GPP")
    summary(modis.lai)
    summary(modis.gpp)
    
    modis.wide <- rbind(modis.lai, modis.gpp)
    summary(modis.wide)
    
    modis <- stack(modis.wide[,which(substr(names(modis.wide), 5, 8)=="Mean")])
    names(modis) <- c("mean", "Site")
    modis$Site <- as.factor(substr(modis$Site, 1, 3)) # Extract the site from the former column names
    modis[,c("Year", "Var")] <- modis.wide[,c("Year", "Var")] # Add in the other important identifiers
    modis <- modis[,c("Site", "Var", "Year", "mean")] # Rearrange to a prettier, more logical order
    modis$sd <- stack(modis.wide[,which(substr(names(modis.wide), 5, 6)=="SD")])[,1]
    modis$min <- stack(modis.wide[,which(substr(names(modis.wide), 5, 7)=="Min")])[,1]
    modis$max <- stack(modis.wide[,which(substr(names(modis.wide), 5, 7)=="Max")])[,1]
    summary(modis)
    
    # Convert GPP from kg/m2/8-days to kg/m2/s
    modis[modis$Var=="GPP",c("mean", "sd", "min", "max")] <- modis[modis$Var=="GPP",c("mean", "sd", "min", "max")] / (8*24*60*60) 
  
    modis$Site <- factor(modis$Site, levels=c("PDL", "PBL", "PUN", "PMB", "PHA", "PHO"))

    png(file.path(fig.dir, "ModisData_Year.png"), height=8, width=10, units="in", res=220)
    ggplot(data=modis) +
      facet_grid(Var~Site, scales="free_y") +
      geom_ribbon(aes(x=Year, ymin=min, ymax=max, fill=Site), alpha=0.5) +
      # geom_ribbon(aes(x=Year, ymin=mean-sd, ymax=max+sd, fill=Site), alpha=0.8) +
      geom_line(aes(x=Year, y=mean, color=Site), size=1.5) +
      geom_point(aes(x=Year, y=mean, color=Site), size=2) +
      scale_x_continuous(expand=c(0,0), breaks=seq(min(modis$Year)+1, max(modis$Year), by=4)) +
      theme_bw() 
    dev.off()
    
    write.csv(modis, file.path(out.dir, "MODIS_year_sites_summary.csv"), row.names=F)
  }
  # ------------------

  
  # ------------------
  #    2.07. Tree Rings - AGB, AGBi (AGB increment)
  # ------------------
  {
    
  }
  # ------------------
  
}
# -------------------------------------------

