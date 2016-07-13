# -------------------------------------------
# Quantifying and comparing periods & magnitude of change through time & space among met drivers, 
# models, & ecosystem variables
#
# Author: Christy Rollinson, crollinson@gmail.com
#
# ------------------
# Objectives:
# ------------------
# A. Model Benchmarking Score Card
#    1. Summary by benchmark dataset of model & benchmark mean & SD, pair-wise bias, and the number 
#       of sites the model was within the observed range/95%CI
# B. Correlating Overall Model Behaviorwith Model Characteristics
#    1. Use Ordination with the benchmarking scores to see if models group and whether model characteristics
#       explain those relationships
# ------------------
#
#
# ------------------
# Workflow
# ------------------
# 1. Set up file structure, libraries, etc.
#
# 2. Extract & Format Benchmarks
#    2.01. Climate
#    2.02. Composition -- STEPPS 1 & 2
#    2.03. Composition -- FIA Biomass
#    2.04. Composition -- SetVeg
#    2.05. Biomass -- NBCD
#    2.06. Biomass -- SetVeg
#    2.07. Biomass -- BabySTEPPS
#    2.08. Carbon Fluxes - Flux Towers
#    2.09. Carbon Fluxes - MODIS
#    2.10. LAI -- MODIS
#         --
#
# 3. Comparisons with Models
#    3.1. Climate
#         -- Region
#         -- HIPS
#    3.2. Composition
#         -- Paleo
#         -- SetVeg
#         -- Modern (FIA)
#    3.3. Biomass
#         -- Paleo
#         -- SetVeg
#         -- Modern - Raster (NBCD)
#         -- Modern - Site
#    3.4. Carbon Fluxes
#         -- Modern - Site
#         -- Modern - Raster
#    3.5. LAI
#         -- Modern
#
# 4. Model behavior (ordination)
# ------------------
# -------------------------------------------
rm(list=ls())

# -------------------------------------------
# 1. Set up file structure, libraries, etc.
# -------------------------------------------
library(car)
library(ggplot2)
library(mgcv)
setwd("~/Dropbox/PalEON_CR/PalEON_MIP_Site/Analyses/Change-and-Stability") # Path to this project github repository: https://github.com/PalEON-Project/MIP-Change-and-Stability.git
# path.gamm.func <- "~/Desktop/R_Functions/"  # Path to github repository of my GAMM helper functions: https://github.com/crollinson/R_Functions.git
inputs    <- "Data/" # Path to my cleaned model output

mip.utils <- "~/Dropbox/PalEON_CR/MIP_Utils/" # Path to PalEON MIP Utility repository: https://github.com/PalEON-Project/MIP_Utils.git

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
# ------------------
#    2.01. Climate
# ------------------
#         --
# ------------------

# ------------------
#    2.02. Composition -- STEPPS 1 & 2
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
      stepps2.pft[[m]][stepps2.pft[[m]]$site==s,"upr"] <- apply(dat.temp, 1, quantile, 0.025)
      stepps2.pft[[m]][stepps2.pft[[m]]$site==s,"lwr"] <- apply(dat.temp, 1, quantile, 0.975)
      
    }
  }
  
  pdf(file.path(fig.dir, "STEPPS2_ModelPFTs.pdf"))
  for(m in names(stepps2.pft)){
    print(
      ggplot(data=stepps2.pft[[m]]) +
        geom_ribbon(aes(x=year, ymin=lwr, ymax=upr, fill=site), alpha=0.4) +
        geom_line(aes(x=year, y=mean, color=site, linetype=pft)) +
        ggtitle(m)
    )
  }
  dev.off()
  
  save(stepps2.pft, file=file.path(out.dir, "STEPPS2_ModelPFTs.RData"))
  # --------
}
# ------------------

# ------------------
#    2.03. Composition -- FIA Biomass
# ------------------
# ------------------

# ------------------
#    2.04. Composition -- SetVeg
# ------------------
# ------------------

# ------------------
#    2.05. Biomass -- NBCD
# ------------------
{
  library(raster); library(rgdal); library(sp)
  
  # Crop the NBCD raster to just paleon to make it easier to deal with
  # nbcd <- raster(file.path(raw.dir, "NBCD_countrywide_biomass_240m_raster", "NBCD_countrywide_biomass_mosaic.tif"))
  # paleon <- raster("~/Dropbox/PalEON_CR/env_regional/env_paleon/domain_mask/paleon_domain.nc")

  # Transform the paleon domain to make life easier
  # paleon2 <- projectRaster(paleon, crs=projection(nbcd), filename=file.path(raw.dir, "NBCD_countrywide_biomass_240m_raster", "paleon_reproject"), overwrite=T)

  # Crop to the paleon domain
  # nbcd2 <- crop(nbcd, paleon2, filename=file.path(raw.dir, "NBCD_countrywide_biomass_240m_raster", "NBCD_paleon_domain"), overwrite=T)
  # plot(nbcd2)

  # Transforming the cropped grid to latlon and from Mg/cell to MgC/Ha and saving it to work with
  # nbcd3 <- projectRaster(nbcd2, crs=projection(paleon), filename=file.path(raw.dir, "NBCD_countrywide_biomass_240m_raster", "NBCD_paleon_latlon"), overwrite=T)
  # nbcd3
  
  nbcd3 <- raster(file.path(raw.dir, "NBCD_countrywide_biomass_240m_raster", "NBCD_paleon_latlon"))
  nbcd3
  plot(nbcd3)
  
  hips.bm <- data.frame(Site=hips$Site)
  # Extracting the NBCD values just in the 0.5 x 0.5-degree cell for each site
  # Convert to MgC/Ha (from Mg Biomass per cell)
  for(i in 1:nrow(hips)){
    print(hips$Site[i])
    site.ext <- extent(c(hips$lonmin[i], hips$lonmax[i], hips$latmin[i], hips$latmax[i]))
    bm.site <- crop(nbcd3, extent(c(hips$lonmin[i], hips$lonmax[i], hips$latmin[i], hips$latmax[i])), filename=file.path(raw.dir, "NBCD_countrywide_biomass_240m_raster", paste0("NBCD_", hips$Site[i])), overwrite=T)

    hips.bm[i,"nbcd.mean"] <- mean(getValues(bm.site)/5.76*0.5)     
    hips.bm[i,"nbcd.sd"  ] <- sd(getValues(bm.site)/5.76*0.5)     
    hips.bm[i,"nbcd.min" ] <- min(getValues(bm.site)/5.76*0.5)     
    hips.bm[i,"nbcd.max" ] <- max(getValues(bm.site)/5.76*0.5)     
           
  } # extent=xmin, xmax, ymin, ymax
  hips.bm
  
  write.csv(hips.bm, file.path(out.dir, "NBCD_sites_summary.csv"), row.names=F)

}
# ------------------

# ------------------
#    2.06. Biomass -- SetVeg
# ------------------
# ------------------

# ------------------
#    2.07. Biomass -- BabySTEPPS
# ------------------
# ------------------

# ------------------
#    2.08. Carbon Fluxes - Flux Towers
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
#    2.09. Carbon Fluxes - MODIS
# ------------------
# ------------------

# ------------------
#    2.10. LAI -- MODIS
# ------------------
# ------------------



# -------------------------------------------

# -------------------------------------------
# 3. Comparisons with Models
# -------------------------------------------
# ------------------
#    3.1. Climate
# ------------------
#         -- Region
#         -- HIPS
# ------------------

# ------------------
#    3.2. Composition
# ------------------
#         -- Paleo
#         -- SetVeg
#         -- Modern (FIA)
# ------------------

# ------------------
#    3.3. Biomass
# ------------------
#         -- Paleo
#         -- SetVeg
#         -- Modern - Raster (NBCD)
#         -- Modern - Site
# ------------------

# ------------------
#    3.4. Carbon Fluxes
# ------------------
#         -- Modern - Site
#         -- Modern - Raster
# ------------------

# ------------------
#    3.5. LAI
# ------------------
#         -- Modern
# ------------------
# -------------------------------------------

# -------------------------------------------
# 4. Model behavior (ordination)
# -------------------------------------------
# -------------------------------------------
