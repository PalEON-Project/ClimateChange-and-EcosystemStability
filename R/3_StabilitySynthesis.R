# -------------------------------------------
# Quantifying and comparing periods & magnitude of change through time & space among met drivers, models, & ecosystem variables
# Author: Christy Rollinson, crollinson@gmail.com
#
# ------------------
# Objectives:
# ------------------
# A. Met Driver Analysis
#    1. Spatially coherent shifts in change -- do all sites show similar timing & magnitude of change in the drivers?
#    2. Met Drivers vs. Paleo Data -- do met drivers show similar timing & rates of change as paleo proxy data?
#
# B. Ecosystem Analysis
#    1. Which ecosytem variables were most/least stable at/before 1850?
#       -- How does this compare to STEPPS output
#    2. Does is ecosystem instability temporally aligned with met instability or is it lagged?
#    2. Which models had greatest/least tendencies towards pre-1850 stability?
#    3. Is modern stability similar or different from what was observed in the past?
#        -- by model/variable
# ------------------
#
#
# ------------------
# Workflow
# ------------------
# 1. Set up file structure, libraries, etc.
#
# 2. Analyzing Met Drivers
#    2.1. Site vs. Regional Trends -- change relative to spinup period
#         -- Change going into 1850
#         -- Change at end of simulations
#    2.2. Site Comparisons -- all broken down into pre & post 1850
#         -- mean rate of unstable periods
#         -- mean duration of unstable periods
#         -- percent time spent in unstable state
#         -- Time of greatest instability
#    2.3. Drivers vs. Paleo-Proxy Comparisions
#         -- times of greatest instability
#         -- rates of unstable period
#         -- duration of unstable periods
#
# 3. Analysing Ecosystem Model Output
#    3.1. Generate a synthesizing figure
#         -- compare models & variables; site probably not important
#    3.2. Calculate average instability statistics by variable & model (pre- & post-1850)
#         -- mean rate of unstable periods
#         -- mean duration of unstable periods
#         -- percent time spent in unstable state
#         -- measure of interannual variability
#              - SD as % of mean?
#    3.3. Compare rates & durations of instability by variable (pre & post-1850)
#         -- are the most variables stable pre-1850 also the most stable post 1900?
#         -- are slow processes/states more or less stable than fast ones?
#    3.4. Compare rates & durations of instability by models
#         -- are the most stable models pre-1850 also the most stable post 1900?
#    3.5. Can we quantify the relationship between interannual variability and ecosystem stability?
#         -- something that applies to both models & variables 
#         -- can we use somethign about interannual variability to predict stability today & in the past?
# ------------------
# -------------------------------------------
rm(list=ls())

# -------------------------------------------
# 1. Set up file structure, libraries, etc.
# -------------------------------------------
library(car)
library(ggplot2); library(scales)
library(mgcv)
# setwd("~/Dropbox/PalEON_CR/PalEON_MIP_Site/Analyses/Change-and-Stability") # Path to this project github repository: https://github.com/PalEON-Project/MIP-Change-and-Stability.git
setwd("~/Desktop/Research/PalEON_CR/PalEON_MIP_Site/Analyses/Change-and-Stability") # Path to this project github repository: https://github.com/PalEON-Project/MIP-Change-and-Stability.git
# path.gamm.func <- "~/Desktop/R_Functions/"  # Path to github repository of my GAMM helper functions: https://github.com/crollinson/R_Functions.git
inputs    <- "Data/" # Path to my cleaned model output

mip.utils <- "~/Dropbox/PalEON_CR/MIP_Utils/" # Path to PalEON MIP Utility repository: https://github.com/PalEON-Project/MIP_Utils.git

out.dir <- "Data/StabilitySynthesis" # Path to where the analysis output should go
fig.dir <- "Figures/StabilitySynthesis" # Path to where figures should go

if(!dir.exists(out.dir)) dir.create(out.dir)
if(!dir.exists(fig.dir)) dir.create(fig.dir)

# Colors used for graphing
model.colors <- read.csv("../../Model.Colors.csv")
model.colors $Model.Order <- recode(model.colors$Model, "'CLM4.5-BGC'='01'; 'CLM4.5-CN'='02'; 'ED2'='03'; 'ED2-LU'='04';  'JULES-STATIC'='05'; 'JULES-TRIFFID'='06'; 'LINKAGES'='07'; 'LPJ-GUESS'='08'; 'LPJ-WSL'='09'; 'SiBCASA'='10'")
levels(model.colors$Model.Order)[1:10] <- c("CLM-BGC", "CLM-CN", "ED2", "ED2-LU", "JULES-STATIC", "JULES-TRIFFID", "LINKAGES", "LPJ-GUESS", "LPJ-WSL", "SiBCASA")
model.colors <- model.colors[order(model.colors$Model.Order),]

model.colors$Model.Family <- ifelse(substr(model.colors$Model,1,3)=="CLM", "CLM", 
                                    ifelse(substr(model.colors$Model, 1, 3)=="ED2", "ED2",
                                           ifelse(substr(model.colors$Model, 1, 3)=="LPJ", "LPJ",
                                                  ifelse(substr(model.colors$Model,1,5)=="JULES", "JULES",
                                                         toupper(model.colors$Model)))))
model.colors$Model.Family <- as.factor(model.colors$Model.Family)
model.colors$shape <- as.numeric(model.colors$Model.Family)
model.colors

# -------------------------------------------


# ------------------
# 2.0 Load Ecosystem Stability Results
# ------------------
{
  # Load in individual data sets
  soilcarb <- read.csv("Data/EcosysChange/StabilityCalcs_SoilCarb_850-1850_1900-2010.csv")
  fcomp <- read.csv("Data/EcosysChange/StabilityCalcs_Fcomp_850-1850_1900-2010.csv")
  agb <- read.csv("Data/EcosysChange/StabilityCalcs_AGB_850-1850_1900-2010.csv")
  lai <- read.csv("Data/EcosysChange/StabilityCalcs_LAI_850-1850_1900-2010.csv")
  nee <- read.csv("Data/EcosysChange/StabilityCalcs_NEE_850-1850_1900-2010.csv")
  gpp <- read.csv("Data/EcosysChange/StabilityCalcs_GPP_850-1850_1900-2010.csv")
  
  # Make one mega-dataframe
  ecosys.change <- rbind(soilcarb, fcomp, agb, lai, nee, gpp)
  
  ecosys.change$var  <- factor(ecosys.change$var , levels=c("GPP", "NEE", "LAI", "AGB", "Fcomp", "SoilCarb"))
  ecosys.change$Site <- factor(ecosys.change$Site, levels=c("PDL", "PBL", "PUN", "PMB", "PHA", "PHO"))
  summary(ecosys.change)
  
  # Change fluxes to per year to make better graphs
  sec2yr <- 1*60*60*24*365 # 1 sec * 60 sec/min * 60 min/day * 365 day/yr
  vars.flux <- c("NEE", "GPP")
  ecosys.change[ecosys.change$var %in% vars.flux, c("Y", "mean", "lwr", "upr", "deriv.mean", "deriv.lwr", "deriv.upr")] <- ecosys.change[ecosys.change$var %in% vars.flux, c("Y", "mean", "lwr", "upr", "deriv.mean", "deriv.lwr", "deriv.upr")]*sec2yr
  
  
  # Graphing the data to make sure it read in okay
  for(v in levels(ecosys.change$var)){
    # print(
    #   ggplot(data=ecosys.change[ecosys.change$var==v,]) + 
    #     facet_wrap(~Model.Name, scales="free_y") +
    #     geom_line(aes(x=Year, y=Y, color=Site), size=0.5, alpha=0.3) +
    #     geom_ribbon(data=ecosys.change[ecosys.change$var==v & ecosys.change$Year<1850,], aes(x=Year, ymin=lwr, ymax=upr, fill=Site), alpha=0.3) +
    #     geom_ribbon(data=ecosys.change[ecosys.change$var==v & ecosys.change$Year>1900,], aes(x=Year, ymin=lwr, ymax=upr, fill=Site), alpha=0.3) +
    #     geom_line(data=ecosys.change[ecosys.change$var==v & ecosys.change$Year<1850,], aes(x=Year, y=mean, color=Site), size=1, alpha=0.2) +
    #     geom_line(data=ecosys.change[ecosys.change$var==v & ecosys.change$Year>1900,], aes(x=Year, y=mean, color=Site), size=1, alpha=0.2) +
    #     geom_line(aes(x=Year, y=mean.sig, color=Site), size=1.5, alpha=1) +
    #     geom_vline(xintercept=1850, linetype="dashed") +
    #     geom_vline(xintercept=1900, linetype="dashed") +
    #     scale_x_continuous(expand=c(0,0), name="Year") +
    #     scale_y_continuous(expand=c(0,0), name="ecosys.change") +
    #     ggtitle(v) +
    #     # scale_color_manual(values=col.model) +
    #     # scale_fill_manual(values=col.model) + 
    #     theme_bw()
    # )
  }
  
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
  ecosys.change.mod <- aggregate(ecosys.change[,c("Y", "Y.anom", "mean.anom", "lwr.anom", "upr.anom", "deriv.mean", "deriv.lwr", "deriv.upr")],
                                 by=ecosys.change[,c("Year", "var", "Model", "Model.Order")],
                                 FUN=mean, na.rm=T)
  ecosys.change.mod$n.sig <- aggregate(ecosys.change[,"sig"],
                                       by=ecosys.change[,c("Year", "var", "Model", "Model.Order")],
                                       FUN=function(x){length(which(x == "*"))})[,"x"]
  ecosys.change.mod$mod.sig <- as.factor(ifelse(ecosys.change.mod$n.sig>3, "*", NA))
  summary(ecosys.change.mod)
  # ------------------
}
# ------------------

# ------------------
# 3. Analyzing sensitivity 
#    4.1 average years per site per model of instability (pre/post); percentage?
#    4.2 max rate of change pre/post
#    4.3 model traits as predictors of pre-settlement stability
# ------------------
{
  summary(ecosys.change)
  # ------------------
  # 4.1 average years per site per model of instability (pre/post); percentage?
  # ------------------
  {
    # Aggregating to get some of the easy stats
    ecosys.stats.site <- merge(data.frame(Site=levels(ecosys.change$Site)), data.frame(var=levels(ecosys.change$var)))
    ecosys.stats.site <- merge(ecosys.stats.site, data.frame(Model=levels(ecosys.change$Model), Model.Order=levels(ecosys.change$Model.Order)))
    summary(ecosys.stats.site)
    
    # Set up a list to paralellize the calculations
    stats.list <- list()
    
    pb <- txtProgressBar(min = 0, max = nrow(ecosys.stats.site), style = 3)
    i=1
    
    # Figuring out site-level stability stats
    for(v in unique(ecosys.stats.site$var)){
      for(s in unique(ecosys.stats.site$Site)){
        # print(paste0(" ---- ", s, " ---- "))
        
        for(m in unique(ecosys.stats.site$Model)){
          setTxtProgressBar(pb, i)
          
          if(is.na(mean(ecosys.change[ecosys.change$var==v & ecosys.change$Site==s  & ecosys.change$Model==m ,"deriv.mean"], na.rm=T))){
            i=i+1
            next
          }
          stats.list[[paste(v, s, m, sep="_")]] <- ecosys.change[ecosys.change$var==v & ecosys.change$Site==s  & ecosys.change$Model==m,]
          i=i+1
         }  # End Model loop
      } # End site loop
      # print(ecosys.stats.site[ecosys.stats.site$var==v,])
    } # End var loop
    
    # Function to find periods of instabilty (this is the long step)
    instability.duration <- function(change.data){
      # Set up vectors for each site/var
      ins1 <- vector() # Vector for pre-1850 instability
      ins2 <- vector() # vector for post-1900 instability
      dur=0 # Start duration count at 0
      
      # ---------------
      # Go through each year to find duration of instability periods because I don't know how else to do it
      # ---------------
      for(y in min(change.data$Year):max(ecosys.change$Year)){ 
        if(!is.na(change.data[change.data$Year==y,"sig"])){
          dur=dur+1 # If this is a period of significant change add it to the year count
          if(y==max(change.data$Year)){ # If this is our last year and it's instable, record the duration
            ins2 <- c(ins2, dur)
          }
        } else { # this year is instable 
          # If we have a logged duration of instability, record it
          if(dur > 0 & y<=1850){
            ins1 <- c(ins1, dur)
          } else if(dur > 0 & y>=1900) {
            ins2 <- c(ins2, dur)
          }
          dur=0 # Set duration back to 0
        }
      } # End  year loop
      # ---------------
      inst.dur <- list(pre.set=ins1, modern=ins2)
      return(inst.dur)
    }

    # Run the caclulation in parallel
    # test <- list()
    # test[["GPP_PDL_clm.bgc"]] <- stats.list[["GPP_PDL_clm.bgc"]]
    # test[["NEE_PDL_clm.bgc"]] <- stats.list[["NEE_PDL_clm.bgc"]]
    # test[["AGB_PDL_clm.bgc"]] <- stats.list[["AGB_PDL_clm.bgc"]]
    # test[["LAI_PDL_clm.bgc"]] <- stats.list[["LAI_PDL_clm.bgc"]]
    # test.out <- instability.duration(test)
    
    library(parallel)
    cores.use <- min(8, length(stats.list))
    dat.out <- mclapply(stats.list, instability.duration, mc.cores=cores.use)
    
    
    # Figuring out site-level stability stats
    i=1
    for(v in unique(ecosys.stats.site$var)){
      for(s in unique(ecosys.stats.site$Site)){
        # print(paste0(" ---- ", s, " ---- "))
        
        for(m in unique(ecosys.stats.site$Model)){
          name.ind <- paste(v, s, m, sep="_")
          setTxtProgressBar(pb, i)
          
          if(is.na(mean(ecosys.change[ecosys.change$var==v & ecosys.change$Site==s  & ecosys.change$Model==m ,"deriv.mean"], na.rm=T))){
            i=i+1
            next
          }
          
          # ---------------
          # Insert data for paleo period
          # ---------------
          if(length(dat.out[[name.ind]]$pre.set)==0) {
            ecosys.stats.site[ecosys.stats.site$Site==s & ecosys.stats.site$var==v & ecosys.stats.site$Model==m, "n.inst.paleo"] <- length(ecosys.change[ecosys.change$Site==s & ecosys.change$var==v & ecosys.change$Model==m & ecosys.change$Year<1850 & !is.na(ecosys.change$sig), "deriv.mean"])
            ecosys.stats.site[ecosys.stats.site$Site==s & ecosys.stats.site$var==v & ecosys.stats.site$Model==m, "rate.mean.paleo"] <- 0
            ecosys.stats.site[ecosys.stats.site$Site==s & ecosys.stats.site$var==v & ecosys.stats.site$Model==m, "rate.sd.paleo"] <- NA
            
            
            # Saving the mean duration
            ecosys.stats.site[ecosys.stats.site$Site==s & ecosys.stats.site$var==v & ecosys.stats.site$Model==m, "dur.mean.paleo"] <- 0
            ecosys.stats.site[ecosys.stats.site$Site==s & ecosys.stats.site$var==v & ecosys.stats.site$Model==m, "dur.sd.paleo"  ] <- NA
            
            # Saving the max significant rate of change & time that occurs
            # max.paleo <- which(ecosys.change$Site==s & ecosys.change$var==v & ecosys.change$Model==m & ecosys.change$Year<1850 & !is.na(ecosys.change$sig) & abs(ecosys.change$deriv.mean)==abs(max(ecosys.change[ecosys.change$Site==s & ecosys.change$var==v & ecosys.change$Model==m & ecosys.change$Year<1850 & !is.na(ecosys.change$sig), "deriv.mean"], na.rm=T)))
            ecosys.stats.site[ecosys.stats.site$Site==s & ecosys.stats.site$var==v & ecosys.stats.site$Model==m, "rate.max.paleo"] <- 0
            ecosys.stats.site[ecosys.stats.site$Site==s & ecosys.stats.site$var==v & ecosys.stats.site$Model==m, "rate.max.paleo.yr"] <- NA
          } else {
            # Finding the number of years in instability
            ecosys.stats.site[ecosys.stats.site$Site==s & ecosys.stats.site$var==v & ecosys.stats.site$Model==m, "n.inst.paleo"] <- length(ecosys.change[ecosys.change$Site==s & ecosys.change$var==v & ecosys.change$Model==m & ecosys.change$Year<1850 & !is.na(ecosys.change$sig), "deriv.mean"])
            ecosys.stats.site[ecosys.stats.site$Site==s & ecosys.stats.site$var==v & ecosys.stats.site$Model==m, "rate.mean.paleo"] <- mean(ecosys.change[ecosys.change$Site==s & ecosys.change$var==v & ecosys.change$Model==m & ecosys.change$Year<1850 & !is.na(ecosys.change$sig), "deriv.mean"])
            ecosys.stats.site[ecosys.stats.site$Site==s & ecosys.stats.site$var==v & ecosys.stats.site$Model==m, "rate.sd.paleo"] <- sd(ecosys.change[ecosys.change$Site==s & ecosys.change$var==v & ecosys.change$Model==m & ecosys.change$Year<1850 & !is.na(ecosys.change$sig), "deriv.mean"])
            
            # Saving the mean duration
            ecosys.stats.site[ecosys.stats.site$Site==s & ecosys.stats.site$var==v & ecosys.stats.site$Model==m, "dur.mean.paleo"] <- mean(dat.out[[name.ind]]$pre.set, na.rm=T)
            ecosys.stats.site[ecosys.stats.site$Site==s & ecosys.stats.site$var==v & ecosys.stats.site$Model==m, "dur.sd.paleo"] <- sd(dat.out[[name.ind]]$pre.set, na.rm=T)
            
            # Saving the max significant rate of change & time that occurs
            max.paleo <- which(ecosys.change$Site==s & ecosys.change$var==v & ecosys.change$Model==m & ecosys.change$Year<1850 & !is.na(ecosys.change$sig) & abs(ecosys.change$deriv.mean)==abs(max(ecosys.change[ecosys.change$Site==s & ecosys.change$var==v & ecosys.change$Model==m & ecosys.change$Year<1850 & !is.na(ecosys.change$sig), "deriv.mean"], na.rm=T)))
            ecosys.stats.site[ecosys.stats.site$Site==s & ecosys.stats.site$var==v & ecosys.stats.site$Model==m, "rate.max.paleo"] <- ecosys.change[max.paleo, "deriv.mean"]
            ecosys.stats.site[ecosys.stats.site$Site==s & ecosys.stats.site$var==v & ecosys.stats.site$Model==m, "rate.max.paleo.yr"] <- ecosys.change[max.paleo, "Year"]
          }
          # ---------------
          
          # ---------------
          # Insert data for Modern period
          # ---------------
          if(length(dat.out[[name.ind]]$modern)==0) {
            ecosys.stats.site[ecosys.stats.site$Site==s & ecosys.stats.site$var==v & ecosys.stats.site$Model==m, "n.inst.modrn"] <- length(ecosys.change[ecosys.change$Site==s & ecosys.change$var==v & ecosys.change$Model==m & ecosys.change$Year>1900 & !is.na(ecosys.change$sig), "deriv.mean"])
            ecosys.stats.site[ecosys.stats.site$Site==s & ecosys.stats.site$var==v & ecosys.stats.site$Model==m, "rate.mean.modrn"] <- 0
            ecosys.stats.site[ecosys.stats.site$Site==s & ecosys.stats.site$var==v & ecosys.stats.site$Model==m, "rate.sd.modrn"] <- NA
            
            # Saving the mean duration
            ecosys.stats.site[ecosys.stats.site$Site==s & ecosys.stats.site$var==v & ecosys.stats.site$Model==m, "dur.mean.modrn"] <- 0
            ecosys.stats.site[ecosys.stats.site$Site==s & ecosys.stats.site$var==v & ecosys.stats.site$Model==m, "dur.sd.modrn"  ] <- NA
            
            # Saving the max significant rate of change & time that occurs
            # max.modrn <- which(ecosys.change$Site==s & ecosys.change$var==v & ecosys.change$Model==m & ecosys.change$Year<1850 & !is.na(ecosys.change$sig) & abs(ecosys.change$deriv.mean)==abs(max(ecosys.change[ecosys.change$Site==s & ecosys.change$var==v & ecosys.change$Model==m & ecosys.change$Year<1850 & !is.na(ecosys.change$sig), "deriv.mean"], na.rm=T)))
            ecosys.stats.site[ecosys.stats.site$Site==s & ecosys.stats.site$var==v & ecosys.stats.site$Model==m, "rate.max.modrn"] <- 0
            ecosys.stats.site[ecosys.stats.site$Site==s & ecosys.stats.site$var==v & ecosys.stats.site$Model==m, "rate.max.modrn.yr"] <- NA
          } else {
            # Finding the number of years in instability
            ecosys.stats.site[ecosys.stats.site$Site==s & ecosys.stats.site$var==v & ecosys.stats.site$Model==m, "n.inst.modrn"] <- length(ecosys.change[ecosys.change$Site==s & ecosys.change$var==v & ecosys.change$Model==m & ecosys.change$Year>1900 & !is.na(ecosys.change$sig), "deriv.mean"])
            ecosys.stats.site[ecosys.stats.site$Site==s & ecosys.stats.site$var==v & ecosys.stats.site$Model==m, "rate.mean.modrn"] <- mean(ecosys.change[ecosys.change$Site==s & ecosys.change$var==v & ecosys.change$Model==m & ecosys.change$Year>1900 & !is.na(ecosys.change$sig), "deriv.mean"])
            ecosys.stats.site[ecosys.stats.site$Site==s & ecosys.stats.site$var==v & ecosys.stats.site$Model==m, "rate.sd.modrn"] <- sd(ecosys.change[ecosys.change$Site==s & ecosys.change$var==v & ecosys.change$Model==m & ecosys.change$Year>1900 & !is.na(ecosys.change$sig), "deriv.mean"])
            
            # Saving the mean duration
            ecosys.stats.site[ecosys.stats.site$Site==s & ecosys.stats.site$var==v & ecosys.stats.site$Model==m, "dur.mean.modrn"] <- mean(dat.out[[name.ind]]$modern, na.rm=T)
            ecosys.stats.site[ecosys.stats.site$Site==s & ecosys.stats.site$var==v & ecosys.stats.site$Model==m, "dur.sd.modrn"] <- sd(dat.out[[name.ind]]$modern, na.rm=T)
            
            # Saving the max significant rate of change & time that occurs
            max.modrn <- which(ecosys.change$Site==s & ecosys.change$var==v & ecosys.change$Model==m & ecosys.change$Year>1900 & !is.na(ecosys.change$sig) & abs(ecosys.change$deriv.mean)==abs(max(ecosys.change[ecosys.change$Site==s & ecosys.change$var==v & ecosys.change$Model==m & ecosys.change$Year>1900 & !is.na(ecosys.change$sig), "deriv.mean"], na.rm=T)))
            ecosys.stats.site[ecosys.stats.site$Site==s & ecosys.stats.site$var==v & ecosys.stats.site$Model==m, "rate.max.modrn"] <- ecosys.change[max.modrn, "deriv.mean"]
            ecosys.stats.site[ecosys.stats.site$Site==s & ecosys.stats.site$var==v & ecosys.stats.site$Model==m, "rate.max.modrn.yr"] <- ecosys.change[max.modrn, "Year"]
          }
          
          
          i=i+1
        }  # End Model loop
      } # End site loop
      # print(ecosys.stats.site[ecosys.stats.site$var==v,])
    } # End var loop
    
    summary(ecosys.stats.site)
    
    write.csv(ecosys.stats.site, file.path(out.dir, "Stability_Models_Stats_All.csv"), row.names=F)
    
    
    
    # Aggregating to the stats to variable level, leveraging the spatial sd
    models.stats.means <- aggregate(ecosys.stats.site[,5:ncol(ecosys.stats.site)],
                                    by=ecosys.stats.site[,c("Model", "Model.Order", "var")],
                                    FUN=mean, na.rm=T)
    summary(models.stats.means)
    
    models.stats.sd <- aggregate(ecosys.stats.site[,5:ncol(ecosys.stats.site)],
                                 by=ecosys.stats.site[,c("Model", "Model.Order","var")],
                                 FUN=sd, na.rm=T)
    summary(models.stats.sd)
 
    # Aggregating to the stats to variable level, leveraging the spatial sd
    ensemble.stats.means <- aggregate(models.stats.means[,4:ncol(models.stats.means)],
                                    by=list(models.stats.means[,c("var")]),
                                    FUN=mean, na.rm=T)
    names(ensemble.stats.means)[1] <- "var"
    ensemble.stats.means
    
    ensemble.stats.sd <- aggregate(models.stats.means[,4:ncol(models.stats.means)],
                                 by=list(models.stats.means[,c("var")]),
                                 FUN=sd, na.rm=T)
    names(ensemble.stats.sd)[1] <- "var"
    ensemble.stats.sd
    
    # Making  a publication-style table summarizing periods of statisticlaly significant change in the ecosys drivers
    # NOTE: Standard Deviation based on spatial variability
    library(stringr)
    ensemble.summary <- data.frame(var=ensemble.stats.means$var, 
                                      n.paleo=paste0(str_pad(round(ensemble.stats.means$n.inst.paleo,0), 3, pad=" "), " +/- ", str_pad(round(ensemble.stats.sd$n.inst.paleo, 0), 3, pad=" ")),
                                      rate.paleo.mean = paste0(format(ensemble.stats.means$rate.mean.paleo,digits=3,scientific=T), " +/- ", format(ensemble.stats.sd$rate.mean.paleo, digits=3, scientific=T)),
                                      rate.paleo.max = paste0(format(ensemble.stats.means$rate.max.paleo,digits=3, scientific=T), " +/- ", format(ensemble.stats.sd$rate.max.paleo, digits=3, scientific=T)),
                                      
                                      n.modern=paste0(str_pad(round(ensemble.stats.means$n.inst.modrn,0), 2, pad=" "), " +/- ", str_pad(round(ensemble.stats.sd$n.inst.modrn, 0), 2, pad=" ")),
                                      rate.modern.mean = paste0(format(ensemble.stats.means$rate.mean.modrn,digits=3,scientific=T), " +/- ", format(ensemble.stats.sd$rate.mean.modrn, digits=3, scientific=T)),
                                      rate.modern.max = paste0(format(ensemble.stats.means$rate.max.modrn,digits=3, scientific=T), " +/- ", format(ensemble.stats.sd$rate.max.modrn, digits=3, scientific=T))
    )
    ensemble.summary$var  <- factor(ensemble.summary$var , levels=c("GPP", "NEE", "LAI", "AGB", "Fcomp", "SoilCarb"))
    ensemble.summary      <- ensemble.summary[order(ensemble.summary$var),]
    ensemble.summary
    
    write.csv(ensemble.summary, file.path(out.dir, "Stability_EnsembleSummary.csv"), row.names=F)
  }
  # ------------------
  
  # ------------------
  # 4.3 model traits as predictors of pre-settlement stability
  # ------------------
  {
    library(nlme)
    
    # Create a list with informaiton about the different models
    summary(ecosys.change.mod)
    unique(ecosys.change.mod$Model)
    model.char <- list()
    model.char[["clm.bgc"      ]] <- list(time.step =   30.0, drivers=c("tair", "precipf", "swdown", "lwdown", "qair", "press", "wind", "CO2"      ), veg.scheme="static", pft.scheme="basic"       )
    model.char[["clm.cn"       ]] <- list(time.step =   30.0, drivers=c("tair", "precipf", "swdown", "lwdown", "qair", "press", "wind", "CO2"      ), veg.scheme="static", pft.scheme="basic"       )
    model.char[["ed2"          ]] <- list(time.step =    7.5, drivers=c("tair", "precipf", "swdown", "lwdown", "qair", "press", "wind", "CO2"      ), veg.scheme="dynamic", pft.scheme="succession" )
    model.char[["ed2.lu"       ]] <- list(time.step =    7.5, drivers=c("tair", "precipf", "swdown", "lwdown", "qair", "press", "wind", "CO2", "LU"), veg.scheme="dynamic", pft.scheme="succession" )
    model.char[["jules.stat"   ]] <- list(time.step =   60.0, drivers=c("tair", "precipf", "swdown", "lwdown", "qair", "press", "wind", "CO2"      ), veg.scheme="static" , pft.scheme="basic"      )
    model.char[["jules.triffid"]] <- list(time.step =   60.0, drivers=c("tair", "precipf", "swdown", "lwdown", "qair", "press", "wind", "CO2"      ), veg.scheme="dynamic", pft.scheme="basic"      )
    model.char[["linkages"     ]] <- list(time.step =  5.4e5, drivers=c("tair", "precipf"                                                          ), veg.scheme="dynamic", pft.scheme="species"    )
    model.char[["lpj.guess"    ]] <- list(time.step = 1440.0, drivers=c("tair", "precipf", "swdown",                                    "CO2",  "N"), veg.scheme="dynamic", pft.scheme="bioclimatic")
    model.char[["lpj.wsl"      ]] <- list(time.step = 1440.0, drivers=c("tair", "precipf", "swdown", "lwdown",                          "CO2"      ), veg.scheme="dynamic", pft.scheme="bioclimatic")
    model.char[["sibcasa"      ]] <- list(time.step =     NA, drivers=c("tair", "precipf", "swdown", "lwdown", "qair", "press", "wind", "CO2"      ), veg.scheme="static" , pft.scheme="ecosystem"  )
    
    # Load the stability Stats
    ecosys.stats.site <- read.csv(file.path(out.dir, "Stability_Models_Stats_All.csv"))
    summary(ecosys.stats.site)
    
    # Attaching model information to the stability calcs
    for(m in unique(ecosys.stats.site$Model)){
      ecosys.stats.site[ecosys.stats.site$Model==m, "timestep"  ] <- model.char[[m]]$time.step
      ecosys.stats.site[ecosys.stats.site$Model==m, "veg.scheme"] <- model.char[[m]]$veg.scheme
      ecosys.stats.site[ecosys.stats.site$Model==m, "pft.scheme"] <- model.char[[m]]$pft.scheme
      ecosys.stats.site[ecosys.stats.site$Model==m, "drivers"   ] <- length(model.char[[m]]$drivers)
    }
    ecosys.stats.site$veg.scheme <- as.factor(ecosys.stats.site$veg.scheme)
    ecosys.stats.site$pft.scheme <- as.factor(ecosys.stats.site$pft.scheme)
    summary(ecosys.stats.site)
    
    # Using a mixed model to see if certain characteristics predict the mean rate of paleo change
    # NOTE: -1 assumes default is no change (rate = 0)
    paleo.rate.mean <- lme(abs(rate.mean.paleo) ~ timestep + veg.scheme + (pft.scheme-1) + (var -1) +1, random=list(var=~1, Model=~1, Site=~1), data=ecosys.stats.site[,], na.action=na.omit)
    anova(paleo.rate.mean) # veg scheme, pft scheme
    summary(paleo.rate.mean)
    
    unique(ecosys.stats.site$var)
    paleo.gpp.rate.mean <- lme(abs(rate.mean.paleo) ~ timestep + veg.scheme + (pft.scheme-1) +1, random=list(Model=~1, Site=~1), data=ecosys.stats.site[ecosys.stats.site$var=="GPP",], na.action=na.omit)
    anova(paleo.gpp.rate.mean) # Veg scheme
    summary(paleo.gpp.rate.mean) # Static = less stable
    
    paleo.nee.rate.mean <- lme(abs(rate.mean.paleo) ~ timestep + veg.scheme + (pft.scheme-1) + 1, random=list(Model=~1, Site=~1), data=ecosys.stats.site[ecosys.stats.site$var=="NEE",], na.action=na.omit)
    anova(paleo.nee.rate.mean) # Time step sig
    summary(paleo.nee.rate.mean) # Fast time step = less stable?
    
    paleo.lai.rate.mean <- lme(abs(rate.mean.paleo) ~ timestep + veg.scheme + (pft.scheme-1) +1, random=list(Model=~1, Site=~1), data=ecosys.stats.site[ecosys.stats.site$var=="LAI",], na.action=na.omit)
    anova(paleo.lai.rate.mean) # veg, pft scheme
    summary(paleo.lai.rate.mean) # Static models more stable; succession pft more stable?

    paleo.agb.rate.mean <- lme(abs(rate.mean.paleo) ~ timestep + veg.scheme + (pft.scheme-1) +1, random=list(Model=~1, Site=~1), data=ecosys.stats.site[ecosys.stats.site$var=="AGB",], na.action=na.omit)
    anova(paleo.agb.rate.mean) # time step, veg scheme
    summary(paleo.agb.rate.mean) # Fast time step = less stable?; static veg = more stable
    
    paleo.fcomp.rate.mean <- lme(abs(rate.mean.paleo) ~ timestep + veg.scheme +( pft.scheme-1) +1, random=list(Model=~1, Site=~1), data=ecosys.stats.site[ecosys.stats.site$var=="Fcomp",], na.action=na.omit)
    anova(paleo.fcomp.rate.mean) # time step
    summary(paleo.fcomp.rate.mean) # Fast time step = less stable?; static veg = more stable
    
    paleo.soilc.rate.mean <- lme(abs(rate.mean.paleo) ~ timestep + veg.scheme + (pft.scheme-1) +1, random=list(Model=~1, Site=~1), data=ecosys.stats.site[ecosys.stats.site$var=="SoilCarb",], na.action=na.omit)
    anova(paleo.soilc.rate.mean) # veg scheme, pft scheme
    summary(paleo.fcomp.rate.mean) # static veg = more stable
   
    
    # Using a mixed model to see if certain characteristics predict the number of unstable years
    # NOTE: -1 assumes default is no change (yrs = 0)
    paleo.yrs <- lme(n.inst.paleo ~ timestep + veg.scheme + pft.scheme + var-1, random=list(var=~1, Model=~1, Site=~1), data=ecosys.stats.site[,], na.action=na.omit)
    anova(paleo.yrs) # time step, Veg scheme, pft scheme
    summary(paleo.yrs)
    
    unique(ecosys.stats.site$var)
    paleo.gpp.yrs <- lme(n.inst.paleo ~ timestep + veg.scheme + pft.scheme-1, random=list(Model=~1, Site=~1), data=ecosys.stats.site[ecosys.stats.site$var=="GPP",], na.action=na.omit)
    anova(paleo.gpp.yrs) # time step, veg scheme
    summary(paleo.gpp.yrs) # Static = more stable; slow= more stable
    
    paleo.nee.yrs <- lme(n.inst.paleo ~ timestep + veg.scheme + pft.scheme-1, random=list(Model=~1, Site=~1), data=ecosys.stats.site[ecosys.stats.site$var=="NEE",], na.action=na.omit)
    anova(paleo.nee.yrs) # time step, veg scheme
    summary(paleo.nee.yrs) # Static = more stable; slow= LESS stable
    
    paleo.lai.yrs <- lme(n.inst.paleo ~ timestep + veg.scheme + pft.scheme-1, random=list(Model=~1, Site=~1), data=ecosys.stats.site[ecosys.stats.site$var=="LAI",], na.action=na.omit)
    anova(paleo.lai.yrs) # time, veg, pft scheme, 
    summary(paleo.lai.yrs) # Static = more stable, slow = more stable
    
    paleo.agb.yrs <- lme(n.inst.paleo ~ timestep + veg.scheme + pft.scheme-1, random=list(Model=~1, Site=~1), data=ecosys.stats.site[ecosys.stats.site$var=="AGB",], na.action=na.omit)
    anova(paleo.agb.yrs) # time step, veg scheme
    summary(paleo.agb.yrs) # Static = more stable, slow = more stable
    
    paleo.fcomp.yrs <- lme(n.inst.paleo ~ timestep + veg.scheme + pft.scheme-1, random=list(Model=~1, Site=~1), data=ecosys.stats.site[ecosys.stats.site$var=="Fcomp",], na.action=na.omit)
    anova(paleo.fcomp.yrs) # veg scheme
    summary(paleo.fcomp.yrs) # static = more stable
    
    paleo.soilc.yrs <- lme(n.inst.paleo ~ timestep + veg.scheme + pft.scheme-1, random=list(Model=~1, Site=~1), data=ecosys.stats.site[ecosys.stats.site$var=="SoilCarb",], na.action=na.omit)
    anova(paleo.soilc.yrs) # time, veg, pft scheme, 
    summary(paleo.fcomp.yrs) # static = more stable, slow = LESS stable
    
  }
  # ------------------
}
# ------------------

