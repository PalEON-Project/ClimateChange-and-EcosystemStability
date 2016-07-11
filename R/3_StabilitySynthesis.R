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
library(ggplot2)
library(mgcv)
setwd("~/Dropbox/PalEON_CR/PalEON_MIP_Site/Analyses/Change-and-Stability") # Path to this project github repository: https://github.com/PalEON-Project/MIP-Change-and-Stability.git
# path.gamm.func <- "~/Desktop/R_Functions/"  # Path to github repository of my GAMM helper functions: https://github.com/crollinson/R_Functions.git
inputs    <- "Data/" # Path to my cleaned model output

mip.utils <- "~/Dropbox/PalEON_CR/MIP_Utils/" # Path to PalEON MIP Utility repository: https://github.com/PalEON-Project/MIP_Utils.git

out.dir <- "Data/StabilitySynthesis" # Path to where the analysis output should go
fig.dir <- "Figures/StabilitySynthesis" # Path to where figures should go

if(!dir.exists(out.dir)) dir.create(out.dir)
if(!dir.exists(fig.dir)) dir.create(fig.dir)
# -------------------------------------------

# -------------------------------------------
# 2. Analyzing Met Drivers
# -------------------------------------------
# ------------------
# 2.0 Load met w/ stability analysis; summary figure.
# ------------------
met <- read.csv("Data/Met/StabilityCalcs_SiteMet_850-1850_1900-2010.csv")
met$var <- factor(met$var, levels=c("tair", "precipf", "swdown", "lwdown", "qair", "psurf", "wind"))
met$Site <- factor(met$Site, levels=c("PDL", "PBL", "PUN", "PMB", "PHA", "PHO"))
met[met$var=="tair",c("Y", "mean", "lwr", "upr", "mean.sig")] <- met[met$var=="tair",c("Y", "mean", "lwr", "upr", "mean.sig")]-273.15 # Temp to celcius
met[met$var=="precipf",c("Y", "mean", "lwr", "upr", "deriv.mean", "deriv.lwr", "deriv.upr", "mean.sig")] <- met[met$var=="precipf",c("Y", "mean", "lwr", "upr", "deriv.mean", "deriv.lwr", "deriv.upr", "mean.sig")]*60*60*24*365 # Precip to 
summary(met)

# Graphing the data to make sure it read in okay
for(v in levels(met$var)){
  print(
  ggplot(data=met[met$var==v,]) + 
    facet_grid(var~., scales="free_y") +
    geom_line(aes(x=Year, y=Y, color=Site), size=0.5, alpha=0.3) +
    geom_ribbon(data=met[met$var==v & met$Year<1850,], aes(x=Year, ymin=lwr, ymax=upr, fill=Site), alpha=0.3) +
    geom_ribbon(data=met[met$var==v & met$Year>1900,], aes(x=Year, ymin=lwr, ymax=upr, fill=Site), alpha=0.3) +
    geom_line(data=met[met$var==v & met$Year<1850,], aes(x=Year, y=mean, color=Site), size=1, alpha=0.2) +
    geom_line(data=met[met$var==v & met$Year>1900,], aes(x=Year, y=mean, color=Site), size=1, alpha=0.2) +
    geom_line(aes(x=Year, y=mean.sig, color=Site), size=2, alpha=1) +
    geom_vline(xintercept=1850, linetype="dashed") +
    geom_vline(xintercept=1900, linetype="dashed") +
    scale_x_continuous(expand=c(0,0), name="Year") +
    scale_y_continuous(expand=c(0,0), name="Met") +
    ggtitle(v) +
    # scale_color_manual(values=col.model) +
    # scale_fill_manual(values=col.model) + 
    theme_bw()
  )
}

# Aggregating across sites by year to try and find cohesive patterns to show
ref.window = 850:869 # standardizing everything to the climate of the 20-year spinup
summary(met)

# Making placeholders for our anomaly variables
met[, c("Y.anom", "mean.anom", "lwr.anom", "upr.anom", "mean.sig.anom")] <- NA
for(v in unique(met$var)){
  print(paste0(" **** ", v, " **** "))
  for(s in unique(met.stats.site$Site)){
    # print(paste0(" ---- ", s, " ---- "))
    
    # Find the reference (spinup) state
    ref <- mean(met[met$var==v & met$Site==s & met$Year %in% ref.window, "Y"], na.rm=T)
    met[met$var==v & met$Site==s, c("Y.anom", "mean.anom", "lwr.anom", "upr.anom", "mean.sig.anom")] <- met[met$var==v & met$Site==s, c("Y", "mean", "lwr", "upr", "mean.sig")] - ref
  }
}
summary(met)


# Aggregate met up to the regional average; using range across sites for the upr/lwr
met.yrs <- aggregate(met[,c("Y", "Y.anom", "mean.anom", "lwr.anom", "upr.anom", "deriv.mean", "deriv.lwr", "deriv.upr")],
                     by=met[,c("Year", "var")],
                     FUN=mean, na.rm=T)
met.yrs$n.sig <- aggregate(met[,"sig"],
                           by=met[,c("Year", "var")],
                           FUN=function(x){length(which(x == "*"))})[,3]
summary(met.yrs)

png(file.path(fig.dir, "Stability_Test.png"), height=11, widht=8.5, units="in", res=180)
ggplot(data=met.yrs) +
  facet_grid(var~., scales="free_y") +
  geom_line(aes(x=Year, y=Y.anom), size=0.5, alpha=0.3) +
  geom_ribbon(data=met[met$Year<1850,], aes(x=Year, ymin=lwr.anom, ymax=upr.anom), alpha=0.3) +
  geom_ribbon(data=met[met$Year>1900,], aes(x=Year, ymin=lwr.anom, ymax=upr.anom), alpha=0.3) +
  geom_line(data=met[met$Year<1850,], aes(x=Year, y=mean.anom), size=1, alpha=0.2) +
  geom_line(data=met[met$Year>1900,], aes(x=Year, y=mean.anom), size=1, alpha=0.2) +
  geom_point(aes(x=Year, y=mean.anom, size=as.factor(n.sig)), alpha=1) +
  geom_vline(xintercept=1850, linetype="dashed") +
  geom_vline(xintercept=1900, linetype="dashed") +
  scale_x_continuous(expand=c(0,0), name="Year") +
  scale_y_continuous(expand=c(0,0), name="Met") + 
  scale_size_manual(values=seq(0, 1, length.out=7)) +
  theme_bw()
dev.off()

# ------------------


# ------------------
# 2.1. Site vs. Regional Trends -- change relative to spinup period
# ------------------
#         -- Change going into 1850
#         -- Change at end of simulations
# ------------------

# ------------------
# 2.2. Site Comparisons -- all broken down into pre & post 1850
# ------------------
#         -- mean rate of unstable periods
#         -- mean duration of unstable periods
#         -- percent time spent in unstable state
#         -- Time of greatest instability

# Aggregating to get some of the easy stats
met.stats.site <- merge(data.frame(Site=levels(met$Site)), data.frame(var=levels(met$var)))
summary(met.stats.site)

# Figuring out site-level stability stats
for(v in unique(met.stats.site$var)){
  print(paste0(" **** ", v, " **** "))
  for(s in unique(met.stats.site$Site)){
    print(paste0(" ---- ", s, " ---- "))
    
    # Set up vectors for each site/var
    ins1 <- vector() # Vector for pre-1850 instability
    ins2 <- vector() # vector for post-1900 instability
    dur=0 # Start duration count at 0
  
    # ---------------
    # Go through each year because I don't know how else to do it
    # ---------------
    for(y in min(met$Year):max(met$Year)){ 
      if(!is.na(met[met$var==v & met$Site==s & met$Year==y,"sig"])){
        dur=dur+1 # If this is a period of significant change add it to the year count
        if(y==max(met$Year)){ # If this is our last year and it's instable, record the duration
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
    
    # ---------------
    # Insert data for paleo period
    # ---------------
    if(length(ins1)==0) {
      met.stats.site[met.stats.site$Site==s & met.stats.site$var==v, "n.inst.paleo"] <- length(met[met$Site==s & met$var==v & met$Year<1850 & !is.na(met$sig), "deriv.mean"])
      met.stats.site[met.stats.site$Site==s & met.stats.site$var==v, "rate.mean.paleo"] <- 0
      met.stats.site[met.stats.site$Site==s & met.stats.site$var==v, "rate.sd.paleo"] <- NA

            
      # Saving the mean duration
      met.stats.site[met.stats.site$Site==s & met.stats.site$var==v, "dur.mean.paleo"] <- 0
      met.stats.site[met.stats.site$Site==s & met.stats.site$var==v, "dur.sd.paleo"  ] <- NA

      # Saving the max significant rate of change & time that occurs
      # max.paleo <- which(met$Site==s & met$var==v & met$Year<1850 & !is.na(met$sig) & abs(met$deriv.mean)==abs(max(met[met$Site==s & met$var==v & met$Year<1850 & !is.na(met$sig), "deriv.mean"], na.rm=T)))
      met.stats.site[met.stats.site$Site==s & met.stats.site$var==v, "rate.max.paleo"] <- 0
      met.stats.site[met.stats.site$Site==s & met.stats.site$var==v, "rate.max.paleo.yr"] <- NA
    } else {
      # Finding the number of years in instability
      met.stats.site[met.stats.site$Site==s & met.stats.site$var==v, "n.inst.paleo"] <- length(met[met$Site==s & met$var==v & met$Year<1850 & !is.na(met$sig), "deriv.mean"])
      met.stats.site[met.stats.site$Site==s & met.stats.site$var==v, "rate.mean.paleo"] <- mean(met[met$Site==s & met$var==v & met$Year<1850 & !is.na(met$sig), "deriv.mean"])
      met.stats.site[met.stats.site$Site==s & met.stats.site$var==v, "rate.sd.paleo"] <- sd(met[met$Site==s & met$var==v & met$Year<1850 & !is.na(met$sig), "deriv.mean"])
      
      # Saving the mean duration
      met.stats.site[met.stats.site$Site==s & met.stats.site$var==v, "dur.mean.paleo"] <- mean(ins1, na.rm=T)
      met.stats.site[met.stats.site$Site==s & met.stats.site$var==v, "dur.sd.paleo"] <- sd(ins1, na.rm=T)

      # Saving the max significant rate of change & time that occurs
      max.paleo <- which(met$Site==s & met$var==v & met$Year<1850 & !is.na(met$sig) & abs(met$deriv.mean)==abs(max(met[met$Site==s & met$var==v & met$Year<1850 & !is.na(met$sig), "deriv.mean"], na.rm=T)))
      met.stats.site[met.stats.site$Site==s & met.stats.site$var==v, "rate.max.paleo"] <- met[max.paleo, "deriv.mean"]
      met.stats.site[met.stats.site$Site==s & met.stats.site$var==v, "rate.max.paleo.yr"] <- met[max.paleo, "Year"]
    }
    # ---------------
    
    # ---------------
    # Insert data for Modern period
    # ---------------
    if(length(ins2)==0) {
      met.stats.site[met.stats.site$Site==s & met.stats.site$var==v, "n.inst.modrn"] <- length(met[met$Site==s & met$var==v & met$Year>1900 & !is.na(met$sig), "deriv.mean"])
      met.stats.site[met.stats.site$Site==s & met.stats.site$var==v, "rate.mean.modrn"] <- 0
      met.stats.site[met.stats.site$Site==s & met.stats.site$var==v, "rate.sd.modrn"] <- NA
      
      # Saving the mean duration
      met.stats.site[met.stats.site$Site==s & met.stats.site$var==v, "dur.mean.modrn"] <- 0
      met.stats.site[met.stats.site$Site==s & met.stats.site$var==v, "dur.sd.modrn"  ] <- NA
      
      # Saving the max significant rate of change & time that occurs
      # max.modrn <- which(met$Site==s & met$var==v & met$Year<1850 & !is.na(met$sig) & abs(met$deriv.mean)==abs(max(met[met$Site==s & met$var==v & met$Year<1850 & !is.na(met$sig), "deriv.mean"], na.rm=T)))
      met.stats.site[met.stats.site$Site==s & met.stats.site$var==v, "rate.max.modrn"] <- 0
      met.stats.site[met.stats.site$Site==s & met.stats.site$var==v, "rate.max.modrn.yr"] <- NA
    } else {
      # Finding the number of years in instability
      met.stats.site[met.stats.site$Site==s & met.stats.site$var==v, "n.inst.modrn"] <- length(met[met$Site==s & met$var==v & met$Year>1900 & !is.na(met$sig), "deriv.mean"])
      met.stats.site[met.stats.site$Site==s & met.stats.site$var==v, "rate.mean.modrn"] <- mean(met[met$Site==s & met$var==v & met$Year>1900 & !is.na(met$sig), "deriv.mean"])
      met.stats.site[met.stats.site$Site==s & met.stats.site$var==v, "rate.sd.modrn"] <- sd(met[met$Site==s & met$var==v & met$Year>1900 & !is.na(met$sig), "deriv.mean"])
      
      # Saving the mean duration
      met.stats.site[met.stats.site$Site==s & met.stats.site$var==v, "dur.mean.modrn"] <- mean(ins2, na.rm=T)
      met.stats.site[met.stats.site$Site==s & met.stats.site$var==v, "dur.sd.modrn"] <- sd(ins2, na.rm=T)
      
      # Saving the max significant rate of change & time that occurs
      max.modrn <- which(met$Site==s & met$var==v & met$Year>1900 & !is.na(met$sig) & abs(met$deriv.mean)==abs(max(met[met$Site==s & met$var==v & met$Year>1900 & !is.na(met$sig), "deriv.mean"], na.rm=T)))
      met.stats.site[met.stats.site$Site==s & met.stats.site$var==v, "rate.max.modrn"] <- met[max.modrn, "deriv.mean"]
      met.stats.site[met.stats.site$Site==s & met.stats.site$var==v, "rate.max.modrn.yr"] <- met[max.modrn, "Year"]
    }
    # ---------------
    
  } # End site loop
  print(met.stats.site[met.stats.site$var==v,])
} # End var loop


met.stats.site[met.stats.site$var=="tair",]

# Aggregating to the stats to variable level, leveraging the spatial sd
met.stats.means <- aggregate(met.stats.site[,3:ncol(met.stats.site)], by=list(met.stats.site$var), FUN=mean, na.rm=T)
names(met.stats.means)[1] <- "var"
met.stats.means

met.stats.sd <- aggregate(met.stats.site[,3:ncol(met.stats.site)], by=list(met.stats.site$var), FUN=sd, na.rm=T)
names(met.stats.sd)[1] <- "var"
met.stats.sd

# Making  a publication-style table summarizing periods of statisticlaly significant change in the met drivers
# NOTE: Standard Deviation based on spatial variability
met.summary <- data.frame(var=met.stats.means$var, 
                          n.paleo=paste0(str_pad(round(met.stats.means$n.inst.paleo,0), 3, pad=" "), " +/- ", str_pad(round(met.stats.sd$n.inst.paleo, 0), 3, pad=" ")),
                          rate.paleo.mean = paste0(format(met.stats.means$rate.mean.paleo,digits=3,scientific=T), " +/- ", format(met.stats.sd$rate.mean.paleo, digits=3, scientific=T)),
                          rate.paleo.max = paste0(format(met.stats.means$rate.max.paleo,digits=3, scientific=T), " +/- ", format(met.stats.sd$rate.max.paleo, digits=3, scientific=T)),

                          n.modern=paste0(str_pad(round(met.stats.means$n.inst.modrn,0), 2, pad=" "), " +/- ", str_pad(round(met.stats.sd$n.inst.modrn, 0), 2, pad=" ")),
                          rate.modern.mean = paste0(format(met.stats.means$rate.mean.modrn,digits=3,scientific=T), " +/- ", format(met.stats.sd$rate.mean.modrn, digits=3, scientific=T)),
                          rate.modern.max = paste0(format(met.stats.means$rate.max.modrn,digits=3, scientific=T), " +/- ", format(met.stats.sd$rate.max.modrn, digits=3, scientific=T))
                          )
met.summary$var <- factor(met.summary$var, levels=c("tair", "precipf", "swdown", "lwdown", "qair", "psurf", "wind"))
met.summary <- met.summary[order(met.summary$var),]
met.summary

write.csv(met.summary, file.path(out.dir, "MetSummary_SigChange_AcrossSites.csv"), row.names=F)

# ------------------

# ------------------
# 2.3. Drivers vs. Paleo-Proxy Comparisions
# ------------------
#         -- times of greatest instability
#         -- rates of unstable period
#         -- duration of unstable periods
# ------------------
# -------------------------------------------

# -------------------------------------------
# 3. Analysing Ecosystem Model Output
# -------------------------------------------
# ------------------
#    3.1. Generate a synthesizing figure
# ------------------
#         -- compare models & variables; site probably not important
# ------------------

# ------------------
#    3.2. Calculate average instability statistics by variable & model (pre- & post-1850)
# ------------------
#         -- mean rate of unstable periods
#         -- mean duration of unstable periods
#         -- percent time spent in unstable state
#         -- measure of interannual variability
#              - SD as % of mean?
# ------------------

# ------------------
#    3.3. Compare rates & durations of instability by variable (pre & post-1850)
# ------------------
#         -- are the most variables stable pre-1850 also the most stable post 1900?
#         -- are slow processes/states more or less stable than fast ones?
# ------------------

# ------------------
#    3.4. Compare rates & durations of instability by models
# ------------------
#         -- are the most stable models pre-1850 also the most stable post 1900?
# ------------------

# ------------------
#    3.5. Can we quantify the relationship between interannual variability and ecosystem stability?
# ------------------
#         -- something that applies to both models & variables 
#         -- can we use somethign about interannual variability to predict stability today & in the past?
# ------------------
# -------------------------------------------
