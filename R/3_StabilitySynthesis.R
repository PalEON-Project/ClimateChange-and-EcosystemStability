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
setwd("~/Dropbox/PalEON_CR/PalEON_MIP_Site/Analyses/Change-and-Stability") # Path to this project github repository: https://github.com/PalEON-Project/MIP-Change-and-Stability.git
# path.gamm.func <- "~/Desktop/R_Functions/"  # Path to github repository of my GAMM helper functions: https://github.com/crollinson/R_Functions.git
inputs    <- "Data/" # Path to my cleaned model output

mip.utils <- "~/Dropbox/PalEON_CR/MIP_Utils/" # Path to PalEON MIP Utility repository: https://github.com/PalEON-Project/MIP_Utils.git

out.dir <- "Data/StabilitySynthesis" # Path to where the analysis output should go
fig.dir <- "Figures/StabilitySynthesis" # Path to where figures should go

if(!dir.exists(out.dir)) dir.create(out.dir)
if(!dir.exists(fig.dir)) dir.create(fig.dir)
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
# 3.0. Generate a synthesizing figure
# ------------------
#         -- compare models & variables; site probably not important
{
  # Averaging across models to be able to look at the ensemble variability of 
  # Note: min/max weren't working, so we'll take the 0 and 100 percentiles
  ecosys.change.region <- aggregate(ecosys.change.mod[,c("Y", "Y.anom", "mean.anom", "lwr.anom", "upr.anom", "deriv.mean", "n.sig")],
                                    by=ecosys.change.mod[,c("Year", "var")],
                                    FUN=mean, na.rm=T)
  ecosys.change.region[,c("Y.anom.min", "n.min")] <- aggregate(ecosys.change.mod[,c("Y.anom", "n.sig")],
                                                                          by=ecosys.change.mod[,c("Year", "var")],
                                                                          FUN=quantile, 0, na.rm=T)[,c("Y.anom", "n.sig")]
  ecosys.change.region[,c("Y.anom.max", "n.max")] <- aggregate(ecosys.change.mod[,c("Y.anom", "n.sig")],
                                                                          by=ecosys.change.mod[,c("Year", "var")],
                                                                          FUN=quantile, 1, na.rm=T)[,c("Y.anom", "n.sig")]
  ecosys.change.region$mod.sig <- aggregate(ecosys.change.mod[,"mod.sig"],
                                            by=ecosys.change.mod[,c("Year", "var")],
                                            FUN=function(x){length(which(x == "*"))})[,"x"]
 
  # doing some smoothing to make cleaner graphs
  library(zoo)
  for(v in unique(ecosys.change.region$var)){
      ecosys.change.region[ecosys.change.region$var==v, "Y.anom.10"] <- rollapply(ecosys.change.region[ecosys.change.region$var==v, "Y.anom"], width=10, FUN=mean, fill=NA)
      ecosys.change.region[ecosys.change.region$var==v, "Y.anom.min.10"] <- rollapply(ecosys.change.region[ecosys.change.region$var==v, "Y.anom.min"], width=10, FUN=mean, fill=NA)
      ecosys.change.region[ecosys.change.region$var==v, "Y.anom.max.10"] <- rollapply(ecosys.change.region[ecosys.change.region$var==v, "Y.anom.max"], width=10, FUN=mean, fill=NA)
  }
  summary(ecosys.change.region)
  
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
  
  # ------------------
  # Figuring out the most useful way to graph thigns
  # ------------------
{  # Setting up the panels we'll need to go by variable and allow variable scaling
  # Order: (may need to swithc AGB & Fcomp, but that shouldn't matter)
  # 1. GPP
  # 2. NEE
  # 3. LAI
  # 4. AGB
  # 5. Fcomp
  # 6. Soil Carb
  
  # Graphing 
  ggplot(data=ecosys.change.region[,]) +
    facet_grid(var~., scales="free_y") +
    scale_x_continuous(expand=c(0,0), name="Year") +
    scale_y_continuous(expand=c(0,0), name=paste0("Deviation from Spinup")) +
    # geom_line(aes(x=Year, y=Y.anom), size=0.5, alpha=0.3) +
    geom_ribbon(aes(x=Year, ymin=Y.anom.min.10, ymax=Y.anom.max.10), alpha=0.25) +
    geom_hline(yintercept=0, linetype="dashed") +
    geom_ribbon(data=ecosys.change.region[ ecosys.change.region$Year<1850,], 
                aes(x=Year, ymin=lwr.anom, ymax=upr.anom), alpha=0.3) +
    geom_ribbon(data=ecosys.change.region[ecosys.change.region$Year>1900,], 
                aes(x=Year, ymin=lwr.anom, ymax=upr.anom), alpha=0.3) +
    geom_line(data=ecosys.change.region[ecosys.change.region$Year<1850,], 
              aes(x=Year, y=mean.anom), size=1, alpha=0.2) +
    geom_line(data=ecosys.change.region[ecosys.change.region$Year>1900,], 
              aes(x=Year, y=mean.anom), size=1, alpha=0.2) +
    geom_point(aes(x=Year, y=mean.anom, size=as.factor(mod.sig), color=abs(deriv.mean)), alpha=1) +
    geom_vline(xintercept=1850, linetype="dashed") +
    geom_vline(xintercept=1900, linetype="dashed") +
    scale_size_manual(name="# Models Showing Change", values=seq(0, 4, length.out=11), breaks=seq(0,10, by=1)) +
    scale_color_gradient(low="gray25", high="red", guide="colorbar", name=paste0("Rate of \nChange")) +
    theme_bw() +
    theme(legend.position="top")
  
  
  
  library(gridExtra); library(gtable)
  library(cowplot)
  
  extract.legend.size <- function(dat.plot){
    # Function from: http://stackoverflow.com/questions/16501999/positioning-two-legends-independently-in-a-faceted-ggplot2-plot/17470321#17470321
    g_legend<-function(a.gplot){
      tmp <- ggplot_gtable(ggplot_build(a.gplot))
      leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
      legend <- tmp$grobs[[leg]]
      return(legend)
    }
    
    # Trying to put the size legend above and the color scale on the right:
    # http://stackoverflow.com/questions/13143894/how-do-i-position-two-legends-independently-in-ggplot
    plot.size <- ggplot(data=dat.plot[dat.plot$var=="LAI",]) + # LAI is the only variable to have all 10 models show instability
      facet_grid(var~Site, scales="fixed") +
      geom_point(aes(x=Year, y=mean.anom, size=as.factor(n.sig)), alpha=1) +
      theme_bw() +
      scale_size_manual(values=seq(0, 5, length.out=11), breaks=seq(0,10, by=1)) +
      guides(size=guide_legend(title="# Models Showing Change", nrow=1)) +
      theme(legend.position="top",
            legend.key=element_blank(),
            legend.key.height=unit(2, "lines")
      ) +
      theme(plot.margin=unit(c(1,1,0,1), "lines"))
    
    # extract size legend
    legend.size <- g_legend(plot.size)
    
    return(legend.size)
  }
  extract.legend.color <- function(dat.plot, var.plot, sites){
    # Function from: http://stackoverflow.com/questions/16501999/positioning-two-legends-independently-in-a-faceted-ggplot2-plot/17470321#17470321
    g_legend<-function(a.gplot){
      tmp <- ggplot_gtable(ggplot_build(a.gplot))
      leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
      legend <- tmp$grobs[[leg]]
      return(legend)
    }
    # 
    # Trying to put the size legend above and the color scale on the right:
    # http://stackoverflow.com/questions/13143894/how-do-i-position-two-legends-independently-in-ggplot
    plot.color <- ggplot(data=dat.plot[dat.plot$var==var.plot & dat.plot$Site %in% sites,]) +
      facet_grid(var~Site, scales="fixed") +
      geom_point(aes(x=Year, y=mean.anom, color=abs(deriv.mean)), alpha=1) +
      scale_color_gradient(low="gray25", high="red", guide="colorbar", name=paste0("Rate of \nChange")) +
      guides(colorbar=guide_legend(title.theme=element_text(margin=c(0,5,0,0))))  +
      theme_bw() +
      theme(legend.key.height=unit(0.8, "lines")) +
      theme(axis.text.x=element_blank(),
            axis.title.x=element_blank(),
            axis.text.y=element_text(margin=margin(0,10,0,0)),
            axis.title.y=element_text(margin=margin(0,10,0,0)),
            axis.ticks.length=unit(-0.5, unit="lines")) +
      theme(strip.text.x=element_blank()) +
      theme(plot.margin=unit(c(0,1,0.5,1), "lines"))
    
    # extract size legend
    legend.color <- g_legend(plot.color)
    # legend.color <- heights=
    # Take the legend from p1
    # legend.color <- gtable_filter(ggplot_gtable(ggplot_build(plot.color)), "guide-box") 
    # legGrob <- grobTree(legend.color)
    
    return(legend.color)
  }
  plot.blank <- function(dat.plot){
    plot.blank <- ggplot(data=dat.plot) +
      geom_blank()
    
    return(plot.blank)
  }
  plot.top.site <- function(dat.plot, var.plot, sites){
    plot.top <- ggplot(data=dat.plot[dat.plot$var==var.plot & dat.plot$Site %in% sites,]) +
      theme_bw() +
      facet_grid(var~Site, scales="fixed") +
      # geom_line(aes(x=Year, y=Y.anom), size=0.5, alpha=0.3) +
      geom_ribbon(aes(x=Year, ymin=Y.anom.min.10, ymax=Y.anom.max.10), alpha=0.25) +
      geom_ribbon(data=dat.plot[dat.plot$var==var.plot & dat.plot$Site %in% sites & dat.plot$Year<1850,], 
                  aes(x=Year, ymin=lwr.anom, ymax=upr.anom), alpha=0.3) +
      geom_ribbon(data=dat.plot[dat.plot$var==var.plot & dat.plot$Site %in% sites & dat.plot$Year>1900,], 
                  aes(x=Year, ymin=lwr.anom, ymax=upr.anom), alpha=0.3) +
      geom_line(data=dat.plot[dat.plot$var==var.plot & dat.plot$Site %in% sites & dat.plot$Year<1850,], 
                aes(x=Year, y=mean.anom), size=1, alpha=0.2) +
      geom_line(data=dat.plot[dat.plot$var==var.plot & dat.plot$Site %in% sites & dat.plot$Year>1900,], 
                aes(x=Year, y=mean.anom), size=1, alpha=0.2) +
      geom_point(aes(x=Year, y=mean.anom, size=as.factor(n.sig), color=abs(deriv.mean)), alpha=1) +
      geom_vline(xintercept=1850, linetype="dashed") +
      geom_vline(xintercept=1900, linetype="dashed") +
      scale_x_continuous(expand=c(0,0), name="Year") +
      scale_y_continuous(expand=c(0,0), name=paste0("Deviation from Spinup")) +
      scale_size_manual(values=seq(0, 5, length.out=11), breaks=seq(0,10, by=1)) +
      scale_color_gradient(low="gray25", high="red", guide="colorbar", name=paste0("Rate of \nChange")) +
      guides(size=F) +
      theme(legend.key.height=unit(0.8, "lines")) +
      theme(axis.text.x=element_blank(),
            axis.title.x=element_blank(),
            axis.text.y=element_text(margin=margin(0,10,0,0)),
            axis.title.y=element_text(margin=margin(0,10,0,0), size=rel(0.8)),
            axis.ticks.length=unit(-0.5, unit="lines")) +
      theme(plot.margin=unit(c(0.1,1,0.5,1), "lines"))
    
    return(plot.top)
  }
  plot.middle.site <- function(dat.plot, var.plot, sites){
    
    plot.mid <- ggplot(data=dat.plot[dat.plot$var==var.plot & dat.plot$Site %in% sites,]) +
      facet_grid(var~Site, scales="fixed") +
      # geom_line(aes(x=Year, y=Y.anom), size=0.5, alpha=0.3) +
      geom_ribbon(aes(x=Year, ymin=Y.anom.min.10, ymax=Y.anom.max.10), alpha=0.25) +
      geom_ribbon(data=dat.plot[dat.plot$var==var.plot & dat.plot$Site %in% sites & dat.plot$Year<1850,], 
                  aes(x=Year, ymin=lwr.anom, ymax=upr.anom), alpha=0.3) +
      geom_ribbon(data=dat.plot[dat.plot$var==var.plot & dat.plot$Site %in% sites & dat.plot$Year>1900,], 
                  aes(x=Year, ymin=lwr.anom, ymax=upr.anom), alpha=0.3) +
      geom_line(data=dat.plot[dat.plot$var==var.plot & dat.plot$Site %in% sites & dat.plot$Year<1850,], 
                aes(x=Year, y=mean.anom), size=1, alpha=0.2) +
      geom_line(data=dat.plot[dat.plot$var==var.plot & dat.plot$Site %in% sites & dat.plot$Year>1900,], 
                aes(x=Year, y=mean.anom), size=1, alpha=0.2) +
      geom_point(aes(x=Year, y=mean.anom, size=as.factor(n.sig), color=abs(deriv.mean)), alpha=1) +
      geom_vline(xintercept=1850, linetype="dashed") +
      geom_vline(xintercept=1900, linetype="dashed") +
      scale_x_continuous(expand=c(0,0), name="Year") +
      scale_y_continuous(expand=c(0,0), name=paste0("Deviation from Spinup")) +
      scale_size_manual(values=seq(0, 5, length.out=11), breaks=seq(0,10, by=1)) +
      scale_color_gradient(low="gray25", high="red", guide="colorbar", name=paste0("Rate of \nChange")) +
      guides(size=F) +
      theme_bw() +
      theme(legend.key.height=unit(0.8, "lines")) +
      theme(axis.text.x=element_blank(),
            axis.title.x=element_blank(),
            axis.text.y=element_text(margin=margin(0,10,0,0)),
            axis.title.y=element_text(margin=margin(0,10,0,0), size=rel(0.8)),
            axis.ticks.length=unit(-0.5, unit="lines")) +
      theme(strip.text.x=element_blank()) +
      theme(plot.margin=unit(c(0,1,0.5,1), "lines"))
    
    return(plot.mid)
  }
  plot.bottom.site <- function(dat.plot, var.plot, sites){
    
    plot.mid <- ggplot(data=dat.plot[dat.plot$var==var.plot & dat.plot$Site %in% sites,]) +
      facet_grid(var~Site, scales="fixed") +
      # geom_line(aes(x=Year, y=Y.anom), size=0.5, alpha=0.3) +
      geom_ribbon(aes(x=Year, ymin=Y.anom.min.10, ymax=Y.anom.max.10), alpha=0.25) +
      geom_ribbon(data=dat.plot[dat.plot$var==var.plot & dat.plot$Site %in% sites & dat.plot$Year<1850,], 
                  aes(x=Year, ymin=lwr.anom, ymax=upr.anom), alpha=0.3) +
      geom_ribbon(data=dat.plot[dat.plot$var==var.plot & dat.plot$Site %in% sites & dat.plot$Year>1900,], 
                  aes(x=Year, ymin=lwr.anom, ymax=upr.anom), alpha=0.3) +
      geom_line(data=dat.plot[dat.plot$var==var.plot & dat.plot$Site %in% sites & dat.plot$Year<1850,], 
                aes(x=Year, y=mean.anom), size=1, alpha=0.2) +
      geom_line(data=dat.plot[dat.plot$var==var.plot & dat.plot$Site %in% sites & dat.plot$Year>1900,], 
                aes(x=Year, y=mean.anom), size=1, alpha=0.2) +
      geom_point(aes(x=Year, y=mean.anom, size=as.factor(n.sig), color=abs(deriv.mean)), alpha=1) +
      geom_vline(xintercept=1850, linetype="dashed") +
      geom_vline(xintercept=1900, linetype="dashed") +
      scale_x_continuous(expand=c(0,0), name="Year") +
      scale_y_continuous(expand=c(0,0), name=paste0("Deviation from Spinup")) +
      scale_size_manual(values=seq(0, 5, length.out=11), breaks=seq(0,10, by=1)) +
      scale_color_gradient(low="gray25", high="red", guide="colorbar", name=paste0("Rate of \nChange")) +
      guides(size=F) +
      theme_bw() +
      theme(legend.key.height=unit(0.8, "lines")) +
      theme(axis.text.x=element_text(margin=margin(10,0,0,0)),
            axis.title.x=element_text(margin=margin(10,0,0,0)),
            axis.text.y=element_text(margin=margin(0,10,0,0)),
            axis.title.y=element_text(margin=margin(0,10,0,0), size=rel(0.8)),
            axis.ticks.length=unit(-0.5, unit="lines")) +
      theme(strip.text.x=element_blank()) +
      theme(plot.margin=unit(c(0,1,1,1), "lines"))
    
    return(plot.mid)
  }
  
  legend.size <- extract.legend.size(dat.plot=ecosys.change.site)
  # plot.blank1 <- plot.blank(dat.plot=ecosys.change.site)
  
  for(s in unique(ecosys.change.site$Site)){
    print(paste0(" ** ", s))
    
    plot.gpp <- plot.top.site(dat.plot=ecosys.change.site, var="GPP", sites=s)
    lgnd.gpp <- extract.legend.color(dat.plot=ecosys.change.site, var="GPP", sites=s)
    plot.nee <- plot.middle.site(dat.plot=ecosys.change.site, var="NEE", sites=s)
    lgnd.nee <- extract.legend.color(dat.plot=ecosys.change.site, var="NEE", sites=s)
    plot.lai <- plot.middle.site(dat.plot=ecosys.change.site, var="LAI", sites=s)
    lgnd.lai <- extract.legend.color(dat.plot=ecosys.change.site, var="LAI", sites=s)
    plot.agb <- plot.middle.site(dat.plot=ecosys.change.site, var="AGB", sites=s)
    lgnd.agb <- extract.legend.color(dat.plot=ecosys.change.site, var="AGB", sites=s)
    plot.fcomp <- plot.middle.site(dat.plot=ecosys.change.site, var="Fcomp", sites=s)
    lgnd.fcomp <- extract.legend.color(dat.plot=ecosys.change.site, var="Fcomp", sites=s)
    plot.soilc <- plot.bottom.site(dat.plot=ecosys.change.site, var="SoilCarb", sites=s)
    lgnd.soilc <- extract.legend.color(dat.plot=ecosys.change.site, var="SoilCarb", sites=s)
    
    # 
    # plot.soilc2 <- ggplotGrob(plot.soilc )
    # plot.soilc2$grobs
    # # plot.soilc2$grobs[[9]]$heights = unit.c(unit(0.1, "null"), unit(20, "npc"), unit(0.1, "null"))
    # # plot.soilc2$heights[[6]] <- unit(0.3, "null")
    # # plot.soilc2$widths[[6]] <- unit(10, "null")
    # plot(plot.soilc2)
    if(s == "PHA"){
      mar.plot <- list(gpp=theme(plot.margin=unit(c(0.1,1,0.4,1), "lines")),
                       nee=theme(plot.margin=unit(c(0.1,1.3,0.4,0.75), "lines")),
                       lai=theme(plot.margin=unit(c(0.1,0.95,0.4,0.75), "lines")),
                       agb=theme(plot.margin=unit(c(0.1,1.5,0.4,0.9), "lines")),
                       fcomp=theme(plot.margin=unit(c(0.1,1.2,0.4,0.7), "lines")),
                       soilc=theme(plot.margin=unit(c(0.1,1.5,0.5,0.85), "lines")))
    } else if(s=="PHO"){
      mar.plot <- list(gpp=theme(plot.margin=unit(c(0.1,1.3,0.4,1), "lines")),
                       nee=theme(plot.margin=unit(c(0.1,0.9,0.4,0.75), "lines")),
                       lai=theme(plot.margin=unit(c(0.1,0.9,0.4,1.3), "lines")),
                       agb=theme(plot.margin=unit(c(0.1,1.5,0.4,0.9), "lines")),
                       fcomp=theme(plot.margin=unit(c(0.1,1.2,0.4,0.7), "lines")),
                       soilc=theme(plot.margin=unit(c(0.1,1.5,0.5,0.85), "lines")))
    } else if(s=="PUN"){
      mar.plot <- list(gpp=theme(plot.margin=unit(c(0.1,1,0.4,1), "lines")),
                       nee=theme(plot.margin=unit(c(0.1,0.65,0.4,0.75), "lines")),
                       lai=theme(plot.margin=unit(c(0.1,1.4,0.4,1.3), "lines")),
                       agb=theme(plot.margin=unit(c(0.1,1.65,0.4,0.9), "lines")),
                       fcomp=theme(plot.margin=unit(c(0.1,1.05,0.4,0.7), "lines")),
                       soilc=theme(plot.margin=unit(c(0.1,1.70,0.5,0.85), "lines")))
    } else if(s=="PMB"){
      mar.plot <- list(gpp=theme(plot.margin=unit(c(0.1,0.95,0.4,1), "lines")),
                       nee=theme(plot.margin=unit(c(0.1,0.9,0.4,0.8), "lines")),
                       lai=theme(plot.margin=unit(c(0.1,1,0.4,1.35), "lines")),
                       agb=theme(plot.margin=unit(c(0.1,1.25,0.4,0.9), "lines")),
                       fcomp=theme(plot.margin=unit(c(0.1,1.0,0.4,0.7), "lines")),
                       soilc=theme(plot.margin=unit(c(0.1,1.25,0.5,0.85), "lines")))
    } else if(s=="PBL"){
      mar.plot <- list(gpp=theme(plot.margin=unit(c(0.1,0.9,0.4,1.05), "lines")),
                       nee=theme(plot.margin=unit(c(0.1,0.9,0.4,0.8), "lines")),
                       lai=theme(plot.margin=unit(c(0.1,0.95,0.4,1.3), "lines")),
                       agb=theme(plot.margin=unit(c(0.1,1.2,0.4,1.25), "lines")),
                       fcomp=theme(plot.margin=unit(c(0.1,0.6,0.4,0.9), "lines")),
                       soilc=theme(plot.margin=unit(c(0.1,1.2,0.5,0.8), "lines")))
    } else if(s=="PDL"){
      mar.plot <- list(gpp=theme(plot.margin=unit(c(0.1,0.7,0.4,0.8), "lines")),
                       nee=theme(plot.margin=unit(c(0.1,0.7,0.4,0.6), "lines")),
                       lai=theme(plot.margin=unit(c(0.1,0.7,0.4,1.3), "lines")),
                       agb=theme(plot.margin=unit(c(0.1,0.95,0.4,0.9), "lines")),
                       fcomp=theme(plot.margin=unit(c(0.1,0.3,0.4,0.7), "lines")),
                       soilc=theme(plot.margin=unit(c(0.1,0.9,0.5,0.85), "lines")))
    }
    
    
    png(file.path(fig.dir, paste0("Stability_Site_",s,".png")), height=11.5, width=8.5, units="in", res=180)
    print(
      plot_grid(legend.size, 
                plot.gpp + mar.plot$gpp, 
                plot.nee + mar.plot$nee, 
                plot.lai + mar.plot$lai, 
                plot.agb + mar.plot$agb, 
                plot.fcomp + mar.plot$fcomp, 
                plot.soilc + mar.plot$soilc, 
                ncol=1, nrow=7, rel_heights = c(0.05, rep(0.9/5, 5), 0.9/5+0.05))
    )
    dev.off()
  }
  
  
  # set up plots grid
  plot.gpp2 <- plot.top.site(dat.plot=ecosys.change.site, var="GPP", sites=unique(ecosys.change.site$Site))
  # lgnd.gpp <- extract.legend.color(dat.plot=ecosys.change.site, var="GPP", sites="PHA")
  plot.nee2 <- plot.middle.site(dat.plot=ecosys.change.site, var="NEE", sites=unique(ecosys.change.site$Site))
  # lgnd.nee <- extract.legend.color(dat.plot=ecosys.change.site, var="NEE", sites="PHA")
  plot.lai2 <- plot.middle.site(dat.plot=ecosys.change.site, var="LAI", sites=unique(ecosys.change.site$Site))
  # lgnd.lai <- extract.legend.color(dat.plot=ecosys.change.site, var="LAI", sites="PHA")
  plot.agb2 <- plot.middle.site(dat.plot=ecosys.change.site, var="AGB", sites=unique(ecosys.change.site$Site))
  # lgnd.agb <- extract.legend.color(dat.plot=ecosys.change.site, var="AGB", sites="PHA")
  plot.fcomp2 <- plot.middle.site(dat.plot=ecosys.change.site, var="Fcomp", sites=unique(ecosys.change.site$Site))
  # lgnd.fcomp <- extract.legend.color(dat.plot=ecosys.change.site, var="Fcomp", sites="PHA")
  plot.soilc2 <- plot.bottom.site(dat.plot=ecosys.change.site, var="SoilCarb", sites=unique(ecosys.change.site$Site))
  # lgnd.soilc <- extract.legend.color(dat.plot=ecosys.change.site, var="SoilCarb", sites="PHA")
  png(file.path(fig.dir, paste0("Stability_Site_All.png")), height=11.5, width=15.5, units="in", res=220)
  plot_grid(legend.size, plot.gpp2, plot.nee2, plot.lai2, plot.agb2, plot.fcomp2, plot.soilc2, ncol=1, nrow=7, rel_heights = c(0.03, rep(0.9/5, 5), 0.9/5+0.07))
  dev.off()
  
  # Graphing by variable  
  for(v in unique(ecosys.change.mod$var)){
    paste0("** ",v)
    png(file.path(fig.dir, paste0("Stability_Var_",v,".png")), height=8.5, width=8.5, units="in", res=180)
    print(
      ggplot(data=ecosys.change.mod[ecosys.change.mod$var==v,]) +
        facet_wrap(~Model.Order, scales="fixed") +
        geom_line(aes(x=Year, y=Y.anom), size=0.5, alpha=0.3) +
        geom_ribbon(data=ecosys.change.mod[ecosys.change.mod$var==v & ecosys.change.mod$Year<1850,], aes(x=Year, ymin=lwr.anom, ymax=upr.anom), alpha=0.3) +
        geom_ribbon(data=ecosys.change.mod[ecosys.change.mod$var==v & ecosys.change.mod$Year>1900,], aes(x=Year, ymin=lwr.anom, ymax=upr.anom), alpha=0.3) +
        geom_line(data=ecosys.change.mod[ecosys.change.mod$var==v & ecosys.change.mod$Year<1850,], aes(x=Year, y=mean.anom), size=1, alpha=0.2) +
        geom_line(data=ecosys.change.mod[ecosys.change.mod$var==v & ecosys.change.mod$Year>1900,], aes(x=Year, y=mean.anom), size=1, alpha=0.2) +
        geom_point(aes(x=Year, y=mean.anom, size=as.factor(n.sig), color=abs(deriv.mean)), alpha=1) +
        geom_vline(xintercept=1850, linetype="dashed") +
        geom_vline(xintercept=1900, linetype="dashed") +
        scale_x_continuous(expand=c(0,0), name="Year") +
        scale_y_continuous(expand=c(0,0), name=paste0(v, " Deviation from Spinup")) +
        scale_size_manual(values=seq(0, 3, length.out=length(unique(ecosys.change.mod$n.sig)))) +
        scale_color_gradient(low="gray25", high="red", guide="colorbar", name=paste0("Rate of \nChange\n(", v,"/yr)")) +
        guides(size=guide_legend(title="# Sites \nShowing \nChange")) +
        ggtitle(v) +
        theme_bw() +
        theme(plot.title=element_text(face="bold", size=rel(2)))
    )
    dev.off()
    
    png(file.path(fig.dir, paste0("Stability_Var_",v,"_FreeY.png")), height=8.5, width=8.5, units="in", res=180)
    print(
      ggplot(data=ecosys.change.mod[ecosys.change.mod$var==v,]) +
        facet_wrap(~Model.Order, scales="free_y") +
        geom_line(aes(x=Year, y=Y.anom), size=0.5, alpha=0.3) +
        geom_ribbon(data=ecosys.change.mod[ecosys.change.mod$var==v & ecosys.change.mod$Year<1850,], aes(x=Year, ymin=lwr.anom, ymax=upr.anom), alpha=0.3) +
        geom_ribbon(data=ecosys.change.mod[ecosys.change.mod$var==v & ecosys.change.mod$Year>1900,], aes(x=Year, ymin=lwr.anom, ymax=upr.anom), alpha=0.3) +
        geom_line(data=ecosys.change.mod[ecosys.change.mod$var==v & ecosys.change.mod$Year<1850,], aes(x=Year, y=mean.anom), size=1, alpha=0.2) +
        geom_line(data=ecosys.change.mod[ecosys.change.mod$var==v & ecosys.change.mod$Year>1900,], aes(x=Year, y=mean.anom), size=1, alpha=0.2) +
        geom_point(aes(x=Year, y=mean.anom, size=as.factor(n.sig), color=abs(deriv.mean)), alpha=1) +
        geom_vline(xintercept=1850, linetype="dashed") +
        geom_vline(xintercept=1900, linetype="dashed") +
        scale_x_continuous(expand=c(0,0), name="Year") +
        scale_y_continuous(expand=c(0,0), name=paste0(v, " Deviation from Spinup")) +
        scale_size_manual(values=seq(0, 3, length.out=length(unique(ecosys.change.mod$n.sig)))) +
        scale_color_gradient(low="gray25", high="red", guide="colorbar", name=paste0("Rate of \nChange\n(", v,"/yr)")) +
        guides(size=guide_legend(title="# Sites \nShowing \nChange")) +
        ggtitle(v) +
        theme_bw() +
        theme(plot.title=element_text(face="bold", size=rel(2)))  
    )
    dev.off()
  }
  
  # Graphing by Model
  for(m in unique(ecosys.change.mod$Model.Order)){
    paste0("** ",m)
    png(file.path(fig.dir, paste0("Stability_Model_",m,".png")), height=8.5, width=8.5, units="in", res=220)
    print(
      ggplot(data=ecosys.change.mod[ecosys.change.mod$Model.Order==m,]) +
        # facet_wrap(~var, scales="free_y") +
        facet_grid(var~., scales="free_y") +
        geom_line(aes(x=Year, y=Y.anom), size=0.5, alpha=0.3) +
        geom_ribbon(data=ecosys.change.mod[ecosys.change.mod$Model.Order==m & ecosys.change.mod$Year<1850,], aes(x=Year, ymin=lwr.anom, ymax=upr.anom), alpha=0.3) +
        geom_ribbon(data=ecosys.change.mod[ecosys.change.mod$Model.Order==m & ecosys.change.mod$Year>1900,], aes(x=Year, ymin=lwr.anom, ymax=upr.anom), alpha=0.3) +
        geom_line(data=ecosys.change.mod[ecosys.change.mod$Model.Order==m & ecosys.change.mod$Year<1850,], aes(x=Year, y=mean.anom), size=1, alpha=0.2) +
        geom_line(data=ecosys.change.mod[ecosys.change.mod$Model.Order==m & ecosys.change.mod$Year>1900,], aes(x=Year, y=mean.anom), size=1, alpha=0.2) +
        geom_point(aes(x=Year, y=mean.anom, size=as.factor(n.sig), color=abs(deriv.mean)), alpha=1) +
        geom_vline(xintercept=1850, linetype="dashed") +
        geom_vline(xintercept=1900, linetype="dashed") +
        scale_x_continuous(expand=c(0,0), name="Year") +
        scale_y_continuous(expand=c(0,0), name=paste0(m, " Deviation from Spinup")) +
        scale_size_manual(values=seq(0, 3, length.out=length(unique(ecosys.change.mod$n.sig)))) +
        scale_color_gradient(low="gray25", high="red", guide="colorbar", name=paste0("Rate of \nChange\n(", m,"/yr)")) +
        guides(size=guide_legend(title="# Sites \nShowing \nChange")) +
        ggtitle(m) +
        theme_bw() +
        theme(plot.title=element_text(face="bold", size=rel(2)))
    )
    dev.off()
  }
  
  for(v in unique(ecosys.change.site$var)){
    paste0("** ",v)
    png(file.path(fig.dir, paste0("Stability_Var_",v,"_Sites.png")), height=8.5, width=8.5, units="in", res=180)
    print(
      ggplot(data=ecosys.change.site[ecosys.change.site$var==v,]) +
        facet_wrap(~Site, scales="fixed") +
        geom_line(aes(x=Year, y=Y.anom), size=0.5, alpha=0.3) +
        geom_ribbon(data=ecosys.change.site[ecosys.change.site$var==v & ecosys.change.site$Year<1850,], aes(x=Year, ymin=lwr.anom, ymax=upr.anom), alpha=0.3) +
        geom_ribbon(data=ecosys.change.site[ecosys.change.site$var==v & ecosys.change.site$Year>1900,], aes(x=Year, ymin=lwr.anom, ymax=upr.anom), alpha=0.3) +
        geom_line(data=ecosys.change.site[ecosys.change.site$var==v & ecosys.change.site$Year<1850,], aes(x=Year, y=mean.anom), size=1, alpha=0.2) +
        geom_line(data=ecosys.change.site[ecosys.change.site$var==v & ecosys.change.site$Year>1900,], aes(x=Year, y=mean.anom), size=1, alpha=0.2) +
        geom_point(aes(x=Year, y=mean.anom, size=as.factor(n.sig), color=abs(deriv.mean)), alpha=1) +
        geom_vline(xintercept=1850, linetype="dashed") +
        geom_vline(xintercept=1900, linetype="dashed") +
        scale_x_continuous(expand=c(0,0), name="Year") +
        scale_y_continuous(expand=c(0,0), name=paste0(v, " Deviation from Spinup")) +
        scale_size_manual(values=seq(0, 3, length.out=length(unique(ecosys.change.site$n.sig)))) +
        scale_color_gradient(low="gray25", high="red", guide="colorbar", name=paste0("Mean \nRate of \nChange")) +
        guides(size=guide_legend(title="# Models \nShowing \nChange")) +
        ggtitle(v) +
        theme_bw() +
        theme(plot.title=element_text(face="bold", size=rel(2)))
    )
    dev.off()
  }
  
  # Showing all variables by site
  png(file.path(fig.dir, paste0("Stability_Site_AllSites.png")), height=8.5, width=11.5, units="in", res=220)
  print(
    ggplot(data=ecosys.change.site[,]) +
      facet_grid(var~Site, scales="free_y") +
      geom_line(aes(x=Year, y=Y.anom), size=0.5, alpha=0.3) +
      geom_ribbon(data=ecosys.change.site[ecosys.change.site$Site==i & ecosys.change.site$Year<1850,], aes(x=Year, ymin=lwr.anom, ymax=upr.anom), alpha=0.3) +
      geom_ribbon(data=ecosys.change.site[ecosys.change.site$Site==i & ecosys.change.site$Year>1900,], aes(x=Year, ymin=lwr.anom, ymax=upr.anom), alpha=0.3) +
      geom_line(data=ecosys.change.site[ecosys.change.site$Site==i & ecosys.change.site$Year<1850,], aes(x=Year, y=mean.anom), size=1, alpha=0.2) +
      geom_line(data=ecosys.change.site[ecosys.change.site$Site==i & ecosys.change.site$Year>1900,], aes(x=Year, y=mean.anom), size=1, alpha=0.2) +
      geom_point(aes(x=Year, y=mean.anom, size=as.factor(n.sig), color=abs(deriv.mean)), alpha=1) +
      geom_vline(xintercept=1850, linetype="dashed") +
      geom_vline(xintercept=1900, linetype="dashed") +
      scale_x_continuous(expand=c(0,0), name="Year") +
      scale_y_continuous(expand=c(0,0), name=paste0(" Deviation from Spinup")) +
      scale_size_manual(values=seq(0, 3, length.out=length(unique(ecosys.change.site$n.sig)))) +
      scale_color_gradient(low="gray25", high="red", guide="colorbar", name=paste0("Mean \nRate of \nChange")) +
      guides(size=guide_legend(title="# Models \nShowing \nChange")) +
      # ggtitle(i) +
      theme_bw() +
      theme(plot.title=element_text(face="bold", size=rel(2)))
  )
  dev.off()
  
  
  for(i in unique(ecosys.change.site$Site)){
    paste0("** ",i)
    png(file.path(fig.dir, paste0("Stability_Site_",i,".png")), height=8.5, width=8.5, units="in", res=180)
    print(
      ggplot(data=ecosys.change.site[ecosys.change.site$Site==i,]) +
        facet_grid(var~Site, scales="free_y") +
        geom_line(aes(x=Year, y=Y.anom), size=0.5, alpha=0.3) +
        geom_ribbon(data=ecosys.change.site[ecosys.change.site$Site==i & ecosys.change.site$Year<1850,], aes(x=Year, ymin=lwr.anom, ymax=upr.anom), alpha=0.3) +
        geom_ribbon(data=ecosys.change.site[ecosys.change.site$Site==i & ecosys.change.site$Year>1900,], aes(x=Year, ymin=lwr.anom, ymax=upr.anom), alpha=0.3) +
        geom_line(data=ecosys.change.site[ecosys.change.site$Site==i & ecosys.change.site$Year<1850,], aes(x=Year, y=mean.anom), size=1, alpha=0.2) +
        geom_line(data=ecosys.change.site[ecosys.change.site$Site==i & ecosys.change.site$Year>1900,], aes(x=Year, y=mean.anom), size=1, alpha=0.2) +
        geom_point(aes(x=Year, y=mean.anom, size=as.factor(n.sig), color=abs(deriv.mean)), alpha=1) +
        geom_vline(xintercept=1850, linetype="dashed") +
        geom_vline(xintercept=1900, linetype="dashed") +
        scale_x_continuous(expand=c(0,0), name="Year") +
        scale_y_continuous(expand=c(0,0), name=paste0(" Deviation from Spinup")) +
        scale_size_manual(values=seq(0, 3, length.out=length(unique(ecosys.change.site$n.sig)))) +
        scale_color_gradient(low="gray25", high="red", guide="colorbar", name=paste0("Mean \nRate of \nChange")) +
        guides(size=guide_legend(title="# Models \nShowing \nChange")) +
        ggtitle(i) +
        theme_bw() +
        theme(plot.title=element_text(face="bold", size=rel(2)))
    )
    dev.off()
  }
  }
  # ------------------
  
}
# ------------------

# ------------------
# 4. Analyzing sensitivity 
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
    
    # Figuring out site-level stability stats
    for(v in unique(ecosys.stats.site$var)){
      print(paste0(" **** ", v, " **** "))
      for(s in unique(ecosys.stats.site$Site)){
        print(paste0(" ---- ", s, " ---- "))
        
        for(m in unique(ecosys.stats.site$Model)){
          # Set up vectors for each site/var
          ins1 <- vector() # Vector for pre-1850 instability
          ins2 <- vector() # vector for post-1900 instability
          dur=0 # Start duration count at 0
          
          # ---------------
          # Go through each year to find duration of instability periods because I don't know how else to do it
          # ---------------
          for(y in min(ecosys.change$Year):max(ecosys.change$Year)){ 
            if(!is.na(ecosys.change[ecosys.change$var==v & ecosys.change$Site==s  & ecosys.change$Model==m & ecosys.change$Year==y,"sig"])){
              dur=dur+1 # If this is a period of significant change add it to the year count
              if(y==max(ecosys.change$Year)){ # If this is our last year and it's instable, record the duration
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
            ecosys.stats.site[ecosys.stats.site$Site==s & ecosys.stats.site$var==v & ecosys.stats.site$Model==m, "dur.mean.paleo"] <- mean(ins1, na.rm=T)
            ecosys.stats.site[ecosys.stats.site$Site==s & ecosys.stats.site$var==v & ecosys.stats.site$Model==m, "dur.sd.paleo"] <- sd(ins1, na.rm=T)
            
            # Saving the max significant rate of change & time that occurs
            max.paleo <- which(ecosys.change$Site==s & ecosys.change$var==v & ecosys.change$Model==m & ecosys.change$Year<1850 & !is.na(ecosys.change$sig) & abs(ecosys.change$deriv.mean)==abs(max(ecosys.change[ecosys.change$Site==s & ecosys.change$var==v & ecosys.change$Model==m & ecosys.change$Year<1850 & !is.na(ecosys.change$sig), "deriv.mean"], na.rm=T)))
            ecosys.stats.site[ecosys.stats.site$Site==s & ecosys.stats.site$var==v & ecosys.stats.site$Model==m, "rate.max.paleo"] <- ecosys.change[max.paleo, "deriv.mean"]
            ecosys.stats.site[ecosys.stats.site$Site==s & ecosys.stats.site$var==v & ecosys.stats.site$Model==m, "rate.max.paleo.yr"] <- ecosys.change[max.paleo, "Year"]
          }
          # ---------------
          
          # ---------------
          # Insert data for Modern period
          # ---------------
          if(length(ins2)==0) {
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
            ecosys.stats.site[ecosys.stats.site$Site==s & ecosys.stats.site$var==v & ecosys.stats.site$Model==m, "dur.mean.modrn"] <- mean(ins2, na.rm=T)
            ecosys.stats.site[ecosys.stats.site$Site==s & ecosys.stats.site$var==v & ecosys.stats.site$Model==m, "dur.sd.modrn"] <- sd(ins2, na.rm=T)
            
            # Saving the max significant rate of change & time that occurs
            max.modrn <- which(ecosys.change$Site==s & ecosys.change$var==v & ecosys.change$Model==m & ecosys.change$Year>1900 & !is.na(ecosys.change$sig) & abs(ecosys.change$deriv.mean)==abs(max(ecosys.change[ecosys.change$Site==s & ecosys.change$var==v & ecosys.change$Model==m & ecosys.change$Year>1900 & !is.na(ecosys.change$sig), "deriv.mean"], na.rm=T)))
            ecosys.stats.site[ecosys.stats.site$Site==s & ecosys.stats.site$var==v & ecosys.stats.site$Model==m, "rate.max.modrn"] <- ecosys.change[max.modrn, "deriv.mean"]
            ecosys.stats.site[ecosys.stats.site$Site==s & ecosys.stats.site$var==v & ecosys.stats.site$Model==m, "rate.max.modrn.yr"] <- ecosys.change[max.modrn, "Year"]
          }
          # ---------------
        }  # End Model loop
      } # End site loop
      print(ecosys.stats.site[ecosys.stats.site$var==v,])
    } # End var loop
    
    
    ecosys.stats.site[ecosys.stats.site$var=="tair",]
    
    # Aggregating to the stats to variable level, leveraging the spatial sd
    ecosys.stats.means <- aggregate(ecosys.stats.site[,3:ncol(ecosys.stats.site)], by=list(ecosys.stats.site$var), FUN=mean, na.rm=T)
    names(ecosys.stats.means)[1] <- "var"
    ecosys.stats.means
    
    ecosys.stats.sd <- aggregate(ecosys.stats.site[,3:ncol(ecosys.stats.site)], by=list(ecosys.stats.site$var), FUN=sd, na.rm=T)
    names(ecosys.stats.sd)[1] <- "var"
    ecosys.stats.sd
    
    # Making  a publication-style table summarizing periods of statisticlaly significant change in the ecosys drivers
    # NOTE: Standard Deviation based on spatial variability
    ecosys.stab.summary <- data.frame(var=ecosys.stats.means$var, 
                                      n.paleo=paste0(str_pad(round(ecosys.stats.means$n.inst.paleo,0), 3, pad=" "), " +/- ", str_pad(round(ecosys.stats.sd$n.inst.paleo, 0), 3, pad=" ")),
                                      rate.paleo.mean = paste0(format(ecosys.stats.means$rate.mean.paleo,digits=3,scientific=T), " +/- ", format(ecosys.stats.sd$rate.mean.paleo, digits=3, scientific=T)),
                                      rate.paleo.max = paste0(format(ecosys.stats.means$rate.max.paleo,digits=3, scientific=T), " +/- ", format(ecosys.stats.sd$rate.max.paleo, digits=3, scientific=T)),
                                      
                                      n.modern=paste0(str_pad(round(ecosys.stats.means$n.inst.modrn,0), 2, pad=" "), " +/- ", str_pad(round(ecosys.stats.sd$n.inst.modrn, 0), 2, pad=" ")),
                                      rate.modern.mean = paste0(format(ecosys.stats.means$rate.mean.modrn,digits=3,scientific=T), " +/- ", format(ecosys.stats.sd$rate.mean.modrn, digits=3, scientific=T)),
                                      rate.modern.max = paste0(format(ecosys.stats.means$rate.max.modrn,digits=3, scientific=T), " +/- ", format(ecosys.stats.sd$rate.max.modrn, digits=3, scientific=T))
    )
    ecosys.stab.summary$var  <- factor(ecosys.stab.summary$var , levels=c("GPP", "NEE", "LAI", "AGB", "Fcomp", "SoilCarb"))
    ecosys.stab.summary$Site <- factor(ecosys.stab.summary$Site, levels=c("PDL", "PBL", "PUN", "PMB", "PHA", "PHO"))
    ecosys.stab.summary      <- ecosys.stab.summary[order(ecosys.stab.summary$var),]
    ecosys.stab.summary
    
    write.csv(ecosys.stab.summary, file.path(out.dir, "EcosysSummary_SigChange_AcrossSites.csv"), row.names=F)
  }
  # ------------------
  
  # ------------------
  # 4.2 max rate of change pre/post
  # ------------------
  # ------------------
  
  # ------------------
  # 4.3 model traits as predictors of pre-settlement stability
  # ------------------
  # ------------------
}
# ------------------

