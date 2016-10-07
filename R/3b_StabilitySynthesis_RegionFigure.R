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
#    2.1. Generate a synthesizing figure
# ------------------
#         -- compare models & variables; site probably not important
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
}
# ------------------

# ------------------
# Actually making the figure
# ------------------
{
  # --------
  # setting up the graphing functions
  # --------
  {
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
      plot.size <- ggplot(data=dat.plot[dat.plot$var=="AGB",]) + # LAI is the only variable to have all 10 models show instability
        facet_grid(var~., scales="fixed") +
        geom_point(aes(x=Year, y=mean.anom, size=as.factor(mod.sig)), alpha=1) +
        theme_bw() +
        scale_size_manual(values=seq(0, 4, length.out=11), breaks=seq(0,10, by=1)) +
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
    extract.legend.color <- function(dat.plot, var.plot){
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
      plot.color <- ggplot(data=dat.plot[dat.plot$var==var.plot,]) +
        facet_grid(var~., scales="fixed") +
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
      plot.top <- ggplot(data=dat.plot[dat.plot$var==var.plot,]) +
        theme_bw() +
        facet_grid(var~., scales="fixed") +
        # geom_line(aes(x=Year, y=Y.anom), size=0.5, alpha=0.3) +
        geom_ribbon(aes(x=Year, ymin=Y.anom.min.10, ymax=Y.anom.max.10), alpha=0.25) +
        geom_ribbon(data=dat.plot[dat.plot$var==var.plot & dat.plot$Year<1850,], 
                    aes(x=Year, ymin=lwr.anom, ymax=upr.anom), alpha=0.3) +
        geom_ribbon(data=dat.plot[dat.plot$var==var.plot & dat.plot$Year>1900,], 
                    aes(x=Year, ymin=lwr.anom, ymax=upr.anom), alpha=0.3) +
        geom_line(data=dat.plot[dat.plot$var==var.plot & dat.plot$Year<1850,], 
                  aes(x=Year, y=mean.anom), size=1, alpha=0.2) +
        geom_line(data=dat.plot[dat.plot$var==var.plot & dat.plot$Year>1900,], 
                  aes(x=Year, y=mean.anom), size=1, alpha=0.2) +
        geom_point(aes(x=Year, y=mean.anom, size=as.factor(mod.sig), color=abs(deriv.mean)), alpha=1) +
        geom_vline(xintercept=1850, linetype="dashed") +
        geom_vline(xintercept=1900, linetype="dashed") +
        scale_x_continuous(expand=c(0,0), name="Year") +
        scale_y_continuous(expand=c(0,0), name=paste0("Deviation from Spinup")) +
        scale_size_manual(values=seq(0, 4, length.out=11), breaks=seq(0,10, by=1)) +
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
      
      plot.mid <- ggplot(data=dat.plot[dat.plot$var==var.plot,]) +
        facet_grid(var~., scales="fixed") +
        # geom_line(aes(x=Year, y=Y.anom), size=0.5, alpha=0.3) +
        geom_ribbon(aes(x=Year, ymin=Y.anom.min.10, ymax=Y.anom.max.10), alpha=0.25) +
        geom_ribbon(data=dat.plot[dat.plot$var==var.plot & dat.plot$Year<1850,], 
                    aes(x=Year, ymin=lwr.anom, ymax=upr.anom), alpha=0.3) +
        geom_ribbon(data=dat.plot[dat.plot$var==var.plot & dat.plot$Year>1900,], 
                    aes(x=Year, ymin=lwr.anom, ymax=upr.anom), alpha=0.3) +
        geom_line(data=dat.plot[dat.plot$var==var.plot & dat.plot$Year<1850,], 
                  aes(x=Year, y=mean.anom), size=1, alpha=0.2) +
        geom_line(data=dat.plot[dat.plot$var==var.plot & dat.plot$Year>1900,], 
                  aes(x=Year, y=mean.anom), size=1, alpha=0.2) +
        geom_point(aes(x=Year, y=mean.anom, size=as.factor(mod.sig), color=abs(deriv.mean)), alpha=1) +
        geom_vline(xintercept=1850, linetype="dashed") +
        geom_vline(xintercept=1900, linetype="dashed") +
        scale_x_continuous(expand=c(0,0), name="Year") +
        scale_y_continuous(expand=c(0,0), name=paste0("Deviation from Spinup")) +
        scale_size_manual(values=seq(0, 4, length.out=11), breaks=seq(0,10, by=1)) +
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
      
      plot.bottom <- ggplot(data=dat.plot[dat.plot$var==var.plot,]) +
        facet_grid(var~., scales="fixed") +
        # geom_line(aes(x=Year, y=Y.anom), size=0.5, alpha=0.3) +
        geom_ribbon(aes(x=Year, ymin=Y.anom.min.10, ymax=Y.anom.max.10), alpha=0.25) +
        geom_ribbon(data=dat.plot[dat.plot$var==var.plot & dat.plot$Year<1850,], 
                    aes(x=Year, ymin=lwr.anom, ymax=upr.anom), alpha=0.3) +
        geom_ribbon(data=dat.plot[dat.plot$var==var.plot & dat.plot$Year>1900,], 
                    aes(x=Year, ymin=lwr.anom, ymax=upr.anom), alpha=0.3) +
        geom_line(data=dat.plot[dat.plot$var==var.plot & dat.plot$Year<1850,], 
                  aes(x=Year, y=mean.anom), size=1, alpha=0.2) +
        geom_line(data=dat.plot[dat.plot$var==var.plot & dat.plot$Year>1900,], 
                  aes(x=Year, y=mean.anom), size=1, alpha=0.2) +
        geom_point(aes(x=Year, y=mean.anom, size=as.factor(mod.sig), color=abs(deriv.mean)), alpha=1) +
        geom_vline(xintercept=1850, linetype="dashed") +
        geom_vline(xintercept=1900, linetype="dashed") +
        scale_x_continuous(expand=c(0,0), name="Year") +
        scale_y_continuous(expand=c(0,0), name=paste0("Deviation from Spinup")) +
        scale_size_manual(values=seq(0, 4, length.out=11), breaks=seq(0,10, by=1)) +
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
      
      return(plot.bottom)
    }
  }
  # --------
  
  # --------
  # Doing the graphs
  # --------
  legend.size <- extract.legend.size(dat.plot=ecosys.change.region)
  plot.gpp <- plot.top.site(dat.plot=ecosys.change.region, var.plot="GPP")
  plot.nee <- plot.middle.site(dat.plot=ecosys.change.region, var.plot="NEE")
  plot.lai <- plot.middle.site(dat.plot=ecosys.change.region, var.plot="LAI")
  plot.agb <- plot.middle.site(dat.plot=ecosys.change.region, var.plot="AGB")
  plot.fcomp <- plot.middle.site(dat.plot=ecosys.change.region, var.plot="Fcomp")
  plot.soilc <- plot.bottom.site(dat.plot=ecosys.change.region, var.plot="SoilCarb")
  
  png(file.path(fig.dir, paste0("Stability_Region.png")), height=11.5, width=8.5, units="in", res=180)
  print(
    plot_grid(legend.size, 
              plot.gpp   + theme(plot.margin=unit(c(0.1,1.0,0.4,1.0), "lines")), 
              plot.nee   + theme(plot.margin=unit(c(0.1,1.0,0.4,0.75), "lines")), 
              plot.lai   + theme(plot.margin=unit(c(0.1,1.0,0.4,1.5), "lines")), 
              plot.agb   + theme(plot.margin=unit(c(0.1,1.25,0.4,1.1), "lines")), 
              plot.fcomp + theme(plot.margin=unit(c(0.1,0.9,0.4,0.6), "lines")), 
              plot.soilc + theme(plot.margin=unit(c(0.1,1.2,0.5,0.8), "lines")), 
              ncol=1, nrow=7, rel_heights = c(0.05, rep(0.9/5, 5), 0.9/5+0.05))
  )
  dev.off()
  
  # --------
  
}
# ------------------
