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

# -------------------------------------------
# 2. Graphing Stability Met Drivers
# -------------------------------------------
# ------------------
# 2.0 Met stability analysis & exploratory figures
# ------------------
{
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
  # standardizing everything to the climate of the 20-year spinup since that's what should be stable
  ref.window = 850:869 
  summary(met)
  
  # Making placeholders for our anomaly variables
  met[, c("Y.anom", "mean.anom", "lwr.anom", "upr.anom", "mean.sig.anom")] <- NA
  for(v in unique(met$var)){
    print(paste0(" **** ", v, " **** "))
    for(s in unique(met$Site)){
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
                             FUN=function(x){length(which(x == "*"))})[,"x"]
  summary(met.yrs)
  
  # Add more user-friendly names to figures
  met.yrs$var.name <- recode(met.yrs$var, "'tair'='Temperature'; 'precipf'='Precipitation'; 'swdown'='Shortwave Rad.'; 'lwdown'='Longwave Rad.'; 'qair'='Humidity'; 'psurf'='Pressure'; 'wind'='Wind'")
  met.yrs$var.name <- factor(met.yrs$var.name, levels=c("Temperature", "Precipitation", "Shortwave Rad.", "Longwave Rad.", "Humidity", "Pressure", "Wind"))
  summary(met.yrs)
  
  # Applying a decadal smoothing just to see what our gam looks like relative to that
  library(zoo)
  for(v in unique(met.yrs$var)){
    met.yrs[met.yrs$var==v, "Y.anom.10"] <- rollapply(met.yrs[met.yrs$var==v, "Y.anom"], width=10, FUN=mean, fill=NA)
  }
  summary(met.yrs)
  
  
  summary(met)
  # Doing some recoding to make pretier graphs
  met$Stability <- as.factor(ifelse(is.na(met$sig), "Stable", "Change"))
  met$Stability <- factor(met$Stability, levels=c("Stable", "Change"))
  
  extract.legend.stability <- function(dat.plot, var.plot){
    # Function from: http://stackoverflow.com/questions/16501999/positioning-two-legends-independently-in-a-faceted-ggplot2-plot/17470321#17470321
    g_legend<-function(a.gplot){
      tmp <- ggplot_gtable(ggplot_build(a.gplot))
      leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
      legend <- tmp$grobs[[leg]]
      return(legend)
    }
    
    # Trying to put the size legend above and the color scale on the right:
    # http://stackoverflow.com/questions/13143894/how-do-i-position-two-legends-independently-in-ggplot
    plot.size <- ggplot(data=dat.plot[dat.plot$var==var.plot,]) + # LAI is the only variable to have all 10 models show instability
      facet_grid(var~Site, scales="fixed") +
      geom_point(aes(x=Year, y=mean.anom, size=Stability), alpha=1) +
      theme_bw() +
      scale_size_manual(values=c(1, 5)) +
      guides(size=guide_legend(title="Climate Stability", nrow=1)) +
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
  
  
  plot.top.raw <- function(dat.plot, var.plot, sites){
    plot.top <- ggplot(data=dat.plot[dat.plot$var==var.plot & dat.plot$Site %in% sites,]) +
      theme_bw() +
      facet_grid(var~Site, scales="fixed") +
      geom_line(aes(x=Year, y=Y), size=0.5, alpha=0.3) +
      geom_hline(data=dat.plot[dat.plot$var==var.plot & dat.plot$Site %in% sites & dat.plot$Year<1850,], aes(yintercept=mean(Y)), linetype="dashed") +
      # geom_ribbon(aes(x=Year, ymin=Y.anom.min.10, ymax=Y.anom.max.10), alpha=0.25) +
      geom_ribbon(data=dat.plot[dat.plot$var==var.plot & dat.plot$Site %in% sites & dat.plot$Year<1850,], 
                  aes(x=Year, ymin=lwr, ymax=upr), alpha=0.3) +
      geom_ribbon(data=dat.plot[dat.plot$var==var.plot & dat.plot$Site %in% sites & dat.plot$Year>1900,], 
                  aes(x=Year, ymin=lwr, ymax=upr), alpha=0.3) +
      geom_point(data=dat.plot[dat.plot$var==var.plot & dat.plot$Site %in% sites,], 
                 aes(x=Year, y=mean, color=abs(deriv.mean), size=Stability), alpha=1) +
      # geom_point(data=dat.plot[dat.plot$var==var.plot & dat.plot$Site %in% sites & dat.plot$Year>1900,], 
      # aes(x=Year, y=mean, color=abs(deriv.mean)), size=1, alpha=1) +
      # geom_point(aes(x=Year, y=mean, size=as.factor(sig), color=abs(deriv.mean)), alpha=1) +
      geom_vline(xintercept=1850, linetype="dashed") +
      geom_vline(xintercept=1900, linetype="dashed") +
      scale_x_continuous(expand=c(0,0), name="Year") +
      scale_y_continuous(expand=c(0,0), name=paste0(var.plot)) +
      scale_size_manual(values=c(1, 5)) +
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
  plot.middle.raw <- function(dat.plot, var.plot, sites){
    
    plot.mid <- ggplot(data=dat.plot[dat.plot$var==var.plot & dat.plot$Site %in% sites,]) +
      theme_bw() +
      facet_grid(var~Site, scales="fixed") +
      geom_line(aes(x=Year, y=Y), size=0.5, alpha=0.3) +
      geom_hline(data=dat.plot[dat.plot$var==var.plot & dat.plot$Site %in% sites & dat.plot$Year<1850,], aes(yintercept=mean(Y)), linetype="dashed") +
      # geom_ribbon(aes(x=Year, ymin=Y.anom.min.10, ymax=Y.anom.max.10), alpha=0.25) +
      geom_ribbon(data=dat.plot[dat.plot$var==var.plot & dat.plot$Site %in% sites & dat.plot$Year<1850,], 
                  aes(x=Year, ymin=lwr, ymax=upr), alpha=0.3) +
      geom_ribbon(data=dat.plot[dat.plot$var==var.plot & dat.plot$Site %in% sites & dat.plot$Year>1900,], 
                  aes(x=Year, ymin=lwr, ymax=upr), alpha=0.3) +
      geom_point(data=dat.plot[dat.plot$var==var.plot & dat.plot$Site %in% sites,], 
                 aes(x=Year, y=mean, color=abs(deriv.mean), size=Stability), alpha=1) +
      # geom_point(data=dat.plot[dat.plot$var==var.plot & dat.plot$Site %in% sites & dat.plot$Year>1900,], 
      # aes(x=Year, y=mean, color=abs(deriv.mean)), size=1, alpha=1) +
      # geom_point(aes(x=Year, y=mean, size=as.factor(sig), color=abs(deriv.mean)), alpha=1) +
      geom_vline(xintercept=1850, linetype="dashed") +
      geom_vline(xintercept=1900, linetype="dashed") +
      scale_x_continuous(expand=c(0,0), name="Year") +
      scale_y_continuous(expand=c(0,0), name=paste0(var.plot)) +
      scale_size_manual(values=c(1, 5)) +
      scale_color_gradient(low="gray25", high="red", guide="colorbar", name=paste0("Rate of \nChange")) +
      guides(size=F) +
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
  plot.bottom.raw <- function(dat.plot, var.plot, sites){
    
    plot.bottom <- ggplot(data=dat.plot[dat.plot$var==var.plot & dat.plot$Site %in% sites,]) +
      theme_bw() +
      facet_grid(var~Site, scales="fixed") +
      geom_line(aes(x=Year, y=Y), size=0.5, alpha=0.3) +
      geom_hline(data=dat.plot[dat.plot$var==var.plot & dat.plot$Site %in% sites & dat.plot$Year<1850,], aes(yintercept=mean(Y)), linetype="dashed") +
      # geom_ribbon(aes(x=Year, ymin=Y.anom.min.10, ymax=Y.anom.max.10), alpha=0.25) +
      geom_ribbon(data=dat.plot[dat.plot$var==var.plot & dat.plot$Site %in% sites & dat.plot$Year<1850,], 
                  aes(x=Year, ymin=lwr, ymax=upr), alpha=0.3) +
      geom_ribbon(data=dat.plot[dat.plot$var==var.plot & dat.plot$Site %in% sites & dat.plot$Year>1900,], 
                  aes(x=Year, ymin=lwr, ymax=upr), alpha=0.3) +
      geom_point(data=dat.plot[dat.plot$var==var.plot & dat.plot$Site %in% sites,], 
                 aes(x=Year, y=mean, color=abs(deriv.mean), size=Stability), alpha=1) +
      # geom_point(data=dat.plot[dat.plot$var==var.plot & dat.plot$Site %in% sites & dat.plot$Year>1900,], 
      # aes(x=Year, y=mean, color=abs(deriv.mean)), size=1, alpha=1) +
      # geom_point(aes(x=Year, y=mean, size=as.factor(sig), color=abs(deriv.mean)), alpha=1) +
      geom_vline(xintercept=1850, linetype="dashed") +
      geom_vline(xintercept=1900, linetype="dashed") +
      scale_x_continuous(expand=c(0,0), name="Year") +
      scale_y_continuous(expand=c(0,0), name=paste0(var.plot)) +
      scale_size_manual(values=c(1, 5)) +
      scale_color_gradient(low="gray25", high="red", guide="colorbar", name=paste0("Rate of \nChange")) +
      guides(size=F) +
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
  
  legend.size <- extract.legend.stability(dat.plot=met, var.plot="tair")
  # plot.blank1 <- plot.blank(dat.plot=ecosys.change.site)
  
  for(s in unique(ecosys.change.site$Site)){
    print(paste0(" ** ", s))
    
    # Putting all the vars together
    plot.tair <- plot.top.raw(dat.plot=met, var="tair", sites=s)
    plot.precipf <- plot.middle.raw(dat.plot=met, var="precipf", sites=s)
    plot.swdown <- plot.middle.raw(dat.plot=met, var="swdown", sites=s)
    plot.lwdown <- plot.middle.raw(dat.plot=met, var="lwdown", sites=s)
    plot.qair <- plot.middle.raw(dat.plot=met, var="qair", sites=s)
    plot.psurf <- plot.middle.raw(dat.plot=met, var="psurf", sites=s)
    plot.wind <- plot.bottom.raw(dat.plot=met, var="wind", sites=s)
    
    # A second graph of lwdown so it can also be the bottom panel
    plot.lwdown2 <- plot.bottom.raw(dat.plot=met, var="lwdown", sites=s)
    # 
    # plot.soilc2 <- ggplotGrob(plot.soilc )
    # plot.soilc2$grobs
    # # plot.soilc2$grobs[[9]]$heights = unit.c(unit(0.1, "null"), unit(20, "npc"), unit(0.1, "null"))
    # # plot.soilc2$heights[[6]] <- unit(0.3, "null")
    # # plot.soilc2$widths[[6]] <- unit(10, "null")
    # plot(plot.soilc2)
    # if(s == "PHA"){
    mar.plot <- list(tair   =theme(plot.margin=unit(c(0.1, 1.0, 0.4, 2.2), "lines")),
                     precipf=theme(plot.margin=unit(c(0.1, 1.0, 0.4, 1.1), "lines")),
                     swdown =theme(plot.margin=unit(c(0.1, 1.0, 0.4, 1.5), "lines")),
                     lwdown =theme(plot.margin=unit(c(0.1, 0.8, 0.4, 1.5), "lines")),
                     qair   =theme(plot.margin=unit(c(0.1, 0.75,0.4, 0.6), "lines")),
                     psurf  =theme(plot.margin=unit(c(0.1, 1.0, 0.4, 0.7), "lines")),
                     wind   =theme(plot.margin=unit(c(0.1, 0.35, 0.5, 1.55), "lines")),
                     lwdown2=theme(plot.margin=unit(c(0.1, 0.8, 0.5, 1.5), "lines"))
    )
    # }
    
    
    
    png(file.path(fig.dir, paste0("Stability_Met_Site_",s,".png")), height=11.5, width=8.5, units="in", res=180)
    print(
      plot_grid(legend.size, 
                plot.tair + mar.plot$tair, 
                plot.precipf + mar.plot$precipf, 
                plot.swdown + mar.plot$swdown, 
                plot.lwdown + mar.plot$lwdown, 
                plot.qair + mar.plot$qair, 
                plot.psurf + mar.plot$psurf, 
                plot.wind + mar.plot$wind, 
                ncol=1, nrow=8, rel_heights = c(0.05, rep(0.9/6, 6), 0.9/6+0.05))
    )
    dev.off()
    
    png(file.path(fig.dir, paste0("Stability_Met_Site_",s,"_ClimRad.png")), height=8.5, width=8.5, units="in", res=180)
    print(
      plot_grid(legend.size, 
                plot.tair + mar.plot$tair, 
                plot.precipf + mar.plot$precipf, 
                plot.swdown + mar.plot$swdown, 
                plot.lwdown2 + mar.plot$lwdown2, 
                ncol=1, nrow=5, rel_heights = c(0.3, 1.1, 1, 1, 1.3))
    )
    dev.off()
    
  }
  
  
  
  png(file.path(fig.dir, "Stability_Met_Region_4Vars.png"), height=11, width=8.5, units="in", res=180)
  print(
    ggplot(data=met.yrs[met.yrs$var %in% c("tair", "precipf", "swdown", "lwdown"),]) +
      facet_grid(var.name~., scales="free_y") +
      # facet_wrap(~var, scales="free_y") +
      geom_line(aes(x=Year, y=Y.anom), size=0.5, alpha=0.3) +
      geom_ribbon(data=met.yrs[met.yrs$var %in% c("tair", "precipf", "swdown", "lwdown") & met.yrs$Year<1850,], aes(x=Year, ymin=lwr.anom, ymax=upr.anom), alpha=0.3) +
      geom_ribbon(data=met.yrs[met.yrs$var %in% c("tair", "precipf", "swdown", "lwdown") & met.yrs$Year>1900,], aes(x=Year, ymin=lwr.anom, ymax=upr.anom), alpha=0.3) +
      geom_line(data=met.yrs[met.yrs$var %in% c("tair", "precipf", "swdown", "lwdown") & met.yrs$Year<1850,], aes(x=Year, y=mean.anom), size=1, alpha=0.2) +
      geom_line(data=met.yrs[met.yrs$var %in% c("tair", "precipf", "swdown", "lwdown") & met.yrs$Year>1900,], aes(x=Year, y=mean.anom), size=1, alpha=0.2) +
      geom_point(aes(x=Year, y=mean.anom, size=as.factor(n.sig)), alpha=1) +
      geom_vline(xintercept=1850, linetype="dashed") +
      geom_vline(xintercept=1900, linetype="dashed") +
      scale_x_continuous(expand=c(0,0), name="Year") +
      scale_y_continuous(expand=c(0,0), name="Deviation from Spinup State") + 
      scale_size_manual(values=seq(0, 3, length.out=length(unique(met.yrs$n.sig)))) +
      guides(size=guide_legend(nrow=1, title="# Sites Showing Change")) +
      theme_bw() +
      theme(legend.position="top")
  )
  dev.off()
  
  
  
  png(file.path(fig.dir, "Stability_Met_Region_AllVars.png"), width=11, height=8.5, units="in", res=180)
  print(
    ggplot(data=met.yrs) +
      facet_grid(var~., scales="free_y") +
      # facet_wrap(~var, scales="free_y") +
      geom_line(aes(x=Year, y=Y.anom), size=0.5, alpha=0.3) +
      geom_ribbon(data=met.yrs[met.yrs$Year<1850,], aes(x=Year, ymin=lwr.anom, ymax=upr.anom), alpha=0.3) +
      geom_ribbon(data=met.yrs[met.yrs$Year>1900,], aes(x=Year, ymin=lwr.anom, ymax=upr.anom), alpha=0.3) +
      geom_line(data=met.yrs[met.yrs$Year<1850,], aes(x=Year, y=mean.anom), size=1, alpha=0.2) +
      geom_line(data=met.yrs[met.yrs$Year>1900,], aes(x=Year, y=mean.anom), size=1, alpha=0.2) +
      geom_point(aes(x=Year, y=mean.anom, size=as.factor(n.sig)), alpha=1) +
      geom_vline(xintercept=1850, linetype="dashed") +
      geom_vline(xintercept=1900, linetype="dashed") +
      scale_x_continuous(expand=c(0,0), name="Year") +
      scale_y_continuous(expand=c(0,0), name="met.yrs") + 
      scale_size_manual(values=seq(0, 3, length.out=length(unique(met.yrs$n.sig)))) +
      guides(size=guide_legend(nrow=1, title="# Sites Showing Change")) +
      theme_bw() +
      theme(legend.position="top")
  )
  dev.off()
  
  
  
  
  
}
# ------------------


# ------------------
# 2.2. Site Comparisons -- all broken down into pre & post 1850
# ------------------
#         -- mean rate of unstable periods
#         -- mean duration of unstable periods
#         -- percent time spent in unstable state
#         -- Time of greatest instability
{
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
}
# ------------------
