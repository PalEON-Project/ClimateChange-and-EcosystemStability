library(ggplot2); library(gridExtra); library(grid); library(scales)
setwd("~/Dropbox/PalEON_CR/PalEON_MIP_Site/Analyses/Change-and-Stability")

# -------------------------------------------
# Figure 1: Regional Met Change
# -------------------------------------------
{

fig.out <- "Figures/Met"
met.region <- read.csv("Data/Met/Met_Region_Annual.csv")
paleon.states <- map_data("state")

hips <- data.frame(Site = c( "PHA",  "PHO",  "PUN",  "PBL",  "PDL",  "PMB"),
                   lat  = c( 42.54,  45.25,  46.22,  46.28,  47.17,  43.61),
                   lon  = c(-72.18, -68.73, -89.53, -94.58, -95.17, -82.83),
                   lat2 = c( 42.75,  45.25,  46.25,  46.25,  47.25,  43.75),
                   lon2 = c(-72.25, -68.75, -89.75, -94.75, -95.25, -82.75))

tair.mean <- ggplot(data=met.region[met.region$Time == "Modern",]) +
  facet_grid(Time ~ .) +
  geom_raster(aes(x=lon, y=lat, fill=tair)) +
  geom_path(data=paleon.states, aes(x=long, y=lat, group=group), size=0.2) +
  geom_point(data=hips, aes(x=lon, y=lat), color="gray50", size=3) +
  scale_x_continuous(limits=range(met.region$lon), expand=c(0,0), breaks=seq(-95, -60, by=10), name="Longitude") +
  scale_y_continuous(limits=range(met.region$lat), expand=c(0,0), breaks=seq( 37.5,  50, by= 4), name="Latitude") +
  scale_fill_gradient(low="blue4", high="red2", name="Tair (C)") +
  coord_fixed(ratio=1) +
  theme(legend.position=c(0.85, 0.3),
        legend.background=element_blank(),
        legend.key.height=unit(0.3, "lines"),
        legend.text=element_text(size=unit(8, "points")),
        legend.title=element_text(size=unit(8, "points")),
        panel.background=element_blank(),
        panel.border=element_rect(size=0.25, fill=NA),
        axis.ticks.length=unit(-0.3, "lines"),
        axis.text.x=element_blank(),
        axis.title.x=element_blank(),
        axis.text.y=element_text(margin=unit(c(0.5, 1.5, 0.5, 0.5), "lines")),
        strip.text=element_blank(),
        plot.margin=unit(c(0.5, 0.5, 0.5, 1), "lines"))

tair.diff1 <- ggplot(data=met.region[met.region$Time == "diff.Settlement",]) +
  facet_grid(Time ~ .) +
  geom_raster(aes(x=lon, y=lat, fill=tair)) +
  geom_path(data=paleon.states, aes(x=long, y=lat, group=group), size=0.2) +
  geom_point(data=hips, aes(x=lon, y=lat), color="gray50", size=3) +
  scale_x_continuous(limits=range(met.region$lon), expand=c(0,0), breaks=seq(-95, -60, by=10), name="Longitude") +
  scale_y_continuous(limits=range(met.region$lat), expand=c(0,0), breaks=seq( 37.5,  50, by= 4), name="Latitude") +
  scale_fill_gradient2(low="blue4", high="red2", name="Tair (C)", midpoint=0, limits=range(met.region[substr(met.region$Time,1,4)=="diff", "tair"])) +
  coord_fixed(ratio=1) +
  theme(legend.position=c(0.85, 0.3),
        legend.background=element_blank(),
        legend.key.height=unit(0.3, "lines"),
        legend.text=element_text(size=unit(6, "points")),
        legend.title=element_text(size=unit(8, "points")),
        panel.background=element_blank(),
        panel.border=element_rect(size=0.25, fill=NA),
        axis.ticks.length=unit(-0.3, "lines"),
        axis.text.x=element_blank(),
        axis.title.x=element_blank(),
        axis.text.y=element_text(margin=unit(c(0.5, 1.5, 0.5, 0.5), "lines")),
        strip.text=element_blank(),
        plot.margin=unit(c(0.5, 0.5, 0.5, 1), "lines"))

tair.diff2 <- ggplot(data=met.region[met.region$Time == "diff.Spinup",]) +
  facet_grid(Time ~ .) +
  geom_raster(aes(x=lon, y=lat, fill=tair)) +
  geom_path(data=paleon.states, aes(x=long, y=lat, group=group), size=0.2) +
  geom_point(data=hips, aes(x=lon, y=lat), color="gray50", size=3) +
  scale_x_continuous(limits=range(met.region$lon), expand=c(0,0), breaks=seq(-95, -60, by=10), name="Longitude") +
  scale_y_continuous(limits=range(met.region$lat), expand=c(0,0), breaks=seq( 37.5,  50, by= 4), name="Latitude") +
  scale_fill_gradient2(low="blue4", high="red2", name="Tair (C)", midpoint=0, limits=range(met.region[substr(met.region$Time,1,4)=="diff", "tair"])) +
  coord_fixed(ratio=1) +
  theme(legend.position=c(0.85, 0.3),
        legend.background=element_blank(),
        legend.key.height=unit(0.3, "lines"),
        legend.text=element_text(size=unit(6, "points")),
        legend.title=element_text(size=unit(8, "points")),
        panel.background=element_blank(),
        panel.border=element_rect(size=0.25, fill=NA),
        axis.ticks.length=unit(-0.3, "lines"),
        axis.text.x=element_text(margin=unit(c(1.5, 0.5, 0.5, 0.5), "lines")),
        axis.text.y=element_text(margin=unit(c(0.5, 1.5, 0.5, 0.5), "lines")),
        # axis.text.x=element_blank(),
        # axis.title.x=element_blank(),
        # axis.text.y=element_blank(),
        # axis.title.y=element_blank(),
        strip.text=element_blank(),
        plot.margin=unit(c(0.5, 0.5, 1, 1), "lines"))

precipf.mean <- ggplot(data=met.region[met.region$Time == "Modern",]) +
  facet_grid(Time ~ .) +
  geom_raster(aes(x=lon, y=lat, fill=precipf)) +
  geom_path(data=paleon.states, aes(x=long, y=lat, group=group), size=0.2) +
  geom_point(data=hips, aes(x=lon, y=lat), color="gray50", size=3) +
  scale_x_continuous(limits=range(met.region$lon), expand=c(0,0), breaks=seq(-95, -60, by=10), name="Longitude") +
  scale_y_continuous(limits=range(met.region$lat), expand=c(0,0), breaks=seq( 37.5,  50, by= 4), name="Latitude") +
  scale_fill_gradient(low="gray90", high="blue2", name="Precip (mm)") +
  coord_fixed(ratio=1) +
  theme(legend.position=c(0.85, 0.3),
        legend.background=element_blank(),
        legend.key.height=unit(0.3, "lines"),
        legend.text=element_text(size=unit(6, "points")),
        legend.title=element_text(size=unit(8, "points")),
        panel.background=element_blank(),
        panel.border=element_rect(size=0.25, fill=NA),
        axis.ticks.length=unit(-0.3, "lines"),
        axis.text.x=element_blank(),
        axis.title.x=element_blank(),
        axis.text.y=element_blank(),
        axis.title.y=element_blank(),
        # strip.text=element_blank(),
        plot.margin=unit(c(0.5, 1, 0.5, 0.5), "lines"))


precip.diff1 <- ggplot(data=met.region[met.region$Time == "diff.Settlement",]) +
  facet_grid(Time ~ .) +
  geom_raster(aes(x=lon, y=lat, fill=precipf)) +
  geom_path(data=paleon.states, aes(x=long, y=lat, group=group), size=0.2) +
  geom_point(data=hips, aes(x=lon, y=lat), color="gray50", size=3) +
  scale_x_continuous(limits=range(met.region$lon), expand=c(0,0), breaks=seq(-95, -60, by=10), name="Longitude") +
  scale_y_continuous(limits=range(met.region$lat), expand=c(0,0), breaks=seq( 37.5,  50, by= 4), name="Latitude") +
  scale_fill_gradient2(low="red2", high="blue2", name="Precip (mm)", midpoint=0, limits=range(met.region[substr(met.region$Time,1,4)=="diff", "precipf"])) +
  coord_fixed(ratio=1) +
  theme(legend.position=c(0.85, 0.3),
        legend.background=element_blank(),
        legend.key.height=unit(0.3, "lines"),
        legend.text=element_text(size=unit(6, "points")),
        legend.title=element_text(size=unit(8, "points")),
        panel.background=element_blank(),
        panel.border=element_rect(size=0.25, fill=NA),
        axis.ticks.length=unit(-0.3, "lines"),
        axis.text.x=element_blank(),
        axis.title.x=element_blank(),
        axis.text.y=element_blank(),
        axis.title.y=element_blank(),
        # strip.text=element_blank(),
        plot.margin=unit(c(0.5, 1, 0.5, 0.5), "lines"))

precip.diff2 <- ggplot(data=met.region[met.region$Time == "diff.Spinup",]) +
  facet_grid(Time ~ .) +
  geom_raster(aes(x=lon, y=lat, fill=precipf)) +
  geom_path(data=paleon.states, aes(x=long, y=lat, group=group), size=0.2) +
  geom_point(data=hips, aes(x=lon, y=lat), color="gray50", size=3) +
  scale_x_continuous(limits=range(met.region$lon), expand=c(0,0), breaks=seq(-95, -60, by=10), name="Longitude") +
  scale_y_continuous(limits=range(met.region$lat), expand=c(0,0), breaks=seq( 37.5,  50, by= 4), name="Latitude") +
  scale_fill_gradient2(low="red2", high="blue2", name="Precip (mm)", midpoint=0, limits=range(met.region[substr(met.region$Time,1,4)=="diff", "precipf"])) +
  coord_fixed(ratio=1) +
  theme(legend.position=c(0.85, 0.3),
        legend.background=element_blank(),
        legend.key.height=unit(0.3, "lines"),
        legend.text=element_text(size=unit(6, "points")),
        legend.title=element_text(size=unit(8, "points")),
        panel.background=element_blank(),
        panel.border=element_rect(size=0.25, fill=NA),
        axis.ticks.length=unit(-0.3, "lines"),
        axis.text.x=element_text(margin=unit(c(1.5, 0.5, 0.5, 0.5), "lines")),
        axis.text.y=element_blank(),
        axis.title.y=element_blank(),
        # strip.text=element_blank(),
        plot.margin=unit(c(0.5, 1, 1, 0.5), "lines"))


png("Figures/Met/Fig1_MetChange.png", height=5, width=7.5, unit="in", res=220)
# pdf(file.path(fig.out, "MetChange.pdf"), height=5, width=8)
grid.newpage()
pushViewport(viewport(layout=grid.layout(nrow=3,ncol=2, widths=c(1.1, 1.0), heights=c(1.0, 1.0, 1.5))))
print(tair.mean   , vp = viewport(layout.pos.row = 1, layout.pos.col=1))
print(tair.diff1  , vp = viewport(layout.pos.row = 2, layout.pos.col=1))
print(tair.diff2  , vp = viewport(layout.pos.row = 3, layout.pos.col=1))
print(precipf.mean, vp = viewport(layout.pos.row = 1, layout.pos.col=2))
print(precip.diff1, vp = viewport(layout.pos.row = 2, layout.pos.col=2))
print(precip.diff2, vp = viewport(layout.pos.row = 3, layout.pos.col=2))
dev.off()
}
# -------------------------------------------

# -------------------------------------------
# Supplemental Figure 1: Regional Met Change -- ALl drivers
# -------------------------------------------
{
  
  fig.out <- "Figures/Met"
  met.region <- read.csv("Data/Met/Met_Region_Annual.csv")
  paleon.states <- map_data("state")
  
  met.stack <- stack(met.region[,c("tair", "precipf", "swdown", "lwdown", "psurf", "qair", "wind")])
  names(met.stack) <- c("value", "var")
  met.stack[, c("lat", "lon", "Time")] <- met.region[,c("lat", "lon", "Time")]
  
  hips <- data.frame(Site = c( "PHA",  "PHO",  "PUN",  "PBL",  "PDL",  "PMB"),
                     lat  = c( 42.54,  45.25,  46.22,  46.28,  47.17,  43.61),
                     lon  = c(-72.18, -68.73, -89.53, -94.58, -95.17, -82.83),
                     lat2 = c( 42.75,  45.25,  46.25,  46.25,  47.25,  43.75),
                     lon2 = c(-72.25, -68.75, -89.75, -94.75, -95.25, -82.75))
  
  # Themes and margins specifcations 
  {
  theme.base1 <- theme(legend.position=c(0.85, 0.3),
                      legend.background=element_blank(),
                      legend.key.height=unit(0.3, "lines"),
                      legend.text=element_text(size=unit(6, "points")),
                      legend.title=element_text(size=unit(8, "points")),
                      panel.background=element_blank(),
                      panel.border=element_rect(size=0.25, fill=NA),
                      axis.ticks.length=unit(-0.3, "lines"))
  theme.base2 <- theme(legend.position=c(0.95, 0.3),
                      legend.background=element_blank(),
                      legend.key.height=unit(0.3, "lines"),
                      legend.text=element_text(size=unit(6, "points")),
                      legend.title=element_text(size=unit(8, "points")),
                      panel.background=element_blank(),
                      panel.border=element_rect(size=0.25, fill=NA),
                      axis.ticks.length=unit(-0.3, "lines"))
  margin.topleft <- theme(axis.text.x=element_blank(),
                          axis.title.x=element_blank(),
                          axis.text.y=element_text(margin=unit(c(0.5, 1.5, 0.5, 0.5), "lines")),
                          strip.text.y=element_blank(),
                          plot.margin=unit(c(1, 0.5, 0.5, 1), "lines"))
  margin.topmid <- theme(axis.text.x=element_blank(),
                         axis.title.x=element_blank(),
                         axis.text.y=element_blank(),
                         axis.title.y=element_blank(),
                         strip.text.y=element_blank(),
                         plot.margin=unit(c(1, 0.5, 0.5, 0.5), "lines"))
  margin.topright <- theme(axis.text.x=element_blank(),
                           axis.title.x=element_blank(),
                           axis.text.y=element_blank(),
                           axis.title.y=element_blank(),
                           # strip.text=element_blank(),
                           plot.margin=unit(c(1, 1, 0.5, 0.5), "lines"))
  margin.midleft  <- theme(axis.text.x=element_blank(),
                           axis.title.x=element_blank(),
                           axis.text.y=element_text(margin=unit(c(0.5, 1.5, 0.5, 0.5), "lines")),
                           strip.text=element_blank(),
                           plot.margin=unit(c(0.5, 0.5, 0.5, 1), "lines"))
  margin.midmid   <- theme(axis.text.x=element_blank(),
                           axis.title.x=element_blank(),
                           axis.text.y=element_blank(),
                           axis.title.y=element_blank(),
                           strip.text=element_blank(),
                           plot.margin=unit(c(0.5, 0.5, 0.5, 0.5), "lines"))
  margin.midright <- theme(axis.text.x=element_blank(),
                           axis.title.x=element_blank(),
                           axis.text.y=element_blank(),
                           axis.title.y=element_blank(),
                           strip.text.x=element_blank(),
                           plot.margin=unit(c(0.5, 1, 0.5, 0.5), "lines"))
  margin.lowleft <- theme(axis.text.x=element_text(margin=unit(c(1.5, 0.5, 0.5, 0.5), "lines")),
                          # axis.title.x=element_blank(),
                          # axis.title.y=element_blank(),
                          axis.text.y=element_text(margin=unit(c(0.5, 1.5, 0.5, 0.5), "lines")),
                          strip.text=element_blank(),
                          plot.margin=unit(c(0.5, 0.5, 1, 1), "lines"))
  margin.lowmid  <- theme(axis.text.x=element_text(margin=unit(c(1.5, 0.5, 0.5, 0.5), "lines")),
                          # axis.title.x=element_blank(),
                          axis.title.y=element_blank(),
                          axis.text.y=element_blank(),
                          strip.text=element_blank(),
                          plot.margin=unit(c(0.5, 0.5, 1, 0.5), "lines"))
  margin.lowright  <- theme(axis.text.x=element_text(margin=unit(c(1.5, 0.5, 0.5, 0.5), "lines")),
                            # axis.title.x=element_blank(),
                            axis.title.y=element_blank(),
                            axis.text.y=element_blank(),
                            strip.text.x=element_blank(),
                            plot.margin=unit(c(0.5, 1, 1, 0.5), "lines"))
  }

  tair.mean <- { ggplot(data=met.stack[met.stack$Time=="Modern" & met.stack$var=="tair",]) +
    facet_grid(var ~ Time) +
    geom_raster(aes(x=lon, y=lat, fill=value)) +
    geom_path(data=paleon.states, aes(x=long, y=lat, group=group), size=0.2) +
    geom_point(data=hips, aes(x=lon, y=lat), color="gray50", size=3) +
    scale_x_continuous(limits=range(met.region$lon), expand=c(0,0), breaks=seq(-95, -60, by=10), name="Longitude") +
    scale_y_continuous(limits=range(met.region$lat), expand=c(0,0), breaks=seq( 37.5,  50, by= 4), name="Latitude") +
    scale_fill_gradient(low="blue4", high="red2") +
    coord_fixed(ratio=1) +
    theme.base2 + margin.topleft}

  tair.diff1 <- {ggplot(data=met.stack[substr(met.stack$Time,1,4) == "diff" & met.stack$var=="tair",]) +
    facet_grid(var ~ Time) +
    geom_raster(aes(x=lon, y=lat, fill=value)) +
    geom_path(data=paleon.states, aes(x=long, y=lat, group=group), size=0.2) +
    geom_point(data=hips, aes(x=lon, y=lat), color="gray50", size=3) +
    scale_x_continuous(limits=range(met.region$lon), expand=c(0,0), breaks=seq(-95, -60, by=10), name="Longitude") +
    scale_y_continuous(limits=range(met.region$lat), expand=c(0,0), breaks=seq( 37.5,  50, by= 4), name="Latitude") +
    scale_fill_gradient2(low="blue4", high="red2", midpoint=0) +
    coord_fixed(ratio=1) +
    theme.base2 + margin.topright}
  
  precipf.mean <- { ggplot(data=met.stack[met.stack$Time=="Modern" & met.stack$var=="precipf",]) +
      facet_grid(var ~ Time) +
      geom_raster(aes(x=lon, y=lat, fill=value)) +
      geom_path(data=paleon.states, aes(x=long, y=lat, group=group), size=0.2) +
      geom_point(data=hips, aes(x=lon, y=lat), color="gray50", size=3) +
      scale_x_continuous(limits=range(met.region$lon), expand=c(0,0), breaks=seq(-95, -60, by=10), name="Longitude") +
      scale_y_continuous(limits=range(met.region$lat), expand=c(0,0), breaks=seq( 37.5,  50, by= 4), name="Latitude") +
      scale_fill_gradient(low="gray80", high="blue2") +
      coord_fixed(ratio=1) +
      theme.base1 + margin.midleft}
  
  precipf.diff1 <- {ggplot(data=met.stack[substr(met.stack$Time,1,4) == "diff" & met.stack$var=="precipf",]) +
      facet_grid(var ~ Time) +
      geom_raster(aes(x=lon, y=lat, fill=value)) +
      geom_path(data=paleon.states, aes(x=long, y=lat, group=group), size=0.2) +
      geom_point(data=hips, aes(x=lon, y=lat), color="gray50", size=3) +
      scale_x_continuous(limits=range(met.region$lon), expand=c(0,0), breaks=seq(-95, -60, by=10), name="Longitude") +
      scale_y_continuous(limits=range(met.region$lat), expand=c(0,0), breaks=seq( 37.5,  50, by= 4), name="Latitude") +
      scale_fill_gradient2(low="red3", high="blue2", midpoint=0) +
      coord_fixed(ratio=1) +
      theme.base2 + margin.midright}
  
  swdown.mean <- { ggplot(data=met.stack[met.stack$Time=="Modern" & met.stack$var=="swdown",]) +
      facet_grid(var ~ Time) +
      geom_raster(aes(x=lon, y=lat, fill=value)) +
      geom_path(data=paleon.states, aes(x=long, y=lat, group=group), size=0.2) +
      geom_point(data=hips, aes(x=lon, y=lat), color="gray50", size=3) +
      scale_x_continuous(limits=range(met.region$lon), expand=c(0,0), breaks=seq(-95, -60, by=10), name="Longitude") +
      scale_y_continuous(limits=range(met.region$lat), expand=c(0,0), breaks=seq( 37.5,  50, by= 4), name="Latitude") +
      scale_fill_gradient(low="lightsteelblue3", high="darkgoldenrod2") +
      coord_fixed(ratio=1) +
      theme.base1 + margin.midleft}
  
  swdown.diff1 <- {ggplot(data=met.stack[substr(met.stack$Time,1,4) == "diff" & met.stack$var=="swdown",]) +
      facet_grid(var ~ Time) +
      geom_raster(aes(x=lon, y=lat, fill=value)) +
      geom_path(data=paleon.states, aes(x=long, y=lat, group=group), size=0.2) +
      geom_point(data=hips, aes(x=lon, y=lat), color="gray50", size=3) +
      scale_x_continuous(limits=range(met.region$lon), expand=c(0,0), breaks=seq(-95, -60, by=10), name="Longitude") +
      scale_y_continuous(limits=range(met.region$lat), expand=c(0,0), breaks=seq( 37.5,  50, by= 4), name="Latitude") +
      scale_fill_gradient2(low="lightsteelblue4", mid="darkgoldenrod2", midpoint=0) +
      coord_fixed(ratio=1) +
      theme.base2 + margin.midright}
  
  lwdown.mean <- { ggplot(data=met.stack[met.stack$Time=="Modern" & met.stack$var=="lwdown",]) +
      facet_grid(var ~ Time) +
      geom_raster(aes(x=lon, y=lat, fill=value)) +
      geom_path(data=paleon.states, aes(x=long, y=lat, group=group), size=0.2) +
      geom_point(data=hips, aes(x=lon, y=lat), color="gray50", size=3) +
      scale_x_continuous(limits=range(met.region$lon), expand=c(0,0), breaks=seq(-95, -60, by=10), name="Longitude") +
      scale_y_continuous(limits=range(met.region$lat), expand=c(0,0), breaks=seq( 37.5,  50, by= 4), name="Latitude") +
      scale_fill_gradient(low="lightsteelblue3", high="orange3") +
      coord_fixed(ratio=1) +
      theme.base1 + margin.midleft}
  
  lwdown.diff1 <- {ggplot(data=met.stack[substr(met.stack$Time,1,4) == "diff" & met.stack$var=="lwdown",]) +
      facet_grid(var ~ Time) +
      geom_raster(aes(x=lon, y=lat, fill=value)) +
      geom_path(data=paleon.states, aes(x=long, y=lat, group=group), size=0.2) +
      geom_point(data=hips, aes(x=lon, y=lat), color="gray50", size=3) +
      scale_x_continuous(limits=range(met.region$lon), expand=c(0,0), breaks=seq(-95, -60, by=10), name="Longitude") +
      scale_y_continuous(limits=range(met.region$lat), expand=c(0,0), breaks=seq( 37.5,  50, by= 4), name="Latitude") +
      scale_fill_gradient2(low="lightsteelblue4", mid="gray90", high="orange3", midpoint=0) +
      coord_fixed(ratio=1) +
      theme.base2 + margin.midright }

  qair.mean <- { ggplot(data=met.stack[met.stack$Time=="Modern" & met.stack$var=="qair",]) +
      facet_grid(var ~ Time) +
      geom_raster(aes(x=lon, y=lat, fill=value)) +
      geom_path(data=paleon.states, aes(x=long, y=lat, group=group), size=0.2) +
      geom_point(data=hips, aes(x=lon, y=lat), color="gray50", size=3) +
      scale_x_continuous(limits=range(met.region$lon), expand=c(0,0), breaks=seq(-95, -60, by=10), name="Longitude") +
      scale_y_continuous(limits=range(met.region$lat), expand=c(0,0), breaks=seq( 37.5,  50, by= 4), name="Latitude") +
      scale_fill_gradient(low="gray80", high="dodgerblue2") +
      coord_fixed(ratio=1) +
      theme.base1 + margin.midleft}
  
  qair.diff1 <- {ggplot(data=met.stack[substr(met.stack$Time,1,4) == "diff" & met.stack$var=="qair",]) +
      facet_grid(var ~ Time) +
      geom_raster(aes(x=lon, y=lat, fill=value)) +
      geom_path(data=paleon.states, aes(x=long, y=lat, group=group), size=0.2) +
      geom_point(data=hips, aes(x=lon, y=lat), color="gray50", size=3) +
      scale_x_continuous(limits=range(met.region$lon), expand=c(0,0), breaks=seq(-95, -60, by=10), name="Longitude") +
      scale_y_continuous(limits=range(met.region$lat), expand=c(0,0), breaks=seq( 37.5,  50, by= 4), name="Latitude") +
      scale_fill_gradient2(low="red2", mid="gray80",  high="dodgerblue2", midpoint=0) +
      coord_fixed(ratio=1) +
      theme.base2 + margin.midright} 
  
  psurf.mean <- { ggplot(data=met.stack[met.stack$Time=="Modern" & met.stack$var=="psurf",]) +
      facet_grid(var ~ Time) +
      geom_raster(aes(x=lon, y=lat, fill=value)) +
      geom_path(data=paleon.states, aes(x=long, y=lat, group=group), size=0.2) +
      geom_point(data=hips, aes(x=lon, y=lat), color="gray50", size=3) +
      scale_x_continuous(limits=range(met.region$lon), expand=c(0,0), breaks=seq(-95, -60, by=10), name="Longitude") +
      scale_y_continuous(limits=range(met.region$lat), expand=c(0,0), breaks=seq( 37.5,  50, by= 4), name="Latitude") +
      scale_fill_gradient(low="gray80", high="goldenrod3") +
      coord_fixed(ratio=1) +
      theme.base1 + margin.midleft}
  
  psurf.diff1 <- {ggplot(data=met.stack[substr(met.stack$Time,1,4) == "diff" & met.stack$var=="psurf",]) +
      facet_grid(var ~ Time) +
      geom_raster(aes(x=lon, y=lat, fill=value)) +
      geom_path(data=paleon.states, aes(x=long, y=lat, group=group), size=0.2) +
      geom_point(data=hips, aes(x=lon, y=lat), color="gray50", size=3) +
      scale_x_continuous(limits=range(met.region$lon), expand=c(0,0), breaks=seq(-95, -60, by=10), name="Longitude") +
      scale_y_continuous(limits=range(met.region$lat), expand=c(0,0), breaks=seq( 37.5,  50, by= 4), name="Latitude") +
      scale_fill_gradient2(low="darkslategray4", mid="gray80",  high="goldenrod3", midpoint=0) +
      coord_fixed(ratio=1) +
      theme.base2 + margin.midright} 
  
  wind.mean <- { ggplot(data=met.stack[met.stack$Time=="Modern" & met.stack$var=="wind",]) +
      facet_grid(var ~ Time) +
      geom_raster(aes(x=lon, y=lat, fill=value)) +
      geom_path(data=paleon.states, aes(x=long, y=lat, group=group), size=0.2) +
      geom_point(data=hips, aes(x=lon, y=lat), color="gray50", size=3) +
      scale_x_continuous(limits=range(met.region$lon), expand=c(0,0), breaks=seq(-95, -60, by=10), name="Longitude") +
      scale_y_continuous(limits=range(met.region$lat), expand=c(0,0), breaks=seq( 37.5,  50, by= 4), name="Latitude") +
      scale_fill_gradient(low="gray80", high="darkslategray4") +
      coord_fixed(ratio=1) +
      theme.base1 + margin.lowleft}
  
  wind.diff1 <- {ggplot(data=met.stack[substr(met.stack$Time,1,4) == "diff" & met.stack$var=="wind",]) +
      facet_grid(var ~ Time) +
      geom_raster(aes(x=lon, y=lat, fill=value)) +
      geom_path(data=paleon.states, aes(x=long, y=lat, group=group), size=0.2) +
      geom_point(data=hips, aes(x=lon, y=lat), color="gray50", size=3) +
      scale_x_continuous(limits=range(met.region$lon), expand=c(0,0), breaks=seq(-95, -60, by=10), name="Longitude") +
      scale_y_continuous(limits=range(met.region$lat), expand=c(0,0), breaks=seq( 37.5,  50, by= 4), name="Latitude") +
      scale_fill_gradient2(low="red", mid="gray80",  high="darkslategray4", midpoint=0) +
      coord_fixed(ratio=1) +
      theme.base2 + margin.lowright} 
  
  
  png("Figures/Met/SuppFig1_MetChange_All.png", height=9, width=8, unit="in", res=220)
  # pdf(file.path(fig.out, "MetChange.pdf"), height=5, width=8)
  grid.newpage()
  pushViewport(viewport(layout=grid.layout(nrow=7,ncol=2, widths=c(1, 1.65), heights=c(1.4, 1.0, 1.0, 1.0, 1.0, 1.0, 1.75))))
  print(tair.mean   , vp = viewport(layout.pos.row = 1, layout.pos.col=1))
  print(tair.diff1  , vp = viewport(layout.pos.row = 1, layout.pos.col=2))
  print(precipf.mean, vp = viewport(layout.pos.row = 2, layout.pos.col=1))
  print(precipf.diff1, vp = viewport(layout.pos.row = 2, layout.pos.col=2))
  print(swdown.mean , vp = viewport(layout.pos.row = 3, layout.pos.col=1))
  print(swdown.diff1, vp = viewport(layout.pos.row = 3, layout.pos.col=2))
  print(lwdown.mean , vp = viewport(layout.pos.row = 4, layout.pos.col=1))
  print(lwdown.diff1, vp = viewport(layout.pos.row = 4, layout.pos.col=2))
  print(qair.mean   , vp = viewport(layout.pos.row = 5, layout.pos.col=1))
  print(qair.diff1  , vp = viewport(layout.pos.row = 5, layout.pos.col=2))
  print(psurf.mean  , vp = viewport(layout.pos.row = 6, layout.pos.col=1))
  print(psurf.diff1 , vp = viewport(layout.pos.row = 6, layout.pos.col=2))
  print(wind.mean   , vp = viewport(layout.pos.row = 7, layout.pos.col=1))
  print(wind.diff1  , vp = viewport(layout.pos.row = 7, layout.pos.col=2))
  dev.off()
}
# -------------------------------------------
