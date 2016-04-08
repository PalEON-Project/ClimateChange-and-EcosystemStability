# -------------------------------------------
# Assessing the scales and magnitude of change in ecosystem state

# 1. Stastical detection of significant change
#    - Use TP regression splines
# 2. Comparison with paleoclimate reconstructions
# -------------------------------------------
rm(list=ls())


# -------------------------------------------
# Load libraries; set file paths
# -------------------------------------------
library(car)
library(ggplot2)
library(mgcv)
setwd("~/Desktop/Research/PalEON_CR/PalEON_MIP_Site/Analyses/Change-and-Stability")

inputs    <- "~/Desktop/Research/PalEON_CR/PalEON_MIP_Site/phase1a_output_variables"
dat.out <- "Data/EcosysChange"
fig.out <- "Figures/EcosysChange"

if(!dir.exists(dat.out)) dir.create(dat.out)
if(!dir.exists(fig.out)) dir.create(fig.out)
# -------------------------------------------


# -------------------------------------------
# Reading in raw ecosystem model output
# -------------------------------------------
ecosys <- read.csv(file.path(inputs, "PalEON_MIP_Yearly.csv"))
ecosys$Model.Order <- recode(ecosys$Model, "'clm.bgc'='01'; 'clm.cn'='02'; 'ed2'='03'; 'ed2.lu'='04';  'jules.stat'='05'; 'jules.triffid'='06'; 'linkages'='07'; 'lpj.guess'='08'; 'lpj.wsl'='09'; 'sibcasa'='10'")
levels(ecosys$Model.Order) <- c("CLM-BGC", "CLM-CN", "ED2", "ED2-LU", "JULES-STATIC", "JULES-TRIFFID", "LINKAGES", "LPJ-GUESS", "LPJ-WSL", "SiBCASA")
summary(ecosys)

# Colors used for graphing
model.colors <- read.csv("../../Model.Colors.csv")
model.colors $Model.Order <- recode(model.colors$Model, "'CLM4.5-BGC'='01'; 'CLM4.5-CN'='02'; 'ED2'='03'; 'ED2-LU'='04';  'JULES-STATIC'='05'; 'JULES-TRIFFID'='06'; 'LINKAGES'='07'; 'LPJ-GUESS'='08'; 'LPJ-WSL'='09'; 'SiBCASA'='10'")
levels(model.colors$Model.Order)[1:10] <- c("CLM-BGC", "CLM-CN", "ED2", "ED2-LU", "JULES-STATIC", "JULES-TRIFFID", "LINKAGES", "LPJ-GUESS", "LPJ-WSL", "SiBCASA")
model.colors <- model.colors[order(model.colors$Model.Order),]
model.colors
# -------------------------------------------

# -------------------------------------------
# Detecting change in the models through time
# -------------------------------------------
# Exploratory Graphs
{
col.model <- paste(model.colors[model.colors$Model.Order %in% unique(ecosys$Model.Order),"color"])

png(file.path(fig.out, "NEE_Sites_Smooth.png"), height=11, width=8.5, units="in", res=180)
ggplot(data=ecosys) +
  facet_wrap(~Site) +
  geom_line(aes(x=Year, y=NEE, color=Model), size=0.25, alpha=0.3) +
  geom_smooth(aes(x=Year, y=NEE, color=Model), se=F, size=1.5) +
  scale_y_continuous(expand=c(0,0)) +
  scale_x_continuous(expand=c(0,0)) +
  scale_color_manual(values=col.model) +
  theme_bw()
dev.off()

png(file.path(fig.out, "NPP_Sites_Smooth.png"), height=11, width=8.5, units="in", res=180)
ggplot(data=ecosys) +
  facet_wrap(~Site) +
  geom_line(aes(x=Year, y=NPP, color=Model), size=0.25, alpha=0.3) +
  geom_smooth(aes(x=Year, y=NPP, color=Model), se=F, size=1.5) +
  scale_y_continuous(expand=c(0,0)) +
  scale_x_continuous(expand=c(0,0)) +
  scale_color_manual(values=col.model) +
  theme_bw()
dev.off()

png(file.path(fig.out, "AGB_Sites_Smooth.png"), height=11, width=8.5, units="in", res=180)
ggplot(data=ecosys) +
  facet_wrap(~Site) +
  geom_line(aes(x=Year, y=AGB, color=Model), size=0.25, alpha=0.3) +
  geom_smooth(aes(x=Year, y=AGB, color=Model), se=F, size=1.5) +
  scale_y_continuous(expand=c(0,0)) +
  scale_x_continuous(expand=c(0,0)) +
  scale_color_manual(values=col.model) +
  theme_bw()
dev.off()

png(file.path(fig.out, "Evergreen_Sites_Smooth.png"), height=11, width=8.5, units="in", res=180)
ggplot(data=ecosys) +
  facet_wrap(~Site) +
  geom_line(aes(x=Year, y=Evergreen, color=Model), size=0.25, alpha=0.3) +
  geom_smooth(aes(x=Year, y=Evergreen, color=Model), se=F, size=1.5) +
  scale_y_continuous(expand=c(0,0)) +
  scale_x_continuous(expand=c(0,0)) +
  scale_color_manual(values=col.model) +
  theme_bw()
dev.off()

png(file.path(fig.out, "SoilCarb_Sites_Smooth.png"), height=11, width=8.5, units="in", res=180)
ggplot(data=ecosys) +
  facet_wrap(~Site) +
  geom_line(aes(x=Year, y=SoilCarb, color=Model), size=0.25, alpha=0.3) +
  geom_smooth(aes(x=Year, y=SoilCarb, color=Model), se=F, size=1.5) +
  scale_y_continuous(expand=c(0,0)) +
  scale_x_continuous(expand=c(0,0)) +
  scale_color_manual(values=col.model) +
  theme_bw()
dev.off()
}

# Fit some simple gams with temporal thin-plate regression splines on year by model across sites
gam.nee <- gam(NEE ~ s(Year,by=Model) + Model*Site, data=ecosys)
gam.npp <- gam(NPP ~ s(Year,by=Model) + Model*Site, data=ecosys)
gam.agb <- gam(AGB ~ s(Year,by=Model) + Model*Site, data=ecosys)
gam.evg <- gam(Evergreen ~ s(Year,by=Model) + Model*Site, data=ecosys)
gam.sc  <- gam(SoilCarb ~ s(Year,by=Model) + Model*Site, data=ecosys)

# Doing a fancy prediction to get full 95% CIs
source("~/Desktop/Research/R_Functions/Calculate_GAMM_Posteriors.R")
nee.predict <- post.distns(model.gam=gam.nee, model.name="NEE"      , n=100, newdata=ecosys, vars="Year", terms=T)
npp.predict <- post.distns(model.gam=gam.npp, model.name="NPP"      , n=100, newdata=ecosys, vars="Year", terms=T)
agb.predict <- post.distns(model.gam=gam.agb, model.name="AGB"      , n=100, newdata=ecosys, vars="Year", terms=T)
evg.predict <- post.distns(model.gam=gam.evg, model.name="Evergreen", n=100, newdata=ecosys[!ecosys$Model=="sibcasa",], vars="Year", terms=T)
sc.predict  <- post.distns(model.gam=gam.sc , model.name="SoilCarb" , n=100, newdata=ecosys, vars="Year", terms=T)

# A little bit of re-formatting before we rbind the predictions together
nee.predict$Effect <- as.factor("NEE")
npp.predict$Effect <- as.factor("NPP")
agb.predict$Effect <- as.factor("AGB")
evg.predict$Effect <- as.factor("Evergreen")
sc.predict $Effect <- as.factor("SoilCarb")

nee.predict$Year <- ecosys$Year
npp.predict$Year <- ecosys$Year
agb.predict$Year <- ecosys$Year
evg.predict$Year <- ecosys[!ecosys$Model=="sibcasa","Year"]
sc.predict $Year <- ecosys$Year

nee.predict$Model <- ecosys$Model
npp.predict$Model <- ecosys$Model
agb.predict$Model <- ecosys$Model
evg.predict$Model <- ecosys[!ecosys$Model=="sibcasa","Model"]
sc.predict $Model <- ecosys$Model

nee.predict$Model.Order <- ecosys$Model.Order
npp.predict$Model.Order <- ecosys$Model.Order
agb.predict$Model.Order <- ecosys$Model.Order
evg.predict$Model.Order <- ecosys[!ecosys$Model=="sibcasa","Model.Order"]
sc.predict $Model.Order <- ecosys$Model.Order

ecosys.change <- rbind(nee.predict, npp.predict, agb.predict, evg.predict, sc.predict)
ecosys.change <- aggregate(ecosys.change[,c("mean", "lwr", "upr")], 
                           by=ecosys.change[,c("Effect","Model", "Model.Order", "Year")],
                           FUN=mean)
summary(ecosys.change)

col.model <- paste(model.colors[model.colors$Model.Order %in% unique(ecosys.change$Model.Order),"color"])

png(file.path(fig.out, "Ecosystems_TemporalTrends.png"), height=8, width=9, units="in", res=180)
ggplot(data=ecosys.change) +
  facet_grid(Effect~., scales="free_y") +
  geom_ribbon(aes(x=Year, ymin=lwr, ymax=upr, fill=Model), alpha=0.5) +
  geom_line(aes(x=Year, y=mean, color=Model), size=1) +
  geom_vline(xintercept=1850, linetype="dashed") +
  scale_x_continuous(expand=c(0,0), name="Year") +
  scale_y_continuous(expand=c(0,0), name="Temporal Trend") +
  scale_color_manual(values=col.model) +
  scale_fill_manual(values=col.model) + 
  theme_bw()
dev.off()
# -------------------------------------------
