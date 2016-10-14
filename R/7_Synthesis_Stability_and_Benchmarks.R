# -------------------------------------------
# See if paleo stability predicts modern performance
# Author: Christy Rollinson, crollinson@gmail.com
#
# ------------------
# Objectives:
# ------------------
# See if paleo stability predicts modern performance
# ------------------
#
#
# ------------------
# Workflow
# ------------------
# ------------------
# -------------------------------------------
rm(list=ls())

# -------------------------------------------
# 1. Set up file structure, libraries, etc.
# -------------------------------------------
library(car)
library(ggplot2); library(scales)
library(mgcv)
library(nlme)
setwd("~/Dropbox/PalEON_CR/PalEON_MIP_Site/Analyses/Change-and-Stability") # Path to this project github repository: https://github.com/PalEON-Project/MIP-Change-and-Stability.git
# setwd("~/Desktop/Research/PalEON_CR/PalEON_MIP_Site/Analyses/Change-and-Stability") # Path to this project github repository: https://github.com/PalEON-Project/MIP-Change-and-Stability.git
# path.gamm.func <- "~/Desktop/R_Functions/"  # Path to github repository of my GAMM helper functions: https://github.com/crollinson/R_Functions.git
inputs    <- "Data/" # Path to my cleaned model output

mip.utils <- "~/Dropbox/PalEON_CR/MIP_Utils/" # Path to PalEON MIP Utility repository: https://github.com/PalEON-Project/MIP_Utils.git

out.dir <- "Data/StabilityBenchmarkSynthesis" # Path to where the analysis output should go
fig.dir <- "Figures/StabilityBenchmarkSynthesis" # Path to where figures should go

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
# 2. load the datasets
# ------------------
{
  # Load the stability Stats
  ecosys.stats.site <- read.csv("Data/StabilitySynthesis/Stability_Models_Stats_All.csv")
  summary(ecosys.stats.site)
  
  # Load the benchmark stats
  models.bench <- read.csv("Data/BenchmarkingSynthesis/Benchmarks_Summary_All.csv")
  names(models.bench)[which(names(models.bench)=="Var")] <- "var"
  summary(models.bench)
  
  stability.bench <- merge(models.bench, ecosys.stats.site, all.x=F)
  summary(stability.bench)
}
# ------------------

# ------------------
# Modeling benchmark performance as a function of paleo stabiltiy
# ------------------
unique(stability.bench$var)

mod.all <- lme(abs(bias) ~ (n.inst.paleo + rate.mean.paleo + rate.max.paleo)*(var-1) - n.inst.paleo - rate.mean.paleo - rate.max.paleo, random=list(Dataset=~1, Site=~1, Model=~1), data=stability.bench[,], na.action=na.omit)
anova(mod.all)
summary(mod.all)

mod.agb.nbcdmax <- lme(abs(bias) ~ n.inst.paleo + rate.mean.paleo + rate.max.paleo, random=list(Site=~1, Model=~1), data=stability.bench[stability.bench$Dataset=="NBCD.max",], na.action=na.omit)
summary(mod.agb.nbcdmax)

mod.agb.setveg <- lme(abs(bias) ~ n.inst.paleo + rate.mean.paleo + rate.max.paleo, random=list(Site=~1, Model=~1), data=stability.bench[stability.bench$Dataset=="SetVeg" & stability.bench$var=="AGB",], na.action=na.omit)
summary(mod.agb.setveg)

# agb.nbcd.setveg <- lme()

mod.n <- lme(abs(bias) ~ n.inst.paleo*var - n.inst.paleo - var -1 , random=list(Dataset=~1, Site=~1, Model=~1), data=stability.bench[,], na.action=na.omit)
summary(mod.n)

mod.rate.mean <- lme(abs(bias) ~ rate.mean.paleo*var - rate.mean.paleo - var -1 , random=list(Dataset=~1, Site=~1, Model=~1), data=stability.bench[,], na.action=na.omit)
summary(mod.rate.mean)

mod.rate.max <- lme(abs(bias) ~ rate.max.paleo*var - rate.max.paleo - var -1 , random=list(Dataset=~1, Site=~1, Model=~1), data=stability.bench[,], na.action=na.omit)
summary(mod.rate.max)


mod.gpp <- lme(abs(bias) ~ n.inst.paleo + rate.mean.paleo + rate.max.paleo , random=list(Dataset=~1, Site=~1, Model=~1), data=stability.bench[stability.bench$var=="GPP",], na.action=na.omit)
anova(mod.gpp)
summary(mod.gpp) # n, mean rate; max rate semi-sig

mod.nee <- lme(abs(bias) ~ n.inst.paleo + rate.mean.paleo + rate.max.paleo , random=list(Dataset=~1, Site=~1, Model=~1), data=stability.bench[stability.bench$var=="NEE",], na.action=na.omit)
anova(mod.nee)
summary(mod.nee) # mean rate

mod.lai <- lme(abs(bias) ~ n.inst.paleo + rate.mean.paleo + rate.max.paleo -1, random=list(Dataset=~1, Site=~1, Model=~1), data=stability.bench[stability.bench$var=="LAI",], na.action=na.omit)
anova(mod.lai)
summary(mod.lai) # n, mean rate, max rate

mod.agb <- lme(abs(bias) ~ n.inst.paleo + rate.mean.paleo + rate.max.paleo -1, random=list(Dataset=~1, Benchmark=~1, Site=~1, Model=~1), data=stability.bench[stability.bench$var=="AGB",], na.action=na.omit)
anova(mod.agb)
summary(mod.agb) # max rate

mod.agb.nbcdmax <- lme(abs(bias) ~ n.inst.paleo + rate.mean.paleo + rate.max.paleo, random=list(Site=~1, Model=~1), data=stability.bench[stability.bench$Dataset=="NBCD.max",], na.action=na.omit)
summary(mod.agb.nbcdmax)

mod.agb.setveg <- lme(abs(bias) ~ n.inst.paleo + rate.mean.paleo + rate.max.paleo, random=list(Site=~1, Model=~1), data=stability.bench[stability.bench$Dataset=="SetVeg" & stability.bench$var=="AGB",], na.action=na.omit)
summary(mod.agb.setveg)



mod.fcomp <- lme(abs(bias) ~ n.inst.paleo + rate.mean.paleo + rate.max.paleo -1, random=list(Dataset=~1, Benchmark=~1, Site=~1, Model=~1), data=stability.bench[stability.bench$var=="Fcomp",], na.action=na.omit)
anova(mod.fcomp)
summary(mod.fcomp) # max rate

# ------------------