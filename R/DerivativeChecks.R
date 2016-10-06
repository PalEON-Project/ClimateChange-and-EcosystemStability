# -------------------------------------------
# Derivative calculations from earlier don't line up with the change in temperature observed
# in the same spline... what's up with that?
# 
# Comparing the derivative calculation with first-differences from the sims
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
met.dir <- "Data/Met"

inputs    <- "Data/" # Path to my cleaned model output

path.gamm.func <- "~/Desktop/R_Functions/"  # Path to github repository of my GAMM helper functions: https://github.com/crollinson/R_Functions.git
mip.utils <- "~/Dropbox/Research/PalEON_CR/MIP_Utils/" # Path to PalEON MIP Utility repository: https://github.com/PalEON-Project/MIP_Utils.git

out.dir <- "Data/StabilitySynthesis_DerivChecks" # Path to where the analysis output should go
fig.dir <- "Figures/StabilitySynthesis_DerivChecks" # Path to where figures should go

if(!dir.exists(out.dir)) dir.create(out.dir)
if(!dir.exists(fig.dir)) dir.create(fig.dir)
# -------------------------------------------

# -------------------------------------------
# -------------------------------------------
source("R/0_TimeAnalysis.R")
source(file.path(path.gamm.func, "Calculate_GAMM_Derivs.R"))

# Load the met data
met.sites <- read.csv(file.path(met.dir, "Met_Sites_Annual.csv"))
summary(met.sites)

# subset just the PHA temperature data to test & work with for easiness
dat <- met.sites[met.sites$Site=="PHA" & met.sites$Year<1850,]
summary(dat)

# ------------------
# Running the guts of 0_TimeAnalysis.R so we can isolate problems
# Removed fac.fit to ignore site for now
# ------------------
# Setting up the data frame with the generalized column names
Y="tair"
k.freq=25
# Determine the number of knots to use in the model
k.use = round((max(dat$Year) - min(dat$Year))/k.freq,0)

dat.temp <- data.frame(Year=dat$Year)
dat.temp[,"Y"] <- dat[,Y] - mean(dat[dat$Year<869,Y])
# dat.temp[,"fac.fit"] <- dat[,fac.fit]

# Fit the gam
# mod.gam.raw <- mod.gam
mod.gam <- gam(Y ~ s(Year, k=k.use), data=dat.temp)
summary(mod.gam)
par(mfrow=c(1,1))
plot(mod.gam)
abline(h=0, col="red", lty="dashed")

# Fitting the posteriors; KEEPING THE SIMS!
source(file.path(path.gamm.func, "Calculate_GAMM_Posteriors.R"))
mod.predict <- post.distns(model.gam=mod.gam, model.name=Y, n=100, newdata=dat.temp, vars="Year", terms=T, return.sims=T)
mod.predict$ci$Year <- dat$Year
mod.predict$sims$Year <- dat$Year
# mod.predict[,fac.fit] <- dat.temp$fac.fit
head(mod.predict$ci)

plot(mean ~ Year, data=mod.predict$ci, ylim=range(mod.predict$ci[,c("mean", "lwr", "upr")]), type="l")
lines(lwr ~ Year, data=mod.predict$ci, lty="dashed")
lines(upr ~ Year, data=mod.predict$ci, lty="dashed")
abline(h=mean(mod.predict$ci[mod.predict$ci$Year<=0869, "mean"]), col="red", lty="dashed")
# abline(h=0, col="red", lty="dashed")

# Fitting the derivatives, KEEPING THE SIMS!
source(file.path(path.gamm.func, "Calculate_GAMM_Derivs.R"))
mod.deriv <- calc.derivs(mod.gam, newdata=dat.temp, vars="Year", return.sims=T)
summary(mod.deriv)

summary(mod.deriv$ci)

plot(mean ~ Year, data=mod.deriv$ci, ylim=range(mod.deriv$ci[,c("mean", "lwr", "upr")]), type="l")
lines(lwr ~ Year, data=mod.deriv$ci, lty="dashed")
lines(upr ~ Year, data=mod.deriv$ci, lty="dashed")
abline(h=0, col="red", lty="dashed")


# Calculating the derivative as the first difference from the SIMS posteriors
library(zoo)
sims.out <- mod.predict$sims
# sims.out <- sims.out[order(sims.out$x, decreasing=T),] # Flipping the order so we can get the change from the *previous* timestep

head(sims.out[,1:10])
diff.calc <- sims.out
diff.calc[,which(substr(names(sims.out),1,1)=="X")] <- rbind(NA, apply(sims.out[,which(substr(names(sims.out),1,1)=="X")], 2, FUN=diff, na.pad=T))
head(diff.calc[,1:10])

diff.ci <- diff.calc[,c("Model", "Effect","x")]
diff.ci$Year <- diff.ci$x
diff.ci$mean <- apply(diff.calc[,which(substr(names(diff.calc),1,1)=="X")], 1, FUN=mean, na.rm=T)
diff.ci$lwr <- apply(diff.calc[,which(substr(names(diff.calc),1,1)=="X")], 1, FUN=quantile, 0.025, na.rm=T)
diff.ci$upr <- apply(diff.calc[,which(substr(names(diff.calc),1,1)=="X")], 1, FUN=quantile, 0.975, na.rm=T)
diff.ci$sig <- as.factor(ifelse(diff.ci$lwr*diff.ci$upr>0, "*", NA))
summary(diff.ci)

# diff.ci <- diff.ci[order(diff.ci$x, decreasing=F),] # Flipping the order so we can get the change from the *previous* timestep

plot(mean ~ Year, data=diff.ci, ylim=range(diff.ci[,c("mean", "lwr", "upr")], na.rm=T), type="l")
lines(lwr ~ Year, data=diff.ci, lty="dashed")
lines(upr ~ Year, data=diff.ci, lty="dashed")
abline(h=0, col="red", lty="dashed")


# test <- c(1,2,4,7,13)
# diff(test)

# Comparing signficiant periods of change based on derivative and sim first diffs
par(mfrow=c(2,1))
plot(mean ~ Year, data=mod.deriv$ci, ylim=range(mod.deriv$ci[,c("mean", "lwr", "upr")]), type="l", 
     main="Derivs", ylab="Rate of Change")
lines(lwr ~ Year, data=mod.deriv$ci, lty="dashed")
lines(upr ~ Year, data=mod.deriv$ci, lty="dashed")
points(mean ~ Year, data=mod.deriv$ci[!is.na(mod.deriv$ci$sig),], pch=19, cex=0.5, col="blue")
abline(h=mean(mod.deriv$ci[mod.deriv$ci$Year<=0869, "mean"]), col="red", lty="dashed")
# abline(h=mean(mod.deriv$ci[mod.deriv$ci$Year<=0869, "upr"]), col="red", lty="dashed")
# abline(h=mean(mod.deriv$ci[mod.deriv$ci$Year<=0869, "lwr"]), col="red", lty="dashed")

plot(mean ~ Year, data=diff.ci, ylim=range(diff.ci[,c("mean", "lwr", "upr")],na.rm=T), type="l", 
     main="Diffs", ylab="Rate of Change")
lines(lwr ~ Year, data=diff.ci, lty="dashed")
lines(upr ~ Year, data=diff.ci, lty="dashed")
points(mean ~ Year, data=diff.ci[!is.na(diff.ci$sig),], pch=19, cex=0.5, col="blue")
abline(h=mean(diff.ci[diff.ci$Year<=0869, "mean"],na.rm=T), col="red", lty="dashed")
# abline(h=mean(diff.ci[diff.ci$Year<=0869, "upr"]), col="red", lty="dashed")
# abline(h=mean(diff.ci[diff.ci$Year<=0869, "lwr"]), col="red", lty="dashed")



par(mfrow=c(2,1))
plot(mean ~ Year, data=mod.predict$ci, ylim=range(mod.predict$ci[,c("mean", "lwr", "upr")]), type="l", 
     main="Derivs", ylab="Temp Deviation")
lines(lwr ~ Year, data=mod.predict$ci, lty="dashed")
lines(upr ~ Year, data=mod.predict$ci, lty="dashed")
points(mean ~ Year, data=mod.predict$ci[!is.na(mod.deriv$ci$sig),], pch=19, cex=0.5, col="blue")
abline(h=mean(mod.predict$ci[mod.predict$ci$Year<=0869, "mean"]), col="red", lty="dashed")
# abline(h=mean(mod.predict$ci[mod.predict$ci$Year<=0869, "upr"]), col="red", lty="dashed")
# abline(h=mean(mod.predict$ci[mod.predict$ci$Year<=0869, "lwr"]), col="red", lty="dashed")

plot(mean ~ Year, data=mod.predict$ci, ylim=range(mod.predict$ci[,c("mean", "lwr", "upr")]), type="l", 
     main="Diffs", ylab="Temp Deviation")
lines(lwr ~ Year, data=mod.predict$ci, lty="dashed")
lines(upr ~ Year, data=mod.predict$ci, lty="dashed")
points(mean ~ Year, data=mod.predict$ci[!is.na(diff.ci$sig),], pch=19, cex=0.5, col="blue")
abline(h=mean(mod.predict$ci[mod.predict$ci$Year<=0869, "mean"]), col="red", lty="dashed")
# abline(h=mean(mod.predict$ci[mod.predict$ci$Year<=0869, "upr"]), col="red", lty="dashed")
# abline(h=mean(mod.predict$ci[mod.predict$ci$Year<=0869, "lwr"]), col="red", lty="dashed")

# ------------------

par(mfrow=c(1,1))
plot(mean ~ Year, data=mod.predict$ci, ylim=range(mod.predict$ci[,c("mean", "lwr", "upr")]), type="l", 
     main="Derivs", ylab="Temp Deviation", lwd=2, col="red")
iters <- sample(1:100, 10)
for(i in 1:length(iters)){
  lines(mod.predict$sims[,paste0("X", iters[i])] ~ mod.predict$sims$Year)
}
lines(mean ~ Year, data=mod.predict$ci, lwd=4, col="red")
# -------------------------------------------


