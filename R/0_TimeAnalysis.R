# Run the temporal gams
analyze.time <- function(dat, Y, fac.fit, k.freq, path.gamm.func){
  # dat    = data set to fit the gam to 
  # Y      = what response variable of interest
  # fac.by = what factor we want to use to fit different gams (i.e. sites or models)
  # k.freq = number of years for each knot; even those these aren't necessarily 
  #          spaced evenly, it helps make sure that we give models with different 
  #          lengths similar flexibilities
  # path.gamm.func = file path to my R_Functions github repository
  
  # Setting up the data frame with the generalized column names
  dat.temp <- data.frame(Year=dat$Year)
  dat.temp[,"Y"] <- dat[,Y]
  dat.temp[,"fac.fit"] <- dat[,fac.fit]
  
  # Determine the number of knots to use in the model
  k.use = round((max(dat$Year) - min(dat$Year))/k.freq,0)

  # Fit the gam
  mod.gam <- gam(Y ~ s(Year, by=fac.fit, k=k.use) + fac.fit, data=dat.temp)
  
  # Doing a fancy prediction to get full 95% CIs
  source(file.path(path.gamm.func, "Calculate_GAMM_Posteriors.R"))
  mod.predict <- post.distns(model.gam=mod.gam, model.name=Y, n=100, newdata=dat.temp, vars="Year", terms=F)
  mod.predict$Year <- dat$Year
  mod.predict[,fac.fit] <- dat.temp$fac.fit
  
  # Calculating the derivatives to identify periods of significant change
  source(file.path(path.gamm.func, "Calculate_GAMM_Derivs.R"))
  mod.deriv    <- calc.derivs(mod.gam, newdata=dat.temp, vars="Year")
  names(mod.deriv)[which(names(mod.deriv) %in% c("mean", "lwr", "upr"))] <- paste0("deriv.", c("mean", "lwr", "upr"))
  summary(mod.deriv)
  
  mod.out <- list()
  mod.out[["model"]] <- mod.gam
  mod.out[["out"]] <- cbind(mod.predict, mod.deriv[,!names(mod.deriv) %in% c("Year", "fac.fit")])
    
  
  return(mod.out)
}
