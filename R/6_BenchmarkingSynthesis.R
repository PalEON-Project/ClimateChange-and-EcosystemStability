# -------------------------------------------
# Sythesizing benchmarks from script 5 to come up with an overall assement of 
# model performance
# Author: Christy Rollinson, crollinson@gmail.com
#
# ------------------
# Objectives:
# ------------------
# A. Model Benchmarking Score Card
#    1. Summary by benchmark dataset of model & benchmark mean & SD, pair-wise bias, and the number 
#       of sites the model was within the observed range/95%CI
# B. Correlating Overall Model Behaviorwith Model Characteristics
#    1. Use Ordination with the benchmarking scores to see if models group and whether model characteristics
#       explain those relationships
# ------------------
#
#
# ------------------
# Workflow
# ------------------
# 1. Set up file structure, libraries, etc.
# 2. Load and combine benchmark scores into tables
# 3. Summary by Benchmark (Table + Figure)
# 4. Model Score Cards
# 5. Ordination attempt
# ------------------
# -------------------------------------------
rm(list=ls())

# -------------------------------------------
# 1. Set up file structure, libraries, etc.
# -------------------------------------------
{
  library(car)
  library(ggplot2)
  # library(mgcv)
  # setwd("~/Desktop/Research/PalEON_CR/PalEON_MIP_Site/Analyses/Change-and-Stability") # Path to this project github repository: https://github.com/PalEON-Project/MIP-Change-and-Stability.git
  setwd("~/Dropbox/PalEON_CR/PalEON_MIP_Site/Analyses/Change-and-Stability") # Path to this project github repository: https://github.com/PalEON-Project/MIP-Change-and-Stability.git
  # path.gamm.func <- "~/Desktop/R_Functions/"  # Path to github repository of my GAMM helper functions: https://github.com/crollinson/R_Functions.git
  inputs    <- "Data/" # Path to my cleaned model output
  
  # mip.utils <- "~/Desktop/Research/PalEON_CR/MIP_Utils/" # Path to PalEON MIP Utility repository: https://github.com/PalEON-Project/MIP_Utils.git
  # path.raw <- "~/Desktop/Research/PalEON_CR/PalEON_MIP_Site/phase1a_model_output/" # Path to raw model output
  mip.utils <- "~/Dropbox/PalEON_CR/MIP_Utils/" # Path to PalEON MIP Utility repository: https://github.com/PalEON-Project/MIP_Utils.git
  path.raw <- "~/Dropbox/PalEON_CR/PalEON_MIP_Site/phase1a_model_output/" # Path to raw model output
  
  bench.dir <- "Data/Benchmarking"
  out.dir <- "Data/BenchmarkingSynthesis" # Path to where the analysis output should go
  fig.dir <- "Figures/BenchmarkingSynthesis" # Path to where figures should go
  raw.dir <- "raw_data/Benchmarks/"
  
  if(!dir.exists(out.dir)) dir.create(out.dir)
  if(!dir.exists(fig.dir)) dir.create(fig.dir)
  
  hips <- data.frame(Site = c( "PHA",  "PHO",  "PUN",  "PBL",  "PDL",  "PMB"),
                     lat  = c( 42.54,  45.25,  46.22,  46.28,  47.17,  43.61),
                     lon  = c(-72.18, -68.73, -89.53, -94.58, -95.17, -82.83),
                     latmin = c( 42.50,  45.00,  46.00,  46.00,  47.00,  43.50),
                     latmax = c( 43.00,  45.50,  46.50,  46.50,  47.50,  44.00),
                     lonmin = c(-72.50, -69.00, -90.00, -95.00, -95.50, -83.00),
                     lonmax = c(-72.00, -68.50, -89.50, -94.50, -95.00, -82.50))
  
  
  # Load the baseline ecosys file
  ecosys <- read.csv("Data/PalEON_MIP_Yearly_withFcomp.csv")
  
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
  
  # Getting a list of the models
  model.names <- data.frame(Model=unique(ecosys$Model), Model.Order=unique(ecosys$Model.Order))
  model.names$Model.Name <- recode(model.names$Model, "'ed2'='ED2'; 'ed2.lu'='ED2-LU'; 
                                   'clm.bgc'='CLM45BGC'; 'clm.cn'='CLM45CN'; 
                                   'lpj.wsl'='LPJ-WSL'; 'lpj.guess'='LPJ-GUESS'; 
                                   'jules.stat'='JULES'; 'jules.triffid'='JULES_TRIFFID';
                                   'linkages'='LINKAGES'; 'sibcasa'='SiBCASA'")
  model.names$Version    <- recode(model.names$Model, "'ed2'='v7'; 'ed2.lu'='v8'; 
                                   'clm.bgc'='v5.1'; 'clm.cn'='v3.1'; 
                                   'lpj.wsl'='v6'; 'lpj.guess'='v6'; 
                                   'jules.stat'='v2'; 'jules.triffid'='v1';
                                   'linkages'='v2.1'; 'sibcasa'='v1'")
  model.names$Mod.Type    <- recode(model.names$Model, "'ed2'='ED2'; 'ed2.lu'='ED2'; 
                                    'clm.bgc'='CLM'; 'clm.cn'='CLM'; 
                                    'lpj.wsl'='LPJ.WSL'; 'lpj.guess'='LPJ.GUESS'; 
                                    'jules.stat'='JULES'; 'jules.triffid'='JULES';
                                    'linkages'='LINKAGES'; 'sibcasa'='SIBCASA'")
}
# -------------------------------------------

# -------------------------------------------
# 2. Load and combine benchmark scores into tables
#    Note: Loading just the model performance scores to simplify things
# -------------------------------------------
{
  cols.keep <- c("Dataset", "Var", "Model", "Model.Order", "bench.mean", "bench.sd", "model.mean", "model.sd", "bias", "bias.per", "bias.abs", "bias.sd", "bias.per.sd", "bias.abs.sd")
  
  # 2.01. Composition -- STEPPS 1 (contibutor/contact: Andria)
  comp.stepps <- read.csv(file.path(bench.dir, "Composition_STEPPS2_ModelBenchmarks_Summary.csv"))
  comp.stepps$Dataset <- as.factor("STEPPS")
  comp.stepps$Var <- as.factor("Fcomp")
  names(comp.stepps)[which(names(comp.stepps)=="stepps.mean")] <- "bench.mean"
  names(comp.stepps)[which(names(comp.stepps)=="stepps.sd")]   <- "bench.sd"
  comp.stepps <- comp.stepps[,cols.keep]
  comp.stepps
  
  # 2.02. Composition -- SetVeg (contributor/contact: Simon? Jody?)
  comp.setveg <- read.csv(file.path(bench.dir, "Composition_SetVeg_ModelBenchmarks_Summary.csv"))
  comp.setveg$Dataset <- as.factor("SetVeg")
  comp.setveg$Var <- as.factor("Fcomp")
  names(comp.setveg)[which(names(comp.setveg)=="svcomp.mean")] <- "bench.mean"
  names(comp.setveg)[which(names(comp.setveg)=="svcomp.sd")]   <- "bench.sd"
  comp.setveg <- comp.setveg[,cols.keep]
  comp.setveg
  
  # 2.03. Biomass -- SetVeg  (contributor/contact: Jody? Kelly?)
  biom.setveg <- read.csv(file.path(bench.dir, "Biomass_SetVeg_ModelBenchmarks_Summary.csv"))
  biom.setveg$Dataset <- as.factor("SetVeg")
  biom.setveg$Var <- as.factor("AGB")
  names(biom.setveg)[which(names(biom.setveg)=="svbiom.mean")] <- "bench.mean"
  names(biom.setveg)[which(names(biom.setveg)=="svbiom.sd")]   <- "bench.sd"
  biom.setveg <- biom.setveg[,cols.keep]
  biom.setveg
  
  # 2.04. Biomass -- NBCD (contributor/contact: Christy)
  biom.nbcd <- read.csv(file.path(bench.dir, "Biomass_NBCD_ModelBenchmarks_Summary.csv"))
  biom.nbcd$Dataset <- as.factor(paste("NBCD", biom.nbcd$Ref, sep="."))
  biom.nbcd$Var <- as.factor("AGB")
  names(biom.nbcd)[which(names(biom.nbcd)=="nbcd.mean")] <- "bench.mean"
  names(biom.nbcd)[which(names(biom.nbcd)=="nbcd.sd")]   <- "bench.sd"
  biom.nbcd <- biom.nbcd[,cols.keep]
  biom.nbcd
  
  # 2.05. Carbon Fluxes - Flux Towers (contributor/contact: Dave)
  flux.tower <- read.csv(file.path(bench.dir, "Fluxes_FluxTowers_ModelBenchmarks_Summary.csv"))
  flux.tower$Dataset <- as.factor("Ameriflux")
  names(flux.tower)[which(names(flux.tower)=="Flux")] <- "Var"
  names(flux.tower)[which(names(flux.tower)=="tower.mean")] <- "bench.mean"
  names(flux.tower)[which(names(flux.tower)=="tower.sd")]   <- "bench.sd"
  flux.tower <- flux.tower[,cols.keep]
  flux.tower
  
  # 2.06. MODIS - LAI, GPP (contributor/contact: Bethany)
  modis <- read.csv(file.path(bench.dir, "MODIS_ModelBenchmarks_Summary.csv"))
  modis$Dataset <- as.factor("MODIS")
  names(modis)[which(names(modis)=="Flux")] <- "Var"
  names(modis)[which(names(modis)=="modis.mean")] <- "bench.mean"
  names(modis)[which(names(modis)=="modis.sd")]   <- "bench.sd"
  modis <- modis[,cols.keep]
  modis
  
  
  # # Combine all benchmarks into 1 data frame:
  models.bench <- rbind(comp.stepps, comp.setveg, biom.setveg, biom.nbcd, flux.tower, modis)
  models.bench$Benchmark <- as.factor(paste(models.bench$Dataset, models.bench$Var, sep="-"))
  summary(models.bench)
 
  write.csv(models.bench, file.path(out.dir, "Benchmarks_Summary_All.png"), row.names=F) 
}
# -------------------------------------------



# -------------------------------------------
# 3. Summary by Benchmark (Table + Figure)
# -------------------------------------------
{
  col.model <- model.colors[model.colors$Model.Order %in% unique(models.bench$Model.Order),"color"]
  shape.model <- model.colors[model.colors$Model.Order %in% unique(models.bench$Model.Order),"shape"]
  
  models.bench$Var  <- factor(models.bench$Var , levels=c("GPP", "RE", "NEE", "LAI", "AGB", "Fcomp"))
  summary(models.bench)

  # Making a summary figure
  png(file.path(fig.dir, "Benchmarks_Summary_All.png"), height=6, width=8, units="in", res=220)
  print(
  ggplot(data=models.bench) +
    facet_grid(.~Var, scales="free_x") +
    geom_hline(yintercept=0, linetype="dashed") +
    # geom_hline(yintercept=-100, linetype="dashed", size=0.5, color="gray50") +
    geom_pointrange(aes(x=Dataset,  y=bias.per*100, ymin=(bias.per-bias.per.sd)*100, ymax=(bias.per+bias.per.sd)*100, color=Model.Order, shape=Model.Order), position=position_jitter(0.75), size=0.8, alpha=0.8) +
    scale_color_manual(values=paste(col.model)) +
    scale_shape_manual(values=shape.model+14) +
    scale_y_continuous(name="Percent Bias") +
    theme_bw() +
    theme(axis.text.x=element_text(angle=45, hjust=1)) +
    coord_cartesian(ylim=c(min(models.bench$bias.per-models.bench$bias.per.sd, na.rm=T)*100, quantile(models.bench$bias.per+models.bench$bias.per.sd, 0.95, na.rm=T)*100))
  )
  dev.off()
  }
# -------------------------------------------

# -------------------------------------------
# 4. Model Score Cards
# -------------------------------------------
{
  # Figure out which columns are numeric should be rounded
  cols.num <- vector()
  for(i in 1:ncol(models.bench)){
    if(is.numeric(models.bench[,i])) cols.num <- c(cols.num, i)
  }
  
  for(m in unique(models.bench$Model.Order)){
    dat.now <- models.bench[models.bench$Model.Order==m,]
    
    # Round everything to make cleaner graphs;
    # Note: not doing this to the raw data so we don't create headaches
    dat.now[,cols.num] <- round(dat.now[,cols.num], 1)
    dat.now <- dat.now[order(dat.now$Var),] # Organizing by variable
    
    dat.print <- dat.now[,c("Var", "Dataset")]
    dat.print$Benchmark <- paste0(dat.now$bench.mean, " (", dat.now$bench.sd   , ")")
    dat.print$Model     <- paste0(dat.now$model.mean, " (", dat.now$model.sd   , ")")
    dat.print$Bias      <- paste0(dat.now$bias      , " (", dat.now$bias.sd    , ")")
    dat.print$Bias.Abs  <- paste0(dat.now$bias.abs  , " (", dat.now$bias.abs.sd, ")")
    
    write.csv(dat.print, file.path(out.dir, paste0("ScoreCard_", m, ".csv")), row.names=F)
  }
}
# -------------------------------------------




# -------------------------------------------
# 5. Ordinate models & benchmarks to see if there are clear trends or if it ends up
#    being mostly linear
# -------------------------------------------
{
  # Reshape data to be in Observation (Model) x Variable (Benchmark) format
  # Using the mean absolute bias so we don't have to deal with negative values
  library(reshape2)
  bench.dat <- recast(data=models.bench[,c("Model.Order", "Benchmark", "bias.abs")], Model.Order ~ Benchmark)
  row.names(bench.dat) <- bench.dat$Model.Order
  bench.dat <- bench.dat[,!names(bench.dat)=="Model.Order"]
  bench.dat
  
  # Flipping the way of dealign with the data, just to try; 
  # NOTE: I don't think this is what we want to do, but we'll try it
  bench.dat2 <- recast(data=models.bench[,c("Model.Order", "Benchmark", "bias.abs")], Benchmark ~ Model.Order)
  row.names(bench.dat2) <- bench.dat2$Benchmark
  bench.dat2 <- bench.dat2[,!names(bench.dat2)=="Benchmark"]
  bench.dat2
  
  # Dealing with NAs
  #  -- For now, just making 0 although this is a BAD assumption
  bench.dat[is.na(bench.dat)] <- 0
  bench.dat
  
  bench.dat2[is.na(bench.dat2)] <- 0
  bench.dat2
  
  # Running an NMDS ordination using Vegan
  library(vegan)
  set.seed(1624)
  bench.ord <- metaMDS(bench.dat, dist="bray")
  summary(bench.ord)
  
  bench.plot <- ordiplot(bench.ord, type="none")
  points(bench.plot, "species")
  text(bench.plot, "sites")
  
  # Trying it with the axes flipped
  bench.ord2 <- metaMDS(bench.dat2, dist="bray")
  summary(bench.ord2)
  
  bench.plot2 <- ordiplot(bench.ord2, type="none")
  # points(bench.plot2, "sites")
  text(bench.plot2, "species")
  text(bench.plot2, "sites", col="blue", cex=0.75)

}
# -------------------------------------------
