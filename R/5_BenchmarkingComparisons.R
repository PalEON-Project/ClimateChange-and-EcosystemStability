# -------------------------------------------
# Quantifying and comparing periods & magnitude of change through time & space among met drivers, 
# models, & ecosystem variables
#
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
#
# 2. Model Benchmarking
#    2.01. Composition -- STEPPS 1 (contibutor/contact: Andria)
#    2.02. Composition -- SetVeg (contributor/contact: Simon? Jody?)
#    2.03. Biomass -- SetVeg  (contributor/contact: Jody? Kelly?)
#    2.04. Biomass -- NBCD (contributor/contact: Christy)
#    2.05. Carbon Fluxes - Flux Towers (contributor/contact: Dave)
#    2.06. MODIS - LAI, GPP (contributor/contact: Bethany)
#    2.07. Tree RIngs -- AGB, AGBincrement (contributor/contact: Ross, Alex)
# ------------------
# -------------------------------------------
rm(list=ls())

# -------------------------------------------
# 1. Set up file structure, libraries, etc.
# -------------------------------------------
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


out.dir <- "Data/Benchmarking" # Path to where the analysis output should go
fig.dir <- "Figures/Benchmarking" # Path to where figures should go
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
# -------------------------------------------

# -------------------------------------------
# 2. Extract & Format Benchmarks
# -------------------------------------------
{
  # ------------------
  #    2.01. Composition -- STEPPS 1
  #     -- Finding the dominant PFT by model specifications
  # ------------------
  {
	  # Load STEPPS benchmark & PFT crosswalk
	  #  -- This is just the Fcomp of the dominant PFT
	  load(file.path(out.dir, "STEPPS2_ModelPFTs.RData"))
	  summary(stepps2.pft$ED2)
  
	  # Read in the model-STEPPS crosswalk
	  pft.xwalk <- read.csv("raw_data/STEPPS_conversion.csv")
	  pft.xwalk$taxon <- toupper(pft.xwalk$STEPPS)
	  pft.xwalk$taxon <- recode(pft.xwalk$taxon, "'OTHER CONIFER'='OTHER.CONIFER'; 'OTHER HARDWOOD'='OTHER.HARDWOOD'")
	  summary(pft.xwalk)
  
	  # Load the table that says which PFT
	  pft.mod <- read.csv("raw_data/ModelPFT_Table.csv")
	  pft.mod

	  source(file.path(mip.utils, "Phase1_sites/extract_output_site.R"))
  
	  
	  for(m in unique(ecosys$Model)){
  		if(m == "sibcasa") next # Sibcasa has no composition, so skip
  		mod=ifelse(m=="linkages", "LINKAGES.STEPPS",paste(model.names[model.names$Model==m, "Mod.Type"]))
  		sites.stepps <- unique(stepps2.pft[[mod]]$site)
  		print(paste0(m))
  		model.name = model.names[model.names$Model==m, "Model.Name"]
  		pft.tmp <- extract.paleon.site(model=model.name, 
  									   model.dir=file.path(path.raw, paste(model.name, model.names[model.names$Model==m, "Version"], sep=".")), 
  									   sites=paste(sites.stepps), 
  									   vars="Fcomp")
  
  		if(m %in% c("jules.stat")){
  		  pft.tmp <- extract.paleon.site(model=model.name, 
  										 model.dir=file.path(path.raw, paste(model.name, model.names[model.names$Model==m, "Version"], sep=".")), 
  										 sites=paste(sites.stepps), 
  										 vars="LAI")  
  		  lai.sums <- apply(pft.tmp$LAI, c(1,3), FUN=sum)
  		  pft.tmp$Fcomp <- pft.tmp$LAI
  		  for(i in 1:dim(pft.tmp$LAI)[3]){
  			pft.tmp$Fcomp[,,i] <- pft.tmp$LAI[,,i]/lai.sums[,i]
  		  }
  		}
  	
  		if(m=="lpj.guess") pft.tmp$Fcomp <- pft.tmp$Fcomp[,1:(dim(pft.tmp$Fcomp)[2]-1),]
  	
  		# Aggregate to annual
  		if(dim(pft.tmp$Fcomp)[1]>1161){
  		  yr.rows <- seq(1, dim(pft.tmp$Fcomp)[1], by=12)
  		  out.tmp <- array(dim=c(length(yr.rows), c(dim(pft.tmp$Fcomp)[2:3])))
  		  for(i in 1:length(yr.rows)){
  			out.tmp[i,,] <- apply(pft.tmp$Fcomp[yr.rows[i]:(yr.rows[i]+11),,], c(2,3), FUN=mean)
  		  }
  		  pft.tmp$Fcomp <- out.tmp
  		}
  		dimnames(pft.tmp$Fcomp)[[3]] <- sites.stepps
  		dimnames(pft.tmp$Fcomp)[[1]] <- 849+1:dim(pft.tmp$Fcomp)[1]
  	
  		for(i in 1:dim(pft.tmp$Fcomp)[3]){
  		  site.now <- dimnames(pft.tmp$Fcomp)[[3]][i]
  		  pft.take <- unique(stepps2.pft[[mod]][stepps2.pft[[mod]]$site==site.now, "pft"])
  		  col.max <- which(paste(pft.mod[,mod])==paste(pft.take))
  		  if(length(col.max)==1){
  			fcomp <- pft.tmp$Fcomp[,col.max,i]
  		  } else {
  			fcomp <- apply(pft.tmp$Fcomp[,col.max,i], 1, sum) 
  		  }
  		  pft.max <- data.frame(Model=m,
  								Site=as.factor(dimnames(pft.tmp$Fcomp)[[3]][i]),
  								Year=as.numeric(dimnames(pft.tmp$Fcomp)[[1]]),
  								PFTmax=pft.take,
  								Fcomp=fcomp)
  	  
  		  # # Aggregate to STEPPS time
  	      # stepps2.pft[[mod]][, m] <- NA
  		  for(yr in unique(stepps2.pft[[mod]]$year)){
  			stepps2.pft[[mod]][stepps2.pft[[mod]]$site == site.now & stepps2.pft[[mod]]$year==yr, m] <- mean(pft.max[pft.max$Year>=yr-50 & pft.max$Year<=yr+50,"Fcomp"])
  		  }      
  		} # End site loop
  		# Store some calculations
  		stepps2.pft[[mod]][,paste0(m, ".bias")] <- stepps2.pft[[mod]][,m]-stepps2.pft[[mod]][,"mean"]
  		# stepps2.pft[[mod]][,paste0(m, ".bias.per")] <- (stepps2.pft[[mod]][,m])/stepps2.pft[[mod]][,"mean"]-1
  		stepps2.pft[[mod]][,paste0(m, ".check")] <- as.factor(ifelse(stepps2.pft[[mod]][,m]>=stepps2.pft[[mod]][,"lwr"] & stepps2.pft[[mod]][,m]<=stepps2.pft[[mod]][,"upr"], "*", NA))
	  } # end model loop
  
	  # Change LINKAGES.STEPPS to just LINKAGES
	  names(stepps2.pft)[which(names(stepps2.pft)=="LINKAGES.STEPPS")] <- "LINKAGES"
	  
	  # STEPPS Benchmarking Performance
	  for(i in 1:length(stepps2.pft)){
  		models.pull <- paste(model.names[model.names$Mod.Type==names(stepps2.pft)[i], "Model"])
  		if(length(models.pull)>1){
  		  dat.base <- stack(stepps2.pft[[i]][,models.pull])
  		  names(dat.base) <- c("model.mean", "Model")
  		  dat.base[,c("bias")] <- stack(stepps2.pft[[i]][,paste0(models.pull, ".bias")])[,1]
  		} else{
  		  dat.base <- data.frame(model.mean=stepps2.pft[[i]][,models.pull], Model=models.pull)
  		  dat.base[,c("bias")] <- stepps2.pft[[i]][,paste0(models.pull, ".bias")]
  		}
  		
  		dat.base[,c("year", "site", "pft")] <- stepps2.pft[[i]][,c("year", "site", "pft")]
  		dat.base[,c("stepps.mean", "stepps.lwr", "stepps.upr")] <- stepps2.pft[[i]][,c("mean", "lwr", "upr")]
  	
  		dat.base <- dat.base[,c("site", "year", "Model", "pft", "stepps.mean", "stepps.lwr", "stepps.upr", "model.mean", "bias")]
  		if(i==1){
  		  stepps.bench <- dat.base
  		} else {
  		  stepps.bench <- rbind(stepps.bench, dat.base)
  		}
	  }
	  stepps.bench$bias.per <- stepps.bench$bias/stepps.bench$stepps.mean
	  stepps.bench$bias.abs <- abs(stepps.bench$bias)
	  stepps.bench$CI.check <- as.factor(ifelse(stepps.bench$model.mean>=stepps.bench$stepps.lwr & stepps.bench$model.mean<=stepps.bench$stepps.upr, "*", NA))
	  write.csv(stepps.bench, file.path(out.dir, "Composition_STEPPS2_ModelBenchmarks.csv"), row.names=F)

	  png(file.path(fig.dir, "Composition_STEPPS2_ModelBenchmarks.png"), height=8, width=10, units="in", res=180)
	  ggplot(data=stepps.bench) +
  		facet_wrap(~Model) +
  		geom_ribbon(aes(x=year, ymin=stepps.lwr, ymax=stepps.upr, fill=site), alpha=0.3) +
  		geom_line(aes(x=year, y=stepps.mean, color=site), linetype="dashed") +
  		geom_line(aes(x=year, y=model.mean, color=site), size=2) +
  		theme_bw() +
  		theme(legend.position="top")
	  dev.off()

	  # Lookign at some quick summary stats by model
	  stepps.bench2 <- aggregate(stepps.bench[,c("stepps.mean", "stepps.lwr", "stepps.upr", "model.mean", "bias", "bias.per", "bias.abs")],
								               by=stepps.bench[,c("Model", "site", "pft")],
								               FUN=mean, na.rm=T)
	  stepps.bench2[,c("stepps.sd", "model.sd", "bias.sd", "bias.per.sd", "bias.abs.sd")] <- aggregate(stepps.bench[,c("stepps.mean", "model.mean", "bias", "bias.per", "bias.abs")],
	                                                                                    by=stepps.bench[,c("Model", "site", "pft")],
	                                                                                    FUN=sd, na.rm=T)[,c("stepps.mean", "model.mean", "bias", "bias.per", "bias.abs")]
	  stepps.bench2$CI <- as.factor(ifelse(stepps.bench2$model.mean>=stepps.bench2$stepps.lwr & stepps.bench2$model.mean<=stepps.bench2$stepps.upr, "in", "out"))
	  stepps.bench2 <- merge(stepps.bench2, model.names[,c("Model", "Model.Order")], all.x=T)
	  

	  png(file.path(fig.dir, "Composition_STEPPS2_ModelBenchmarks_SiteMeans.png"), height=8, width=10, units="in", res=180)
	  ggplot(data=stepps.bench2[,]) +
  		facet_wrap(~site) +
  		geom_pointrange(aes(x=Model.Order, y=stepps.mean, ymin=stepps.lwr, ymax=stepps.upr), color="gray") +
  		geom_point(aes(x=Model.Order, y=model.mean, color=Model.Order), size=3) +
  		scale_y_continuous(name="Fcomp STEPPS Dominant PFT") +
  		scale_x_discrete(name="Model") +
  		scale_color_manual(values=paste(model.colors$color), name="Model") + 
  		theme_bw() +
  		theme(axis.text.x=element_blank())
	  dev.off()

	  stepps.performance <- aggregate(stepps.bench2[,c("stepps.mean", "model.mean", "bias", "bias.per", "bias.abs")],
	                                  by=stepps.bench2[,c("Model", "Model.Order")],
	                                  FUN=mean, na.rm=T)
	  stepps.performance[,c("stepps.sd", "model.sd", "bias.sd", "bias.per.sd", "bias.abs.sd")] <- aggregate(stepps.bench2[,c("stepps.mean", "model.mean", "bias", "bias.per", "bias.abs")],
	                                                                                         by=stepps.bench2[,c("Model", "Model.Order")],
	                                                                                         FUN=sd, na.rm=T)[,c("stepps.mean", "model.mean", "bias", "bias.per", "bias.abs")]
	  stepps.performance
	  # On average, JULES-STATIC & CLM-BGC (two static veg models) had the the lowest bias, 
	  #  but CLM-LINKAGES, CLM-CN, and LPJ-WSL (two dynamic veg models) had the most consistent bias
	  write.csv(stepps.performance, file.path(out.dir, "Composition_STEPPS2_ModelBenchmarks_Summary.csv"), row.names=F)
  }
  # ------------------

  # ------------------
  #    2.02. Composition -- SetVeg
  # ------------------
  {
    # Load SetVeg benchmark & PFT crosswalk
    #  -- This is just the Fcomp of the dominant PFT
    setveg.pft <- read.csv(file.path(out.dir, "SetVeg_ModelPFTs.csv"))
    names(setveg.pft) <- c("Mod.Type", "Site", "pft", "svcomp.mean", "svcomp.lwr", "svcomp.upr")
    summary(setveg.pft)
    
    # Merging in the true model information
    setveg.pft$Mod.Type <- as.factor(ifelse(setveg.pft$Mod.Type=="LINKAGES.SETVEG", "LINKAGES", paste(setveg.pft$Mod.Type)))
    setveg.pft <- merge(setveg.pft, model.names[,c("Mod.Type", "Model", "Model.Order")], all.x=T)
    summary(setveg.pft)
    
    # Read in the model-STEPPS crosswalk
    pft.xwalk <- read.csv("raw_data/SetVeg_conversion.csv")
    pft.xwalk$taxon <- toupper(pft.xwalk$SetVeg)
    # pft.xwalk$taxon <- recode(pft.xwalk$taxon, "'OTHER CONIFER'='OTHER.CONIFER'; 'OTHER HARDWOOD'='OTHER.HARDWOOD'")
    summary(pft.xwalk)
    
    # Load the table that says which PFT
    pft.mod <- read.csv("raw_data/ModelPFT_Table.csv")
    pft.mod
    
    
    source(file.path(mip.utils, "Phase1_sites/extract_output_site.R"))
    
    for(m in unique(setveg.pft$Model)){
      if(m == "sibcasa") next # Sibcasa has no composition, so skip
      mod=paste(model.names[model.names$Model==m, "Mod.Type"])
      sites.setveg <- unique(setveg.pft[setveg.pft$Model==m, "Site"])
      print(paste0(m))
      model.name = model.names[model.names$Model==m, "Model.Name"]
      if(m %in% c("jules.stat")){
        pft.tmp <- extract.paleon.site(model=model.name, 
                                       model.dir=file.path(path.raw, paste(model.name, model.names[model.names$Model==m, "Version"], sep=".")), 
                                       sites=paste(sites.setveg), 
                                       vars="LAI")  
        lai.sums <- apply(pft.tmp$LAI, c(1,3), FUN=sum)
        pft.tmp$Fcomp <- pft.tmp$LAI
        for(i in 1:dim(pft.tmp$LAI)[3]){
          pft.tmp$Fcomp[,,i] <- pft.tmp$LAI[,,i]/lai.sums[,i]
        }
      } else {
        pft.tmp <- extract.paleon.site(model=model.name, 
                                       model.dir=file.path(path.raw, paste(model.name, model.names[model.names$Model==m, "Version"], sep=".")), 
                                       sites=paste(sites.setveg), 
                                       vars="Fcomp")
      }
      
      if(m=="lpj.guess") pft.tmp$Fcomp <- pft.tmp$Fcomp[,1:(dim(pft.tmp$Fcomp)[2]-1),]
      
      # Aggregate to annual
      if(dim(pft.tmp$Fcomp)[1]>1161){
        yr.rows <- seq(1, dim(pft.tmp$Fcomp)[1], by=12)
        out.tmp <- array(dim=c(length(yr.rows), c(dim(pft.tmp$Fcomp)[2:3])))
        for(i in 1:length(yr.rows)){
          out.tmp[i,,] <- apply(pft.tmp$Fcomp[yr.rows[i]:(yr.rows[i]+11),,], c(2,3), FUN=mean)
        }
        pft.tmp$Fcomp <- out.tmp
      }
      dimnames(pft.tmp$Fcomp)[[3]] <- sites.setveg
      dimnames(pft.tmp$Fcomp)[[1]] <- 849+1:dim(pft.tmp$Fcomp)[1]
      
      # Extract the fcomp for each site
      for(i in 1:dim(pft.tmp$Fcomp)[3]){
        site.now <- dimnames(pft.tmp$Fcomp)[[3]][i]
        pft.take <- unique(setveg.pft[setveg.pft$Model==m & setveg.pft$Site==site.now, "pft"])
        col.max <- which(paste(pft.mod[,mod])==paste(pft.take))
        if(length(col.max)==1){
          fcomp <- pft.tmp$Fcomp[,col.max,i]
        } else {
          fcomp <- apply(pft.tmp$Fcomp[,col.max,i], 1, sum) 
        }

        # Extract the Settlement Vegetation Era mean (1800-1850)
        setveg.pft[setveg.pft$Model==m & setveg.pft$Site==site.now, "model.mean"] <- round(mean(fcomp[which(names(fcomp)>=1800 & names(fcomp)<=1850)]),6)
        setveg.pft[setveg.pft$Model==m & setveg.pft$Site==site.now, "model.sd"] <- round(sd(fcomp[which(names(fcomp)>=1800 & names(fcomp)<=1850)]),6)
      } # End site loop
    } # end model loop

    # # Store some calculations
    setveg.pft$bias <- setveg.pft$model.mean - setveg.pft$svcomp.mean
    setveg.pft$bias.per <- setveg.pft$bias/setveg.pft$svcomp.mean
    setveg.pft$bias.abs <- abs(setveg.pft$bias)
    setveg.pft$range.check <- as.factor(ifelse(setveg.pft$model.mean + 2*setveg.pft$model.sd>=setveg.pft$svcomp.lwr & setveg.pft$model.mean - 2*setveg.pft$model.sd<=setveg.pft$svcomp.upr,
                                               "*", NA))
    summary(setveg.pft)
    write.csv(setveg.pft, file.path(out.dir, "Composition_SetVeg_ModelBenchmarks.csv"), row.names=F)
    
    # Order the sites from W to E
    setveg.pft$Site <- factor(setveg.pft$Site, levels=c("PDL", "PBL", "PUN", "PMB", "PHA", "PHO"))
    
    # Graphing the model bias in the setveg record
    png(file.path(fig.dir, "Composition_SetVeg_ModelBenchmarks_byModel.png"), height=8, width=10, units="in", res=180)
    ggplot(data=setveg.pft[,]) +
      facet_wrap(~Site) +
      geom_pointrange(aes(x=Model.Order, y=svcomp.mean, ymin=svcomp.lwr, ymax=svcomp.upr), color="black") +
      geom_pointrange(aes(x=Model.Order, y=model.mean, ymin=model.mean-2*model.sd, ymax=model.mean+2*model.sd, color=Model.Order)) +
      scale_y_continuous(name="Fcomp SetVeg Dominant PFT") +
      scale_x_discrete(name="Model") +
      scale_color_manual(values=paste(model.colors$color), name="Model") + 
      theme_bw() +
      theme(axis.text.x=element_text(angle=45, hjust=1))
    dev.off()
    
    png(file.path(fig.dir, "Composition_SetVeg_ModelBenchmarks_bySite.png"), height=8, width=10, units="in", res=180)
    ggplot(data=setveg.pft[,]) +
      facet_wrap(~Model.Order) +
      geom_pointrange(aes(x=Site, y=svcomp.mean, ymin=svcomp.lwr, ymax=svcomp.upr), color="black") +
      geom_pointrange(aes(x=Site, y=model.mean, ymin=model.mean-2*model.sd, ymax=model.mean+2*model.sd, color=Model.Order)) +
      scale_y_continuous(name="Fcomp SetVeg Dominant PFT") +
      scale_x_discrete(name="Model") +
      scale_color_manual(values=paste(model.colors$color), name="Model") + 
      theme_bw() +
      theme(axis.text.x=element_text(angle=45, hjust=1))
    dev.off()
    
    # Getting the average bias for each model
    setveg.mod <- aggregate(setveg.pft[,c("svcomp.mean", "model.mean", "bias", "bias.per", "bias.abs")],
                            by=setveg.pft[,c("Model", "Model.Order")],
                            FUN=mean)
    
    setveg.mod[,c("svcomp.sd", "model.sd", "bias.sd", "bias.per.sd", "bias.abs.sd")] <- aggregate(setveg.pft[,c("svcomp.mean", "model.mean", "bias", "bias.per", "bias.abs")],
                                                                  by=setveg.pft[,c("Model", "Model.Order")],
                                                                  FUN=sd)[,c("svcomp.mean", "model.mean", "bias", "bias.per", "bias.abs")]
    setveg.mod
    write.csv(setveg.mod, file.path(out.dir, "Composition_SetVeg_ModelBenchmarks_Summary.csv"), row.names=F)
    
    # Getting the average bias for each model
    setveg.site <- aggregate(setveg.pft[,c("svcomp.mean", "model.mean", "bias", "bias.per", "bias.abs")],
                             by=list(setveg.pft$Site),
                             FUN=mean)
    names(setveg.site)[1] <- "Site"
    
    setveg.site[,c("svcomp.sd", "model.sd", "bias.sd", "bias.per.sd", "bias.abs.sd")] <- aggregate(setveg.pft[,c("svcomp.mean", "model.mean", "bias", "bias.per", "bias.abs")],
                                                                   by=list(setveg.pft$Site),
                                                                   FUN=sd)[,c("svcomp.mean", "model.mean", "bias", "bias.per", "bias.abs")]
    setveg.site
    write.csv(setveg.site, file.path(out.dir, "Composition_SetVeg_SiteBenchmarks_Summary.csv"), row.names=F)
  }
  # ------------------
  
  # ------------------
  #    2.03. Biomass -- SetVeg
  # ------------------
  {
    # Read in the set veg biomass estimate
	  sv.biom <- read.csv(file.path(out.dir, "SetVeg_Biomass.csv"))
	  names(sv.biom)[4:7] <- c("svbiom.mean", "svbiom.sd", "svbiom.lwr", "svbiom.upr")
	  sv.biom <- sv.biom[complete.cases(sv.biom),] # get rid of sites outside of setveg biomass
	  sv.biom[,c("svbiom.mean", "svbiom.sd", "svbiom.lwr", "svbiom.upr")] <- sv.biom[,c("svbiom.mean", "svbiom.sd", "svbiom.lwr", "svbiom.upr")]*0.1*.5 # Convert Mg HA-1 to KgC m-2 
	  sv.biom
	  
	  mod.svbiom <- aggregate(ecosys[ecosys$Year>=1800 & ecosys$Year<=1850 & ecosys$Site %in% sv.biom$Site, "AGB"], 
	                          by=ecosys[ecosys$Year>=1800 & ecosys$Year<=1850 & ecosys$Site %in% sv.biom$Site, c("Site", "Model", "Model.Order")],
	                          FUN=mean)
	  names(mod.svbiom)[4] <- c("model.mean")
	  mod.svbiom$model.sd <- aggregate(ecosys[ecosys$Year>=1800 & ecosys$Year<=1850 & ecosys$Site %in% sv.biom$Site, "AGB"], 
	                                     by=ecosys[ecosys$Year>=1800 & ecosys$Year<=1850 & ecosys$Site %in% sv.biom$Site, c("Site", "Model", "Model.Order")],
	                                     FUN=sd)[,"x"]
	  summary(mod.svbiom)
	  
	  mod.svbiom <- merge(mod.svbiom, sv.biom)
	  
	  mod.svbiom$bias <- mod.svbiom$model.mean - mod.svbiom$svbiom.mean
	  mod.svbiom$bias.per <- mod.svbiom$bias/mod.svbiom$svbiom.mean
	  mod.svbiom$bias.abs <- abs(mod.svbiom$bias)
	  summary(mod.svbiom)
	  
	  write.csv(mod.svbiom, file.path(out.dir, "Biomass_SetVeg_ModelBenchmarks.csv"), row.names=F)
	  
	  # Order the sites from W to E
	  mod.svbiom$Site <- factor(mod.svbiom$Site, levels=c("PDL", "PBL", "PUN"))
	  
	  # Graphing the model bias in the setveg record
	  png(file.path(fig.dir, "Biomass_SetVeg_ModelBenchmarks_byModel.png"), height=8, width=10, units="in", res=180)
	  ggplot(data=mod.svbiom[,]) +
	    facet_wrap(~Site) +
	    geom_hline(aes(yintercept=svbiom.mean), color="black", size=1.5, linetype="solid") +
	    geom_hline(aes(yintercept=svbiom.mean-2*svbiom.sd), color="gray50", size=1, linetype="solid") +
	    geom_hline(aes(yintercept=svbiom.mean+2*svbiom.sd), color="gray50", size=1, linetype="solid") +
	    geom_pointrange(aes(x=Model.Order, y=model.mean, ymin=model.mean-2*model.sd, ymax=model.mean+2*model.sd, color=Model.Order)) +
	    scale_y_continuous(name="Fcomp SetVeg Dominant PFT") +
	    scale_x_discrete(name="Model") +
	    scale_color_manual(values=paste(model.colors$color), name="Model") + 
	    theme_bw() +
	    theme(axis.text.x=element_text(angle=45, hjust=1))
	  dev.off()
	  
	  png(file.path(fig.dir, "Biomass_SetVeg_ModelBenchmarks_bySite.png"), height=8, width=10, units="in", res=180)
	  ggplot(data=mod.svbiom[,]) +
	    facet_wrap(~Model.Order) +
	    geom_pointrange(aes(x=Site, y=svbiom.mean, ymin=svbiom.mean-2*svbiom.sd, ymax=svbiom.mean+2*svbiom.sd), color="black") +
	    geom_pointrange(aes(x=Site, y=model.mean, ymin=model.mean-2*model.sd, ymax=model.mean+2*model.sd, color=Model.Order)) +
	    scale_y_continuous(name="Fcomp SetVeg Dominant PFT") +
	    scale_x_discrete(name="Model") +
	    scale_color_manual(values=paste(model.colors$color), name="Model") + 
	    theme_bw() +
	    theme(axis.text.x=element_text(angle=45, hjust=1))
	  dev.off()
	  
	  # Getting the average bias for each model
	  setveg.mod <- aggregate(mod.svbiom[,c("svbiom.mean", "model.mean", "bias", "bias.per", "bias.abs")],
	                          by=mod.svbiom[,c("Model", "Model.Order")],
	                          FUN=mean)
	  
	  setveg.mod[,c("svbiom.sd", "model.sd", "bias.sd", "bias.per.sd", "bias.abs.sd")] <- aggregate(mod.svbiom[,c("svbiom.mean", "model.mean", "bias", "bias.per", "bias.abs")],
	                                                                                         by=mod.svbiom[,c("Model", "Model.Order")],
	                                                                                         FUN=sd)[,c("svbiom.mean", "model.mean", "bias", "bias.per", "bias.abs")]
	  setveg.mod
	  write.csv(setveg.mod, file.path(out.dir, "Biomass_SetVeg_ModelBenchmarks_Summary.csv"), row.names=F)
	  
	  # Getting the average bias for each model
	  setveg.site <- aggregate(mod.svbiom[,c("svbiom.mean", "model.mean", "bias", "bias.per", "bias.abs")],
	                           by=list(mod.svbiom$Site),
	                           FUN=mean)
	  names(setveg.site)[1] <- "Site"
	  
	  setveg.site[,c("svbiom.sd", "model.sd", "bias.sd", "bias.per.sd", "bias.abs.sd")] <- aggregate(mod.svbiom[,c("svbiom.mean", "model.mean", "bias", "bias.per", "bias.abs")],
	                                                                 by=list(mod.svbiom$Site),
	                                                                 FUN=sd)[,c("svbiom.mean", "model.mean", "bias", "bias.per", "bias.abs")]
	  setveg.site
	  write.csv(setveg.site, file.path(out.dir, "Biomass_SetVeg_SiteBenchmarks_Summary.csv"), row.names=F)
  }
  # ------------------
        
    
  # ------------------
  #    2.04. Biomass -- NBCD
  # ------------------
	{
	  # Read in the NBCD data
	  bm.nbcd <- read.csv(file.path(out.dir, "NBCD_sites_summary.csv"))
	  summary(bm.nbcd)
  
	  # NBCD combines FIA with Landsat (1999-2002), so lets pull Model AGB estimates from that time period
	  nbcd.bench <- aggregate(ecosys[ecosys$Year>=1999 & ecosys$Year<=2002, "AGB"], 
		  					  by=ecosys[ecosys$Year>=1999 & ecosys$Year<=2002,c("Model", "Model.Order", "Site")], 
			  				  FUN=mean)
	  names(nbcd.bench)[which(names(nbcd.bench)=="x")] <- "model.mean"
	  nbcd.bench[,"model.sd"] <- aggregate(ecosys[ecosys$Year>=1999 & ecosys$Year<=2002, "AGB"], 
				                      			   by=ecosys[ecosys$Year>=1999 & ecosys$Year<=2002,c("Model", "Model.Order", "Site")], 
				                      			   FUN=sd)[,"x"]
	  nbcd.bench[nbcd.bench$Model=="jules.stat",c("model.mean", "model.sd")] <- NA # Putting NAs in Jules because it doesn't do Biomass
	  summary(nbcd.bench)
  
	  # Merging NBCD with the model results
	  nbcd.bench <- merge(nbcd.bench, bm.nbcd)
  
	  # Calculating some statistics
	  # NOTE: Calculating these off of the max value in the grid cell to get the values for the most in-tact forest of the region
	  nbcd.bench1 <- nbcd.bench
	  nbcd.bench1$Ref <- as.factor("mean")
	  nbcd.bench1$bias <- nbcd.bench1$model.mean - nbcd.bench1$nbcd.mean
	  nbcd.bench1$bias.sd <- nbcd.bench1$bias/nbcd.bench1$nbcd.sd
	  nbcd.bench1$bias.per <- nbcd.bench1$bias/nbcd.bench1$nbcd.mean
	  nbcd.bench1$bias.abs <- abs(nbcd.bench1$bias)
	  
	  nbcd.bench$Ref <- as.factor("max")
	  nbcd.bench$bias <- nbcd.bench$model.mean - nbcd.bench$nbcd.max
	  nbcd.bench$bias.sd <- nbcd.bench$bias/nbcd.bench$nbcd.sd
	  nbcd.bench$bias.per <- nbcd.bench$bias/nbcd.bench$nbcd.max
	  nbcd.bench$bias.abs <- abs(nbcd.bench$bias)

	  nbcd.bench <- rbind(nbcd.bench, nbcd.bench1)

	  nbcd.bench$range.check <- as.factor(ifelse(nbcd.bench$model.mean>=nbcd.bench$nbcd.min & nbcd.bench$model.mean<=nbcd.bench$nbcd.max, "*", NA))
	  
	  summary(nbcd.bench)
	  
	  write.csv(nbcd.bench, file.path(out.dir, "Biomass_NBCD_ModelBenchmarks.csv"), row.names=F)
  
	  nbcd.bench[, "Range"] <- as.factor(ifelse(is.na(nbcd.bench$range.check), "out", "in"))
	  
	  # Order the sites from W to E
	  nbcd.bench$Site <- factor(nbcd.bench$Site, levels=c("PDL", "PBL", "PUN", "PMB", "PHA", "PHO"))
	  
	  png(file.path(fig.dir, "Biomass_NBCD_ModelBenchmarks_bySitel.png"), height=8, width=10, units="in", res=180)
	  ggplot(data=nbcd.bench[]) +
  		facet_wrap(~Site) +
  		geom_hline(aes(yintercept=nbcd.mean), color="black", size=1.5, linetype="solid") +
  		geom_hline(aes(yintercept=nbcd.mean-nbcd.sd), color="gray50", size=1, linetype="solid") +
  		geom_hline(aes(yintercept=nbcd.mean+nbcd.sd), color="gray50", size=1, linetype="solid") +
  		geom_hline(aes(yintercept=nbcd.max), color="red", size=1, linetype="dashed") +
  	  geom_hline(aes(yintercept=nbcd.min), color="red", size=1, linetype="dashed") +   
  		geom_pointrange(aes(x=Model.Order, y=model.mean, ymin=model.mean-model.sd, ymax=model.mean+model.sd, color=Model.Order), size=1.5) +
  		scale_y_continuous(name="AGB (kgC m-2)", limits=c(0,max(nbcd.bench$model.mean+nbcd.bench$model.sd, na.rm=T))) +
  		scale_x_discrete(name="Model") +
  		scale_color_manual(values=paste(model.colors$color), name="Model") + 
  		# scale_alpha_manual(values=c(1, 0.5), name="NBCD Range") +
  		theme_bw() +
	    theme(axis.text.x=element_text(angle=45, hjust=1))
	  dev.off()

	  png(file.path(fig.dir, "Biomass_NBCD_ModelBenchmarks_byModel.png"), height=8, width=10, units="in", res=180)
	  ggplot(data=nbcd.bench[]) +
	    facet_wrap(~Model.Order) +
	    # geom_hline(aes(yintercept=nbcd.max), color="black", size=2, linetype="dashed") +
	    # geom_hline(aes(yintercept=nbcd.mean), color="black", size=1.5, linetype="solid") +
	    # geom_hline(aes(yintercept=nbcd.mean-nbcd.sd), color="gray50", size=1, linetype="solid") +
	    # geom_hline(aes(yintercept=nbcd.mean+nbcd.sd), color="gray50", size=1, linetype="solid") +
	    geom_pointrange(aes(x=Site, y=nbcd.mean, ymin=nbcd.mean-nbcd.sd, ymax=nbcd.mean+nbcd.sd), size=1.25) +
	    geom_pointrange(aes(x=Site, y=model.mean, ymin=model.mean-model.sd, ymax=model.mean+model.sd, color=Model.Order), size=1.25) +
	    scale_y_continuous(name="AGB (kgC m-2)") +
	    scale_x_discrete(name="Model") +
	    scale_color_manual(values=paste(model.colors$color), name="Model", guid=guide_legend(ncol=2)) + 
	    theme_bw() +
	    theme(axis.text.x=element_text(angle=45, hjust=1)) +
	    theme(legend.position=c(0.75, 0.15))
	  dev.off()
	  
	  # Lookign at some quick summary stats by model
	  nbcd.performance <- aggregate(nbcd.bench[,c("nbcd.mean", "model.mean", "bias", "bias.per", "bias.abs")],
	                                by=nbcd.bench[,c("Model", "Model.Order", "Ref")],
	                                FUN=mean, na.rm=T)
	  nbcd.performance[,c("nbcd.sd", "model.sd", "bias.sd", "bias.per.sd", "bias.abs.sd")] <- aggregate(nbcd.bench[,c("nbcd.mean", "model.mean", "bias", "bias.per", "bias.abs")],
	                                                                                    by=nbcd.bench[,c("Model", "Model.Order", "Ref")],
	                                                                                    FUN=sd, na.rm=T)[,c("nbcd.mean", "model.mean", "bias", "bias.per", "bias.abs")]
	  nbcd.performance
	  # On average, JULES-TRIFFID had the the lowest bias, but CLM was the most consistent in it's bias (lowest SD.bias)
	  write.csv(nbcd.performance, file.path(out.dir, "Biomass_NBCD_ModelBenchmarks_Summary.csv"), row.names=F)
	}
  # ------------------
  
  # ------------------
  #    2.05. Carbon Fluxes - Flux Towers
  # ------------------
  { 
    # Read in the data
    flux <- read.csv(file.path(out.dir, "FluxData_Year.csv"))
    # names(flux) <- c("Site", "TowerName", "Year", "NEE.tower", "RE.tower", "GPP.tower")
    # flux$NEE.tower <- flux$NEE.tower*-1 # NEE sign from flux tower is backwards? is acutally NEP??
    summary(flux)
    
    # Stacking everything together to make the benchmarking easier
    flux.stack <- stack(flux[,c("NEE", "RE", "GPP")])
    names(flux.stack) <- c("tower.mean", "Flux")
    flux.stack[,c("Site", "TowerName", "Year")] <- flux[,c("SITE", "SiteName", "YEAR")]
    summary(flux.stack)
    
    summary(ecosys)
    ecosys$RE <- ecosys$AutoResp + ecosys$HeteroResp
    ecosys.stack <- stack(ecosys[ecosys$Year>=min(flux.stack$Year), c("GPP", "NEE", "RE")])
    names(ecosys.stack) <- c("model.mean", "Flux")
    ecosys.stack[,c("Model", "Model.Order", "Site", "Year")] <- ecosys[ecosys$Year>=min(flux.stack$Year), c("Model", "Model.Order", "Site", "Year")]
    
    flux.mod <- merge(flux.stack, ecosys.stack, all.x=T, all.y=F)
    flux.mod$model.mean <- flux.mod$model.mean*60*60*24*365.25 # Change to kgC m-2 s-1 to kgC m-2 yr-1
    summary(flux.mod)
    
    flux.mod$bias <- flux.mod$model.mean - flux.mod$tower.mean
    flux.mod$bias.per <- flux.mod$bias/flux.mod$tower.mean
    flux.mod$bias.abs <- abs(flux.mod$bias)
    summary(flux.mod)
    write.csv(flux.mod, file.path(out.dir, "Fluxes_FluxTowers_ModelBenchmarks.csv"), row.names=F)
    
    
    png(file.path(fig.dir, "NEE_FluxTowers_TimeSeries.png"), height=8, width=10, units="in", res=180)
    ggplot(data=flux.mod[flux.mod$Flux=="NEE",]) +
      facet_grid(Site~.) +
      geom_line(aes(x=Year, y=model.mean, color=Model.Order), size=0.75) +
      geom_line(aes(x=Year, y=tower.mean, linetype=TowerName), size=2) +
      geom_hline(yintercept=0, linetype="dashed", color="black") +
      scale_y_continuous(name="Flux kgC m-2 yr-1") +
      scale_x_continuous(name="Year") +
      scale_color_manual(values=paste(model.colors$color), name="Model") + 
      scale_linetype_manual(values=c("solid", "solid", "solid", "dotted")) +
      theme_bw() 
    dev.off()
    
    png(file.path(fig.dir, "RE_FluxTowers_TimeSeries.png"), height=8, width=10, units="in", res=180)
    ggplot(data=flux.mod[flux.mod$Flux=="RE",]) +
      facet_grid(Site~.) +
      geom_line(aes(x=Year, y=model.mean, color=Model.Order), size=0.75) +
      geom_line(aes(x=Year, y=tower.mean, linetype=TowerName), size=2) +
      scale_y_continuous(name="RE kgC m-2 yr-1") +
      scale_x_continuous(name="Year") +
      scale_color_manual(values=paste(model.colors$color), name="Model") + 
      scale_linetype_manual(values=c("solid", "solid", "solid", "dashed")) +
      theme_bw() 
    dev.off()
    
    png(file.path(fig.dir, "GPP_FluxTowers_TimeSeries.png"), height=8, width=10, units="in", res=180)
    ggplot(data=flux.mod[flux.mod$Flux=="GPP",]) +
      facet_grid(Site~.) +
      geom_line(aes(x=Year, y=model.mean, color=Model.Order), size=0.75) +
      geom_line(aes(x=Year, y=tower.mean, linetype=TowerName), size=2) +
      scale_y_continuous(name="GPP kgC m-2 yr-1") +
      scale_x_continuous(name="Year") +
      scale_color_manual(values=paste(model.colors$color), name="Model") + 
      scale_linetype_manual(values=c("solid", "solid", "solid", "dashed")) +
      theme_bw() 
    dev.off()
    
    # Reducing dimensions so that we can get some summary statistics
    # Get rid of multiple towers
    flux.mod2 <- aggregate(flux.mod[,c("tower.mean", "model.mean", "bias", "bias.per", "bias.abs")],
                           by=flux.mod[,c("Flux", "Site", "Model", "Model.Order", "Year")],
                           FUN=mean, na.rm=T)
    summary(flux.mod2)
    
    # Get rid of multiple years
    flux.mod3 <- aggregate(flux.mod2[,c("tower.mean", "model.mean", "bias", "bias.per", "bias.abs")],
                           by=flux.mod2[,c("Flux", "Site", "Model", "Model.Order")],
                           FUN=mean, na.rm=T)
    summary(flux.mod3)
    
    
    # Getting the average bias for each model
    flux.mod.mod <- aggregate(flux.mod3[,c("tower.mean", "model.mean","bias", "bias.per", "bias.abs")],
                              by=flux.mod3[,c("Flux", "Model", "Model.Order")],
                              FUN=mean, na.rm=T)
    flux.mod.mod[,c("tower.sd", "model.sd","bias.sd", "bias.per.sd", "bias.abs.sd")] <- aggregate(flux.mod3[,c("tower.mean", "model.mean","bias", "bias.per", "bias.abs")],
                                                                                           by=flux.mod3[,c("Flux", "Model", "Model.Order")],
                                                                                           FUN=sd, na.rm=T)[,c("tower.mean", "model.mean","bias", "bias.per", "bias.abs")]
    flux.mod.mod[flux.mod.mod$Flux=="NEE",]
    flux.mod.mod[flux.mod.mod$Flux=="GPP",]
    flux.mod.mod[flux.mod.mod$Flux=="RE",]
    write.csv(flux.mod.mod, file.path(out.dir, "Fluxes_FluxTowers_ModelBenchmarks_Summary.csv"), row.names=F)
    
    # Getting the average bias for each site
    flux.mod.site <- aggregate(flux.mod3[,c("tower.mean", "model.mean", "bias", "bias.abs")],
                               by=flux.mod3[,c("Flux", "Site")],
                               FUN=mean, na.rm=T)
    flux.mod.site[,c("tower.sd", "model.sd", "bias.sd", "bias.per.sd", "bias.abs.sd")] <- aggregate(flux.mod3[,c("tower.mean", "model.mean", "bias", "bias.per", "bias.abs")],
                                                                                             by=flux.mod3[,c("Flux", "Site")],
                                                                                             FUN=sd, na.rm=T)[,c("tower.mean", "model.mean", "bias", "bias.per", "bias.abs")]
    flux.mod.site[flux.mod.site$Flux=="NEE",]
    flux.mod.site[flux.mod.site$Flux=="GPP",]
    flux.mod.site[flux.mod.site$Flux=="RE",]
    write.csv(flux.mod.site, file.path(out.dir, "Fluxes_FluxTowers_SiteBenchmarks_Summary.csv"), row.names=F)
  }
  # ------------------
  
  # ------------------
  #    2.06. MODIS -- LAI, GPP
  # ------------------
  {
    modis <- read.csv(file.path(out.dir, "MODIS_year_sites_summary.csv"))
    names(modis) <- c("Site", "Var", "Year", "modis.mean", "modis.sd", "modis.min", "modis.max")
    summary(modis)
    
    summary(ecosys)
    ecosys.stack <- stack(ecosys[ecosys$Year>=min(modis$Year), c("GPP", "LAI")])
    names(ecosys.stack) <- c("model.mean", "Var")
    ecosys.stack[,c("Model", "Model.Order", "Site", "Year")] <- ecosys[ecosys$Year>=min(modis$Year), c("Model", "Model.Order", "Site", "Year")]
    summary(ecosys.stack)

    modis.mod <- merge(modis, ecosys.stack, all.x=T, all.y=F)
    modis.mod[modis.mod$Var=="GPP",c("model.mean", "modis.mean", "modis.sd", "modis.min", "modis.max")] <- modis.mod[modis.mod$Var=="GPP",c("model.mean", "modis.mean", "modis.sd", "modis.min", "modis.max")]*60*60*24*365.25 # Change GPP only to kgC m-2 s-1 to kgC m-2 yr-1
    summary(modis.mod)

    modis.mod$bias <- modis.mod$model.mean - modis.mod$modis.mean
    modis.mod$bias.per <- modis.mod$bias/modis.mod$modis.mean
    modis.mod$bias.abs <- abs(modis.mod$bias)
    summary(modis.mod)
    write.csv(modis.mod, file.path(out.dir, "MODIS_ModelBenchmarks.csv"), row.names=F)

    # Ordering sites W to E
    modis.mod$Site <- factor(modis.mod$Site, levels=c("PDL", "PBL", "PUN", "PMB", "PHA", "PHO"))

    png(file.path(fig.dir, "LAI_MODIS_TimeSeries.png"), height=8, width=10, units="in", res=180)
    ggplot(data=modis.mod[modis.mod$Var=="LAI",]) +
      facet_wrap(~Site) +
      geom_ribbon(aes(x=Year, ymin=modis.min, ymax=modis.max), alpha=0.5) +
      geom_ribbon(aes(x=Year, ymin=modis.mean-modis.sd, ymax=modis.mean+modis.sd), alpha=0.25) +
      geom_line(aes(x=Year, y=modis.mean), size=1) +
      geom_line(aes(x=Year, y=model.mean, color=Model.Order), size=0.5) +
      scale_y_continuous(name="LAI m2 m-2") +
      scale_x_continuous(name="Year", limits=c(2003, 2010), expand=c(0,0)) +
      scale_color_manual(values=paste(model.colors$color), name="Model") +
      theme_bw()
    dev.off()

    png(file.path(fig.dir, "GPP_MODIS_TimeSeries.png"), height=8, width=10, units="in", res=180)
    ggplot(data=modis.mod[modis.mod$Var=="GPP",]) +
      facet_wrap(~Site) +
      geom_ribbon(aes(x=Year, ymin=modis.min, ymax=modis.max), alpha=0.33) +
      geom_ribbon(aes(x=Year, ymin=modis.mean-modis.sd, ymax=modis.mean+modis.sd), alpha=0.5) +
      geom_line(aes(x=Year, y=modis.mean), size=1) +
      geom_line(aes(x=Year, y=model.mean, color=Model.Order), size=0.5) +
      scale_y_continuous(name="GPP kgC m-2 yr-1") +
      scale_x_continuous(name="Year", limits=c(2003, 2010), expand=c(0,0)) +
      scale_color_manual(values=paste(model.colors$color), name="Model") +
      theme_bw()
    dev.off()

    # Reducing dimensions so that we can get some summary statistics
    # Get rid of multiple years
    modis.mod3 <- aggregate(modis.mod[,c("modis.mean", "model.mean", "bias", "bias.per", "bias.abs")],
                           by=modis.mod[,c("Var", "Site", "Model", "Model.Order")],
                           FUN=mean, na.rm=T)
    summary(modis.mod3)


    # Getting the average bias for each model
    modis.mod.mod <- aggregate(modis.mod3[,c("modis.mean", "model.mean", "bias", "bias.per", "bias.abs")],
                               by=modis.mod3[,c("Var", "Model", "Model.Order")],
                               FUN=mean, na.rm=T)
    modis.mod.mod[,c("modis.sd", "model.sd", "bias.sd", "bias.per.sd", "bias.abs.sd")] <- aggregate(modis.mod3[,c("modis.mean", "model.mean", "bias", "bias.per", "bias.abs")],
                                                                                             by=modis.mod3[,c("Var", "Model", "Model.Order")],
                                                                                             FUN=sd, na.rm=T)[,c("modis.mean", "model.mean", "bias", "bias.per", "bias.abs")]
    modis.mod.mod[modis.mod.mod$Var=="LAI",]
    modis.mod.mod[modis.mod.mod$Var=="GPP",]
    write.csv(modis.mod.mod, file.path(out.dir, "MODIS_ModelBenchmarks_Summary.csv"), row.names=F)

    # Getting the average bias for each site
    modis.mod.site <- aggregate(modis.mod3[,c("modis.mean", "model.mean", "bias", "bias.per", "bias.abs")],
                                by=modis.mod3[,c("Var", "Site")],
                                FUN=mean, na.rm=T)
    modis.mod.site[,c("modis.sd", "model.sd", "bias.sd", "bias.per.sd", "bias.abs.sd")] <- aggregate(modis.mod3[,c("modis.mean", "model.mean", "bias", "bias.per", "bias.abs")],
                                                                     by=modis.mod3[,c("Var", "Site")],
                                                                     FUN=sd, na.rm=T)[,c("modis.mean", "model.mean", "bias", "bias.per", "bias.abs")]
    modis.mod.site[modis.mod.site$Var=="LAI",]
    modis.mod.site[modis.mod.site$Var=="GPP",]
    write.csv(modis, file.path(out.dir, "MODIS_ModelBenchmarks_Site.csv"), row.names=F)
  }
  # ------------------

  
  # ------------------
  #    2.07. Tree Rings - AGB, AGBi (AGB increment)
  # ------------------
  {
    
  }
  # ------------------
  
}
# -------------------------------------------
