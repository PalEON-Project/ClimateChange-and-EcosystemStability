# -------------------------------------------
# Assessing the scales and magnitude of change in ecosystem state
# Author: Christy Rollinson, crollinson@gmail.com

# 1. Stastical detection of significant change
#    - Use TP regression splines
# -------------------------------------------
rm(list=ls())


# -------------------------------------------
# Load libraries; set file paths
# -------------------------------------------
library(car)
library(ggplot2)
library(mgcv)
library(parallel)
setwd("~/Dropbox/PalEON_CR/PalEON_MIP_Site/Analyses/Change-and-Stability") # Path to this project github repository: https://github.com/PalEON-Project/MIP-Change-and-Stability.git
path.gamm.func <- "~/Desktop/R_Functions/"  # Path to github repository of my GAMM helper functions: https://github.com/crollinson/R_Functions.git
inputs    <- "~/Dropbox/PalEON_CR/PalEON_MIP_Site/phase1a_output_variables/" # Path to my cleaned model output

mip.utils <- "~/Dropbox/PalEON_CR/MIP_Utils/" # Path to PalEON MIP Utility repository: https://github.com/PalEON-Project/MIP_Utils.git
path.raw <- "~/Dropbox/PalEON_CR/PalEON_MIP_Site/phase1a_model_output/" # Path to raw model output

out.dir <- "Data/EcosysChange" # Path to where the analysis output should go
fig.dir <- "Figures/EcosysChange" # Path to where figures should go

if(!dir.exists(out.dir)) dir.create(out.dir)
if(!dir.exists(fig.dir)) dir.create(fig.dir)
# -------------------------------------------


# -------------------------------------------
# Reading in raw ecosystem model output & Adding the PFT info
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


model.names <- data.frame(Model=unique(ecosys$Model), Model.Order=unique(ecosys$Model.Order))

# ----------------------
# Extracting the dominant PFT at each site for each model
# ----------------------
source(file.path(mip.utils, "Phase1_sites/extract_output_site.R"))

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
sites.all <- c("PHA", "PHO", "PUN", "PBL", "PDL", "PMB")
for(m in unique(ecosys$Model)){
  print(paste0(m))
  if(m=="sibcasa") next
  model.name <- model.names[model.names$Model==m, "Model.Name"]
  pft.tmp <- extract.paleon.site(model=model.name, 
                                 model.dir=file.path(path.raw, paste(model.name, model.names[model.names$Model==m, "Version"], sep=".")), 
                                 sites=sites.all, 
                                 vars="Fcomp")
  if(m %in% c("jules.stat")){
    pft.tmp <- extract.paleon.site(model=model.name, 
                                   model.dir=file.path(path.raw, paste(model.name, model.names[model.names$Model==m, "Version"], sep=".")), 
                                   sites=sites.all, 
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
  dimnames(pft.tmp$Fcomp)[[3]] <- sites.all
  dimnames(pft.tmp$Fcomp)[[1]] <- 849+1:dim(pft.tmp$Fcomp)[1]
  
  for(i in 1:dim(pft.tmp$Fcomp)[3]){
    col.max <- which(apply(pft.tmp$Fcomp[,,i], 2, mean)==max(apply(pft.tmp$Fcomp[,,i], 2, mean)))
    pft.max <- data.frame(Model=m,
                          Site=as.factor(dimnames(pft.tmp$Fcomp)[[3]][i]),
                          Year=as.numeric(dimnames(pft.tmp$Fcomp)[[1]]),
                          PFTmax=as.factor(col.max),
                          Fcomp=pft.tmp$Fcomp[,col.max,i])
    
    if(i==1 & m==unique(ecosys$Model)[1]){
      fcomp <- pft.max
    } else {
      fcomp <- rbind(fcomp, pft.max)
    }
  }
}
ecosys <- merge(ecosys, fcomp, all.x=T, all.y=T)

summary(ecosys)
ggplot(data=ecosys) +
  facet_wrap(~Model) +
  geom_line(aes(x=Year, y=Fcomp, color=Site, linetype=PFTmax)) +
  theme_bw()
# ----------------------

vars <- c("GPP", "NEE", "LAI", "AGB", "SoilCarb", "Fcomp") # Add dominant PFT

# -------------------------------------------


# -------------------------------------------
# Run the gams & all the post-processing -- Full Temporal Extent
# -------------------------------------------
source("R/0_TimeAnalysis.R")
# Just testing on 1 site for now

for(v in vars){
  print(paste0("Processing Var: ", v, " (0850-2010)"))
  
  # Package the data
  dat.list <- list()
  for(m in unique(ecosys[!is.na(ecosys[,v]),"Model"])){
    dat.list[[m]] <- ecosys[ecosys$Model==m, c("Model", "Site", "Year", v)]
  }
  
  # Run the stats
  cores.use <- min(12, length(dat.list))
  dat.out <- mclapply(dat.list, analyze.time, mc.cores=cores.use, Y=v, fac.fit="Site", k.freq=25, path.gamm.func=path.gamm.func)
  
  # Format the output
  for(m in names(dat.out)){
    dat.out[[m]]$out$Model <- as.factor(m)
    dat.out[[m]]$out$var   <- as.factor(v)
    if(m == names(dat.out)[1]){
      dat.out2 <- dat.out[[m]]$out
    } else{
      dat.out2 <- rbind(dat.out2, dat.out[[m]]$out)
    }
  } # End model formatting

  ## Some additional formatting
  dat.out2$mean.sig <- ifelse(dat.out2$sig=="*", dat.out2$mean, NA)
  dat.out2 <- merge(dat.out2, model.names, all.x=T, all.y=F)
  dat.out2$Model.Order <- factor(dat.out2$Model.Order, levels=model.colors[model.colors$Model.Order %in% unique(dat.out2$Model.Order),"Model.Order"])  
  
  # Save the output
  write.csv(dat.out2, file.path(out.dir, paste0("StabilityCalcs_", v, "_850-2010.csv")), row.names=F)
  save(dat.out, file=file.path(out.dir, paste0("StabilityCalcs_", v, "_850-2010.RData")))
  
  # Make and save some figures
  col.model <- paste(model.colors[model.colors$Model.Order %in% unique(dat.out2$Model.Order),"color"])
  
  png(file.path(fig.dir, paste0("Stability_", v, "_850-2010.png")), height=8.5, width=11, units="in", res=180)
  print(
  ggplot(data=dat.out2[,]) + 
    facet_wrap(~Site) +
    geom_line(aes(x=Year, y=Y, color=Model.Order), size=0.5, alpha=0.3) +
    geom_ribbon(aes(x=Year, ymin=lwr, ymax=upr, fill=Model.Order), alpha=0.3) +
    geom_line(aes(x=Year, y=mean, color=Model.Order), size=1, alpha=0.2) +
    geom_line(aes(x=Year, y=mean.sig, color=Model.Order), size=2, alpha=1) +
    geom_vline(xintercept=1850, linetype="dashed") +
    geom_vline(xintercept=1900, linetype="dashed") +
    scale_x_continuous(expand=c(0,0), name="Year") +
    scale_y_continuous(expand=c(0,0), name=v) +
    scale_color_manual(values=col.model) +
    scale_fill_manual(values=col.model) + 
    ggtitle(v) +
    theme_bw()
  )
  dev.off()
  
} # End variable loop
# -------------------------------------------


# -------------------------------------------
# Run the gams & all the post-processing -- pre-1850 plus 1900-2010
# -------------------------------------------
source("R/0_TimeAnalysis.R")
# Just testing on 1 site for now

for(v in vars){
  print(paste0("Processing Var: ", v, " (0850-1850, 1900-2010)"))
  
  # Package the data
  dat.list <- list()
  dat.list2 <- list()
  for(m in unique(ecosys[!is.na(ecosys[,v]),"Model"])){
    dat.list[[m]] <- ecosys[ecosys$Model==m & ecosys$Year<1850, c("Model", "Site", "Year", v)]
    dat.list2[[m]] <- ecosys[ecosys$Model==m & ecosys$Year>1900, c("Model", "Site", "Year", v)]
  }
  
  # Run the stats
  cores.use <- min(12, length(dat.list))
  dat.out  <- mclapply(dat.list, analyze.time, mc.cores=cores.use, Y=v, fac.fit="Site", k.freq=25, path.gamm.func=path.gamm.func)
  dat.out2 <- mclapply(dat.list2, analyze.time, mc.cores=cores.use, Y=v, fac.fit="Site", k.freq=25, path.gamm.func=path.gamm.func)

  # Format the output
  for(m in names(dat.out)){
    dat.out[[m]]$out$Model <- as.factor(m)
    dat.out[[m]]$out$var   <- as.factor(v)
    dat.out2[[m]]$out$Model <- as.factor(m)
    dat.out2[[m]]$out$var   <- as.factor(v)
    if(m == names(dat.out)[1]){
      dat.out3 <- rbind(dat.out[[m]]$out, dat.out2[[m]]$out)
    } else{
      dat.out3 <- rbind(dat.out3, rbind(dat.out[[m]]$out, dat.out2[[m]]$out))
    }
  } # End model formatting
  
  # Adding the non-modeled data
  ecosys2 <- data.frame(Model=ecosys$Model, Site=ecosys$Site, Year=ecosys$Year, Y=ecosys[,v], var=v)
  ecosys2 <- ecosys2[!is.na(ecosys2$Y),]
  dat.out3 <- merge(dat.out3, ecosys2, all.x=T, all.y=T)
  summary(dat.out3)  

  ## Some additional formatting
  dat.out3$mean.sig <- ifelse(dat.out3$sig=="*", dat.out3$mean, NA)
  dat.out3 <- merge(dat.out3, model.names, all.x=T, all.y=F)
  dat.out3$Model.Order <- factor(dat.out3$Model.Order, levels=model.colors[model.colors$Model.Order %in% unique(dat.out3$Model.Order),"Model.Order"])  
  
  # Save the output
  write.csv(dat.out3, file.path(out.dir, paste0("StabilityCalcs_", v, "_850-1850_1900-2010.csv")), row.names=F)
  save(dat.out, dat.out2, file=file.path(out.dir, paste0("StabilityCalcs_", v, "_850-1850_1900-2010.RData")))
  
  # Make and save some figures
  col.model <- paste(model.colors[model.colors$Model.Order %in% unique(dat.out3$Model.Order),"color"])
  
  png(file.path(fig.dir, paste0("Stability_", v, "_850-1850_1900-2010.png")), height=8.5, width=11, units="in", res=180)
  print(
    ggplot(data=dat.out3[,]) + 
      facet_wrap(~Site) +
      geom_line(aes(x=Year, y=Y, color=Model.Order), size=0.5, alpha=0.3) +
      geom_ribbon(data=dat.out3[dat.out3$Year<1850,], aes(x=Year, ymin=lwr, ymax=upr, fill=Model.Order), alpha=0.3) +
      geom_ribbon(data=dat.out3[dat.out3$Year>1900,], aes(x=Year, ymin=lwr, ymax=upr, fill=Model.Order), alpha=0.3) +
      geom_line(data=dat.out3[dat.out3$Year<1850,], aes(x=Year, y=mean, color=Model.Order), size=1, alpha=0.2) +
      geom_line(data=dat.out3[dat.out3$Year>1900,], aes(x=Year, y=mean, color=Model.Order), size=2, alpha=1) +
      geom_line(aes(x=Year, y=mean.sig, color=Model.Order), size=2, alpha=1) +
      geom_vline(xintercept=1850, linetype="dashed") +
      geom_vline(xintercept=1900, linetype="dashed") +
      scale_x_continuous(expand=c(0,0), name="Year") +
      scale_y_continuous(expand=c(0,0), name=v) +
      scale_color_manual(values=col.model) +
      scale_fill_manual(values=col.model) + 
      ggtitle(v) +
      theme_bw()
  )
  dev.off()
  
} # End variable loop
# -------------------------------------------
