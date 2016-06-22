#read data from fluxnet15 dataset

# Location of data on Daves machine
#      D:\GlobalDatasets\FLUXNET06202016

# FLX_US-Ha1_FLUXNET2015_FULLSET_1991-2012_1-1/

folderpath="D:/GlobalDatasets/FLUXNET06202016/"

filepath="D:/GlobalDatasets/FLUXNET06202016/FLX_US-Ha1_FLUXNET2015_FULLSET_1991-2012_1-1/"

# read all the *.csv files in the working directory
#change to reflect file path
tempFilelist = list.files(filepath, pattern="*.csv")
directorylist = list.dirs(folderpath)
substr(tempFilelist,5,10)
substr(tempFilelist,24,38)

setwd(filepath)

tempFilelist

Ha1_annual=read.csv(paste0(filepath,"FLX_US-Ha1_FLUXNET2015_FULLSET_YY_1991-2012_1-1.csv"))
Ha1_Monthly=read.csv(paste0(filepath,"./FLX_US-Ha1_FLUXNET2015_FULLSET_MM_1991-2012_1-1.csv"))
Ha1_daily=read.csv(paste0(filepath,"./FLX_US-Ha1_FLUXNET2015_FULLSET_DD_1991-2012_1-1.csv"))

# Header for fluxnet 2015 FULL dataset
FluxnetFULLheader = read.csv("D:/GlobalDatasets/FLUXNET06202016/FLUXNET2015_FULLSET_variable_list.csv")

  
  #Night time method mean, 5% and 95% vs reference
  
  plot(Ha1_Monthly$GPP_NT_VUT_REF, Ha1_Monthly$GPP_NT_VUT_MEAN, xlim=c(-5,20),xlab="GPP_NT_VUT_REF", ylim=c(-5,20),ylab="Monthly GPP_NT_VUT_Mean+_95%", col=1) 
  par(new=T)
  plot(Ha1_Monthly$GPP_NT_VUT_REF, Ha1_Monthly$GPP_NT_VUT_95, xlim=c(-5,20),xlab="GPP_NT_VUT_REF", ylim=c(-5,20),ylab="Monthly GPP_NT_VUT_Mean+_95%", col=34) 
  par(new=T)
  plot(Ha1_Monthly$GPP_NT_VUT_REF, Ha1_Monthly$GPP_NT_VUT_05, xlim=c(-5,20),xlab="GPP_NT_VUT_REF", ylim=c(-5,20),ylab="Monthly GPP_NT_VUT_Mean+_95%", col=30) 
  
  
  
  plot(Ha1_annual$GPP_NT_VUT_REF, Ha1_annual$GPP_NT_VUT_MEAN, xlim=c(0,2500),xlab="GPP_NT_VUT_REF", ylim=c(0,2500),ylab="Annual GPP_NT_VUT_Mean+_95%", col=1) 
  par(new=T)
  plot(Ha1_annual$GPP_NT_VUT_REF, Ha1_annual$GPP_NT_VUT_95, xlim=c(0,2500),xlab="GPP_NT_VUT_REF", ylim=c(0,2500),ylab="Annual GPP_NT_VUT_Mean+_95%", col=34) 
  par(new=T)
  plot(Ha1_annual$GPP_NT_VUT_REF, Ha1_annual$GPP_NT_VUT_05, xlim=c(0,2500),xlab="GPP_NT_VUT_REF", ylim=c(0,2500),ylab="Annual GPP_NT_VUT_Mean+_95%", col=30) 
  
  #Night time method mean, 5% and 95% vs DAYTIME reference
  
  plot(Ha1_Monthly$GPP_DT_VUT_REF, Ha1_Monthly$GPP_NT_VUT_MEAN, xlim=c(-5,20),xlab="GPP_DT_VUT_REF", ylim=c(-5,20),ylab="Monthly GPP_NT_VUT_Mean+_95%", col=1) 
  par(new=T)
  plot(Ha1_Monthly$GPP_DT_VUT_REF, Ha1_Monthly$GPP_NT_VUT_95, xlim=c(-5,20),xlab="GPP_DT_VUT_REF", ylim=c(-5,20),ylab="Monthly GPP_NT_VUT_Mean+_95%", col=34) 
  par(new=T)
  plot(Ha1_Monthly$GPP_DT_VUT_REF, Ha1_Monthly$GPP_NT_VUT_05, xlim=c(-5,20),xlab="GPP_DT_VUT_REF", ylim=c(-5,20),ylab="Monthly GPP_NT_VUT_Mean+_95%", col=30) 
  
  
  plot(Ha1_annual$GPP_DT_VUT_REF, Ha1_annual$GPP_NT_VUT_MEAN, xlim=c(0,2500),xlab="GPP_NT_VUT_REF", ylim=c(0,2500),ylab="Annual GPP_NT_VUT_Mean+_95%", col=1) 
  par(new=T)
  plot(Ha1_annual$GPP_DT_VUT_REF, Ha1_annual$GPP_NT_VUT_95, xlim=c(0,2500),xlab="GPP_NT_VUT_REF", ylim=c(0,2500),ylab="Annual GPP_NT_VUT_Mean+_95%", col=34) 
  par(new=T)
  plot(Ha1_annual$GPP_DT_VUT_REF, Ha1_annual$GPP_NT_VUT_05, xlim=c(0,2500),xlab="GPP_NT_VUT_REF", ylim=c(0,2500),ylab="Annual GPP_NT_VUT_Mean+_95%", col=30) 
  
  
  #extend for Howland and Sylvania
  
  
getwd()



setwd("D:/Dropbox/rProjectsShare/MIP-Change-and-Stability")