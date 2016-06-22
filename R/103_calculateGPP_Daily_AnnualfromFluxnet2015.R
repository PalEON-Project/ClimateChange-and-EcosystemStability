#read data from fluxnet15 dataset

# Location of data on Daves machine
#      D:\GlobalDatasets\FLUXNET06202016

# FLX_US-Ha1_FLUXNET2015_FULLSET_1991-2012_1-1/

folderpath="D:/GlobalDatasets/FLUXNET06202016/"


# read all the *.csv files in the working directory
#change to reflect file path
tempFilelist = list.files(filepath, pattern="*.csv")
directorylist = list.dirs(folderpath)
substr(tempFilelist,5,10)
substr(tempFilelist,24,38)


filepathHa1="D:/GlobalDatasets/FLUXNET06202016/FLX_US-Ha1_FLUXNET2015_FULLSET_1991-2012_1-1/"
#Harvard
Ha1_annual=read.csv(paste0(filepathHa1,"FLX_US-Ha1_FLUXNET2015_FULLSET_YY_1991-2012_1-1.csv"))
Ha1_Monthly=read.csv(paste0(filepathHa1,"./FLX_US-Ha1_FLUXNET2015_FULLSET_MM_1991-2012_1-1.csv"))
Ha1_daily=read.csv(paste0(filepathHa1,"./FLX_US-Ha1_FLUXNET2015_FULLSET_DD_1991-2012_1-1.csv"))


#Morgan Monroe
filepathMMS="D:/GlobalDatasets/FLUXNET06202016//FLX_US-MMS_FLUXNET2015_FULLSET_1999-2014_1-1"

MMS_annual=read.csv(paste0(filepathMMS,"./FLX_US-MMS_FLUXNET2015_FULLSET_YY_1999-2014_1-1.csv"))
MMS_Monthly=read.csv(paste0(filepathMMS,"./FLX_US-MMS_FLUXNET2015_FULLSET_MM_1999-2014_1-1.csv"))
MMS_daily=read.csv(paste0(filepathMMS,"./FLX_US-MMS_FLUXNET2015_FULLSET_DD_1999-2014_1-1.csv"))

MMS_Meto = read.csv(paste0(filepathMMS,"./FLX_US-MMS_FLUXNET2015_AUXMETEO_1999-2014_1-1.csv"))


#Sylvania
filepathSyv="D:/GlobalDatasets/FLUXNET06202016//FLX_US-Syv_FLUXNET2015_FULLSET_2001-2014_1-1"

Syv_annual=read.csv(paste0(filepathSyv,"./FLX_US-Syv_FLUXNET2015_FULLSET_YY_2001-2014_1-1.csv"))
Syv_Monthly=read.csv(paste0(filepathSyv,"./FLX_US-Syv_FLUXNET2015_FULLSET_MM_2001-2014_1-1.csv"))
Syv_daily=read.csv(paste0(filepathSyv,"./FLX_US-Syv_FLUXNET2015_FULLSET_DD_2001-2014_1-1.csv"))

Syv_Meto = read.csv(paste0(filepathSyv,"./FLX_US-Syv_FLUXNET2015_AUXMETEO_2001-2014_1-1.csv"))



#Willow Creek
filepathWCr="D:/GlobalDatasets/FLUXNET06202016//FLX_US-WCr_FLUXNET2015_FULLSET_1999-2014_1-1"

#FLX_US-WCr_FLUXNET2015_FULLSET_1999-2014_1-1

WCr_annual=read.csv(paste0(filepathWCr,"./FLX_US-WCr_FLUXNET2015_FULLSET_YY_1999-2014_1-1.csv"))
WCr_Monthly=read.csv(paste0(filepathWCr,".//FLX_US-WCr_FLUXNET2015_FULLSET_MM_1999-2014_1-1.csv"))
WCr_daily=read.csv(paste0(filepathWCr,".//FLX_US-WCr_FLUXNET2015_FULLSET_DD_1999-2014_1-1.csv"))

WCr_Meto = read.csv(paste0(filepathWCr,".//FLX_US-WCr_FLUXNET2015_AUXMETEO_1999-2014_1-1.csv"))




# Header for fluxnet 2015 FULL dataset
# combine with any dataset for Units, description etc
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
  
  
  