#Author Dave Moore davidjpmoore@email.arizona.edu
#June 2016
#Purpose: Read data from fluxnet15 dataset


# Location of data on Daves machine
#      D:\GlobalDatasets\FLUXNET06202016

# FLX_US-Ha1_FLUXNET2015_FULLSET_1991-2012_1-1/


# Load required libraries
library(dplyr)
library(ggplot2)
library(tidyr)

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

Ha1_Meteo = read.csv(paste0(filepathMMS,"./FLX_US-Ha1_FLUXNET2015_AUXMETEO_1991-2012_1-1.csv"))

#Morgan Monroe
filepathMMS="D:/GlobalDatasets/FLUXNET06202016//FLX_US-MMS_FLUXNET2015_FULLSET_1999-2014_1-1"

MMS_annual=read.csv(paste0(filepathMMS,"./FLX_US-MMS_FLUXNET2015_FULLSET_YY_1999-2014_1-1.csv"))
MMS_Monthly=read.csv(paste0(filepathMMS,"./FLX_US-MMS_FLUXNET2015_FULLSET_MM_1999-2014_1-1.csv"))
MMS_daily=read.csv(paste0(filepathMMS,"./FLX_US-MMS_FLUXNET2015_FULLSET_DD_1999-2014_1-1.csv"))

MMS_Meteo = read.csv(paste0(filepathMMS,"./FLX_US-MMS_FLUXNET2015_AUXMETEO_1999-2014_1-1.csv"))


#Sylvania
filepathSyv="D:/GlobalDatasets/FLUXNET06202016//FLX_US-Syv_FLUXNET2015_FULLSET_2001-2014_1-1"

Syv_annual=read.csv(paste0(filepathSyv,"./FLX_US-Syv_FLUXNET2015_FULLSET_YY_2001-2014_1-1.csv"))
Syv_Monthly=read.csv(paste0(filepathSyv,"./FLX_US-Syv_FLUXNET2015_FULLSET_MM_2001-2014_1-1.csv"))
Syv_daily=read.csv(paste0(filepathSyv,"./FLX_US-Syv_FLUXNET2015_FULLSET_DD_2001-2014_1-1.csv"))

Syv_Meteo = read.csv(paste0(filepathSyv,"./FLX_US-Syv_FLUXNET2015_AUXMETEO_2001-2014_1-1.csv"))



#Willow Creek
filepathWCr="D:/GlobalDatasets/FLUXNET06202016//FLX_US-WCr_FLUXNET2015_FULLSET_1999-2014_1-1"

#FLX_US-WCr_FLUXNET2015_FULLSET_1999-2014_1-1

WCr_annual=read.csv(paste0(filepathWCr,"./FLX_US-WCr_FLUXNET2015_FULLSET_YY_1999-2014_1-1.csv"))
WCr_Monthly=read.csv(paste0(filepathWCr,".//FLX_US-WCr_FLUXNET2015_FULLSET_MM_1999-2014_1-1.csv"))
WCr_daily=read.csv(paste0(filepathWCr,".//FLX_US-WCr_FLUXNET2015_FULLSET_DD_1999-2014_1-1.csv"))

WCr_Meteo = read.csv(paste0(filepathWCr,".//FLX_US-WCr_FLUXNET2015_AUXMETEO_1999-2014_1-1.csv"))


# Header for fluxnet 2015 FULL dataset
# combine with any dataset for Units, description etc
FluxnetFULLheader = read.csv("D:/GlobalDatasets/FLUXNET06202016/FLUXNET2015_FULLSET_variable_list.csv")



#Annual values
Harvard_AnnualGPP = Ha1_annual %>%
  mutate(SiteName="US-Ha1") %>%
  select(SiteName, TIMESTAMP, GPP_DT_VUT_MEAN, GPP_DT_VUT_05, GPP_DT_VUT_95, GPP_DT_VUT_SE, GPP_DT_CUT_MEAN, GPP_DT_CUT_05, GPP_DT_CUT_95, GPP_DT_CUT_SE, GPP_NT_VUT_MEAN, GPP_NT_VUT_05, GPP_NT_VUT_95, GPP_NT_VUT_SE, GPP_NT_CUT_MEAN, GPP_NT_CUT_05, GPP_NT_CUT_95, GPP_NT_CUT_SE)

MorganMonroe_AnnualGPP = MMS_annual %>%
  mutate(SiteName="US-MMS") %>%
  select(SiteName, TIMESTAMP, GPP_DT_VUT_MEAN, GPP_DT_VUT_05, GPP_DT_VUT_95, GPP_DT_VUT_SE, GPP_DT_CUT_MEAN, GPP_DT_CUT_05, GPP_DT_CUT_95, GPP_DT_CUT_SE, GPP_NT_VUT_MEAN, GPP_NT_VUT_05, GPP_NT_VUT_95, GPP_NT_VUT_SE, GPP_NT_CUT_MEAN, GPP_NT_CUT_05, GPP_NT_CUT_95, GPP_NT_CUT_SE)

WillowCreek_AnnualGPP = WCr_annual %>%
  mutate(SiteName="US-WCr") %>%
  select(SiteName, TIMESTAMP, GPP_DT_VUT_MEAN, GPP_DT_VUT_05, GPP_DT_VUT_95, GPP_DT_VUT_SE, GPP_DT_CUT_MEAN, GPP_DT_CUT_05, GPP_DT_CUT_95, GPP_DT_CUT_SE, GPP_NT_VUT_MEAN, GPP_NT_VUT_05, GPP_NT_VUT_95, GPP_NT_VUT_SE, GPP_NT_CUT_MEAN, GPP_NT_CUT_05, GPP_NT_CUT_95, GPP_NT_CUT_SE)

Sylvania_AnnualGPP = Syv_annual %>%
  mutate(SiteName="US-Syv") %>%
  select(SiteName, TIMESTAMP, GPP_DT_VUT_MEAN, GPP_DT_VUT_05, GPP_DT_VUT_95, GPP_DT_VUT_SE, GPP_DT_CUT_MEAN, GPP_DT_CUT_05, GPP_DT_CUT_95, GPP_DT_CUT_SE, GPP_NT_VUT_MEAN, GPP_NT_VUT_05, GPP_NT_VUT_95, GPP_NT_VUT_SE, GPP_NT_CUT_MEAN, GPP_NT_CUT_05, GPP_NT_CUT_95, GPP_NT_CUT_SE)

 Tempjoin= full_join(Harvard_AnnualGPP, MorganMonroe_AnnualGPP)
 Tempjoin= full_join(Tempjoin, WillowCreek_AnnualGPP)
AnnualGPP_MIPBenchmark = full_join(Tempjoin, Sylvania_AnnualGPP)

save(AnnualGPP_MIPBenchmark, file = "./data/fluxdata/AnnualGPP_MIPBenchmark.FLUXNET.22Jun2016.RData")





#Monthly Values

Harvard_MonthlyGPP = Ha1_Monthly %>%
  mutate(SiteName="US-Ha1") %>%
  select(SiteName, TIMESTAMP, GPP_DT_VUT_MEAN, GPP_DT_VUT_05, GPP_DT_VUT_95, GPP_DT_VUT_SE, GPP_DT_CUT_MEAN, GPP_DT_CUT_05, GPP_DT_CUT_95, GPP_DT_CUT_SE, GPP_NT_VUT_MEAN, GPP_NT_VUT_05, GPP_NT_VUT_95, GPP_NT_VUT_SE, GPP_NT_CUT_MEAN, GPP_NT_CUT_05, GPP_NT_CUT_95, GPP_NT_CUT_SE)

MorganMonroe_MonthlyGPP = MMS_Monthly %>%
  mutate(SiteName="US-MMS") %>%
  select(SiteName, TIMESTAMP, GPP_DT_VUT_MEAN, GPP_DT_VUT_05, GPP_DT_VUT_95, GPP_DT_VUT_SE, GPP_DT_CUT_MEAN, GPP_DT_CUT_05, GPP_DT_CUT_95, GPP_DT_CUT_SE, GPP_NT_VUT_MEAN, GPP_NT_VUT_05, GPP_NT_VUT_95, GPP_NT_VUT_SE, GPP_NT_CUT_MEAN, GPP_NT_CUT_05, GPP_NT_CUT_95, GPP_NT_CUT_SE)

WillowCreek_MonthlyGPP = WCr_Monthly %>%
  mutate(SiteName="US-WCr") %>%
  select(SiteName, TIMESTAMP, GPP_DT_VUT_MEAN, GPP_DT_VUT_05, GPP_DT_VUT_95, GPP_DT_VUT_SE, GPP_DT_CUT_MEAN, GPP_DT_CUT_05, GPP_DT_CUT_95, GPP_DT_CUT_SE, GPP_NT_VUT_MEAN, GPP_NT_VUT_05, GPP_NT_VUT_95, GPP_NT_VUT_SE, GPP_NT_CUT_MEAN, GPP_NT_CUT_05, GPP_NT_CUT_95, GPP_NT_CUT_SE)

Sylvania_MonthlyGPP = Syv_Monthly %>%
  mutate(SiteName="US-Syv") %>%
  select(SiteName, TIMESTAMP, GPP_DT_VUT_MEAN, GPP_DT_VUT_05, GPP_DT_VUT_95, GPP_DT_VUT_SE, GPP_DT_CUT_MEAN, GPP_DT_CUT_05, GPP_DT_CUT_95, GPP_DT_CUT_SE, GPP_NT_VUT_MEAN, GPP_NT_VUT_05, GPP_NT_VUT_95, GPP_NT_VUT_SE, GPP_NT_CUT_MEAN, GPP_NT_CUT_05, GPP_NT_CUT_95, GPP_NT_CUT_SE)

Tempjoin= full_join(Harvard_MonthlyGPP, MorganMonroe_MonthlyGPP)
Tempjoin= full_join(Tempjoin, WillowCreek_MonthlyGPP)
MonthlyGPP_MIPBenchmark = full_join(Tempjoin, Sylvania_MonthlyGPP)

save(MonthlyGPP_MIPBenchmark, file = "./data/fluxdata/MonthlyGPP_MIPBenchmark.FLUXNET.22Jun2016.RData")

###Diagnostic Plots
  
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
  
  
  