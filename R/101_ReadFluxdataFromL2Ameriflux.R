# This Script Reads L2 Ameriflux files from Harvard, Howland, Sylvania and Willow Creek 
# Author: Dave Moore
# Contact: davidjpmoore@email.arizona.edu


# modify working directory and "filepath"

setwd("D:/GlobalDatasets/ameriflux.allsites.L2_data.05Mar2016/")
#note this was updated by Dave Moore - downloaded Gap filled data on June 17th
#link to UA desktop machine

# 
# Howland
# 
filepath="Howland_Forest_Main/gap_filled/"
# read all the *.csv files in the working directory
#change to reflect file path
tempFilelist = list.files(filepath, pattern="*.csv")
# get the AmeriFlux files information (site and years)
files<-substr(tempFilelist,7,14)
# get the AmeriFlux site code
site<-substr(tempFilelist[1],7,9)
# get the data and header from the L2 files
HowlandHo1 = do.call("rbind", lapply(paste0(filepath,tempFilelist), function(x) read.csv(x, skip=20,header=FALSE,stringsAsFactors = FALSE)))
AmFluxheader=read.csv(paste0(filepath,tempFilelist[1]),skip=17, strip.white=TRUE, nrows=1 ,header=FALSE,stringsAsFactors=FALSE)
colnames(HowlandHo1)<-AmFluxheader

# 
# Harvard
# 
#data up to date in March 2016 - still same version on June 17th 
# read all the *.csv files in the working directory
filepath="Harvard_Forest/gap_filled/"
#change to reflect file path
tempFilelist = list.files(filepath, pattern="*.csv")
# get the AmeriFlux files information (site and years)
files<-substr(tempFilelist,7,14)
# get the AmeriFlux site code
site<-substr(tempFilelist[1],7,9)
# get the data and header from the L2 files
HarvardHa1 = do.call("rbind", lapply(paste0(filepath,tempFilelist), function(x) read.csv(x, skip=20,header=FALSE,stringsAsFactors = FALSE)))
AmFluxheader=read.csv(paste0(filepath,tempFilelist[1]),skip=17, strip.white=TRUE, nrows=1 ,header=FALSE,stringsAsFactors=FALSE)
colnames(HarvardHa1)<-AmFluxheader


#Hemlock Tower
#Gap filled data is not available for Harvard Forest

# Willow creek
# 
#data up to date in March 2016 - still same version on June 17th 
# read all the *.csv files in the working directory
filepath="Willow_Creek/gap_filled/"
#change to reflect file path
tempFilelist = list.files(filepath, pattern="*.csv")
# get the AmeriFlux files information (site and years)
files<-substr(tempFilelist,7,14)
# get the AmeriFlux site code
site<-substr(tempFilelist[1],7,9)
# get the data and header from the L2 files
WillowCreekWCr = do.call("rbind", lapply(paste0(filepath,tempFilelist), function(x) read.csv(x, skip=20,header=FALSE,stringsAsFactors = FALSE)))
AmFluxheader=read.csv(paste0(filepath,tempFilelist[1]),skip=17, strip.white=TRUE, nrows=1 ,header=FALSE,stringsAsFactors=FALSE)
colnames(WillowCreekWCr)<-AmFluxheader

# Sylvania
# 
#data up to date in March 2016 - still same version on June 17th 
# read all the *.csv files in the working directory
filepath="Sylvania_Wilderness/gap_filled/"
#change to reflect file path
tempFilelist = list.files(filepath, pattern="*.csv")
# get the AmeriFlux files information (site and years)
files<-substr(tempFilelist,7,14)
# get the AmeriFlux site code
site<-substr(tempFilelist[1],7,9)
# get the data and header from the L2 files
SylvaniaSyv = do.call("rbind", lapply(paste0(filepath,tempFilelist), function(x) read.csv(x, skip=20,header=FALSE,stringsAsFactors = FALSE)))
AmFluxheader=read.csv(paste0(filepath,tempFilelist[1]),skip=17, strip.white=TRUE, nrows=1 ,header=FALSE,stringsAsFactors=FALSE)
colnames(SylvaniaSyv)<-AmFluxheader


setwd("D:/Dropbox/rProjectsShare/MIP-Change-and-Stability/data/fluxdata/")
save(HowlandHo1, file = "HowlandHo1.ameriflux.allsites.L2_data.17Jun2016.RData")
save(HarvardHa1, file = "HarvardHa1.ameriflux.allsites.L2_data.05Mar2016.RData")
save(WillowCreekWCr, file = "WillowCreekWCr.ameriflux.allsites.L2_data.17Jun2016.RData")
save(SylvaniaSyv, file = "SylvaniaSyv.ameriflux.allsites.L2_data.17Jun2016.RData")


