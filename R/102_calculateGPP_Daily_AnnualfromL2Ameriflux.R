# Author: Dave Moore
# Date started: 06-17-2016
# Explore L2 Data and create GPP estimate in gC/m^2/timestep
# 1) How many missing data obs(in gapfilled data)

# Load required libraries
library(dplyr)
library(ggplot2)
library(tidyr)

#load data
#You should have 4 .Rdata files
# mine are in a folder called "./data/fluxdata"
# setwd("D:/Dropbox/rProjectsShare/MIP-Change-and-Stability")

load("./data/fluxdata/HowlandHo1.ameriflux.allsites.L2_data.17Jun2016.RData")
load("./data/fluxdata/HarvardHa1.ameriflux.allsites.L2_data.05Mar2016.RData")
load("./data/fluxdata/WillowCreekWCr.ameriflux.allsites.L2_data.17Jun2016.RData")
load("./data/fluxdata/WillowCreekWCr.ameriflux.allsites.L2_data.17Jun2016.RData")

hist(HarvardHa1$GPP[HarvardHa1$GPP==-9999])
plot(HarvardHa1$GPP[HarvardHa1$GAP==0])

#define constants
rhoH2O = 1000 #density of water (1mm H2O * m2 = 1 kg)

#1000 kg/m^3 
LHVAP = 2502000 #latent heat of evaporation J/kg @ 0 C
SecPerhalfhour = 30*60 #seconds per half hour
SecPerDay = 86400 #seconds per day
#assume no sublimation

#define RHO 1.3 // air density (kg/m^3)
#define CP 1005. // specific heat of air (J/(kg K))
#define GAMMA 66. // psychometric constant (Pa/K)
#define E_STAR_SNOW 0.6 /* approximate saturation vapor pressure at 0 degrees C (kPa) 
# (we assume snow temperature is 0 degrees C or slightly lower) */

# 
# Evaportranspiration from latent heat of evaporation
# 
#Units LE (J s-1 m-2)
#Units E (kg m-2 s-1) 
# Note that 1 kg m-2 of evaporation is the same as a 1 mm of water over a square
# meter (assuming a water density of 1000 kg m-3).

# LE = LHVAP*EVAPOtrans

#molecular weights to convert from moles to g
MolwtC = 12.0107
MolwtO2 =16*2 
MolwtCO2 = 44.01

# VPsat = 0.611.exp[(17.3*T)/(T+237.3)]

# 
# Harvard
# 
HarvardHa1=HarvardHa1 %>%
  filter(GPP>-100) %>% #remove single -9999 code
  #GPP is reported in micromoles CO2 per m^2 per sec
  #we want GPP in gC per m^2 per time step | time step is 30 minutes
  mutate(GPPgC = (GPP/1000000)*(SecPerhalfhour*MolwtC))  #GPP in gC
  

Harvard_daily = group_by(HarvardHa1,YEAR,DOY)
# 

# Harvard_daySum = Harvard_daily %>%
#   summarise( n=n(),ETdaily=sum(EVAPOtrans, na.rm = TRUE), VPDdaily=mean(VPD, na.rm = TRUE),GPPdaily=sum(GPPgC, na.rm = TRUE),Precip=sum(PREC, na.rm = TRUE), LEgapsCT=sum(LEgaps), VPDgapsCT=sum(VPDgaps))

#calculate daily sums
Harvard_daySum = Harvard_daily %>%
  summarise( n=n(), VPDdaily=mean(VPD, na.rm = TRUE),GPPdaily=sum(GPPgC, na.rm = TRUE),Precip=sum(PREC, na.rm = TRUE))

#calculate annual sums
Harvard_annual = group_by(HarvardHa1,YEAR)
Harvard_YearSum = Harvard_annual %>%
  summarise( n=n(), VPDdaily=mean(VPD, na.rm = TRUE),GPPdaily=sum(GPPgC, na.rm = TRUE),Precip=sum(PREC, na.rm = TRUE))

#quick plot Harvard
plot(Harvard_YearSum$YEAR,Harvard_YearSum$GPP)


####
##Howland - main tower
####

HowlandHo1=HowlandHo1 %>%
  filter(GPP>-100) %>% #remove single -9999 code
  #GPP is reported in micromoles CO2 per m^2 per sec
  #we want GPP in gC per m^2 per time step | time step is 30 minutes
  mutate(GPPgC = (GPP/1000000)*(SecPerhalfhour*MolwtC))  #GPP in gC


Howland_daily = group_by(HowlandHo1,YEAR,DOY)
# 

# Howland_daySum = Howland_daily %>%
#   summarise( n=n(),ETdaily=sum(EVAPOtrans, na.rm = TRUE), VPDdaily=mean(VPD, na.rm = TRUE),GPPdaily=sum(GPPgC, na.rm = TRUE),Precip=sum(PREC, na.rm = TRUE), LEgapsCT=sum(LEgaps), VPDgapsCT=sum(VPDgaps))

#calculate daily sums
Howland_daySum = Howland_daily %>%
  summarise( n=n(), VPDdaily=mean(VPD, na.rm = TRUE),GPPdaily=sum(GPPgC, na.rm = TRUE),Precip=sum(PREC, na.rm = TRUE))

#calculate annual sums
Howland_annual = group_by(HowlandHo1,YEAR)
Howland_YearSum = Howland_annual %>%
  summarise( n=n(), VPDdaily=mean(VPD, na.rm = TRUE),GPPdaily=sum(GPPgC, na.rm = TRUE),Precip=sum(PREC, na.rm = TRUE))

#quick plot Howland
plot(Howland_YearSum$YEAR,Howland_YearSum$GPP)



####
## Sylvania
####

WillowCreekWCr=WillowCreekWCr %>%
  filter(GPP>-100) %>% #remove single -9999 code
  #GPP is reported in micromoles CO2 per m^2 per sec
  #we want GPP in gC per m^2 per time step | time step is 30 minutes
  mutate(GPPgC = (GPP/1000000)*(SecPerhalfhour*MolwtC))  #GPP in gC


Sylvania_daily = group_by(WillowCreekWCr,YEAR,DOY)
# 

# Sylvania_daySum = Sylvania_daily %>%
#   summarise( n=n(),ETdaily=sum(EVAPOtrans, na.rm = TRUE), VPDdaily=mean(VPD, na.rm = TRUE),GPPdaily=sum(GPPgC, na.rm = TRUE),Precip=sum(PREC, na.rm = TRUE), LEgapsCT=sum(LEgaps), VPDgapsCT=sum(VPDgaps))

#calculate daily sums
Sylvania_daySum = Sylvania_daily %>%
  summarise( n=n(), VPDdaily=mean(VPD, na.rm = TRUE),GPPdaily=sum(GPPgC, na.rm = TRUE),Precip=sum(PREC, na.rm = TRUE))

#calculate annual sums
Sylvania_annual = group_by(WillowCreekWCr,YEAR)
Sylvania_YearSum = Sylvania_annual %>%
  summarise( n=n(), VPDdaily=mean(VPD, na.rm = TRUE),GPPdaily=sum(GPPgC, na.rm = TRUE),Precip=sum(PREC, na.rm = TRUE))

#quick plot Howland
plot(Sylvania_YearSum$YEAR,Sylvania_YearSum$GPP)


####
## WillowCreek
####

WillowCreekWCr=WillowCreekWCr %>%
  filter(GPP>-100) %>% #remove single -9999 code
  #GPP is reported in micromoles CO2 per m^2 per sec
  #we want GPP in gC per m^2 per time step | time step is 30 minutes
  mutate(GPPgC = (GPP/1000000)*(SecPerhalfhour*MolwtC))  #GPP in gC


WillowCreek_daily = group_by(WillowCreekWCr,YEAR,DOY)
# 

# Sylvania_daySum = Sylvania_daily %>%
#   summarise( n=n(),ETdaily=sum(EVAPOtrans, na.rm = TRUE), VPDdaily=mean(VPD, na.rm = TRUE),GPPdaily=sum(GPPgC, na.rm = TRUE),Precip=sum(PREC, na.rm = TRUE), LEgapsCT=sum(LEgaps), VPDgapsCT=sum(VPDgaps))

#calculate daily sums
WillowCreek_daySum = WillowCreek_daily %>%
  summarise( n=n(), VPDdaily=mean(VPD, na.rm = TRUE),GPPdaily=sum(GPPgC, na.rm = TRUE),Precip=sum(PREC, na.rm = TRUE))

#calculate annual sums
WillowCreek_annual = group_by(WillowCreekWCr,YEAR)
WillowCreek_YearSum = WillowCreek_annual %>%
  summarise( n=n(), VPDdaily=mean(VPD, na.rm = TRUE),GPPdaily=sum(GPPgC, na.rm = TRUE),Precip=sum(PREC, na.rm = TRUE))

#quick plot Howland
plot(WillowCreek_YearSum$YEAR,WillowCreek_YearSum$GPP)

