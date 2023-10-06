################################################################################### 
# Climate sensitivity functions for bees in three ecosystems at the Sevilleta National Wildlife Refuge
# Power analysis to assess the probability of detecting trends for bee populations in which aridity did not predict abundance

# Heat and desiccation tolerances predict bee abundance under climate change
# Melanie R. Kazenel, Karen W. Wright, Terry Griswold, Kenneth D. Whitney, and Jennifer A. Rudgers

# Date: 2023-08-29
# Corresponding author's email: melanie.kazenel@gmail.com
################################################################################### 


# Load required packages
library(dplyr)
library(reshape2)
library(car)
library(MuMIn)
library(ggplot2)
library(nlme)
library(piecewiseSEM)
library(nlmeU)
library(tidyr)
library(stringr)



##### Format data #####

# Read in data frame of square-root transformed bee abundance data and SPEI data
max_abund_wide2<-read.csv("bee_wide_year_2002-2019_no2016or2017_sqrtmaxabund_2023-08-29.csv")

# Read in list of which CSF models ran
which_ran<-read.csv("CSFs_which_models_run_no2016or2017_2023-08-29.csv")
which_ran<-which_ran[-1,] # remove blank row

# Get summary of how many models ran (skip_to_next = FALSE) vs. did not (skip_to_next = TRUE) for each species x ecosystem combination
which_ran_summary<-which_ran%>%group_by(code,ecosystem,skip_to_next)%>%summarise(count=n())

# Create separate data frames of just the species for which all models ran in a given ecosystem
which_ran_subset<-subset(which_ran_summary,skip_to_next==FALSE & count==37)
which_ran_C<-subset(which_ran_subset,ecosystem=="C")
which_ran_G<-subset(which_ran_subset,ecosystem=="G")
which_ran_B<-subset(which_ran_subset,ecosystem=="B")

# For each ecosystem, select columns of interest, and exclude bee species for which CSFs did not run
creo_prelim<-subset(max_abund_wide2,ecosystem=="C")
creo_env<-creo_prelim[,c("ecosystem","transect","year","station_year", "spring6SPEI","monsoon6SPEI","spring6SPEI_prioryear", "monsoon6SPEI_prioryear")]
creo_data<-creo_prelim[,which_ran_C$code]
creo<-cbind(creo_env,creo_data)

black_prelim <- subset(max_abund_wide2, ecosystem == "G")
black_env<-black_prelim[,c("ecosystem","transect","year","station_year", "spring6SPEI","monsoon6SPEI","spring6SPEI_prioryear", "monsoon6SPEI_prioryear")]
black_data<-black_prelim[,which_ran_G$code]
black<-cbind(black_env,black_data)

blue_prelim <- subset(max_abund_wide2, ecosystem == "B" & year != 2002 & year != 2003)
blue_env<-blue_prelim[,c("ecosystem","transect","year","station_year", "spring6SPEI","monsoon6SPEI","spring6SPEI_prioryear", "monsoon6SPEI_prioryear")]
blue_data<-blue_prelim[,which_ran_B$code]
blue<-cbind(blue_env,blue_data)

# Create matrix to hold statistical results
beeCSF_output<-matrix(nrow=1,ncol=182,byrow=TRUE,dimnames=list(c("row1"),c("code","ecosystem","m1_spring6mo_numDF",	"m2_spring6mo_numDF",	"m3_spring6mo_numDF",	"m1_spring6mo_AR1_numDF",	"m1_spring6mo_AR2_numDF",	"m2_spring6mo_AR1_numDF",	"m2_spring6mo_AR2_numDF",	"m3_spring6mo_AR1_numDF",	"m3_spring6mo_AR2_numDF",	"m1_monsoon6mo_numDF",	"m2_monsoon6mo_numDF",	"m3_monsoon6mo_numDF",	"m1_monsoon6mo_AR1_numDF",	"m1_monsoon6mo_AR2_numDF",	"m2_monsoon6mo_AR1_numDF",	"m2_monsoon6mo_AR2_numDF",	"m3_monsoon6mo_AR1_numDF",	"m3_monsoon6mo_AR2_numDF",	"m1_spring6mo_denDF",	"m2_spring6mo_denDF",	"m3_spring6mo_denDF",	"m1_spring6mo_AR1_denDF",	"m1_spring6mo_AR2_denDF",	"m2_spring6mo_AR1_denDF",	"m2_spring6mo_AR2_denDF",	"m3_spring6mo_AR1_denDF",	"m3_spring6mo_AR2_denDF",	"m1_monsoon6mo_denDF",	"m2_monsoon6mo_denDF",	"m3_monsoon6mo_denDF",	"m1_monsoon6mo_AR1_denDF",	"m1_monsoon6mo_AR2_denDF",	"m2_monsoon6mo_AR1_denDF",	"m2_monsoon6mo_AR2_denDF",	"m3_monsoon6mo_AR1_denDF",	"m3_monsoon6mo_AR2_denDF",	"m1_spring6mo_F",	"m2_spring6mo_F",	"m3_spring6mo_F",	"m1_spring6mo_AR1_F",	"m1_spring6mo_AR2_F",	"m2_spring6mo_AR1_F",	"m2_spring6mo_AR2_F",	"m3_spring6mo_AR1_F",	"m3_spring6mo_AR2_F",	"m1_monsoon6mo_F",	"m2_monsoon6mo_F",	"m3_monsoon6mo_F",	"m1_monsoon6mo_AR1_F",	"m1_monsoon6mo_AR2_F",	"m2_monsoon6mo_AR1_F",	"m2_monsoon6mo_AR2_F",	"m3_monsoon6mo_AR1_F",	"m3_monsoon6mo_AR2_F",	"m1_spring6mo_nc",	"m2_spring6mo_nc",	"m3_spring6mo_nc",	"m1_spring6mo_AR1_nc",	"m1_spring6mo_AR2_nc",	"m2_spring6mo_AR1_nc",	"m2_spring6mo_AR2_nc",	"m3_spring6mo_AR1_nc",	"m3_spring6mo_AR2_nc",	"m1_monsoon6mo_nc",	"m2_monsoon6mo_nc",	"m3_monsoon6mo_nc",	"m1_monsoon6mo_AR1_nc",	"m1_monsoon6mo_AR2_nc",	"m2_monsoon6mo_AR1_nc",	"m2_monsoon6mo_AR2_nc",	"m3_monsoon6mo_AR1_nc",	"m3_monsoon6mo_AR2_nc",	"m1_spring6mo_power",	"m2_spring6mo_power",	"m3_spring6mo_power",	"m1_spring6mo_AR1_power",	"m1_spring6mo_AR2_power",	"m2_spring6mo_AR1_power",	"m2_spring6mo_AR2_power",	"m3_spring6mo_AR1_power",	"m3_spring6mo_AR2_power",	"m1_monsoon6mo_power",	"m2_monsoon6mo_power",	"m3_monsoon6mo_power",	"m1_monsoon6mo_AR1_power",	"m1_monsoon6mo_AR2_power",	"m2_monsoon6mo_AR1_power",	"m2_monsoon6mo_AR2_power",	"m3_monsoon6mo_AR1_power",	"m3_monsoon6mo_AR2_power", "m1_spring6mo_lag_numDF",	"m2_spring6mo_lag_numDF","m3_spring6mo_lag_numDF","m1_spring6mo_lag_AR1_numDF",
"m1_spring6mo_lag_AR2_numDF","m2_spring6mo_lag_AR1_numDF","m2_spring6mo_lag_AR2_numDF","m3_spring6mo_lag_AR1_numDF","m3_spring6mo_lag_AR2_numDF","m1_monsoon6mo_lag_numDF","m2_monsoon6mo_lag_numDF","m3_monsoon6mo_lag_numDF","m1_monsoon6mo_lag_AR1_numDF","m1_monsoon6mo_lag_AR2_numDF","m2_monsoon6mo_lag_AR1_numDF", "m2_monsoon6mo_lag_AR2_numDF","m3_monsoon6mo_lag_AR1_numDF","m3_monsoon6mo_lag_AR2_numDF","m1_spring6mo_lag_denDF","m2_spring6mo_lag_denDF","m3_spring6mo_lag_denDF","m1_spring6mo_lag_AR1_denDF","m1_spring6mo_lag_AR2_denDF","m2_spring6mo_lag_AR1_denDF","m2_spring6mo_lag_AR2_denDF","m3_spring6mo_lag_AR1_denDF","m3_spring6mo_lag_AR2_denDF","m1_monsoon6mo_lag_denDF", "m2_monsoon6mo_lag_denDF","m3_monsoon6mo_lag_denDF","m1_monsoon6mo_lag_AR1_denDF","m1_monsoon6mo_lag_AR2_denDF", "m2_monsoon6mo_lag_AR1_denDF","m2_monsoon6mo_lag_AR2_denDF","m3_monsoon6mo_lag_AR1_denDF", "m3_monsoon6mo_lag_AR2_denDF","m1_spring6mo_lag_F","m2_spring6mo_lag_F","m3_spring6mo_lag_F","m1_spring6mo_lag_AR1_F","m1_spring6mo_lag_AR2_F","m2_spring6mo_lag_AR1_F","m2_spring6mo_lag_AR2_F","m3_spring6mo_lag_AR1_F", "m3_spring6mo_lag_AR2_F","m1_monsoon6mo_lag_F","m2_monsoon6mo_lag_F","m3_monsoon6mo_lag_F","m1_monsoon6mo_lag_AR1_F","m1_monsoon6mo_lag_AR2_F","m2_monsoon6mo_lag_AR1_F","m2_monsoon6mo_lag_AR2_F","m3_monsoon6mo_lag_AR1_F",	"m3_monsoon6mo_lag_AR2_F","m1_spring6mo_lag_nc","m2_spring6mo_lag_nc","m3_spring6mo_lag_nc","m1_spring6mo_lag_AR1_nc",	"m1_spring6mo_lag_AR2_nc","m2_spring6mo_lag_AR1_nc","m2_spring6mo_lag_AR2_nc","m3_spring6mo_lag_AR1_nc", "m3_spring6mo_lag_AR2_nc","m1_monsoon6mo_lag_nc","m2_monsoon6mo_lag_nc","m3_monsoon6mo_lag_nc",	"m1_monsoon6mo_lag_AR1_nc","m1_monsoon6mo_lag_AR2_nc","m2_monsoon6mo_lag_AR1_nc","m2_monsoon6mo_lag_AR2_nc",	"m3_monsoon6mo_lag_AR1_nc","m3_monsoon6mo_lag_AR2_nc","m1_spring6mo_lag_power","m2_spring6mo_lag_power", "m3_spring6mo_lag_power","m1_spring6mo_lag_AR1_power","m1_spring6mo_lag_AR2_power","m2_spring6mo_lag_AR1_power",	"m2_spring6mo_lag_AR2_power","m3_spring6mo_lag_AR1_power","m3_spring6mo_lag_AR2_power","m1_monsoon6mo_lag_power", "m2_monsoon6mo_lag_power","m3_monsoon6mo_lag_power","m1_monsoon6mo_lag_AR1_power","m1_monsoon6mo_lag_AR2_power",	"m2_monsoon6mo_lag_AR1_power","m2_monsoon6mo_lag_AR2_power","m3_monsoon6mo_lag_AR1_power",	"m3_monsoon6mo_lag_AR2_power")))


##### Chihuahuan Desert Shrubland #####

# Create a data frame of just the bee abundance matrix (descriptor variables removed)
speciesMatrix <- creo[,9:229]

# Create a vector of species codes
speciesCodes <- colnames(speciesMatrix)


### Loop through each column of speciesMatrix (each species), running CSF models, running a power analysis for each model, and putting the results in the beeCSF_output matrix

for (i in 1:length(speciesMatrix[1,])) {
  
  # save the species code for column i
  speciesCode <- speciesCodes[i]
  
  # create an object with the name of the ecosystem type
  ecosystem <-"C"
  
  # run mixed effects models
  
  # null model
  m_null<-lme(formula(paste(speciesCode, "~1")), random=~1|transect, data=creo, method="ML")
  
  # models including current year's climate
  
  # spring - current year
  m1_spring6mo<-lme(formula(paste(speciesCode, "~spring6SPEI")), random=~1|transect, data=creo, method="ML")
  m2_spring6mo<-update(m1_spring6mo, .~. +I(spring6SPEI^2))
  m3_spring6mo<-update(m2_spring6mo, .~. +I(spring6SPEI^3))
  
  m1_spring6mo_AR1 <- update(m1_spring6mo,.~.,correlation=corAR1(form=~year|transect))
  m1_spring6mo_AR2 <- update(m1_spring6mo,.~.,correlation=corARMA(form=~year|transect,p=2))
  
  m2_spring6mo_AR1 <- update(m2_spring6mo,.~.,correlation=corAR1(form=~year|transect))
  m2_spring6mo_AR2 <- update(m2_spring6mo,.~.,correlation=corARMA(form=~year|transect,p=2))
  
  m3_spring6mo_AR1 <- update(m3_spring6mo,.~.,correlation=corAR1(form=~year|transect))
  m3_spring6mo_AR2 <- update(m3_spring6mo,.~.,correlation=corARMA(form=~year|transect,p=2))
  
  # monsoon - current year
  m1_monsoon6mo<-lme(formula(paste(speciesCode, "~monsoon6SPEI")), random=~1|transect, data=creo, method="ML")
  m2_monsoon6mo<-update(m1_monsoon6mo, .~. +I(monsoon6SPEI^2))
  m3_monsoon6mo<-update(m2_monsoon6mo, .~. +I(monsoon6SPEI^3))
  
  m1_monsoon6mo_AR1 <- update(m1_monsoon6mo,.~.,correlation=corAR1(form=~year|transect))
  m1_monsoon6mo_AR2 <- update(m1_monsoon6mo,.~.,correlation=corARMA(form=~year|transect,p=2))
  
  m2_monsoon6mo_AR1 <- update(m2_monsoon6mo,.~.,correlation=corAR1(form=~year|transect))
  m2_monsoon6mo_AR2 <- update(m2_monsoon6mo,.~.,correlation=corARMA(form=~year|transect,p=2))
  
  m3_monsoon6mo_AR1 <- update(m3_monsoon6mo,.~.,correlation=corAR1(form=~year|transect))
  m3_monsoon6mo_AR2 <- update(m3_monsoon6mo,.~.,correlation=corARMA(form=~year|transect,p=2))
  
  # models including previous year's climate
  
  # spring - previous year
  m1_spring6mo_lag<-lme(formula(paste(speciesCode, "~spring6SPEI_prioryear")), random=~1|transect, data=creo, method="ML")
  m2_spring6mo_lag<-update(m1_spring6mo_lag, .~. +I(spring6SPEI_prioryear^2))
  m3_spring6mo_lag<-update(m2_spring6mo_lag, .~. +I(spring6SPEI_prioryear^3))
  
  m1_spring6mo_lag_AR1 <- update(m1_spring6mo_lag,.~.,correlation=corAR1(form=~year|transect))
  m1_spring6mo_lag_AR2 <- update(m1_spring6mo_lag,.~.,correlation=corARMA(form=~year|transect,p=2))
  
  m2_spring6mo_lag_AR1 <- update(m2_spring6mo_lag,.~.,correlation=corAR1(form=~year|transect))
  m2_spring6mo_lag_AR2 <- update(m2_spring6mo_lag,.~.,correlation=corARMA(form=~year|transect,p=2))
  
  m3_spring6mo_lag_AR1 <- update(m3_spring6mo_lag,.~.,correlation=corAR1(form=~year|transect))
  m3_spring6mo_lag_AR2 <- update(m3_spring6mo_lag,.~.,correlation=corARMA(form=~year|transect,p=2))
  
  # monsoon - previous year
  m1_monsoon6mo_lag<-lme(formula(paste(speciesCode, "~monsoon6SPEI_prioryear")), random=~1|transect, data=creo, method="ML")
  m2_monsoon6mo_lag<-update(m1_monsoon6mo_lag, .~. +I(monsoon6SPEI_prioryear^2))
  m3_monsoon6mo_lag<-update(m2_monsoon6mo_lag, .~. +I(monsoon6SPEI_prioryear^3))
  
  m1_monsoon6mo_lag_AR1 <- update(m1_monsoon6mo_lag,.~.,correlation=corAR1(form=~year|transect))
  m1_monsoon6mo_lag_AR2 <- update(m1_monsoon6mo_lag,.~.,correlation=corARMA(form=~year|transect,p=2))
  
  m2_monsoon6mo_lag_AR1 <- update(m2_monsoon6mo_lag,.~.,correlation=corAR1(form=~year|transect))
  m2_monsoon6mo_lag_AR2 <- update(m2_monsoon6mo_lag,.~.,correlation=corARMA(form=~year|transect,p=2))
  
  m3_monsoon6mo_lag_AR1 <- update(m3_monsoon6mo_lag,.~.,correlation=corAR1(form=~year|transect))
  m3_monsoon6mo_lag_AR2 <- update(m3_monsoon6mo_lag,.~.,correlation=corARMA(form=~year|transect,p=2))
  
  # run power analyses for each model, creating objects containing the statistical output
  
  # current year models
  m1_spring6mo_numDF<-Pwr(m1_spring6mo)[2,1]
  m2_spring6mo_numDF<-Pwr(m2_spring6mo)[2,1]
  m3_spring6mo_numDF<-Pwr(m3_spring6mo)[2,1]
  m1_spring6mo_AR1_numDF<-Pwr(m1_spring6mo_AR1)[2,1]
  m1_spring6mo_AR2_numDF<-Pwr(m1_spring6mo_AR2)[2,1]
  m2_spring6mo_AR1_numDF<-Pwr(m2_spring6mo_AR1)[2,1]
  m2_spring6mo_AR2_numDF<-Pwr(m2_spring6mo_AR2)[2,1]
  m3_spring6mo_AR1_numDF<-Pwr(m3_spring6mo_AR1)[2,1]
  m3_spring6mo_AR2_numDF<-Pwr(m3_spring6mo_AR2)[2,1]
  m1_monsoon6mo_numDF<-Pwr(m1_monsoon6mo)[2,1]
  m2_monsoon6mo_numDF<-Pwr(m2_monsoon6mo)[2,1]
  m3_monsoon6mo_numDF<-Pwr(m3_monsoon6mo)[2,1]
  m1_monsoon6mo_AR1_numDF<-Pwr(m1_monsoon6mo_AR1)[2,1]
  m1_monsoon6mo_AR2_numDF<-Pwr(m1_monsoon6mo_AR2)[2,1]
  m2_monsoon6mo_AR1_numDF<-Pwr(m2_monsoon6mo_AR1)[2,1]
  m2_monsoon6mo_AR2_numDF<-Pwr(m2_monsoon6mo_AR2)[2,1]
  m3_monsoon6mo_AR1_numDF<-Pwr(m3_monsoon6mo_AR1)[2,1]
  m3_monsoon6mo_AR2_numDF<-Pwr(m3_monsoon6mo_AR2)[2,1]
  m1_spring6mo_denDF<-Pwr(m1_spring6mo)[2,2]
  m2_spring6mo_denDF<-Pwr(m2_spring6mo)[2,2]
  m3_spring6mo_denDF<-Pwr(m3_spring6mo)[2,2]
  m1_spring6mo_AR1_denDF<-Pwr(m1_spring6mo_AR1)[2,2]
  m1_spring6mo_AR2_denDF<-Pwr(m1_spring6mo_AR2)[2,2]
  m2_spring6mo_AR1_denDF<-Pwr(m2_spring6mo_AR1)[2,2]
  m2_spring6mo_AR2_denDF<-Pwr(m2_spring6mo_AR2)[2,2]
  m3_spring6mo_AR1_denDF<-Pwr(m3_spring6mo_AR1)[2,2]
  m3_spring6mo_AR2_denDF<-Pwr(m3_spring6mo_AR2)[2,2]
  m1_monsoon6mo_denDF<-Pwr(m1_monsoon6mo)[2,2]
  m2_monsoon6mo_denDF<-Pwr(m2_monsoon6mo)[2,2]
  m3_monsoon6mo_denDF<-Pwr(m3_monsoon6mo)[2,2]
  m1_monsoon6mo_AR1_denDF<-Pwr(m1_monsoon6mo_AR1)[2,2]
  m1_monsoon6mo_AR2_denDF<-Pwr(m1_monsoon6mo_AR2)[2,2]
  m2_monsoon6mo_AR1_denDF<-Pwr(m2_monsoon6mo_AR1)[2,2]
  m2_monsoon6mo_AR2_denDF<-Pwr(m2_monsoon6mo_AR2)[2,2]
  m3_monsoon6mo_AR1_denDF<-Pwr(m3_monsoon6mo_AR1)[2,2]
  m3_monsoon6mo_AR2_denDF<-Pwr(m3_monsoon6mo_AR2)[2,2]
  m1_spring6mo_F<-Pwr(m1_spring6mo)[2,3]
  m2_spring6mo_F<-Pwr(m2_spring6mo)[2,3]
  m3_spring6mo_F<-Pwr(m3_spring6mo)[2,3]
  m1_spring6mo_AR1_F<-Pwr(m1_spring6mo_AR1)[2,3]
  m1_spring6mo_AR2_F<-Pwr(m1_spring6mo_AR2)[2,3]
  m2_spring6mo_AR1_F<-Pwr(m2_spring6mo_AR1)[2,3]
  m2_spring6mo_AR2_F<-Pwr(m2_spring6mo_AR2)[2,3]
  m3_spring6mo_AR1_F<-Pwr(m3_spring6mo_AR1)[2,3]
  m3_spring6mo_AR2_F<-Pwr(m3_spring6mo_AR2)[2,3]
  m1_monsoon6mo_F<-Pwr(m1_monsoon6mo)[2,3]
  m2_monsoon6mo_F<-Pwr(m2_monsoon6mo)[2,3]
  m3_monsoon6mo_F<-Pwr(m3_monsoon6mo)[2,3]
  m1_monsoon6mo_AR1_F<-Pwr(m1_monsoon6mo_AR1)[2,3]
  m1_monsoon6mo_AR2_F<-Pwr(m1_monsoon6mo_AR2)[2,3]
  m2_monsoon6mo_AR1_F<-Pwr(m2_monsoon6mo_AR1)[2,3]
  m2_monsoon6mo_AR2_F<-Pwr(m2_monsoon6mo_AR2)[2,3]
  m3_monsoon6mo_AR1_F<-Pwr(m3_monsoon6mo_AR1)[2,3]
  m3_monsoon6mo_AR2_F<-Pwr(m3_monsoon6mo_AR2)[2,3]
  m1_spring6mo_nc<-Pwr(m1_spring6mo)[2,4]
  m2_spring6mo_nc<-Pwr(m2_spring6mo)[2,4]
  m3_spring6mo_nc<-Pwr(m3_spring6mo)[2,4]
  m1_spring6mo_AR1_nc<-Pwr(m1_spring6mo_AR1)[2,4]
  m1_spring6mo_AR2_nc<-Pwr(m1_spring6mo_AR2)[2,4]
  m2_spring6mo_AR1_nc<-Pwr(m2_spring6mo_AR1)[2,4]
  m2_spring6mo_AR2_nc<-Pwr(m2_spring6mo_AR2)[2,4]
  m3_spring6mo_AR1_nc<-Pwr(m3_spring6mo_AR1)[2,4]
  m3_spring6mo_AR2_nc<-Pwr(m3_spring6mo_AR2)[2,4]
  m1_monsoon6mo_nc<-Pwr(m1_monsoon6mo)[2,4]
  m2_monsoon6mo_nc<-Pwr(m2_monsoon6mo)[2,4]
  m3_monsoon6mo_nc<-Pwr(m3_monsoon6mo)[2,4]
  m1_monsoon6mo_AR1_nc<-Pwr(m1_monsoon6mo_AR1)[2,4]
  m1_monsoon6mo_AR2_nc<-Pwr(m1_monsoon6mo_AR2)[2,4]
  m2_monsoon6mo_AR1_nc<-Pwr(m2_monsoon6mo_AR1)[2,4]
  m2_monsoon6mo_AR2_nc<-Pwr(m2_monsoon6mo_AR2)[2,4]
  m3_monsoon6mo_AR1_nc<-Pwr(m3_monsoon6mo_AR1)[2,4]
  m3_monsoon6mo_AR2_nc<-Pwr(m3_monsoon6mo_AR2)[2,4]
  m1_spring6mo_power<-Pwr(m1_spring6mo)[2,5]
  m2_spring6mo_power<-Pwr(m2_spring6mo)[2,5]
  m3_spring6mo_power<-Pwr(m3_spring6mo)[2,5]
  m1_spring6mo_AR1_power<-Pwr(m1_spring6mo_AR1)[2,5]
  m1_spring6mo_AR2_power<-Pwr(m1_spring6mo_AR2)[2,5]
  m2_spring6mo_AR1_power<-Pwr(m2_spring6mo_AR1)[2,5]
  m2_spring6mo_AR2_power<-Pwr(m2_spring6mo_AR2)[2,5]
  m3_spring6mo_AR1_power<-Pwr(m3_spring6mo_AR1)[2,5]
  m3_spring6mo_AR2_power<-Pwr(m3_spring6mo_AR2)[2,5]
  m1_monsoon6mo_power<-Pwr(m1_monsoon6mo)[2,5]
  m2_monsoon6mo_power<-Pwr(m2_monsoon6mo)[2,5]
  m3_monsoon6mo_power<-Pwr(m3_monsoon6mo)[2,5]
  m1_monsoon6mo_AR1_power<-Pwr(m1_monsoon6mo_AR1)[2,5]
  m1_monsoon6mo_AR2_power<-Pwr(m1_monsoon6mo_AR2)[2,5]
  m2_monsoon6mo_AR1_power<-Pwr(m2_monsoon6mo_AR1)[2,5]
  m2_monsoon6mo_AR2_power<-Pwr(m2_monsoon6mo_AR2)[2,5]
  m3_monsoon6mo_AR1_power<-Pwr(m3_monsoon6mo_AR1)[2,5]
  m3_monsoon6mo_AR2_power<-Pwr(m3_monsoon6mo_AR2)[2,5]
  
  # previous year models
  m1_spring6mo_lag_numDF<-Pwr(m1_spring6mo_lag)[2,1]
  m2_spring6mo_lag_numDF<-Pwr(m2_spring6mo_lag)[2,1]
  m3_spring6mo_lag_numDF<-Pwr(m3_spring6mo_lag)[2,1]
  m1_spring6mo_lag_AR1_numDF<-Pwr(m1_spring6mo_lag_AR1)[2,1]
  m1_spring6mo_lag_AR2_numDF<-Pwr(m1_spring6mo_lag_AR2)[2,1]
  m2_spring6mo_lag_AR1_numDF<-Pwr(m2_spring6mo_lag_AR1)[2,1]
  m2_spring6mo_lag_AR2_numDF<-Pwr(m2_spring6mo_lag_AR2)[2,1]
  m3_spring6mo_lag_AR1_numDF<-Pwr(m3_spring6mo_lag_AR1)[2,1]
  m3_spring6mo_lag_AR2_numDF<-Pwr(m3_spring6mo_lag_AR2)[2,1]
  m1_monsoon6mo_lag_numDF<-Pwr(m1_monsoon6mo_lag)[2,1]
  m2_monsoon6mo_lag_numDF<-Pwr(m2_monsoon6mo_lag)[2,1]
  m3_monsoon6mo_lag_numDF<-Pwr(m3_monsoon6mo_lag)[2,1]
  m1_monsoon6mo_lag_AR1_numDF<-Pwr(m1_monsoon6mo_lag_AR1)[2,1]
  m1_monsoon6mo_lag_AR2_numDF<-Pwr(m1_monsoon6mo_lag_AR2)[2,1]
  m2_monsoon6mo_lag_AR1_numDF<-Pwr(m2_monsoon6mo_lag_AR1)[2,1]
  m2_monsoon6mo_lag_AR2_numDF<-Pwr(m2_monsoon6mo_lag_AR2)[2,1]
  m3_monsoon6mo_lag_AR1_numDF<-Pwr(m3_monsoon6mo_lag_AR1)[2,1]
  m3_monsoon6mo_lag_AR2_numDF<-Pwr(m3_monsoon6mo_lag_AR2)[2,1]
  m1_spring6mo_lag_denDF<-Pwr(m1_spring6mo_lag)[2,2]
  m2_spring6mo_lag_denDF<-Pwr(m2_spring6mo_lag)[2,2]
  m3_spring6mo_lag_denDF<-Pwr(m3_spring6mo_lag)[2,2]
  m1_spring6mo_lag_AR1_denDF<-Pwr(m1_spring6mo_lag_AR1)[2,2]
  m1_spring6mo_lag_AR2_denDF<-Pwr(m1_spring6mo_lag_AR2)[2,2]
  m2_spring6mo_lag_AR1_denDF<-Pwr(m2_spring6mo_lag_AR1)[2,2]
  m2_spring6mo_lag_AR2_denDF<-Pwr(m2_spring6mo_lag_AR2)[2,2]
  m3_spring6mo_lag_AR1_denDF<-Pwr(m3_spring6mo_lag_AR1)[2,2]
  m3_spring6mo_lag_AR2_denDF<-Pwr(m3_spring6mo_lag_AR2)[2,2]
  m1_monsoon6mo_lag_denDF<-Pwr(m1_monsoon6mo_lag)[2,2]
  m2_monsoon6mo_lag_denDF<-Pwr(m2_monsoon6mo_lag)[2,2]
  m3_monsoon6mo_lag_denDF<-Pwr(m3_monsoon6mo_lag)[2,2]
  m1_monsoon6mo_lag_AR1_denDF<-Pwr(m1_monsoon6mo_lag_AR1)[2,2]
  m1_monsoon6mo_lag_AR2_denDF<-Pwr(m1_monsoon6mo_lag_AR2)[2,2]
  m2_monsoon6mo_lag_AR1_denDF<-Pwr(m2_monsoon6mo_lag_AR1)[2,2]
  m2_monsoon6mo_lag_AR2_denDF<-Pwr(m2_monsoon6mo_lag_AR2)[2,2]
  m3_monsoon6mo_lag_AR1_denDF<-Pwr(m3_monsoon6mo_lag_AR1)[2,2]
  m3_monsoon6mo_lag_AR2_denDF<-Pwr(m3_monsoon6mo_lag_AR2)[2,2]
  m1_spring6mo_lag_F<-Pwr(m1_spring6mo_lag)[2,3]
  m2_spring6mo_lag_F<-Pwr(m2_spring6mo_lag)[2,3]
  m3_spring6mo_lag_F<-Pwr(m3_spring6mo_lag)[2,3]
  m1_spring6mo_lag_AR1_F<-Pwr(m1_spring6mo_lag_AR1)[2,3]
  m1_spring6mo_lag_AR2_F<-Pwr(m1_spring6mo_lag_AR2)[2,3]
  m2_spring6mo_lag_AR1_F<-Pwr(m2_spring6mo_lag_AR1)[2,3]
  m2_spring6mo_lag_AR2_F<-Pwr(m2_spring6mo_lag_AR2)[2,3]
  m3_spring6mo_lag_AR1_F<-Pwr(m3_spring6mo_lag_AR1)[2,3]
  m3_spring6mo_lag_AR2_F<-Pwr(m3_spring6mo_lag_AR2)[2,3]
  m1_monsoon6mo_lag_F<-Pwr(m1_monsoon6mo_lag)[2,3]
  m2_monsoon6mo_lag_F<-Pwr(m2_monsoon6mo_lag)[2,3]
  m3_monsoon6mo_lag_F<-Pwr(m3_monsoon6mo_lag)[2,3]
  m1_monsoon6mo_lag_AR1_F<-Pwr(m1_monsoon6mo_lag_AR1)[2,3]
  m1_monsoon6mo_lag_AR2_F<-Pwr(m1_monsoon6mo_lag_AR2)[2,3]
  m2_monsoon6mo_lag_AR1_F<-Pwr(m2_monsoon6mo_lag_AR1)[2,3]
  m2_monsoon6mo_lag_AR2_F<-Pwr(m2_monsoon6mo_lag_AR2)[2,3]
  m3_monsoon6mo_lag_AR1_F<-Pwr(m3_monsoon6mo_lag_AR1)[2,3]
  m3_monsoon6mo_lag_AR2_F<-Pwr(m3_monsoon6mo_lag_AR2)[2,3]
  m1_spring6mo_lag_nc<-Pwr(m1_spring6mo_lag)[2,4]
  m2_spring6mo_lag_nc<-Pwr(m2_spring6mo_lag)[2,4]
  m3_spring6mo_lag_nc<-Pwr(m3_spring6mo_lag)[2,4]
  m1_spring6mo_lag_AR1_nc<-Pwr(m1_spring6mo_lag_AR1)[2,4]
  m1_spring6mo_lag_AR2_nc<-Pwr(m1_spring6mo_lag_AR2)[2,4]
  m2_spring6mo_lag_AR1_nc<-Pwr(m2_spring6mo_lag_AR1)[2,4]
  m2_spring6mo_lag_AR2_nc<-Pwr(m2_spring6mo_lag_AR2)[2,4]
  m3_spring6mo_lag_AR1_nc<-Pwr(m3_spring6mo_lag_AR1)[2,4]
  m3_spring6mo_lag_AR2_nc<-Pwr(m3_spring6mo_lag_AR2)[2,4]
  m1_monsoon6mo_lag_nc<-Pwr(m1_monsoon6mo_lag)[2,4]
  m2_monsoon6mo_lag_nc<-Pwr(m2_monsoon6mo_lag)[2,4]
  m3_monsoon6mo_lag_nc<-Pwr(m3_monsoon6mo_lag)[2,4]
  m1_monsoon6mo_lag_AR1_nc<-Pwr(m1_monsoon6mo_lag_AR1)[2,4]
  m1_monsoon6mo_lag_AR2_nc<-Pwr(m1_monsoon6mo_lag_AR2)[2,4]
  m2_monsoon6mo_lag_AR1_nc<-Pwr(m2_monsoon6mo_lag_AR1)[2,4]
  m2_monsoon6mo_lag_AR2_nc<-Pwr(m2_monsoon6mo_lag_AR2)[2,4]
  m3_monsoon6mo_lag_AR1_nc<-Pwr(m3_monsoon6mo_lag_AR1)[2,4]
  m3_monsoon6mo_lag_AR2_nc<-Pwr(m3_monsoon6mo_lag_AR2)[2,4]
  m1_spring6mo_lag_power<-Pwr(m1_spring6mo_lag)[2,5]
  m2_spring6mo_lag_power<-Pwr(m2_spring6mo_lag)[2,5]
  m3_spring6mo_lag_power<-Pwr(m3_spring6mo_lag)[2,5]
  m1_spring6mo_lag_AR1_power<-Pwr(m1_spring6mo_lag_AR1)[2,5]
  m1_spring6mo_lag_AR2_power<-Pwr(m1_spring6mo_lag_AR2)[2,5]
  m2_spring6mo_lag_AR1_power<-Pwr(m2_spring6mo_lag_AR1)[2,5]
  m2_spring6mo_lag_AR2_power<-Pwr(m2_spring6mo_lag_AR2)[2,5]
  m3_spring6mo_lag_AR1_power<-Pwr(m3_spring6mo_lag_AR1)[2,5]
  m3_spring6mo_lag_AR2_power<-Pwr(m3_spring6mo_lag_AR2)[2,5]
  m1_monsoon6mo_lag_power<-Pwr(m1_monsoon6mo_lag)[2,5]
  m2_monsoon6mo_lag_power<-Pwr(m2_monsoon6mo_lag)[2,5]
  m3_monsoon6mo_lag_power<-Pwr(m3_monsoon6mo_lag)[2,5]
  m1_monsoon6mo_lag_AR1_power<-Pwr(m1_monsoon6mo_lag_AR1)[2,5]
  m1_monsoon6mo_lag_AR2_power<-Pwr(m1_monsoon6mo_lag_AR2)[2,5]
  m2_monsoon6mo_lag_AR1_power<-Pwr(m2_monsoon6mo_lag_AR1)[2,5]
  m2_monsoon6mo_lag_AR2_power<-Pwr(m2_monsoon6mo_lag_AR2)[2,5]
  m3_monsoon6mo_lag_AR1_power<-Pwr(m3_monsoon6mo_lag_AR1)[2,5]
  m3_monsoon6mo_lag_AR2_power<-Pwr(m3_monsoon6mo_lag_AR2)[2,5]
  
  # bind all target output values together
  output_id<-cbind(speciesCode, ecosystem, m1_spring6mo_numDF,	m2_spring6mo_numDF,	m3_spring6mo_numDF,	m1_spring6mo_AR1_numDF,	m1_spring6mo_AR2_numDF,	m2_spring6mo_AR1_numDF,	m2_spring6mo_AR2_numDF,	m3_spring6mo_AR1_numDF,	m3_spring6mo_AR2_numDF,	m1_monsoon6mo_numDF,	m2_monsoon6mo_numDF,	m3_monsoon6mo_numDF,	m1_monsoon6mo_AR1_numDF,	m1_monsoon6mo_AR2_numDF,	m2_monsoon6mo_AR1_numDF,	m2_monsoon6mo_AR2_numDF,	m3_monsoon6mo_AR1_numDF,	m3_monsoon6mo_AR2_numDF,	m1_spring6mo_denDF,	m2_spring6mo_denDF,	m3_spring6mo_denDF,	m1_spring6mo_AR1_denDF,	m1_spring6mo_AR2_denDF,	m2_spring6mo_AR1_denDF,	m2_spring6mo_AR2_denDF,	m3_spring6mo_AR1_denDF,	m3_spring6mo_AR2_denDF,	m1_monsoon6mo_denDF,	m2_monsoon6mo_denDF,	m3_monsoon6mo_denDF,	m1_monsoon6mo_AR1_denDF,	m1_monsoon6mo_AR2_denDF,	m2_monsoon6mo_AR1_denDF,	m2_monsoon6mo_AR2_denDF,	m3_monsoon6mo_AR1_denDF,	m3_monsoon6mo_AR2_denDF,	m1_spring6mo_F,	m2_spring6mo_F,	m3_spring6mo_F,	m1_spring6mo_AR1_F,	m1_spring6mo_AR2_F,	m2_spring6mo_AR1_F,	m2_spring6mo_AR2_F,	m3_spring6mo_AR1_F,	m3_spring6mo_AR2_F,	m1_monsoon6mo_F,	m2_monsoon6mo_F,	m3_monsoon6mo_F,	m1_monsoon6mo_AR1_F,	m1_monsoon6mo_AR2_F,	m2_monsoon6mo_AR1_F,	m2_monsoon6mo_AR2_F,	m3_monsoon6mo_AR1_F,	m3_monsoon6mo_AR2_F,	m1_spring6mo_nc,	m2_spring6mo_nc,	m3_spring6mo_nc,	m1_spring6mo_AR1_nc,	m1_spring6mo_AR2_nc,	m2_spring6mo_AR1_nc,	m2_spring6mo_AR2_nc,	m3_spring6mo_AR1_nc,	m3_spring6mo_AR2_nc,	m1_monsoon6mo_nc,	m2_monsoon6mo_nc,	m3_monsoon6mo_nc,	m1_monsoon6mo_AR1_nc,	m1_monsoon6mo_AR2_nc,	m2_monsoon6mo_AR1_nc,	m2_monsoon6mo_AR2_nc,	m3_monsoon6mo_AR1_nc,	m3_monsoon6mo_AR2_nc,	m1_spring6mo_power,	m2_spring6mo_power,	m3_spring6mo_power,	m1_spring6mo_AR1_power,	m1_spring6mo_AR2_power,	m2_spring6mo_AR1_power,	m2_spring6mo_AR2_power,	m3_spring6mo_AR1_power,	m3_spring6mo_AR2_power,	m1_monsoon6mo_power,	m2_monsoon6mo_power,	m3_monsoon6mo_power,	m1_monsoon6mo_AR1_power,	m1_monsoon6mo_AR2_power,	m2_monsoon6mo_AR1_power,	m2_monsoon6mo_AR2_power,	m3_monsoon6mo_AR1_power,	m3_monsoon6mo_AR2_power, m1_spring6mo_lag_numDF,m2_spring6mo_lag_numDF,m3_spring6mo_lag_numDF,m1_spring6mo_lag_AR1_numDF,	m1_spring6mo_lag_AR2_numDF,m2_spring6mo_lag_AR1_numDF,m2_spring6mo_lag_AR2_numDF,m3_spring6mo_lag_AR1_numDF,m3_spring6mo_lag_AR2_numDF,m1_monsoon6mo_lag_numDF,m2_monsoon6mo_lag_numDF,m3_monsoon6mo_lag_numDF,m1_monsoon6mo_lag_AR1_numDF,m1_monsoon6mo_lag_AR2_numDF,	m2_monsoon6mo_lag_AR1_numDF,m2_monsoon6mo_lag_AR2_numDF,m3_monsoon6mo_lag_AR1_numDF,m3_monsoon6mo_lag_AR2_numDF, m1_spring6mo_lag_denDF,m2_spring6mo_lag_denDF,m3_spring6mo_lag_denDF,m1_spring6mo_lag_AR1_denDF,m1_spring6mo_lag_AR2_denDF,m2_spring6mo_lag_AR1_denDF,m2_spring6mo_lag_AR2_denDF,
m3_spring6mo_lag_AR1_denDF,m3_spring6mo_lag_AR2_denDF,m1_monsoon6mo_lag_denDF,
m2_monsoon6mo_lag_denDF,m3_monsoon6mo_lag_denDF,m1_monsoon6mo_lag_AR1_denDF,	
m1_monsoon6mo_lag_AR2_denDF,m2_monsoon6mo_lag_AR1_denDF,m2_monsoon6mo_lag_AR2_denDF,
m3_monsoon6mo_lag_AR1_denDF,	m3_monsoon6mo_lag_AR2_denDF,	m1_spring6mo_lag_F,	m2_spring6mo_lag_F,
m3_spring6mo_lag_F,	m1_spring6mo_lag_AR1_F,	m1_spring6mo_lag_AR2_F,	m2_spring6mo_lag_AR1_F,
m2_spring6mo_lag_AR2_F,	m3_spring6mo_lag_AR1_F,	m3_spring6mo_lag_AR2_F,	m1_monsoon6mo_lag_F,
m2_monsoon6mo_lag_F,	m3_monsoon6mo_lag_F,	m1_monsoon6mo_lag_AR1_F,	m1_monsoon6mo_lag_AR2_F,
m2_monsoon6mo_lag_AR1_F,m2_monsoon6mo_lag_AR2_F,m3_monsoon6mo_lag_AR1_F,	m3_monsoon6mo_lag_AR2_F,
m1_spring6mo_lag_nc,	m2_spring6mo_lag_nc,	m3_spring6mo_lag_nc,	m1_spring6mo_lag_AR1_nc,
m1_spring6mo_lag_AR2_nc,m2_spring6mo_lag_AR1_nc,	m2_spring6mo_lag_AR2_nc,	m3_spring6mo_lag_AR1_nc,
m3_spring6mo_lag_AR2_nc,	m1_monsoon6mo_lag_nc,	m2_monsoon6mo_lag_nc,	m3_monsoon6mo_lag_nc,
m1_monsoon6mo_lag_AR1_nc,	m1_monsoon6mo_lag_AR2_nc,	m2_monsoon6mo_lag_AR1_nc,
m2_monsoon6mo_lag_AR2_nc,	m3_monsoon6mo_lag_AR1_nc,m3_monsoon6mo_lag_AR2_nc,
m1_spring6mo_lag_power,	m2_spring6mo_lag_power,m3_spring6mo_lag_power,m1_spring6mo_lag_AR1_power,
m1_spring6mo_lag_AR2_power,	m2_spring6mo_lag_AR1_power,	m2_spring6mo_lag_AR2_power,
m3_spring6mo_lag_AR1_power,	m3_spring6mo_lag_AR2_power,	m1_monsoon6mo_lag_power,
m2_monsoon6mo_lag_power,	m3_monsoon6mo_lag_power,	m1_monsoon6mo_lag_AR1_power, m1_monsoon6mo_lag_AR2_power, m2_monsoon6mo_lag_AR1_power, m2_monsoon6mo_lag_AR2_power, m3_monsoon6mo_lag_AR1_power,	m3_monsoon6mo_lag_AR2_power)
  
  # append results for each bee species to the output dataset
  beeCSF_output<-rbind(beeCSF_output,output_id)
  
}



##### Chihuahuan Desert Grassland #####

# Create a data frame of just the bee abundance matrix (descriptor variables removed)
speciesMatrix <- black[,9:224] 

# Create a vector of species codes
speciesCodes <- colnames(speciesMatrix)


### Loop through each column of speciesMatrix (each species), running CSF models, running a power analysis for each model, and putting the results in the beeCSF_output matrix

for (i in 1:length(speciesMatrix[1,])) {
  
  # save the species code for column i
  speciesCode <- speciesCodes[i]
  
  # create an object with the name of the ecosystem type
  ecosystem <-"G"
  
  # run mixed effects models
  
  # null model
  m_null<-lme(formula(paste(speciesCode, "~1")), random=~1|transect, data=black, method="ML")
  
  # models including current year's climate
  
  # spring - current year
  m1_spring6mo<-lme(formula(paste(speciesCode, "~spring6SPEI")), random=~1|transect, data=black, method="ML")
  m2_spring6mo<-update(m1_spring6mo, .~. +I(spring6SPEI^2))
  m3_spring6mo<-update(m2_spring6mo, .~. +I(spring6SPEI^3))
  
  m1_spring6mo_AR1 <- update(m1_spring6mo,.~.,correlation=corAR1(form=~year|transect))
  m1_spring6mo_AR2 <- update(m1_spring6mo,.~.,correlation=corARMA(form=~year|transect,p=2))
  
  m2_spring6mo_AR1 <- update(m2_spring6mo,.~.,correlation=corAR1(form=~year|transect))
  m2_spring6mo_AR2 <- update(m2_spring6mo,.~.,correlation=corARMA(form=~year|transect,p=2))
  
  m3_spring6mo_AR1 <- update(m3_spring6mo,.~.,correlation=corAR1(form=~year|transect))
  m3_spring6mo_AR2 <- update(m3_spring6mo,.~.,correlation=corARMA(form=~year|transect,p=2))
  
  # monsoon - current year
  m1_monsoon6mo<-lme(formula(paste(speciesCode, "~monsoon6SPEI")), random=~1|transect, data=black, method="ML")
  m2_monsoon6mo<-update(m1_monsoon6mo, .~. +I(monsoon6SPEI^2))
  m3_monsoon6mo<-update(m2_monsoon6mo, .~. +I(monsoon6SPEI^3))
  
  m1_monsoon6mo_AR1 <- update(m1_monsoon6mo,.~.,correlation=corAR1(form=~year|transect))
  m1_monsoon6mo_AR2 <- update(m1_monsoon6mo,.~.,correlation=corARMA(form=~year|transect,p=2))
  
  m2_monsoon6mo_AR1 <- update(m2_monsoon6mo,.~.,correlation=corAR1(form=~year|transect))
  m2_monsoon6mo_AR2 <- update(m2_monsoon6mo,.~.,correlation=corARMA(form=~year|transect,p=2))
  
  m3_monsoon6mo_AR1 <- update(m3_monsoon6mo,.~.,correlation=corAR1(form=~year|transect))
  m3_monsoon6mo_AR2 <- update(m3_monsoon6mo,.~.,correlation=corARMA(form=~year|transect,p=2))
  
  # models including previous year's climate
  
  # spring - previous year
  m1_spring6mo_lag<-lme(formula(paste(speciesCode, "~spring6SPEI_prioryear")), random=~1|transect, data=black, method="ML")
  m2_spring6mo_lag<-update(m1_spring6mo_lag, .~. +I(spring6SPEI_prioryear^2))
  m3_spring6mo_lag<-update(m2_spring6mo_lag, .~. +I(spring6SPEI_prioryear^3))
  
  m1_spring6mo_lag_AR1 <- update(m1_spring6mo_lag,.~.,correlation=corAR1(form=~year|transect))
  m1_spring6mo_lag_AR2 <- update(m1_spring6mo_lag,.~.,correlation=corARMA(form=~year|transect,p=2))
  
  m2_spring6mo_lag_AR1 <- update(m2_spring6mo_lag,.~.,correlation=corAR1(form=~year|transect))
  m2_spring6mo_lag_AR2 <- update(m2_spring6mo_lag,.~.,correlation=corARMA(form=~year|transect,p=2))
  
  m3_spring6mo_lag_AR1 <- update(m3_spring6mo_lag,.~.,correlation=corAR1(form=~year|transect))
  m3_spring6mo_lag_AR2 <- update(m3_spring6mo_lag,.~.,correlation=corARMA(form=~year|transect,p=2))
  
  # monsoon - previous year
  m1_monsoon6mo_lag<-lme(formula(paste(speciesCode, "~monsoon6SPEI_prioryear")), random=~1|transect, data=black, method="ML")
  m2_monsoon6mo_lag<-update(m1_monsoon6mo_lag, .~. +I(monsoon6SPEI_prioryear^2))
  m3_monsoon6mo_lag<-update(m2_monsoon6mo_lag, .~. +I(monsoon6SPEI_prioryear^3))
  
  m1_monsoon6mo_lag_AR1 <- update(m1_monsoon6mo_lag,.~.,correlation=corAR1(form=~year|transect))
  m1_monsoon6mo_lag_AR2 <- update(m1_monsoon6mo_lag,.~.,correlation=corARMA(form=~year|transect,p=2))
  
  m2_monsoon6mo_lag_AR1 <- update(m2_monsoon6mo_lag,.~.,correlation=corAR1(form=~year|transect))
  m2_monsoon6mo_lag_AR2 <- update(m2_monsoon6mo_lag,.~.,correlation=corARMA(form=~year|transect,p=2))
  
  m3_monsoon6mo_lag_AR1 <- update(m3_monsoon6mo_lag,.~.,correlation=corAR1(form=~year|transect))
  m3_monsoon6mo_lag_AR2 <- update(m3_monsoon6mo_lag,.~.,correlation=corARMA(form=~year|transect,p=2))
  
  # run power analyses for each model, creating objects containing the statistical output
  
  # current year models
  m1_spring6mo_numDF<-Pwr(m1_spring6mo)[2,1]
  m2_spring6mo_numDF<-Pwr(m2_spring6mo)[2,1]
  m3_spring6mo_numDF<-Pwr(m3_spring6mo)[2,1]
  m1_spring6mo_AR1_numDF<-Pwr(m1_spring6mo_AR1)[2,1]
  m1_spring6mo_AR2_numDF<-Pwr(m1_spring6mo_AR2)[2,1]
  m2_spring6mo_AR1_numDF<-Pwr(m2_spring6mo_AR1)[2,1]
  m2_spring6mo_AR2_numDF<-Pwr(m2_spring6mo_AR2)[2,1]
  m3_spring6mo_AR1_numDF<-Pwr(m3_spring6mo_AR1)[2,1]
  m3_spring6mo_AR2_numDF<-Pwr(m3_spring6mo_AR2)[2,1]
  m1_monsoon6mo_numDF<-Pwr(m1_monsoon6mo)[2,1]
  m2_monsoon6mo_numDF<-Pwr(m2_monsoon6mo)[2,1]
  m3_monsoon6mo_numDF<-Pwr(m3_monsoon6mo)[2,1]
  m1_monsoon6mo_AR1_numDF<-Pwr(m1_monsoon6mo_AR1)[2,1]
  m1_monsoon6mo_AR2_numDF<-Pwr(m1_monsoon6mo_AR2)[2,1]
  m2_monsoon6mo_AR1_numDF<-Pwr(m2_monsoon6mo_AR1)[2,1]
  m2_monsoon6mo_AR2_numDF<-Pwr(m2_monsoon6mo_AR2)[2,1]
  m3_monsoon6mo_AR1_numDF<-Pwr(m3_monsoon6mo_AR1)[2,1]
  m3_monsoon6mo_AR2_numDF<-Pwr(m3_monsoon6mo_AR2)[2,1]
  m1_spring6mo_denDF<-Pwr(m1_spring6mo)[2,2]
  m2_spring6mo_denDF<-Pwr(m2_spring6mo)[2,2]
  m3_spring6mo_denDF<-Pwr(m3_spring6mo)[2,2]
  m1_spring6mo_AR1_denDF<-Pwr(m1_spring6mo_AR1)[2,2]
  m1_spring6mo_AR2_denDF<-Pwr(m1_spring6mo_AR2)[2,2]
  m2_spring6mo_AR1_denDF<-Pwr(m2_spring6mo_AR1)[2,2]
  m2_spring6mo_AR2_denDF<-Pwr(m2_spring6mo_AR2)[2,2]
  m3_spring6mo_AR1_denDF<-Pwr(m3_spring6mo_AR1)[2,2]
  m3_spring6mo_AR2_denDF<-Pwr(m3_spring6mo_AR2)[2,2]
  m1_monsoon6mo_denDF<-Pwr(m1_monsoon6mo)[2,2]
  m2_monsoon6mo_denDF<-Pwr(m2_monsoon6mo)[2,2]
  m3_monsoon6mo_denDF<-Pwr(m3_monsoon6mo)[2,2]
  m1_monsoon6mo_AR1_denDF<-Pwr(m1_monsoon6mo_AR1)[2,2]
  m1_monsoon6mo_AR2_denDF<-Pwr(m1_monsoon6mo_AR2)[2,2]
  m2_monsoon6mo_AR1_denDF<-Pwr(m2_monsoon6mo_AR1)[2,2]
  m2_monsoon6mo_AR2_denDF<-Pwr(m2_monsoon6mo_AR2)[2,2]
  m3_monsoon6mo_AR1_denDF<-Pwr(m3_monsoon6mo_AR1)[2,2]
  m3_monsoon6mo_AR2_denDF<-Pwr(m3_monsoon6mo_AR2)[2,2]
  m1_spring6mo_F<-Pwr(m1_spring6mo)[2,3]
  m2_spring6mo_F<-Pwr(m2_spring6mo)[2,3]
  m3_spring6mo_F<-Pwr(m3_spring6mo)[2,3]
  m1_spring6mo_AR1_F<-Pwr(m1_spring6mo_AR1)[2,3]
  m1_spring6mo_AR2_F<-Pwr(m1_spring6mo_AR2)[2,3]
  m2_spring6mo_AR1_F<-Pwr(m2_spring6mo_AR1)[2,3]
  m2_spring6mo_AR2_F<-Pwr(m2_spring6mo_AR2)[2,3]
  m3_spring6mo_AR1_F<-Pwr(m3_spring6mo_AR1)[2,3]
  m3_spring6mo_AR2_F<-Pwr(m3_spring6mo_AR2)[2,3]
  m1_monsoon6mo_F<-Pwr(m1_monsoon6mo)[2,3]
  m2_monsoon6mo_F<-Pwr(m2_monsoon6mo)[2,3]
  m3_monsoon6mo_F<-Pwr(m3_monsoon6mo)[2,3]
  m1_monsoon6mo_AR1_F<-Pwr(m1_monsoon6mo_AR1)[2,3]
  m1_monsoon6mo_AR2_F<-Pwr(m1_monsoon6mo_AR2)[2,3]
  m2_monsoon6mo_AR1_F<-Pwr(m2_monsoon6mo_AR1)[2,3]
  m2_monsoon6mo_AR2_F<-Pwr(m2_monsoon6mo_AR2)[2,3]
  m3_monsoon6mo_AR1_F<-Pwr(m3_monsoon6mo_AR1)[2,3]
  m3_monsoon6mo_AR2_F<-Pwr(m3_monsoon6mo_AR2)[2,3]
  m1_spring6mo_nc<-Pwr(m1_spring6mo)[2,4]
  m2_spring6mo_nc<-Pwr(m2_spring6mo)[2,4]
  m3_spring6mo_nc<-Pwr(m3_spring6mo)[2,4]
  m1_spring6mo_AR1_nc<-Pwr(m1_spring6mo_AR1)[2,4]
  m1_spring6mo_AR2_nc<-Pwr(m1_spring6mo_AR2)[2,4]
  m2_spring6mo_AR1_nc<-Pwr(m2_spring6mo_AR1)[2,4]
  m2_spring6mo_AR2_nc<-Pwr(m2_spring6mo_AR2)[2,4]
  m3_spring6mo_AR1_nc<-Pwr(m3_spring6mo_AR1)[2,4]
  m3_spring6mo_AR2_nc<-Pwr(m3_spring6mo_AR2)[2,4]
  m1_monsoon6mo_nc<-Pwr(m1_monsoon6mo)[2,4]
  m2_monsoon6mo_nc<-Pwr(m2_monsoon6mo)[2,4]
  m3_monsoon6mo_nc<-Pwr(m3_monsoon6mo)[2,4]
  m1_monsoon6mo_AR1_nc<-Pwr(m1_monsoon6mo_AR1)[2,4]
  m1_monsoon6mo_AR2_nc<-Pwr(m1_monsoon6mo_AR2)[2,4]
  m2_monsoon6mo_AR1_nc<-Pwr(m2_monsoon6mo_AR1)[2,4]
  m2_monsoon6mo_AR2_nc<-Pwr(m2_monsoon6mo_AR2)[2,4]
  m3_monsoon6mo_AR1_nc<-Pwr(m3_monsoon6mo_AR1)[2,4]
  m3_monsoon6mo_AR2_nc<-Pwr(m3_monsoon6mo_AR2)[2,4]
  m1_spring6mo_power<-Pwr(m1_spring6mo)[2,5]
  m2_spring6mo_power<-Pwr(m2_spring6mo)[2,5]
  m3_spring6mo_power<-Pwr(m3_spring6mo)[2,5]
  m1_spring6mo_AR1_power<-Pwr(m1_spring6mo_AR1)[2,5]
  m1_spring6mo_AR2_power<-Pwr(m1_spring6mo_AR2)[2,5]
  m2_spring6mo_AR1_power<-Pwr(m2_spring6mo_AR1)[2,5]
  m2_spring6mo_AR2_power<-Pwr(m2_spring6mo_AR2)[2,5]
  m3_spring6mo_AR1_power<-Pwr(m3_spring6mo_AR1)[2,5]
  m3_spring6mo_AR2_power<-Pwr(m3_spring6mo_AR2)[2,5]
  m1_monsoon6mo_power<-Pwr(m1_monsoon6mo)[2,5]
  m2_monsoon6mo_power<-Pwr(m2_monsoon6mo)[2,5]
  m3_monsoon6mo_power<-Pwr(m3_monsoon6mo)[2,5]
  m1_monsoon6mo_AR1_power<-Pwr(m1_monsoon6mo_AR1)[2,5]
  m1_monsoon6mo_AR2_power<-Pwr(m1_monsoon6mo_AR2)[2,5]
  m2_monsoon6mo_AR1_power<-Pwr(m2_monsoon6mo_AR1)[2,5]
  m2_monsoon6mo_AR2_power<-Pwr(m2_monsoon6mo_AR2)[2,5]
  m3_monsoon6mo_AR1_power<-Pwr(m3_monsoon6mo_AR1)[2,5]
  m3_monsoon6mo_AR2_power<-Pwr(m3_monsoon6mo_AR2)[2,5]
  
  # previous year models
  m1_spring6mo_lag_numDF<-Pwr(m1_spring6mo_lag)[2,1]
  m2_spring6mo_lag_numDF<-Pwr(m2_spring6mo_lag)[2,1]
  m3_spring6mo_lag_numDF<-Pwr(m3_spring6mo_lag)[2,1]
  m1_spring6mo_lag_AR1_numDF<-Pwr(m1_spring6mo_lag_AR1)[2,1]
  m1_spring6mo_lag_AR2_numDF<-Pwr(m1_spring6mo_lag_AR2)[2,1]
  m2_spring6mo_lag_AR1_numDF<-Pwr(m2_spring6mo_lag_AR1)[2,1]
  m2_spring6mo_lag_AR2_numDF<-Pwr(m2_spring6mo_lag_AR2)[2,1]
  m3_spring6mo_lag_AR1_numDF<-Pwr(m3_spring6mo_lag_AR1)[2,1]
  m3_spring6mo_lag_AR2_numDF<-Pwr(m3_spring6mo_lag_AR2)[2,1]
  m1_monsoon6mo_lag_numDF<-Pwr(m1_monsoon6mo_lag)[2,1]
  m2_monsoon6mo_lag_numDF<-Pwr(m2_monsoon6mo_lag)[2,1]
  m3_monsoon6mo_lag_numDF<-Pwr(m3_monsoon6mo_lag)[2,1]
  m1_monsoon6mo_lag_AR1_numDF<-Pwr(m1_monsoon6mo_lag_AR1)[2,1]
  m1_monsoon6mo_lag_AR2_numDF<-Pwr(m1_monsoon6mo_lag_AR2)[2,1]
  m2_monsoon6mo_lag_AR1_numDF<-Pwr(m2_monsoon6mo_lag_AR1)[2,1]
  m2_monsoon6mo_lag_AR2_numDF<-Pwr(m2_monsoon6mo_lag_AR2)[2,1]
  m3_monsoon6mo_lag_AR1_numDF<-Pwr(m3_monsoon6mo_lag_AR1)[2,1]
  m3_monsoon6mo_lag_AR2_numDF<-Pwr(m3_monsoon6mo_lag_AR2)[2,1]
  m1_spring6mo_lag_denDF<-Pwr(m1_spring6mo_lag)[2,2]
  m2_spring6mo_lag_denDF<-Pwr(m2_spring6mo_lag)[2,2]
  m3_spring6mo_lag_denDF<-Pwr(m3_spring6mo_lag)[2,2]
  m1_spring6mo_lag_AR1_denDF<-Pwr(m1_spring6mo_lag_AR1)[2,2]
  m1_spring6mo_lag_AR2_denDF<-Pwr(m1_spring6mo_lag_AR2)[2,2]
  m2_spring6mo_lag_AR1_denDF<-Pwr(m2_spring6mo_lag_AR1)[2,2]
  m2_spring6mo_lag_AR2_denDF<-Pwr(m2_spring6mo_lag_AR2)[2,2]
  m3_spring6mo_lag_AR1_denDF<-Pwr(m3_spring6mo_lag_AR1)[2,2]
  m3_spring6mo_lag_AR2_denDF<-Pwr(m3_spring6mo_lag_AR2)[2,2]
  m1_monsoon6mo_lag_denDF<-Pwr(m1_monsoon6mo_lag)[2,2]
  m2_monsoon6mo_lag_denDF<-Pwr(m2_monsoon6mo_lag)[2,2]
  m3_monsoon6mo_lag_denDF<-Pwr(m3_monsoon6mo_lag)[2,2]
  m1_monsoon6mo_lag_AR1_denDF<-Pwr(m1_monsoon6mo_lag_AR1)[2,2]
  m1_monsoon6mo_lag_AR2_denDF<-Pwr(m1_monsoon6mo_lag_AR2)[2,2]
  m2_monsoon6mo_lag_AR1_denDF<-Pwr(m2_monsoon6mo_lag_AR1)[2,2]
  m2_monsoon6mo_lag_AR2_denDF<-Pwr(m2_monsoon6mo_lag_AR2)[2,2]
  m3_monsoon6mo_lag_AR1_denDF<-Pwr(m3_monsoon6mo_lag_AR1)[2,2]
  m3_monsoon6mo_lag_AR2_denDF<-Pwr(m3_monsoon6mo_lag_AR2)[2,2]
  m1_spring6mo_lag_F<-Pwr(m1_spring6mo_lag)[2,3]
  m2_spring6mo_lag_F<-Pwr(m2_spring6mo_lag)[2,3]
  m3_spring6mo_lag_F<-Pwr(m3_spring6mo_lag)[2,3]
  m1_spring6mo_lag_AR1_F<-Pwr(m1_spring6mo_lag_AR1)[2,3]
  m1_spring6mo_lag_AR2_F<-Pwr(m1_spring6mo_lag_AR2)[2,3]
  m2_spring6mo_lag_AR1_F<-Pwr(m2_spring6mo_lag_AR1)[2,3]
  m2_spring6mo_lag_AR2_F<-Pwr(m2_spring6mo_lag_AR2)[2,3]
  m3_spring6mo_lag_AR1_F<-Pwr(m3_spring6mo_lag_AR1)[2,3]
  m3_spring6mo_lag_AR2_F<-Pwr(m3_spring6mo_lag_AR2)[2,3]
  m1_monsoon6mo_lag_F<-Pwr(m1_monsoon6mo_lag)[2,3]
  m2_monsoon6mo_lag_F<-Pwr(m2_monsoon6mo_lag)[2,3]
  m3_monsoon6mo_lag_F<-Pwr(m3_monsoon6mo_lag)[2,3]
  m1_monsoon6mo_lag_AR1_F<-Pwr(m1_monsoon6mo_lag_AR1)[2,3]
  m1_monsoon6mo_lag_AR2_F<-Pwr(m1_monsoon6mo_lag_AR2)[2,3]
  m2_monsoon6mo_lag_AR1_F<-Pwr(m2_monsoon6mo_lag_AR1)[2,3]
  m2_monsoon6mo_lag_AR2_F<-Pwr(m2_monsoon6mo_lag_AR2)[2,3]
  m3_monsoon6mo_lag_AR1_F<-Pwr(m3_monsoon6mo_lag_AR1)[2,3]
  m3_monsoon6mo_lag_AR2_F<-Pwr(m3_monsoon6mo_lag_AR2)[2,3]
  m1_spring6mo_lag_nc<-Pwr(m1_spring6mo_lag)[2,4]
  m2_spring6mo_lag_nc<-Pwr(m2_spring6mo_lag)[2,4]
  m3_spring6mo_lag_nc<-Pwr(m3_spring6mo_lag)[2,4]
  m1_spring6mo_lag_AR1_nc<-Pwr(m1_spring6mo_lag_AR1)[2,4]
  m1_spring6mo_lag_AR2_nc<-Pwr(m1_spring6mo_lag_AR2)[2,4]
  m2_spring6mo_lag_AR1_nc<-Pwr(m2_spring6mo_lag_AR1)[2,4]
  m2_spring6mo_lag_AR2_nc<-Pwr(m2_spring6mo_lag_AR2)[2,4]
  m3_spring6mo_lag_AR1_nc<-Pwr(m3_spring6mo_lag_AR1)[2,4]
  m3_spring6mo_lag_AR2_nc<-Pwr(m3_spring6mo_lag_AR2)[2,4]
  m1_monsoon6mo_lag_nc<-Pwr(m1_monsoon6mo_lag)[2,4]
  m2_monsoon6mo_lag_nc<-Pwr(m2_monsoon6mo_lag)[2,4]
  m3_monsoon6mo_lag_nc<-Pwr(m3_monsoon6mo_lag)[2,4]
  m1_monsoon6mo_lag_AR1_nc<-Pwr(m1_monsoon6mo_lag_AR1)[2,4]
  m1_monsoon6mo_lag_AR2_nc<-Pwr(m1_monsoon6mo_lag_AR2)[2,4]
  m2_monsoon6mo_lag_AR1_nc<-Pwr(m2_monsoon6mo_lag_AR1)[2,4]
  m2_monsoon6mo_lag_AR2_nc<-Pwr(m2_monsoon6mo_lag_AR2)[2,4]
  m3_monsoon6mo_lag_AR1_nc<-Pwr(m3_monsoon6mo_lag_AR1)[2,4]
  m3_monsoon6mo_lag_AR2_nc<-Pwr(m3_monsoon6mo_lag_AR2)[2,4]
  m1_spring6mo_lag_power<-Pwr(m1_spring6mo_lag)[2,5]
  m2_spring6mo_lag_power<-Pwr(m2_spring6mo_lag)[2,5]
  m3_spring6mo_lag_power<-Pwr(m3_spring6mo_lag)[2,5]
  m1_spring6mo_lag_AR1_power<-Pwr(m1_spring6mo_lag_AR1)[2,5]
  m1_spring6mo_lag_AR2_power<-Pwr(m1_spring6mo_lag_AR2)[2,5]
  m2_spring6mo_lag_AR1_power<-Pwr(m2_spring6mo_lag_AR1)[2,5]
  m2_spring6mo_lag_AR2_power<-Pwr(m2_spring6mo_lag_AR2)[2,5]
  m3_spring6mo_lag_AR1_power<-Pwr(m3_spring6mo_lag_AR1)[2,5]
  m3_spring6mo_lag_AR2_power<-Pwr(m3_spring6mo_lag_AR2)[2,5]
  m1_monsoon6mo_lag_power<-Pwr(m1_monsoon6mo_lag)[2,5]
  m2_monsoon6mo_lag_power<-Pwr(m2_monsoon6mo_lag)[2,5]
  m3_monsoon6mo_lag_power<-Pwr(m3_monsoon6mo_lag)[2,5]
  m1_monsoon6mo_lag_AR1_power<-Pwr(m1_monsoon6mo_lag_AR1)[2,5]
  m1_monsoon6mo_lag_AR2_power<-Pwr(m1_monsoon6mo_lag_AR2)[2,5]
  m2_monsoon6mo_lag_AR1_power<-Pwr(m2_monsoon6mo_lag_AR1)[2,5]
  m2_monsoon6mo_lag_AR2_power<-Pwr(m2_monsoon6mo_lag_AR2)[2,5]
  m3_monsoon6mo_lag_AR1_power<-Pwr(m3_monsoon6mo_lag_AR1)[2,5]
  m3_monsoon6mo_lag_AR2_power<-Pwr(m3_monsoon6mo_lag_AR2)[2,5]
  
  # bind all target output values together
  output_id<-cbind(speciesCode, ecosystem, m1_spring6mo_numDF,	m2_spring6mo_numDF,	m3_spring6mo_numDF,	m1_spring6mo_AR1_numDF,	m1_spring6mo_AR2_numDF,	m2_spring6mo_AR1_numDF,	m2_spring6mo_AR2_numDF,	m3_spring6mo_AR1_numDF,	m3_spring6mo_AR2_numDF,	m1_monsoon6mo_numDF,	m2_monsoon6mo_numDF,	m3_monsoon6mo_numDF,	m1_monsoon6mo_AR1_numDF,	m1_monsoon6mo_AR2_numDF,	m2_monsoon6mo_AR1_numDF,	m2_monsoon6mo_AR2_numDF,	m3_monsoon6mo_AR1_numDF,	m3_monsoon6mo_AR2_numDF,	m1_spring6mo_denDF,	m2_spring6mo_denDF,	m3_spring6mo_denDF,	m1_spring6mo_AR1_denDF,	m1_spring6mo_AR2_denDF,	m2_spring6mo_AR1_denDF,	m2_spring6mo_AR2_denDF,	m3_spring6mo_AR1_denDF,	m3_spring6mo_AR2_denDF,	m1_monsoon6mo_denDF,	m2_monsoon6mo_denDF,	m3_monsoon6mo_denDF,	m1_monsoon6mo_AR1_denDF,	m1_monsoon6mo_AR2_denDF,	m2_monsoon6mo_AR1_denDF,	m2_monsoon6mo_AR2_denDF,	m3_monsoon6mo_AR1_denDF,	m3_monsoon6mo_AR2_denDF,	m1_spring6mo_F,	m2_spring6mo_F,	m3_spring6mo_F,	m1_spring6mo_AR1_F,	m1_spring6mo_AR2_F,	m2_spring6mo_AR1_F,	m2_spring6mo_AR2_F,	m3_spring6mo_AR1_F,	m3_spring6mo_AR2_F,	m1_monsoon6mo_F,	m2_monsoon6mo_F,	m3_monsoon6mo_F,	m1_monsoon6mo_AR1_F,	m1_monsoon6mo_AR2_F,	m2_monsoon6mo_AR1_F,	m2_monsoon6mo_AR2_F,	m3_monsoon6mo_AR1_F,	m3_monsoon6mo_AR2_F,	m1_spring6mo_nc,	m2_spring6mo_nc,	m3_spring6mo_nc,	m1_spring6mo_AR1_nc,	m1_spring6mo_AR2_nc,	m2_spring6mo_AR1_nc,	m2_spring6mo_AR2_nc,	m3_spring6mo_AR1_nc,	m3_spring6mo_AR2_nc,	m1_monsoon6mo_nc,	m2_monsoon6mo_nc,	m3_monsoon6mo_nc,	m1_monsoon6mo_AR1_nc,	m1_monsoon6mo_AR2_nc,	m2_monsoon6mo_AR1_nc,	m2_monsoon6mo_AR2_nc,	m3_monsoon6mo_AR1_nc,	m3_monsoon6mo_AR2_nc,	m1_spring6mo_power,	m2_spring6mo_power,	m3_spring6mo_power,	m1_spring6mo_AR1_power,	m1_spring6mo_AR2_power,	m2_spring6mo_AR1_power,	m2_spring6mo_AR2_power,	m3_spring6mo_AR1_power,	m3_spring6mo_AR2_power,	m1_monsoon6mo_power,	m2_monsoon6mo_power,	m3_monsoon6mo_power,	m1_monsoon6mo_AR1_power,	m1_monsoon6mo_AR2_power,	m2_monsoon6mo_AR1_power,	m2_monsoon6mo_AR2_power,	m3_monsoon6mo_AR1_power,	m3_monsoon6mo_AR2_power, m1_spring6mo_lag_numDF,m2_spring6mo_lag_numDF,m3_spring6mo_lag_numDF,m1_spring6mo_lag_AR1_numDF,	m1_spring6mo_lag_AR2_numDF,m2_spring6mo_lag_AR1_numDF,m2_spring6mo_lag_AR2_numDF,m3_spring6mo_lag_AR1_numDF,m3_spring6mo_lag_AR2_numDF,m1_monsoon6mo_lag_numDF,m2_monsoon6mo_lag_numDF,m3_monsoon6mo_lag_numDF,m1_monsoon6mo_lag_AR1_numDF,m1_monsoon6mo_lag_AR2_numDF,	m2_monsoon6mo_lag_AR1_numDF,m2_monsoon6mo_lag_AR2_numDF,m3_monsoon6mo_lag_AR1_numDF,m3_monsoon6mo_lag_AR2_numDF, m1_spring6mo_lag_denDF,m2_spring6mo_lag_denDF,m3_spring6mo_lag_denDF,m1_spring6mo_lag_AR1_denDF,m1_spring6mo_lag_AR2_denDF,m2_spring6mo_lag_AR1_denDF,m2_spring6mo_lag_AR2_denDF,
                   m3_spring6mo_lag_AR1_denDF,m3_spring6mo_lag_AR2_denDF,m1_monsoon6mo_lag_denDF,
                   m2_monsoon6mo_lag_denDF,m3_monsoon6mo_lag_denDF,m1_monsoon6mo_lag_AR1_denDF,	
                   m1_monsoon6mo_lag_AR2_denDF,m2_monsoon6mo_lag_AR1_denDF,m2_monsoon6mo_lag_AR2_denDF,
                   m3_monsoon6mo_lag_AR1_denDF,	m3_monsoon6mo_lag_AR2_denDF,	m1_spring6mo_lag_F,	m2_spring6mo_lag_F,
                   m3_spring6mo_lag_F,	m1_spring6mo_lag_AR1_F,	m1_spring6mo_lag_AR2_F,	m2_spring6mo_lag_AR1_F,
                   m2_spring6mo_lag_AR2_F,	m3_spring6mo_lag_AR1_F,	m3_spring6mo_lag_AR2_F,	m1_monsoon6mo_lag_F,
                   m2_monsoon6mo_lag_F,	m3_monsoon6mo_lag_F,	m1_monsoon6mo_lag_AR1_F,	m1_monsoon6mo_lag_AR2_F,
                   m2_monsoon6mo_lag_AR1_F,m2_monsoon6mo_lag_AR2_F,m3_monsoon6mo_lag_AR1_F,	m3_monsoon6mo_lag_AR2_F,
                   m1_spring6mo_lag_nc,	m2_spring6mo_lag_nc,	m3_spring6mo_lag_nc,	m1_spring6mo_lag_AR1_nc,
                   m1_spring6mo_lag_AR2_nc,m2_spring6mo_lag_AR1_nc,	m2_spring6mo_lag_AR2_nc,	m3_spring6mo_lag_AR1_nc,
                   m3_spring6mo_lag_AR2_nc,	m1_monsoon6mo_lag_nc,	m2_monsoon6mo_lag_nc,	m3_monsoon6mo_lag_nc,
                   m1_monsoon6mo_lag_AR1_nc,	m1_monsoon6mo_lag_AR2_nc,	m2_monsoon6mo_lag_AR1_nc,
                   m2_monsoon6mo_lag_AR2_nc,	m3_monsoon6mo_lag_AR1_nc,m3_monsoon6mo_lag_AR2_nc,
                   m1_spring6mo_lag_power,	m2_spring6mo_lag_power,m3_spring6mo_lag_power,m1_spring6mo_lag_AR1_power,
                   m1_spring6mo_lag_AR2_power,	m2_spring6mo_lag_AR1_power,	m2_spring6mo_lag_AR2_power,
                   m3_spring6mo_lag_AR1_power,	m3_spring6mo_lag_AR2_power,	m1_monsoon6mo_lag_power,
                   m2_monsoon6mo_lag_power,	m3_monsoon6mo_lag_power,	m1_monsoon6mo_lag_AR1_power, m1_monsoon6mo_lag_AR2_power, m2_monsoon6mo_lag_AR1_power, m2_monsoon6mo_lag_AR2_power, m3_monsoon6mo_lag_AR1_power,	m3_monsoon6mo_lag_AR2_power)
  
  # append results for each bee species to the output dataset
  beeCSF_output<-rbind(beeCSF_output,output_id)
  
}



##### Plains Grassland #####

# Create a data frame of just the bee abundance matrix (descriptor variables removed)
speciesMatrix <- blue[,9:232] 

# Create a vector of species codes
speciesCodes <- colnames(speciesMatrix)


### Loop through each column of speciesMatrix (each species), running CSF models, running a power analysis for each model, and putting the results in the beeCSF_output matrix

for (i in 1:length(speciesMatrix[1,])) {
  
  # save the species code for column i
  speciesCode <- speciesCodes[i]
  
  # create an object with the name of the ecosystem type
  ecosystem <-"B"
  
  # run mixed effects models
  
  # null model
  m_null<-lme(formula(paste(speciesCode, "~1")), random=~1|transect, data=blue, method="ML")
  
  # models including current year's climate
  
  # spring - current year
  m1_spring6mo<-lme(formula(paste(speciesCode, "~spring6SPEI")), random=~1|transect, data=blue, method="ML")
  m2_spring6mo<-update(m1_spring6mo, .~. +I(spring6SPEI^2))
  m3_spring6mo<-update(m2_spring6mo, .~. +I(spring6SPEI^3))
  
  m1_spring6mo_AR1 <- update(m1_spring6mo,.~.,correlation=corAR1(form=~year|transect))
  m1_spring6mo_AR2 <- update(m1_spring6mo,.~.,correlation=corARMA(form=~year|transect,p=2))
  
  m2_spring6mo_AR1 <- update(m2_spring6mo,.~.,correlation=corAR1(form=~year|transect))
  m2_spring6mo_AR2 <- update(m2_spring6mo,.~.,correlation=corARMA(form=~year|transect,p=2))
  
  m3_spring6mo_AR1 <- update(m3_spring6mo,.~.,correlation=corAR1(form=~year|transect))
  m3_spring6mo_AR2 <- update(m3_spring6mo,.~.,correlation=corARMA(form=~year|transect,p=2))
  
  # monsoon - current year
  m1_monsoon6mo<-lme(formula(paste(speciesCode, "~monsoon6SPEI")), random=~1|transect, data=blue, method="ML")
  m2_monsoon6mo<-update(m1_monsoon6mo, .~. +I(monsoon6SPEI^2))
  m3_monsoon6mo<-update(m2_monsoon6mo, .~. +I(monsoon6SPEI^3))
  
  m1_monsoon6mo_AR1 <- update(m1_monsoon6mo,.~.,correlation=corAR1(form=~year|transect))
  m1_monsoon6mo_AR2 <- update(m1_monsoon6mo,.~.,correlation=corARMA(form=~year|transect,p=2))
  
  m2_monsoon6mo_AR1 <- update(m2_monsoon6mo,.~.,correlation=corAR1(form=~year|transect))
  m2_monsoon6mo_AR2 <- update(m2_monsoon6mo,.~.,correlation=corARMA(form=~year|transect,p=2))
  
  m3_monsoon6mo_AR1 <- update(m3_monsoon6mo,.~.,correlation=corAR1(form=~year|transect))
  m3_monsoon6mo_AR2 <- update(m3_monsoon6mo,.~.,correlation=corARMA(form=~year|transect,p=2))
  
  # models including previous year's climate
  
  # spring - previous year
  m1_spring6mo_lag<-lme(formula(paste(speciesCode, "~spring6SPEI_prioryear")), random=~1|transect, data=blue, method="ML")
  m2_spring6mo_lag<-update(m1_spring6mo_lag, .~. +I(spring6SPEI_prioryear^2))
  m3_spring6mo_lag<-update(m2_spring6mo_lag, .~. +I(spring6SPEI_prioryear^3))
  
  m1_spring6mo_lag_AR1 <- update(m1_spring6mo_lag,.~.,correlation=corAR1(form=~year|transect))
  m1_spring6mo_lag_AR2 <- update(m1_spring6mo_lag,.~.,correlation=corARMA(form=~year|transect,p=2))
  
  m2_spring6mo_lag_AR1 <- update(m2_spring6mo_lag,.~.,correlation=corAR1(form=~year|transect))
  m2_spring6mo_lag_AR2 <- update(m2_spring6mo_lag,.~.,correlation=corARMA(form=~year|transect,p=2))
  
  m3_spring6mo_lag_AR1 <- update(m3_spring6mo_lag,.~.,correlation=corAR1(form=~year|transect))
  m3_spring6mo_lag_AR2 <- update(m3_spring6mo_lag,.~.,correlation=corARMA(form=~year|transect,p=2))
  
  # monsoon - previous year
  m1_monsoon6mo_lag<-lme(formula(paste(speciesCode, "~monsoon6SPEI_prioryear")), random=~1|transect, data=blue, method="ML")
  m2_monsoon6mo_lag<-update(m1_monsoon6mo_lag, .~. +I(monsoon6SPEI_prioryear^2))
  m3_monsoon6mo_lag<-update(m2_monsoon6mo_lag, .~. +I(monsoon6SPEI_prioryear^3))
  
  m1_monsoon6mo_lag_AR1 <- update(m1_monsoon6mo_lag,.~.,correlation=corAR1(form=~year|transect))
  m1_monsoon6mo_lag_AR2 <- update(m1_monsoon6mo_lag,.~.,correlation=corARMA(form=~year|transect,p=2))
  
  m2_monsoon6mo_lag_AR1 <- update(m2_monsoon6mo_lag,.~.,correlation=corAR1(form=~year|transect))
  m2_monsoon6mo_lag_AR2 <- update(m2_monsoon6mo_lag,.~.,correlation=corARMA(form=~year|transect,p=2))
  
  m3_monsoon6mo_lag_AR1 <- update(m3_monsoon6mo_lag,.~.,correlation=corAR1(form=~year|transect))
  m3_monsoon6mo_lag_AR2 <- update(m3_monsoon6mo_lag,.~.,correlation=corARMA(form=~year|transect,p=2))
  
  # run power analyses for each model, creating objects containing the statistical output
  
  # current year models
  m1_spring6mo_numDF<-Pwr(m1_spring6mo)[2,1]
  m2_spring6mo_numDF<-Pwr(m2_spring6mo)[2,1]
  m3_spring6mo_numDF<-Pwr(m3_spring6mo)[2,1]
  m1_spring6mo_AR1_numDF<-Pwr(m1_spring6mo_AR1)[2,1]
  m1_spring6mo_AR2_numDF<-Pwr(m1_spring6mo_AR2)[2,1]
  m2_spring6mo_AR1_numDF<-Pwr(m2_spring6mo_AR1)[2,1]
  m2_spring6mo_AR2_numDF<-Pwr(m2_spring6mo_AR2)[2,1]
  m3_spring6mo_AR1_numDF<-Pwr(m3_spring6mo_AR1)[2,1]
  m3_spring6mo_AR2_numDF<-Pwr(m3_spring6mo_AR2)[2,1]
  m1_monsoon6mo_numDF<-Pwr(m1_monsoon6mo)[2,1]
  m2_monsoon6mo_numDF<-Pwr(m2_monsoon6mo)[2,1]
  m3_monsoon6mo_numDF<-Pwr(m3_monsoon6mo)[2,1]
  m1_monsoon6mo_AR1_numDF<-Pwr(m1_monsoon6mo_AR1)[2,1]
  m1_monsoon6mo_AR2_numDF<-Pwr(m1_monsoon6mo_AR2)[2,1]
  m2_monsoon6mo_AR1_numDF<-Pwr(m2_monsoon6mo_AR1)[2,1]
  m2_monsoon6mo_AR2_numDF<-Pwr(m2_monsoon6mo_AR2)[2,1]
  m3_monsoon6mo_AR1_numDF<-Pwr(m3_monsoon6mo_AR1)[2,1]
  m3_monsoon6mo_AR2_numDF<-Pwr(m3_monsoon6mo_AR2)[2,1]
  m1_spring6mo_denDF<-Pwr(m1_spring6mo)[2,2]
  m2_spring6mo_denDF<-Pwr(m2_spring6mo)[2,2]
  m3_spring6mo_denDF<-Pwr(m3_spring6mo)[2,2]
  m1_spring6mo_AR1_denDF<-Pwr(m1_spring6mo_AR1)[2,2]
  m1_spring6mo_AR2_denDF<-Pwr(m1_spring6mo_AR2)[2,2]
  m2_spring6mo_AR1_denDF<-Pwr(m2_spring6mo_AR1)[2,2]
  m2_spring6mo_AR2_denDF<-Pwr(m2_spring6mo_AR2)[2,2]
  m3_spring6mo_AR1_denDF<-Pwr(m3_spring6mo_AR1)[2,2]
  m3_spring6mo_AR2_denDF<-Pwr(m3_spring6mo_AR2)[2,2]
  m1_monsoon6mo_denDF<-Pwr(m1_monsoon6mo)[2,2]
  m2_monsoon6mo_denDF<-Pwr(m2_monsoon6mo)[2,2]
  m3_monsoon6mo_denDF<-Pwr(m3_monsoon6mo)[2,2]
  m1_monsoon6mo_AR1_denDF<-Pwr(m1_monsoon6mo_AR1)[2,2]
  m1_monsoon6mo_AR2_denDF<-Pwr(m1_monsoon6mo_AR2)[2,2]
  m2_monsoon6mo_AR1_denDF<-Pwr(m2_monsoon6mo_AR1)[2,2]
  m2_monsoon6mo_AR2_denDF<-Pwr(m2_monsoon6mo_AR2)[2,2]
  m3_monsoon6mo_AR1_denDF<-Pwr(m3_monsoon6mo_AR1)[2,2]
  m3_monsoon6mo_AR2_denDF<-Pwr(m3_monsoon6mo_AR2)[2,2]
  m1_spring6mo_F<-Pwr(m1_spring6mo)[2,3]
  m2_spring6mo_F<-Pwr(m2_spring6mo)[2,3]
  m3_spring6mo_F<-Pwr(m3_spring6mo)[2,3]
  m1_spring6mo_AR1_F<-Pwr(m1_spring6mo_AR1)[2,3]
  m1_spring6mo_AR2_F<-Pwr(m1_spring6mo_AR2)[2,3]
  m2_spring6mo_AR1_F<-Pwr(m2_spring6mo_AR1)[2,3]
  m2_spring6mo_AR2_F<-Pwr(m2_spring6mo_AR2)[2,3]
  m3_spring6mo_AR1_F<-Pwr(m3_spring6mo_AR1)[2,3]
  m3_spring6mo_AR2_F<-Pwr(m3_spring6mo_AR2)[2,3]
  m1_monsoon6mo_F<-Pwr(m1_monsoon6mo)[2,3]
  m2_monsoon6mo_F<-Pwr(m2_monsoon6mo)[2,3]
  m3_monsoon6mo_F<-Pwr(m3_monsoon6mo)[2,3]
  m1_monsoon6mo_AR1_F<-Pwr(m1_monsoon6mo_AR1)[2,3]
  m1_monsoon6mo_AR2_F<-Pwr(m1_monsoon6mo_AR2)[2,3]
  m2_monsoon6mo_AR1_F<-Pwr(m2_monsoon6mo_AR1)[2,3]
  m2_monsoon6mo_AR2_F<-Pwr(m2_monsoon6mo_AR2)[2,3]
  m3_monsoon6mo_AR1_F<-Pwr(m3_monsoon6mo_AR1)[2,3]
  m3_monsoon6mo_AR2_F<-Pwr(m3_monsoon6mo_AR2)[2,3]
  m1_spring6mo_nc<-Pwr(m1_spring6mo)[2,4]
  m2_spring6mo_nc<-Pwr(m2_spring6mo)[2,4]
  m3_spring6mo_nc<-Pwr(m3_spring6mo)[2,4]
  m1_spring6mo_AR1_nc<-Pwr(m1_spring6mo_AR1)[2,4]
  m1_spring6mo_AR2_nc<-Pwr(m1_spring6mo_AR2)[2,4]
  m2_spring6mo_AR1_nc<-Pwr(m2_spring6mo_AR1)[2,4]
  m2_spring6mo_AR2_nc<-Pwr(m2_spring6mo_AR2)[2,4]
  m3_spring6mo_AR1_nc<-Pwr(m3_spring6mo_AR1)[2,4]
  m3_spring6mo_AR2_nc<-Pwr(m3_spring6mo_AR2)[2,4]
  m1_monsoon6mo_nc<-Pwr(m1_monsoon6mo)[2,4]
  m2_monsoon6mo_nc<-Pwr(m2_monsoon6mo)[2,4]
  m3_monsoon6mo_nc<-Pwr(m3_monsoon6mo)[2,4]
  m1_monsoon6mo_AR1_nc<-Pwr(m1_monsoon6mo_AR1)[2,4]
  m1_monsoon6mo_AR2_nc<-Pwr(m1_monsoon6mo_AR2)[2,4]
  m2_monsoon6mo_AR1_nc<-Pwr(m2_monsoon6mo_AR1)[2,4]
  m2_monsoon6mo_AR2_nc<-Pwr(m2_monsoon6mo_AR2)[2,4]
  m3_monsoon6mo_AR1_nc<-Pwr(m3_monsoon6mo_AR1)[2,4]
  m3_monsoon6mo_AR2_nc<-Pwr(m3_monsoon6mo_AR2)[2,4]
  m1_spring6mo_power<-Pwr(m1_spring6mo)[2,5]
  m2_spring6mo_power<-Pwr(m2_spring6mo)[2,5]
  m3_spring6mo_power<-Pwr(m3_spring6mo)[2,5]
  m1_spring6mo_AR1_power<-Pwr(m1_spring6mo_AR1)[2,5]
  m1_spring6mo_AR2_power<-Pwr(m1_spring6mo_AR2)[2,5]
  m2_spring6mo_AR1_power<-Pwr(m2_spring6mo_AR1)[2,5]
  m2_spring6mo_AR2_power<-Pwr(m2_spring6mo_AR2)[2,5]
  m3_spring6mo_AR1_power<-Pwr(m3_spring6mo_AR1)[2,5]
  m3_spring6mo_AR2_power<-Pwr(m3_spring6mo_AR2)[2,5]
  m1_monsoon6mo_power<-Pwr(m1_monsoon6mo)[2,5]
  m2_monsoon6mo_power<-Pwr(m2_monsoon6mo)[2,5]
  m3_monsoon6mo_power<-Pwr(m3_monsoon6mo)[2,5]
  m1_monsoon6mo_AR1_power<-Pwr(m1_monsoon6mo_AR1)[2,5]
  m1_monsoon6mo_AR2_power<-Pwr(m1_monsoon6mo_AR2)[2,5]
  m2_monsoon6mo_AR1_power<-Pwr(m2_monsoon6mo_AR1)[2,5]
  m2_monsoon6mo_AR2_power<-Pwr(m2_monsoon6mo_AR2)[2,5]
  m3_monsoon6mo_AR1_power<-Pwr(m3_monsoon6mo_AR1)[2,5]
  m3_monsoon6mo_AR2_power<-Pwr(m3_monsoon6mo_AR2)[2,5]
  
  # previous year models
  m1_spring6mo_lag_numDF<-Pwr(m1_spring6mo_lag)[2,1]
  m2_spring6mo_lag_numDF<-Pwr(m2_spring6mo_lag)[2,1]
  m3_spring6mo_lag_numDF<-Pwr(m3_spring6mo_lag)[2,1]
  m1_spring6mo_lag_AR1_numDF<-Pwr(m1_spring6mo_lag_AR1)[2,1]
  m1_spring6mo_lag_AR2_numDF<-Pwr(m1_spring6mo_lag_AR2)[2,1]
  m2_spring6mo_lag_AR1_numDF<-Pwr(m2_spring6mo_lag_AR1)[2,1]
  m2_spring6mo_lag_AR2_numDF<-Pwr(m2_spring6mo_lag_AR2)[2,1]
  m3_spring6mo_lag_AR1_numDF<-Pwr(m3_spring6mo_lag_AR1)[2,1]
  m3_spring6mo_lag_AR2_numDF<-Pwr(m3_spring6mo_lag_AR2)[2,1]
  m1_monsoon6mo_lag_numDF<-Pwr(m1_monsoon6mo_lag)[2,1]
  m2_monsoon6mo_lag_numDF<-Pwr(m2_monsoon6mo_lag)[2,1]
  m3_monsoon6mo_lag_numDF<-Pwr(m3_monsoon6mo_lag)[2,1]
  m1_monsoon6mo_lag_AR1_numDF<-Pwr(m1_monsoon6mo_lag_AR1)[2,1]
  m1_monsoon6mo_lag_AR2_numDF<-Pwr(m1_monsoon6mo_lag_AR2)[2,1]
  m2_monsoon6mo_lag_AR1_numDF<-Pwr(m2_monsoon6mo_lag_AR1)[2,1]
  m2_monsoon6mo_lag_AR2_numDF<-Pwr(m2_monsoon6mo_lag_AR2)[2,1]
  m3_monsoon6mo_lag_AR1_numDF<-Pwr(m3_monsoon6mo_lag_AR1)[2,1]
  m3_monsoon6mo_lag_AR2_numDF<-Pwr(m3_monsoon6mo_lag_AR2)[2,1]
  m1_spring6mo_lag_denDF<-Pwr(m1_spring6mo_lag)[2,2]
  m2_spring6mo_lag_denDF<-Pwr(m2_spring6mo_lag)[2,2]
  m3_spring6mo_lag_denDF<-Pwr(m3_spring6mo_lag)[2,2]
  m1_spring6mo_lag_AR1_denDF<-Pwr(m1_spring6mo_lag_AR1)[2,2]
  m1_spring6mo_lag_AR2_denDF<-Pwr(m1_spring6mo_lag_AR2)[2,2]
  m2_spring6mo_lag_AR1_denDF<-Pwr(m2_spring6mo_lag_AR1)[2,2]
  m2_spring6mo_lag_AR2_denDF<-Pwr(m2_spring6mo_lag_AR2)[2,2]
  m3_spring6mo_lag_AR1_denDF<-Pwr(m3_spring6mo_lag_AR1)[2,2]
  m3_spring6mo_lag_AR2_denDF<-Pwr(m3_spring6mo_lag_AR2)[2,2]
  m1_monsoon6mo_lag_denDF<-Pwr(m1_monsoon6mo_lag)[2,2]
  m2_monsoon6mo_lag_denDF<-Pwr(m2_monsoon6mo_lag)[2,2]
  m3_monsoon6mo_lag_denDF<-Pwr(m3_monsoon6mo_lag)[2,2]
  m1_monsoon6mo_lag_AR1_denDF<-Pwr(m1_monsoon6mo_lag_AR1)[2,2]
  m1_monsoon6mo_lag_AR2_denDF<-Pwr(m1_monsoon6mo_lag_AR2)[2,2]
  m2_monsoon6mo_lag_AR1_denDF<-Pwr(m2_monsoon6mo_lag_AR1)[2,2]
  m2_monsoon6mo_lag_AR2_denDF<-Pwr(m2_monsoon6mo_lag_AR2)[2,2]
  m3_monsoon6mo_lag_AR1_denDF<-Pwr(m3_monsoon6mo_lag_AR1)[2,2]
  m3_monsoon6mo_lag_AR2_denDF<-Pwr(m3_monsoon6mo_lag_AR2)[2,2]
  m1_spring6mo_lag_F<-Pwr(m1_spring6mo_lag)[2,3]
  m2_spring6mo_lag_F<-Pwr(m2_spring6mo_lag)[2,3]
  m3_spring6mo_lag_F<-Pwr(m3_spring6mo_lag)[2,3]
  m1_spring6mo_lag_AR1_F<-Pwr(m1_spring6mo_lag_AR1)[2,3]
  m1_spring6mo_lag_AR2_F<-Pwr(m1_spring6mo_lag_AR2)[2,3]
  m2_spring6mo_lag_AR1_F<-Pwr(m2_spring6mo_lag_AR1)[2,3]
  m2_spring6mo_lag_AR2_F<-Pwr(m2_spring6mo_lag_AR2)[2,3]
  m3_spring6mo_lag_AR1_F<-Pwr(m3_spring6mo_lag_AR1)[2,3]
  m3_spring6mo_lag_AR2_F<-Pwr(m3_spring6mo_lag_AR2)[2,3]
  m1_monsoon6mo_lag_F<-Pwr(m1_monsoon6mo_lag)[2,3]
  m2_monsoon6mo_lag_F<-Pwr(m2_monsoon6mo_lag)[2,3]
  m3_monsoon6mo_lag_F<-Pwr(m3_monsoon6mo_lag)[2,3]
  m1_monsoon6mo_lag_AR1_F<-Pwr(m1_monsoon6mo_lag_AR1)[2,3]
  m1_monsoon6mo_lag_AR2_F<-Pwr(m1_monsoon6mo_lag_AR2)[2,3]
  m2_monsoon6mo_lag_AR1_F<-Pwr(m2_monsoon6mo_lag_AR1)[2,3]
  m2_monsoon6mo_lag_AR2_F<-Pwr(m2_monsoon6mo_lag_AR2)[2,3]
  m3_monsoon6mo_lag_AR1_F<-Pwr(m3_monsoon6mo_lag_AR1)[2,3]
  m3_monsoon6mo_lag_AR2_F<-Pwr(m3_monsoon6mo_lag_AR2)[2,3]
  m1_spring6mo_lag_nc<-Pwr(m1_spring6mo_lag)[2,4]
  m2_spring6mo_lag_nc<-Pwr(m2_spring6mo_lag)[2,4]
  m3_spring6mo_lag_nc<-Pwr(m3_spring6mo_lag)[2,4]
  m1_spring6mo_lag_AR1_nc<-Pwr(m1_spring6mo_lag_AR1)[2,4]
  m1_spring6mo_lag_AR2_nc<-Pwr(m1_spring6mo_lag_AR2)[2,4]
  m2_spring6mo_lag_AR1_nc<-Pwr(m2_spring6mo_lag_AR1)[2,4]
  m2_spring6mo_lag_AR2_nc<-Pwr(m2_spring6mo_lag_AR2)[2,4]
  m3_spring6mo_lag_AR1_nc<-Pwr(m3_spring6mo_lag_AR1)[2,4]
  m3_spring6mo_lag_AR2_nc<-Pwr(m3_spring6mo_lag_AR2)[2,4]
  m1_monsoon6mo_lag_nc<-Pwr(m1_monsoon6mo_lag)[2,4]
  m2_monsoon6mo_lag_nc<-Pwr(m2_monsoon6mo_lag)[2,4]
  m3_monsoon6mo_lag_nc<-Pwr(m3_monsoon6mo_lag)[2,4]
  m1_monsoon6mo_lag_AR1_nc<-Pwr(m1_monsoon6mo_lag_AR1)[2,4]
  m1_monsoon6mo_lag_AR2_nc<-Pwr(m1_monsoon6mo_lag_AR2)[2,4]
  m2_monsoon6mo_lag_AR1_nc<-Pwr(m2_monsoon6mo_lag_AR1)[2,4]
  m2_monsoon6mo_lag_AR2_nc<-Pwr(m2_monsoon6mo_lag_AR2)[2,4]
  m3_monsoon6mo_lag_AR1_nc<-Pwr(m3_monsoon6mo_lag_AR1)[2,4]
  m3_monsoon6mo_lag_AR2_nc<-Pwr(m3_monsoon6mo_lag_AR2)[2,4]
  m1_spring6mo_lag_power<-Pwr(m1_spring6mo_lag)[2,5]
  m2_spring6mo_lag_power<-Pwr(m2_spring6mo_lag)[2,5]
  m3_spring6mo_lag_power<-Pwr(m3_spring6mo_lag)[2,5]
  m1_spring6mo_lag_AR1_power<-Pwr(m1_spring6mo_lag_AR1)[2,5]
  m1_spring6mo_lag_AR2_power<-Pwr(m1_spring6mo_lag_AR2)[2,5]
  m2_spring6mo_lag_AR1_power<-Pwr(m2_spring6mo_lag_AR1)[2,5]
  m2_spring6mo_lag_AR2_power<-Pwr(m2_spring6mo_lag_AR2)[2,5]
  m3_spring6mo_lag_AR1_power<-Pwr(m3_spring6mo_lag_AR1)[2,5]
  m3_spring6mo_lag_AR2_power<-Pwr(m3_spring6mo_lag_AR2)[2,5]
  m1_monsoon6mo_lag_power<-Pwr(m1_monsoon6mo_lag)[2,5]
  m2_monsoon6mo_lag_power<-Pwr(m2_monsoon6mo_lag)[2,5]
  m3_monsoon6mo_lag_power<-Pwr(m3_monsoon6mo_lag)[2,5]
  m1_monsoon6mo_lag_AR1_power<-Pwr(m1_monsoon6mo_lag_AR1)[2,5]
  m1_monsoon6mo_lag_AR2_power<-Pwr(m1_monsoon6mo_lag_AR2)[2,5]
  m2_monsoon6mo_lag_AR1_power<-Pwr(m2_monsoon6mo_lag_AR1)[2,5]
  m2_monsoon6mo_lag_AR2_power<-Pwr(m2_monsoon6mo_lag_AR2)[2,5]
  m3_monsoon6mo_lag_AR1_power<-Pwr(m3_monsoon6mo_lag_AR1)[2,5]
  m3_monsoon6mo_lag_AR2_power<-Pwr(m3_monsoon6mo_lag_AR2)[2,5]
  
  # bind all target output values together
  output_id<-cbind(speciesCode, ecosystem, m1_spring6mo_numDF,	m2_spring6mo_numDF,	m3_spring6mo_numDF,	m1_spring6mo_AR1_numDF,	m1_spring6mo_AR2_numDF,	m2_spring6mo_AR1_numDF,	m2_spring6mo_AR2_numDF,	m3_spring6mo_AR1_numDF,	m3_spring6mo_AR2_numDF,	m1_monsoon6mo_numDF,	m2_monsoon6mo_numDF,	m3_monsoon6mo_numDF,	m1_monsoon6mo_AR1_numDF,	m1_monsoon6mo_AR2_numDF,	m2_monsoon6mo_AR1_numDF,	m2_monsoon6mo_AR2_numDF,	m3_monsoon6mo_AR1_numDF,	m3_monsoon6mo_AR2_numDF,	m1_spring6mo_denDF,	m2_spring6mo_denDF,	m3_spring6mo_denDF,	m1_spring6mo_AR1_denDF,	m1_spring6mo_AR2_denDF,	m2_spring6mo_AR1_denDF,	m2_spring6mo_AR2_denDF,	m3_spring6mo_AR1_denDF,	m3_spring6mo_AR2_denDF,	m1_monsoon6mo_denDF,	m2_monsoon6mo_denDF,	m3_monsoon6mo_denDF,	m1_monsoon6mo_AR1_denDF,	m1_monsoon6mo_AR2_denDF,	m2_monsoon6mo_AR1_denDF,	m2_monsoon6mo_AR2_denDF,	m3_monsoon6mo_AR1_denDF,	m3_monsoon6mo_AR2_denDF,	m1_spring6mo_F,	m2_spring6mo_F,	m3_spring6mo_F,	m1_spring6mo_AR1_F,	m1_spring6mo_AR2_F,	m2_spring6mo_AR1_F,	m2_spring6mo_AR2_F,	m3_spring6mo_AR1_F,	m3_spring6mo_AR2_F,	m1_monsoon6mo_F,	m2_monsoon6mo_F,	m3_monsoon6mo_F,	m1_monsoon6mo_AR1_F,	m1_monsoon6mo_AR2_F,	m2_monsoon6mo_AR1_F,	m2_monsoon6mo_AR2_F,	m3_monsoon6mo_AR1_F,	m3_monsoon6mo_AR2_F,	m1_spring6mo_nc,	m2_spring6mo_nc,	m3_spring6mo_nc,	m1_spring6mo_AR1_nc,	m1_spring6mo_AR2_nc,	m2_spring6mo_AR1_nc,	m2_spring6mo_AR2_nc,	m3_spring6mo_AR1_nc,	m3_spring6mo_AR2_nc,	m1_monsoon6mo_nc,	m2_monsoon6mo_nc,	m3_monsoon6mo_nc,	m1_monsoon6mo_AR1_nc,	m1_monsoon6mo_AR2_nc,	m2_monsoon6mo_AR1_nc,	m2_monsoon6mo_AR2_nc,	m3_monsoon6mo_AR1_nc,	m3_monsoon6mo_AR2_nc,	m1_spring6mo_power,	m2_spring6mo_power,	m3_spring6mo_power,	m1_spring6mo_AR1_power,	m1_spring6mo_AR2_power,	m2_spring6mo_AR1_power,	m2_spring6mo_AR2_power,	m3_spring6mo_AR1_power,	m3_spring6mo_AR2_power,	m1_monsoon6mo_power,	m2_monsoon6mo_power,	m3_monsoon6mo_power,	m1_monsoon6mo_AR1_power,	m1_monsoon6mo_AR2_power,	m2_monsoon6mo_AR1_power,	m2_monsoon6mo_AR2_power,	m3_monsoon6mo_AR1_power,	m3_monsoon6mo_AR2_power, m1_spring6mo_lag_numDF,m2_spring6mo_lag_numDF,m3_spring6mo_lag_numDF,m1_spring6mo_lag_AR1_numDF,	m1_spring6mo_lag_AR2_numDF,m2_spring6mo_lag_AR1_numDF,m2_spring6mo_lag_AR2_numDF,m3_spring6mo_lag_AR1_numDF,m3_spring6mo_lag_AR2_numDF,m1_monsoon6mo_lag_numDF,m2_monsoon6mo_lag_numDF,m3_monsoon6mo_lag_numDF,m1_monsoon6mo_lag_AR1_numDF,m1_monsoon6mo_lag_AR2_numDF,	m2_monsoon6mo_lag_AR1_numDF,m2_monsoon6mo_lag_AR2_numDF,m3_monsoon6mo_lag_AR1_numDF,m3_monsoon6mo_lag_AR2_numDF, m1_spring6mo_lag_denDF,m2_spring6mo_lag_denDF,m3_spring6mo_lag_denDF,m1_spring6mo_lag_AR1_denDF,m1_spring6mo_lag_AR2_denDF,m2_spring6mo_lag_AR1_denDF,m2_spring6mo_lag_AR2_denDF,
                   m3_spring6mo_lag_AR1_denDF,m3_spring6mo_lag_AR2_denDF,m1_monsoon6mo_lag_denDF,
                   m2_monsoon6mo_lag_denDF,m3_monsoon6mo_lag_denDF,m1_monsoon6mo_lag_AR1_denDF,	
                   m1_monsoon6mo_lag_AR2_denDF,m2_monsoon6mo_lag_AR1_denDF,m2_monsoon6mo_lag_AR2_denDF,
                   m3_monsoon6mo_lag_AR1_denDF,	m3_monsoon6mo_lag_AR2_denDF,	m1_spring6mo_lag_F,	m2_spring6mo_lag_F,
                   m3_spring6mo_lag_F,	m1_spring6mo_lag_AR1_F,	m1_spring6mo_lag_AR2_F,	m2_spring6mo_lag_AR1_F,
                   m2_spring6mo_lag_AR2_F,	m3_spring6mo_lag_AR1_F,	m3_spring6mo_lag_AR2_F,	m1_monsoon6mo_lag_F,
                   m2_monsoon6mo_lag_F,	m3_monsoon6mo_lag_F,	m1_monsoon6mo_lag_AR1_F,	m1_monsoon6mo_lag_AR2_F,
                   m2_monsoon6mo_lag_AR1_F,m2_monsoon6mo_lag_AR2_F,m3_monsoon6mo_lag_AR1_F,	m3_monsoon6mo_lag_AR2_F,
                   m1_spring6mo_lag_nc,	m2_spring6mo_lag_nc,	m3_spring6mo_lag_nc,	m1_spring6mo_lag_AR1_nc,
                   m1_spring6mo_lag_AR2_nc,m2_spring6mo_lag_AR1_nc,	m2_spring6mo_lag_AR2_nc,	m3_spring6mo_lag_AR1_nc,
                   m3_spring6mo_lag_AR2_nc,	m1_monsoon6mo_lag_nc,	m2_monsoon6mo_lag_nc,	m3_monsoon6mo_lag_nc,
                   m1_monsoon6mo_lag_AR1_nc,	m1_monsoon6mo_lag_AR2_nc,	m2_monsoon6mo_lag_AR1_nc,
                   m2_monsoon6mo_lag_AR2_nc,	m3_monsoon6mo_lag_AR1_nc,m3_monsoon6mo_lag_AR2_nc,
                   m1_spring6mo_lag_power,	m2_spring6mo_lag_power,m3_spring6mo_lag_power,m1_spring6mo_lag_AR1_power,
                   m1_spring6mo_lag_AR2_power,	m2_spring6mo_lag_AR1_power,	m2_spring6mo_lag_AR2_power,
                   m3_spring6mo_lag_AR1_power,	m3_spring6mo_lag_AR2_power,	m1_monsoon6mo_lag_power,
                   m2_monsoon6mo_lag_power,	m3_monsoon6mo_lag_power,	m1_monsoon6mo_lag_AR1_power, m1_monsoon6mo_lag_AR2_power, m2_monsoon6mo_lag_AR1_power, m2_monsoon6mo_lag_AR2_power, m3_monsoon6mo_lag_AR1_power,	m3_monsoon6mo_lag_AR2_power)
  
  # append results for each bee species to the output dataset
  beeCSF_output<-rbind(beeCSF_output,output_id)
  
}

# Write .csv file of output
beeCSF_output <- data.frame(beeCSF_output) # make matrix into data frame
beeCSF_output2 <- beeCSF_output[-1,] # remove blank row

#write.csv(beeCSF_output2, "CSF_power_analysis_results_2023-08-29.csv",row.names=FALSE)



##### Summary of Power Analysis Results #####

# read in CSF model results and power analysis results, and join the two
CSFs <- read.csv("bee_CSFs_2023-08-29.csv")
power<-read.csv("CSF_power_analysis_results_2023-08-29.csv")
CSFs<-left_join(CSFs,power,by=c("code","ecosystem"))

# subset to just include bee populations for which aridity did not predict abundance (i.e., the null model was best)
null <- subset(CSFs, dAICc_null == 0)

# create a new data frame of desired columns
null2<-null[,c(1,2,547:726)]

# convert data frame from wide to long form
null2_long<-pivot_longer(data=null2,cols=3:182,names_to="name",values_to="value")

# create separate columns for model name and statistic name, and then remove no-longer-needed "name" column
null2_long$stat<-word(null2_long$name, -1, sep = "_")
null2_long$model<-word(null2_long$name, -5, -2, sep = "_")
null2_long$name<-NULL

# convert back to wide form
wide<-pivot_wider(data=null2_long,names_from="stat",values_from="value")

# calculate the maximum power value obtained for a given population
summary<-wide %>% group_by(code,ecosystem) %>% summarise(power=max(power))

# join information to the summary data frame related to the model with the greatest power for each population
summary<-left_join(summary,wide,by=c("code","ecosystem","power"))

# read in the bee species list, and join bee species information to the summary data frame
specieslist<-read.csv("SEVBeeSpeciesList2002-2019_revised2023-07-19.csv")
summary<-left_join(summary,specieslist[,1:5],by="code")

# recode the ecosystem variable
summary$ecosystem<-as.factor(summary$ecosystem)
levels(summary$ecosystem) <- c("Plains Grassland","Desert Shrubland","Desert Grassland")

# reorder columns
summary<-summary[,c(9:12,2,4:8,3)]

# save as .csv file
#write.csv(summary,"power_analysis_table_2023-08-29.csv",row.names=FALSE)
