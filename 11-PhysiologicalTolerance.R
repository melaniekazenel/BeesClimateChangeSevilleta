################################################################################### 
# How do thermal and desiccation tolerances relate to climate sensitivity, abundance change over time, and body mass?

# Heat and desiccation tolerances predict bee abundance under climate change
# Melanie R. Kazenel, Karen W. Wright, Terry Griswold, Kenneth D. Whitney, and Jennifer A. Rudgers

# Date: 2023-12-01
# Corresponding author's email: melanie.kazenel@gmail.com
################################################################################### 

# Load required packages
library(ggplot2)
library(dplyr)
library(Rfast)
library(tidyr)
library(ape)
library(phyr)
library(phylosignal)
library(ggtext)
library(rr2)
library(phylobase)
library(Rfast)
library(patchwork)


##### Initial data formatting #####

# read in data from thermal tolerance (CTMax) trials
ctmax <- read.csv("ctmax_forpub.csv")

# z-score thermal and VPD tolerance data
ctmax$end_temp_z <- (ctmax$end_temp - mean(ctmax$end_temp))/sd(ctmax$end_temp)
ctmax$end_vpd_z <- (ctmax$end_vpd - mean(ctmax$end_vpd, na.rm=TRUE))/sd(ctmax$end_vpd,na.rm=TRUE)

# read in desiccation tolerance data
destol <- read.csv("destol_forpub.csv")

# for bees for which loss of responsiveness was not observed during trials, change end time of the trial to midpoint time between when the bee was last checked and when it was discovered unresponsive
destol$duration_hr<-ifelse(destol$end_observed==1, destol$duration_hr, destol$duration_hr_adjusted)
head(destol[,c("duration_hr","duration_hr_adjusted")])

# z-score desiccation tolerance metrics
destol$duration_hr_z <- (destol$duration_hr - mean(destol$duration_hr))/sd(destol$duration_hr)
destol$cwc_z <- (destol$cwc - mean(destol$cwc, na.rm=TRUE))/sd(destol$cwc,na.rm=TRUE)


##### Calculate sample sizes from thermal and desiccation tolerance trials #####

# desiccation tolerance
summary<-destol %>%group_by(genus) %>% summarise(count_destol=n())
# thermal tolerance
summary2<-ctmax %>%group_by(genus) %>% summarise(count_ctmax=n())
# combine data frames
summary3<-left_join(summary2,summary,by="genus")


##### Calculate genus-level mean thermal and desiccation tolerance metrics, and merge data frames #####

# Thermal tolerance:
# calculate genus-level mean thermal tolerance metrics
means_ctmax <- ctmax %>% group_by(family,genus) %>% summarise(mean_end_temp=mean(end_temp),se_end_temp=sd(end_temp)/sqrt(n()),mean_end_temp_z=mean(end_temp_z),se_end_temp_z=sd(end_temp_z)/sqrt(n()),mean_end_vpd=mean(end_vpd,na.rm=TRUE),se_end_vpd=sd(end_vpd,na.rm=TRUE)/sqrt(n()),mean_end_vpd_z=mean(end_vpd_z,na.rm=TRUE),se_end_vpd_z=sd(end_vpd_z,na.rm=TRUE)/sqrt(n()),count_ctmax=n())
# for Centris, replace NaN values with NAs
means_ctmax[5,7:10]<-NA
# replace NAs with 0s for SEs for Megachile and Coelioxys
means_ctmax[c(13,14),c("se_end_temp","se_end_temp_z","se_end_vpd","se_end_vpd_z")]<-0

# Desiccation tolerance:
# calculate genus-level mean desiccation tolerance metrics
means_destol <- destol %>% group_by(family,genus) %>% summarise(mean_duration_hr=mean(duration_hr),se_duration_hr=sd(duration_hr)/sqrt(n()),mean_duration_hr_z=mean(duration_hr_z),se_duration_hr_z=sd(duration_hr_z)/sqrt(n()),mean_cwc=mean(cwc,na.rm=TRUE),se_cwc=sd(cwc,na.rm=TRUE)/sqrt(n()), mean_cwc_z=mean(cwc_z,na.rm=TRUE),se_cwc_z=sd(cwc_z,na.rm=TRUE)/sqrt(n()),count_destol=n())
# replace NAs with 0s for SEs for Ashmeadiella
means_destol[13,c("se_duration_hr","se_duration_hr_z","se_cwc","se_cwc_z")]<-0

# Join thermal and desiccation tolerance data frames
means2 <- left_join(means_ctmax,means_destol, by=c("family","genus"))

# Calculate combined thermal + desiccation tolerance metric, and its standard error
means2$combined_tol_duration_hr<-means2$mean_end_temp_z+means2$mean_duration_hr_z
means2$combined_tol_duration_hr_se<-means2$se_end_temp_z+means2$se_duration_hr_z


##### Calculate thermal and desiccation tolerance summary statistics #####

# thermal tolerance summary
summary(means2$mean_end_temp)
# view values in sorted order for each genus
ctmaxmeans<-means2[,c("genus","mean_end_temp","count_ctmax")]
ctmaxmeans<-ctmaxmeans[order(ctmaxmeans$mean_end_temp),]

# desiccation tolerance summary
summary(means2$mean_duration_hr)
# view values in sorted order for each genus
destolmeans<-means2[,c("genus","mean_duration_hr","count_destol")]
destolmeans<-destolmeans[order(destolmeans$mean_duration_hr),]


##### Calculate difference in physiological tolerance between most and least tolerant genera #####

# CTMax
summary<-ctmax %>% group_by(genus) %>% summarise(mean=mean(end_temp))
summary(summary$mean)
(max(summary$mean)-min(summary$mean))/min(summary$mean)

# VPD
summary<-ctmax %>% group_by(genus) %>% summarise(mean=mean(end_vpd))
summary(summary$mean)
(max(summary$mean,na.rm=TRUE)-min(summary$mean,na.rm=TRUE))/min(summary$mean,na.rm=TRUE)

# Desiccation tolerance
summary<-destol %>% group_by(genus) %>% summarise(mean=mean(duration_hr))
summary(summary$mean)
(max(summary$mean)-min(summary$mean))/min(summary$mean)


##### Format linear CSF parameter estimate data, and merge with physiological tolerance data #####

# Read in model output from climate sensitivity functions
CSFs <- read.csv("bee_CSFs_2023-12-01.csv")

### Extract the linear parameter estimate from the best model for each population in cases where the linear model was best

# create a separate data frame for each model variant, containing data for the populations for which that model variant was superior, and format the data

monsoon6mo <- subset(CSFs, dAICc_m1_monsoon6mo == 0)
monsoon6mo <- monsoon6mo[,c("code","ecosystem","dAICc_m1_monsoon6mo","ParamEst_lin_m1_monsoon6mo")]
monsoon6mo$best_model <- "m1_monsoon6mo"
names(monsoon6mo)[3] <- "dAICc"
names(monsoon6mo)[4] <- "ParamEstimate"

monsoon6mo_AR1 <- subset(CSFs, dAICc_m1_monsoon6mo_AR1 == 0)
monsoon6mo_AR1 <- monsoon6mo_AR1[,c("code","ecosystem","dAICc_m1_monsoon6mo_AR1","ParamEst_lin_m1_monsoon6mo_AR1")]
monsoon6mo_AR1$best_model <- "m1_monsoon6mo_AR1"
names(monsoon6mo_AR1)[3] <- "dAICc"
names(monsoon6mo_AR1)[4] <- "ParamEstimate"

monsoon6mo_AR2 <- subset(CSFs, dAICc_m1_monsoon6mo_AR2 == 0)
monsoon6mo_AR2 <- monsoon6mo_AR2[,c("code","ecosystem","dAICc_m1_monsoon6mo_AR2","ParamEst_lin_m1_monsoon6mo_AR2")]
monsoon6mo_AR2$best_model <- "m1_monsoon6mo_AR2"
names(monsoon6mo_AR2)[3] <- "dAICc"
names(monsoon6mo_AR2)[4] <- "ParamEstimate"

monsoon6mo_lag <- subset(CSFs, dAICc_m1_monsoon6mo_lag == 0)
monsoon6mo_lag <- monsoon6mo_lag[,c("code","ecosystem","dAICc_m1_monsoon6mo_lag","ParamEst_lin_m1_monsoon6mo_lag")]
monsoon6mo_lag$best_model <- "m1_monsoon6mo_lag"
names(monsoon6mo_lag)[3] <- "dAICc"
names(monsoon6mo_lag)[4] <- "ParamEstimate"

monsoon6mo_lag_AR1 <- subset(CSFs, dAICc_m1_monsoon6mo_lag_AR1 == 0)
monsoon6mo_lag_AR1 <- monsoon6mo_lag_AR1[,c("code","ecosystem","dAICc_m1_monsoon6mo_lag_AR1","ParamEst_lin_m1_monsoon6mo_lag_AR1")]
monsoon6mo_lag_AR1$best_model <- "m1_monsoon6mo_lag_AR1"
names(monsoon6mo_lag_AR1)[3] <- "dAICc"
names(monsoon6mo_lag_AR1)[4] <- "ParamEstimate"

monsoon6mo_lag_AR2 <- subset(CSFs, dAICc_m1_monsoon6mo_lag_AR2 == 0)
monsoon6mo_lag_AR2 <- monsoon6mo_lag_AR2[,c("code","ecosystem","dAICc_m1_monsoon6mo_lag_AR2","ParamEst_lin_m1_monsoon6mo_lag_AR2")]
monsoon6mo_lag_AR2$best_model <- "m1_monsoon6mo_lag_AR2"
names(monsoon6mo_lag_AR2)[3] <- "dAICc"
names(monsoon6mo_lag_AR2)[4] <- "ParamEstimate"

spring6mo <- subset(CSFs, dAICc_m1_spring6mo == 0)
spring6mo <- spring6mo[,c("code","ecosystem","dAICc_m1_spring6mo","ParamEst_lin_m1_spring6mo")]
spring6mo$best_model <- "m1_spring6mo"
names(spring6mo)[3] <- "dAICc"
names(spring6mo)[4] <- "ParamEstimate"

spring6mo_AR1 <- subset(CSFs, dAICc_m1_spring6mo_AR1 == 0)
spring6mo_AR1 <- spring6mo_AR1[,c("code","ecosystem","dAICc_m1_spring6mo_AR1","ParamEst_lin_m1_spring6mo_AR1")]
spring6mo_AR1$best_model <- "m1_spring6mo_AR1"
names(spring6mo_AR1)[3] <- "dAICc"
names(spring6mo_AR1)[4] <- "ParamEstimate"

spring6mo_AR2 <- subset(CSFs, dAICc_m1_spring6mo_AR2 == 0)
spring6mo_AR2 <- spring6mo_AR2[,c("code","ecosystem","dAICc_m1_spring6mo_AR2","ParamEst_lin_m1_spring6mo_AR2")]
spring6mo_AR2$best_model <- "m1_spring6mo_AR2"
names(spring6mo_AR2)[3] <- "dAICc"
names(spring6mo_AR2)[4] <- "ParamEstimate"

spring6mo_lag <- subset(CSFs, dAICc_m1_spring6mo_lag == 0)
spring6mo_lag <- spring6mo_lag[,c("code","ecosystem","dAICc_m1_spring6mo_lag","ParamEst_lin_m1_spring6mo_lag")]
spring6mo_lag$best_model <- "m1_spring6mo_lag"
names(spring6mo_lag)[3] <- "dAICc"
names(spring6mo_lag)[4] <- "ParamEstimate"

spring6mo_lag_AR1 <- subset(CSFs, dAICc_m1_spring6mo_lag_AR1 == 0)
spring6mo_lag_AR1 <- spring6mo_lag_AR1[,c("code","ecosystem","dAICc_m1_spring6mo_lag_AR1","ParamEst_lin_m1_spring6mo_lag_AR1")]
spring6mo_lag_AR1$best_model <- "m1_spring6mo_lag_AR1"
names(spring6mo_lag_AR1)[3] <- "dAICc"
names(spring6mo_lag_AR1)[4] <- "ParamEstimate"

spring6mo_lag_AR2 <- subset(CSFs, dAICc_m1_spring6mo_lag_AR2 == 0)
spring6mo_lag_AR2 <- spring6mo_lag_AR2[,c("code","ecosystem","dAICc_m1_spring6mo_lag_AR2","ParamEst_lin_m1_spring6mo_lag_AR2")]
spring6mo_lag_AR2$best_model <- "m1_spring6mo_lag_AR2"
names(spring6mo_lag_AR2)[3] <- "dAICc"
names(spring6mo_lag_AR2)[4] <- "ParamEstimate"

# bind together the separate data frames for the different model variants, and add a "model type" column indicating that the linear relationship was best
paramest_linear <- bind_rows(monsoon6mo,monsoon6mo_AR1,monsoon6mo_AR2,monsoon6mo_lag,monsoon6mo_lag_AR1,monsoon6mo_lag_AR2,spring6mo,spring6mo_AR1,spring6mo_AR2,spring6mo_lag,spring6mo_lag_AR1,spring6mo_lag_AR2)

paramest_linear$model_type <- "m1"

# positive the data so that positive values indicate an increase in abundance under increasing aridity
paramest_linear$ParamEstimate<-paramest_linear$ParamEstimate*(-1)


### Extract the linear parameter estimate from the best model for each population in cases where the quadratic model was best

# create a separate data frame for each model variant, containing data for the populations for which that model variant was superior, and format the data

monsoon6mo <- subset(CSFs, dAICc_m2_monsoon6mo == 0)
monsoon6mo <- monsoon6mo[,c("code","ecosystem","dAICc_m2_monsoon6mo","ParamEst_lin_m2_monsoon6mo")]
monsoon6mo$best_model <- "m2_monsoon6mo"
names(monsoon6mo)[3] <- "dAICc"
names(monsoon6mo)[4] <- "ParamEstimate"

monsoon6mo_AR1 <- subset(CSFs, dAICc_m2_monsoon6mo_AR1 == 0)
monsoon6mo_AR1 <- monsoon6mo_AR1[,c("code","ecosystem","dAICc_m2_monsoon6mo_AR1","ParamEst_lin_m2_monsoon6mo_AR1")]
monsoon6mo_AR1$best_model <- "m2_monsoon6mo_AR1"
names(monsoon6mo_AR1)[3] <- "dAICc"
names(monsoon6mo_AR1)[4] <- "ParamEstimate"

monsoon6mo_AR2 <- subset(CSFs, dAICc_m2_monsoon6mo_AR2 == 0)
monsoon6mo_AR2 <- monsoon6mo_AR2[,c("code","ecosystem","dAICc_m2_monsoon6mo_AR2","ParamEst_lin_m2_monsoon6mo_AR2")]
monsoon6mo_AR2$best_model <- "m2_monsoon6mo_AR2"
names(monsoon6mo_AR2)[3] <- "dAICc"
names(monsoon6mo_AR2)[4] <- "ParamEstimate"

monsoon6mo_lag <- subset(CSFs, dAICc_m2_monsoon6mo_lag == 0)
monsoon6mo_lag <- monsoon6mo_lag[,c("code","ecosystem","dAICc_m2_monsoon6mo_lag","ParamEst_lin_m2_monsoon6mo_lag")]
monsoon6mo_lag$best_model <- "m2_monsoon6mo_lag"
names(monsoon6mo_lag)[3] <- "dAICc"
names(monsoon6mo_lag)[4] <- "ParamEstimate"

monsoon6mo_lag_AR1 <- subset(CSFs, dAICc_m2_monsoon6mo_lag_AR1 == 0)
monsoon6mo_lag_AR1 <- monsoon6mo_lag_AR1[,c("code","ecosystem","dAICc_m2_monsoon6mo_lag_AR1","ParamEst_lin_m2_monsoon6mo_lag_AR1")]
monsoon6mo_lag_AR1$best_model <- "m2_monsoon6mo_lag_AR1"
names(monsoon6mo_lag_AR1)[3] <- "dAICc"
names(monsoon6mo_lag_AR1)[4] <- "ParamEstimate"

monsoon6mo_lag_AR2 <- subset(CSFs, dAICc_m2_monsoon6mo_lag_AR2 == 0)
monsoon6mo_lag_AR2 <- monsoon6mo_lag_AR2[,c("code","ecosystem","dAICc_m2_monsoon6mo_lag_AR2","ParamEst_lin_m2_monsoon6mo_lag_AR2")]
monsoon6mo_lag_AR2$best_model <- "m2_monsoon6mo_lag_AR2"
names(monsoon6mo_lag_AR2)[3] <- "dAICc"
names(monsoon6mo_lag_AR2)[4] <- "ParamEstimate"

spring6mo <- subset(CSFs, dAICc_m2_spring6mo == 0)
spring6mo <- spring6mo[,c("code","ecosystem","dAICc_m2_spring6mo","ParamEst_lin_m2_spring6mo")]
spring6mo$best_model <- "m2_spring6mo"
names(spring6mo)[3] <- "dAICc"
names(spring6mo)[4] <- "ParamEstimate"

spring6mo_AR1 <- subset(CSFs, dAICc_m2_spring6mo_AR1 == 0)
spring6mo_AR1 <- spring6mo_AR1[,c("code","ecosystem","dAICc_m2_spring6mo_AR1","ParamEst_lin_m2_spring6mo_AR1")]
spring6mo_AR1$best_model <- "m2_spring6mo_AR1"
names(spring6mo_AR1)[3] <- "dAICc"
names(spring6mo_AR1)[4] <- "ParamEstimate"

spring6mo_AR2 <- subset(CSFs, dAICc_m2_spring6mo_AR2 == 0)
spring6mo_AR2 <- spring6mo_AR2[,c("code","ecosystem","dAICc_m2_spring6mo_AR2","ParamEst_lin_m2_spring6mo_AR2")]
spring6mo_AR2$best_model <- "m2_spring6mo_AR2"
names(spring6mo_AR2)[3] <- "dAICc"
names(spring6mo_AR2)[4] <- "ParamEstimate"

spring6mo_lag <- subset(CSFs, dAICc_m2_spring6mo_lag == 0)
spring6mo_lag <- spring6mo_lag[,c("code","ecosystem","dAICc_m2_spring6mo_lag","ParamEst_lin_m2_spring6mo_lag")]
spring6mo_lag$best_model <- "m2_spring6mo_lag"
names(spring6mo_lag)[3] <- "dAICc"
names(spring6mo_lag)[4] <- "ParamEstimate"

spring6mo_lag_AR1 <- subset(CSFs, dAICc_m2_spring6mo_lag_AR1 == 0)
spring6mo_lag_AR1 <- spring6mo_lag_AR1[,c("code","ecosystem","dAICc_m2_spring6mo_lag_AR1","ParamEst_lin_m2_spring6mo_lag_AR1")]
spring6mo_lag_AR1$best_model <- "m2_spring6mo_lag_AR1"
names(spring6mo_lag_AR1)[3] <- "dAICc"
names(spring6mo_lag_AR1)[4] <- "ParamEstimate"

spring6mo_lag_AR2 <- subset(CSFs, dAICc_m2_spring6mo_lag_AR2 == 0)
spring6mo_lag_AR2 <- spring6mo_lag_AR2[,c("code","ecosystem","dAICc_m2_spring6mo_lag_AR2","ParamEst_lin_m2_spring6mo_lag_AR2")]
spring6mo_lag_AR2$best_model <- "m2_spring6mo_lag_AR2"
names(spring6mo_lag_AR2)[3] <- "dAICc"
names(spring6mo_lag_AR2)[4] <- "ParamEstimate"

# bind together the separate data frames for the different model variants, and add a "model type" column indicating that the quadratic relationship was best
paramest_quadratic <- bind_rows(monsoon6mo,monsoon6mo_AR1,monsoon6mo_AR2,monsoon6mo_lag,monsoon6mo_lag_AR1,monsoon6mo_lag_AR2,spring6mo,spring6mo_AR1,spring6mo_AR2,spring6mo_lag,spring6mo_lag_AR1,spring6mo_lag_AR2)

paramest_quadratic$model_type <- "m2"

# Positive the data so that positive values indicate an increase in abundance under increasing aridity
paramest_quadratic$ParamEstimate<-paramest_quadratic$ParamEstimate*(-1)


### Extract the linear parameter estimate from the best model for each population in cases where the cubic model was best

# create a separate data frame for each model variant, containing data for the populations for which that model variant was superior, and format the data

monsoon6mo <- subset(CSFs, dAICc_m3_monsoon6mo == 0)
monsoon6mo <- monsoon6mo[,c("code","ecosystem","dAICc_m3_monsoon6mo","ParamEst_lin_m3_monsoon6mo")]
monsoon6mo$best_model <- "m3_monsoon6mo"
names(monsoon6mo)[3] <- "dAICc"
names(monsoon6mo)[4] <- "ParamEstimate"

monsoon6mo_AR1 <- subset(CSFs, dAICc_m3_monsoon6mo_AR1 == 0)
monsoon6mo_AR1 <- monsoon6mo_AR1[,c("code","ecosystem","dAICc_m3_monsoon6mo_AR1","ParamEst_lin_m3_monsoon6mo_AR1")]
monsoon6mo_AR1$best_model <- "m3_monsoon6mo_AR1"
names(monsoon6mo_AR1)[3] <- "dAICc"
names(monsoon6mo_AR1)[4] <- "ParamEstimate"

monsoon6mo_AR2 <- subset(CSFs, dAICc_m3_monsoon6mo_AR2 == 0)
monsoon6mo_AR2 <- monsoon6mo_AR2[,c("code","ecosystem","dAICc_m3_monsoon6mo_AR2","ParamEst_lin_m3_monsoon6mo_AR2")]
monsoon6mo_AR2$best_model <- "m3_monsoon6mo_AR2"
names(monsoon6mo_AR2)[3] <- "dAICc"
names(monsoon6mo_AR2)[4] <- "ParamEstimate"

monsoon6mo_lag <- subset(CSFs, dAICc_m3_monsoon6mo_lag == 0)
monsoon6mo_lag <- monsoon6mo_lag[,c("code","ecosystem","dAICc_m3_monsoon6mo_lag","ParamEst_lin_m3_monsoon6mo_lag")]
monsoon6mo_lag$best_model <- "m3_monsoon6mo_lag"
names(monsoon6mo_lag)[3] <- "dAICc"
names(monsoon6mo_lag)[4] <- "ParamEstimate"

monsoon6mo_lag_AR1 <- subset(CSFs, dAICc_m3_monsoon6mo_lag_AR1 == 0)
monsoon6mo_lag_AR1 <- monsoon6mo_lag_AR1[,c("code","ecosystem","dAICc_m3_monsoon6mo_lag_AR1","ParamEst_lin_m3_monsoon6mo_lag_AR1")]
monsoon6mo_lag_AR1$best_model <- "m3_monsoon6mo_lag_AR1"
names(monsoon6mo_lag_AR1)[3] <- "dAICc"
names(monsoon6mo_lag_AR1)[4] <- "ParamEstimate"

monsoon6mo_lag_AR2 <- subset(CSFs, dAICc_m3_monsoon6mo_lag_AR2 == 0)
monsoon6mo_lag_AR2 <- monsoon6mo_lag_AR2[,c("code","ecosystem","dAICc_m3_monsoon6mo_lag_AR2","ParamEst_lin_m3_monsoon6mo_lag_AR2")]
monsoon6mo_lag_AR2$best_model <- "m3_monsoon6mo_lag_AR2"
names(monsoon6mo_lag_AR2)[3] <- "dAICc"
names(monsoon6mo_lag_AR2)[4] <- "ParamEstimate"

spring6mo <- subset(CSFs, dAICc_m3_spring6mo == 0)
spring6mo <- spring6mo[,c("code","ecosystem","dAICc_m3_spring6mo","ParamEst_lin_m3_spring6mo")]
spring6mo$best_model <- "m3_spring6mo"
names(spring6mo)[3] <- "dAICc"
names(spring6mo)[4] <- "ParamEstimate"

spring6mo_AR1 <- subset(CSFs, dAICc_m3_spring6mo_AR1 == 0)
spring6mo_AR1 <- spring6mo_AR1[,c("code","ecosystem","dAICc_m3_spring6mo_AR1","ParamEst_lin_m3_spring6mo_AR1")]
spring6mo_AR1$best_model <- "m3_spring6mo_AR1"
names(spring6mo_AR1)[3] <- "dAICc"
names(spring6mo_AR1)[4] <- "ParamEstimate"

spring6mo_AR2 <- subset(CSFs, dAICc_m3_spring6mo_AR2 == 0)
spring6mo_AR2 <- spring6mo_AR2[,c("code","ecosystem","dAICc_m3_spring6mo_AR2","ParamEst_lin_m3_spring6mo_AR2")]
spring6mo_AR2$best_model <- "m3_spring6mo_AR2"
names(spring6mo_AR2)[3] <- "dAICc"
names(spring6mo_AR2)[4] <- "ParamEstimate"

spring6mo_lag <- subset(CSFs, dAICc_m3_spring6mo_lag == 0)
spring6mo_lag <- spring6mo_lag[,c("code","ecosystem","dAICc_m3_spring6mo_lag","ParamEst_lin_m3_spring6mo_lag")]
spring6mo_lag$best_model <- "m3_spring6mo_lag"
names(spring6mo_lag)[3] <- "dAICc"
names(spring6mo_lag)[4] <- "ParamEstimate"

spring6mo_lag_AR1 <- subset(CSFs, dAICc_m3_spring6mo_lag_AR1 == 0)
spring6mo_lag_AR1 <- spring6mo_lag_AR1[,c("code","ecosystem","dAICc_m3_spring6mo_lag_AR1","ParamEst_lin_m3_spring6mo_lag_AR1")]
spring6mo_lag_AR1$best_model <- "m3_spring6mo_lag_AR1"
names(spring6mo_lag_AR1)[3] <- "dAICc"
names(spring6mo_lag_AR1)[4] <- "ParamEstimate"

spring6mo_lag_AR2 <- subset(CSFs, dAICc_m3_spring6mo_lag_AR2 == 0)
spring6mo_lag_AR2 <- spring6mo_lag_AR2[,c("code","ecosystem","dAICc_m3_spring6mo_lag_AR2","ParamEst_lin_m3_spring6mo_lag_AR2")]
spring6mo_lag_AR2$best_model <- "m3_spring6mo_lag_AR2"
names(spring6mo_lag_AR2)[3] <- "dAICc"
names(spring6mo_lag_AR2)[4] <- "ParamEstimate"

# bind together the separate data frames for the different model variants, and add a "model type" column indicating that the cubic relationship was best
paramest_cubic <- bind_rows(monsoon6mo,monsoon6mo_AR1,monsoon6mo_AR2,monsoon6mo_lag,monsoon6mo_lag_AR1,monsoon6mo_lag_AR2,spring6mo,spring6mo_AR1,spring6mo_AR2,spring6mo_lag,spring6mo_lag_AR1,spring6mo_lag_AR2)

paramest_cubic$model_type <- "m3"

# positive the data so that positive values indicate an increase in abundance under increasing aridity
paramest_cubic$ParamEstimate<-paramest_cubic$ParamEstimate*(-1)


# Combine data frames across model types 
paramest_all <- bind_rows(paramest_linear,paramest_quadratic,paramest_cubic)

# read in species list
specieslist<-read.csv("SEVBeeSpeciesList2002-2019_revised2023-07-19.csv")

# add species data to parameter estimate data frame
CSFs_paramest<-left_join(paramest_all,specieslist,by="code")

# read in list of species from thermal and desiccation tolerance trials
list_tol<-read.csv("speciessummaryctmaxdestol_2023-08-29.csv")

# filter CSF parameter estimate data frame to just contain species for which we measured thermal and/or desiccation tolerance
CSFs_paramest <- filter(CSFs_paramest,code %in% list_tol$code)

# for each genus, calculate mean and se of the linear parameter estimate
lin_param_summary <- CSFs_paramest %>% group_by(family, genus) %>% summarise(lin_param_mean=mean(ParamEstimate),lin_param_se=sd(ParamEstimate)/sqrt(n()))

# combine physiological tolerance and parameter estimate data frames
fulldata<-left_join(means2, lin_param_summary, by=c("family","genus"))



##### Format quadratic CSF parameter estimate data, and merge with physiological tolerance data #####

# Read in model output from climate sensitivity functions
CSFs <- read.csv("bee_CSFs_2023-12-01.csv")

### Create a data frame containing the quadratic parameter estimate from the best model for each population in cases where the linear model was best (quadratic parameter estimate = 0 in these cases)

# create a separate data frame for each model variant, containing data for the populations for which that model variant was superior, and format the data

monsoon6mo <- subset(CSFs, dAICc_m1_monsoon6mo == 0)
monsoon6mo <- monsoon6mo[,c("code","ecosystem","dAICc_m1_monsoon6mo")]
monsoon6mo$best_model <- "m1_monsoon6mo"
names(monsoon6mo)[3] <- "dAICc"
monsoon6mo$ParamEstimate<-0

monsoon6mo_AR1 <- subset(CSFs, dAICc_m1_monsoon6mo_AR1 == 0)
monsoon6mo_AR1 <- monsoon6mo_AR1[,c("code","ecosystem","dAICc_m1_monsoon6mo_AR1")]
monsoon6mo_AR1$best_model <- "m1_monsoon6mo_AR1"
names(monsoon6mo_AR1)[3] <- "dAICc"
monsoon6mo_AR1$ParamEstimate<-0

monsoon6mo_AR2 <- subset(CSFs, dAICc_m1_monsoon6mo_AR2 == 0)
monsoon6mo_AR2 <- monsoon6mo_AR2[,c("code","ecosystem","dAICc_m1_monsoon6mo_AR2")]
monsoon6mo_AR2$best_model <- "m1_monsoon6mo_AR2"
names(monsoon6mo_AR2)[3] <- "dAICc"
monsoon6mo_AR2$ParamEstimate<-0

monsoon6mo_lag <- subset(CSFs, dAICc_m1_monsoon6mo_lag == 0)
monsoon6mo_lag <- monsoon6mo_lag[,c("code","ecosystem","dAICc_m1_monsoon6mo_lag")]
monsoon6mo_lag$best_model <- "m1_monsoon6mo_lag"
names(monsoon6mo_lag)[3] <- "dAICc"
monsoon6mo_lag$ParamEstimate<-0

monsoon6mo_lag_AR1 <- subset(CSFs, dAICc_m1_monsoon6mo_lag_AR1 == 0)
monsoon6mo_lag_AR1 <- monsoon6mo_lag_AR1[,c("code","ecosystem","dAICc_m1_monsoon6mo_lag_AR1")]
monsoon6mo_lag_AR1$best_model <- "m1_monsoon6mo_lag_AR1"
names(monsoon6mo_lag_AR1)[3] <- "dAICc"
monsoon6mo_lag_AR1$ParamEstimate<-0

monsoon6mo_lag_AR2 <- subset(CSFs, dAICc_m1_monsoon6mo_lag_AR2 == 0)
monsoon6mo_lag_AR2 <- monsoon6mo_lag_AR2[,c("code","ecosystem","dAICc_m1_monsoon6mo_lag_AR2")]
monsoon6mo_lag_AR2$best_model <- "m1_monsoon6mo_lag_AR2"
names(monsoon6mo_lag_AR2)[3] <- "dAICc"
monsoon6mo_lag_AR2$ParamEstimate<-0

spring6mo <- subset(CSFs, dAICc_m1_spring6mo == 0)
spring6mo <- spring6mo[,c("code","ecosystem","dAICc_m1_spring6mo")]
spring6mo$best_model <- "m1_spring6mo"
names(spring6mo)[3] <- "dAICc"
spring6mo$ParamEstimate<-0

spring6mo_AR1 <- subset(CSFs, dAICc_m1_spring6mo_AR1 == 0)
spring6mo_AR1 <- spring6mo_AR1[,c("code","ecosystem","dAICc_m1_spring6mo_AR1")]
spring6mo_AR1$best_model <- "m1_spring6mo_AR1"
names(spring6mo_AR1)[3] <- "dAICc"
spring6mo_AR1$ParamEstimate<-0

spring6mo_AR2 <- subset(CSFs, dAICc_m1_spring6mo_AR2 == 0)
spring6mo_AR2 <- spring6mo_AR2[,c("code","ecosystem","dAICc_m1_spring6mo_AR2")]
spring6mo_AR2$best_model <- "m1_spring6mo_AR2"
names(spring6mo_AR2)[3] <- "dAICc"
spring6mo_AR2$ParamEstimate<-0

spring6mo_lag <- subset(CSFs, dAICc_m1_spring6mo_lag == 0)
spring6mo_lag <- spring6mo_lag[,c("code","ecosystem","dAICc_m1_spring6mo_lag")]
spring6mo_lag$best_model <- "m1_spring6mo_lag"
names(spring6mo_lag)[3] <- "dAICc"
spring6mo_lag$ParamEstimate<-0

spring6mo_lag_AR1 <- subset(CSFs, dAICc_m1_spring6mo_lag_AR1 == 0)
spring6mo_lag_AR1 <- spring6mo_lag_AR1[,c("code","ecosystem","dAICc_m1_spring6mo_lag_AR1")]
spring6mo_lag_AR1$best_model <- "m1_spring6mo_lag_AR1"
names(spring6mo_lag_AR1)[3] <- "dAICc"
spring6mo_lag_AR1$ParamEstimate<-0

spring6mo_lag_AR2 <- subset(CSFs, dAICc_m1_spring6mo_lag_AR2 == 0)
spring6mo_lag_AR2 <- spring6mo_lag_AR2[,c("code","ecosystem","dAICc_m1_spring6mo_lag_AR2")]
spring6mo_lag_AR2$best_model <- "m1_spring6mo_lag_AR2"
names(spring6mo_lag_AR2)[3] <- "dAICc"
spring6mo_lag_AR2$ParamEstimate<-0

# bind together the separate data frames for the different model variants, and add a "model type" column indicating that the linear relationship was best
paramest_linear <- bind_rows(monsoon6mo,monsoon6mo_AR1,monsoon6mo_AR2,monsoon6mo_lag,monsoon6mo_lag_AR1,monsoon6mo_lag_AR2,spring6mo,spring6mo_AR1,spring6mo_AR2,spring6mo_lag,spring6mo_lag_AR1,spring6mo_lag_AR2)

paramest_linear$model_type <- "m1"


### Extract the quadratic parameter estimate from the best model for each population in cases where the quadratic model was best

# create a separate data frame for each model variant, containing data for the populations for which that model variant was superior, and format the data

monsoon6mo <- subset(CSFs, dAICc_m2_monsoon6mo == 0)
monsoon6mo <- monsoon6mo[,c("code","ecosystem","dAICc_m2_monsoon6mo","ParamEst_quad_m2_monsoon6mo")]
monsoon6mo$best_model <- "m2_monsoon6mo"
names(monsoon6mo)[3] <- "dAICc"
names(monsoon6mo)[4] <- "ParamEstimate"

monsoon6mo_AR1 <- subset(CSFs, dAICc_m2_monsoon6mo_AR1 == 0)
monsoon6mo_AR1 <- monsoon6mo_AR1[,c("code","ecosystem","dAICc_m2_monsoon6mo_AR1","ParamEst_quad_m2_monsoon6mo_AR1")]
monsoon6mo_AR1$best_model <- "m2_monsoon6mo_AR1"
names(monsoon6mo_AR1)[3] <- "dAICc"
names(monsoon6mo_AR1)[4] <- "ParamEstimate"

monsoon6mo_AR2 <- subset(CSFs, dAICc_m2_monsoon6mo_AR2 == 0)
monsoon6mo_AR2 <- monsoon6mo_AR2[,c("code","ecosystem","dAICc_m2_monsoon6mo_AR2","ParamEst_quad_m2_monsoon6mo_AR2")]
monsoon6mo_AR2$best_model <- "m2_monsoon6mo_AR2"
names(monsoon6mo_AR2)[3] <- "dAICc"
names(monsoon6mo_AR2)[4] <- "ParamEstimate"

monsoon6mo_lag <- subset(CSFs, dAICc_m2_monsoon6mo_lag == 0)
monsoon6mo_lag <- monsoon6mo_lag[,c("code","ecosystem","dAICc_m2_monsoon6mo_lag","ParamEst_quad_m2_monsoon6mo_lag")]
monsoon6mo_lag$best_model <- "m2_monsoon6mo_lag"
names(monsoon6mo_lag)[3] <- "dAICc"
names(monsoon6mo_lag)[4] <- "ParamEstimate"

monsoon6mo_lag_AR1 <- subset(CSFs, dAICc_m2_monsoon6mo_lag_AR1 == 0)
monsoon6mo_lag_AR1 <- monsoon6mo_lag_AR1[,c("code","ecosystem","dAICc_m2_monsoon6mo_lag_AR1","ParamEst_quad_m2_monsoon6mo_lag_AR1")]
monsoon6mo_lag_AR1$best_model <- "m2_monsoon6mo_lag_AR1"
names(monsoon6mo_lag_AR1)[3] <- "dAICc"
names(monsoon6mo_lag_AR1)[4] <- "ParamEstimate"

monsoon6mo_lag_AR2 <- subset(CSFs, dAICc_m2_monsoon6mo_lag_AR2 == 0)
monsoon6mo_lag_AR2 <- monsoon6mo_lag_AR2[,c("code","ecosystem","dAICc_m2_monsoon6mo_lag_AR2","ParamEst_quad_m2_monsoon6mo_lag_AR2")]
monsoon6mo_lag_AR2$best_model <- "m2_monsoon6mo_lag_AR2"
names(monsoon6mo_lag_AR2)[3] <- "dAICc"
names(monsoon6mo_lag_AR2)[4] <- "ParamEstimate"

spring6mo <- subset(CSFs, dAICc_m2_spring6mo == 0)
spring6mo <- spring6mo[,c("code","ecosystem","dAICc_m2_spring6mo","ParamEst_quad_m2_spring6mo")]
spring6mo$best_model <- "m2_spring6mo"
names(spring6mo)[3] <- "dAICc"
names(spring6mo)[4] <- "ParamEstimate"

spring6mo_AR1 <- subset(CSFs, dAICc_m2_spring6mo_AR1 == 0)
spring6mo_AR1 <- spring6mo_AR1[,c("code","ecosystem","dAICc_m2_spring6mo_AR1","ParamEst_quad_m2_spring6mo_AR1")]
spring6mo_AR1$best_model <- "m2_spring6mo_AR1"
names(spring6mo_AR1)[3] <- "dAICc"
names(spring6mo_AR1)[4] <- "ParamEstimate"

spring6mo_AR2 <- subset(CSFs, dAICc_m2_spring6mo_AR2 == 0)
spring6mo_AR2 <- spring6mo_AR2[,c("code","ecosystem","dAICc_m2_spring6mo_AR2","ParamEst_quad_m2_spring6mo_AR2")]
spring6mo_AR2$best_model <- "m2_spring6mo_AR2"
names(spring6mo_AR2)[3] <- "dAICc"
names(spring6mo_AR2)[4] <- "ParamEstimate"

spring6mo_lag <- subset(CSFs, dAICc_m2_spring6mo_lag == 0)
spring6mo_lag <- spring6mo_lag[,c("code","ecosystem","dAICc_m2_spring6mo_lag","ParamEst_quad_m2_spring6mo_lag")]
spring6mo_lag$best_model <- "m2_spring6mo_lag"
names(spring6mo_lag)[3] <- "dAICc"
names(spring6mo_lag)[4] <- "ParamEstimate"

spring6mo_lag_AR1 <- subset(CSFs, dAICc_m2_spring6mo_lag_AR1 == 0)
spring6mo_lag_AR1 <- spring6mo_lag_AR1[,c("code","ecosystem","dAICc_m2_spring6mo_lag_AR1","ParamEst_quad_m2_spring6mo_lag_AR1")]
spring6mo_lag_AR1$best_model <- "m2_spring6mo_lag_AR1"
names(spring6mo_lag_AR1)[3] <- "dAICc"
names(spring6mo_lag_AR1)[4] <- "ParamEstimate"

spring6mo_lag_AR2 <- subset(CSFs, dAICc_m2_spring6mo_lag_AR2 == 0)
spring6mo_lag_AR2 <- spring6mo_lag_AR2[,c("code","ecosystem","dAICc_m2_spring6mo_lag_AR2","ParamEst_quad_m2_spring6mo_lag_AR2")]
spring6mo_lag_AR2$best_model <- "m2_spring6mo_lag_AR2"
names(spring6mo_lag_AR2)[3] <- "dAICc"
names(spring6mo_lag_AR2)[4] <- "ParamEstimate"

# bind together the separate data frames for the different model variants, and add a "model type" column indicating that the quadratic relationship was best
paramest_quadratic <- bind_rows(monsoon6mo,monsoon6mo_AR1,monsoon6mo_AR2,monsoon6mo_lag,monsoon6mo_lag_AR1,monsoon6mo_lag_AR2,spring6mo,spring6mo_AR1,spring6mo_AR2,spring6mo_lag,spring6mo_lag_AR1,spring6mo_lag_AR2)

paramest_quadratic$model_type <- "m2"


### Extract the quadratic parameter estimate from the best model for each population in cases where the cubic model was best

# create a separate data frame for each model variant, containing data for the populations for which that model variant was superior, and format the data

monsoon6mo <- subset(CSFs, dAICc_m3_monsoon6mo == 0)
monsoon6mo <- monsoon6mo[,c("code","ecosystem","dAICc_m3_monsoon6mo","ParamEst_quad_m3_monsoon6mo")]
monsoon6mo$best_model <- "m3_monsoon6mo"
names(monsoon6mo)[3] <- "dAICc"
names(monsoon6mo)[4] <- "ParamEstimate"

monsoon6mo_AR1 <- subset(CSFs, dAICc_m3_monsoon6mo_AR1 == 0)
monsoon6mo_AR1 <- monsoon6mo_AR1[,c("code","ecosystem","dAICc_m3_monsoon6mo_AR1","ParamEst_quad_m3_monsoon6mo_AR1")]
monsoon6mo_AR1$best_model <- "m3_monsoon6mo_AR1"
names(monsoon6mo_AR1)[3] <- "dAICc"
names(monsoon6mo_AR1)[4] <- "ParamEstimate"

monsoon6mo_AR2 <- subset(CSFs, dAICc_m3_monsoon6mo_AR2 == 0)
monsoon6mo_AR2 <- monsoon6mo_AR2[,c("code","ecosystem","dAICc_m3_monsoon6mo_AR2","ParamEst_quad_m3_monsoon6mo_AR2")]
monsoon6mo_AR2$best_model <- "m3_monsoon6mo_AR2"
names(monsoon6mo_AR2)[3] <- "dAICc"
names(monsoon6mo_AR2)[4] <- "ParamEstimate"

monsoon6mo_lag <- subset(CSFs, dAICc_m3_monsoon6mo_lag == 0)
monsoon6mo_lag <- monsoon6mo_lag[,c("code","ecosystem","dAICc_m3_monsoon6mo_lag","ParamEst_quad_m3_monsoon6mo_lag")]
monsoon6mo_lag$best_model <- "m3_monsoon6mo_lag"
names(monsoon6mo_lag)[3] <- "dAICc"
names(monsoon6mo_lag)[4] <- "ParamEstimate"

monsoon6mo_lag_AR1 <- subset(CSFs, dAICc_m3_monsoon6mo_lag_AR1 == 0)
monsoon6mo_lag_AR1 <- monsoon6mo_lag_AR1[,c("code","ecosystem","dAICc_m3_monsoon6mo_lag_AR1","ParamEst_quad_m3_monsoon6mo_lag_AR1")]
monsoon6mo_lag_AR1$best_model <- "m3_monsoon6mo_lag_AR1"
names(monsoon6mo_lag_AR1)[3] <- "dAICc"
names(monsoon6mo_lag_AR1)[4] <- "ParamEstimate"

monsoon6mo_lag_AR2 <- subset(CSFs, dAICc_m3_monsoon6mo_lag_AR2 == 0)
monsoon6mo_lag_AR2 <- monsoon6mo_lag_AR2[,c("code","ecosystem","dAICc_m3_monsoon6mo_lag_AR2","ParamEst_quad_m3_monsoon6mo_lag_AR2")]
monsoon6mo_lag_AR2$best_model <- "m3_monsoon6mo_lag_AR2"
names(monsoon6mo_lag_AR2)[3] <- "dAICc"
names(monsoon6mo_lag_AR2)[4] <- "ParamEstimate"

spring6mo <- subset(CSFs, dAICc_m3_spring6mo == 0)
spring6mo <- spring6mo[,c("code","ecosystem","dAICc_m3_spring6mo","ParamEst_quad_m3_spring6mo")]
spring6mo$best_model <- "m3_spring6mo"
names(spring6mo)[3] <- "dAICc"
names(spring6mo)[4] <- "ParamEstimate"

spring6mo_AR1 <- subset(CSFs, dAICc_m3_spring6mo_AR1 == 0)
spring6mo_AR1 <- spring6mo_AR1[,c("code","ecosystem","dAICc_m3_spring6mo_AR1","ParamEst_quad_m3_spring6mo_AR1")]
spring6mo_AR1$best_model <- "m3_spring6mo_AR1"
names(spring6mo_AR1)[3] <- "dAICc"
names(spring6mo_AR1)[4] <- "ParamEstimate"

spring6mo_AR2 <- subset(CSFs, dAICc_m3_spring6mo_AR2 == 0)
spring6mo_AR2 <- spring6mo_AR2[,c("code","ecosystem","dAICc_m3_spring6mo_AR2","ParamEst_quad_m3_spring6mo_AR2")]
spring6mo_AR2$best_model <- "m3_spring6mo_AR2"
names(spring6mo_AR2)[3] <- "dAICc"
names(spring6mo_AR2)[4] <- "ParamEstimate"

spring6mo_lag <- subset(CSFs, dAICc_m3_spring6mo_lag == 0)
spring6mo_lag <- spring6mo_lag[,c("code","ecosystem","dAICc_m3_spring6mo_lag","ParamEst_quad_m3_spring6mo_lag")]
spring6mo_lag$best_model <- "m3_spring6mo_lag"
names(spring6mo_lag)[3] <- "dAICc"
names(spring6mo_lag)[4] <- "ParamEstimate"

spring6mo_lag_AR1 <- subset(CSFs, dAICc_m3_spring6mo_lag_AR1 == 0)
spring6mo_lag_AR1 <- spring6mo_lag_AR1[,c("code","ecosystem","dAICc_m3_spring6mo_lag_AR1","ParamEst_quad_m3_spring6mo_lag_AR1")]
spring6mo_lag_AR1$best_model <- "m3_spring6mo_lag_AR1"
names(spring6mo_lag_AR1)[3] <- "dAICc"
names(spring6mo_lag_AR1)[4] <- "ParamEstimate"

spring6mo_lag_AR2 <- subset(CSFs, dAICc_m3_spring6mo_lag_AR2 == 0)
spring6mo_lag_AR2 <- spring6mo_lag_AR2[,c("code","ecosystem","dAICc_m3_spring6mo_lag_AR2","ParamEst_quad_m3_spring6mo_lag_AR2")]
spring6mo_lag_AR2$best_model <- "m3_spring6mo_lag_AR2"
names(spring6mo_lag_AR2)[3] <- "dAICc"
names(spring6mo_lag_AR2)[4] <- "ParamEstimate"

# bind together the separate data frames for the different model variants, and add a "model type" column indicating that the cubic relationship was best
paramest_cubic <- bind_rows(monsoon6mo,monsoon6mo_AR1,monsoon6mo_AR2,monsoon6mo_lag,monsoon6mo_lag_AR1,monsoon6mo_lag_AR2,spring6mo,spring6mo_AR1,spring6mo_AR2,spring6mo_lag,spring6mo_lag_AR1,spring6mo_lag_AR2)

paramest_cubic$model_type <- "m3"

# Combine data frames across model types 
paramest_all <- bind_rows(paramest_linear,paramest_quadratic,paramest_cubic)

# read in species list
specieslist<-read.csv("SEVBeeSpeciesList2002-2019_revised2023-07-19.csv")

# add species data to parameter estimate data frame
CSFs_paramest<-left_join(paramest_all,specieslist,by="code")

# read in list of species from thermal and desiccation tolerance trials
list_tol<-read.csv("speciessummaryctmaxdestol_2023-08-29.csv")

# filter CSF parameter estimate data frame to just contain species for which we measured thermal and/or desiccation tolerance
CSFs_paramest <- filter(CSFs_paramest,code %in% list_tol$code)

# for each genus, calculate mean and se of the quadratic parameter estimate
quad_param_summary <- CSFs_paramest %>% group_by(family, genus) %>% summarise(quad_param_mean=mean(ParamEstimate),quad_param_se=sd(ParamEstimate)/sqrt(n()))

# combine physiological tolerance and parameter estimate data frames
fulldata<-left_join(fulldata, quad_param_summary, by=c("family","genus"))


##### Add data on body mass and change in abundance over time #####

## Body mass

# read in mass data for each bee species
mass<-read.csv("SEVBeeBodyMassData_2023-08-29_forpub.csv")

# calculate average mass for each species
body_mass<-mass %>% group_by(code,family,genus,species) %>% summarise(mean_mass_mg=mean(mass_mg),se_mass_mg=sd(mass_mg)/sqrt(n()))

# create dataset of genus-level body mass averages
body_mass_genus <- body_mass %>% group_by(family,genus) %>% summarise(mean_mass=mean(mean_mass_mg, na.rm=TRUE))

# z-score mass data
body_mass_genus$mean_mass_z <- (body_mass_genus$mean_mass - mean(body_mass_genus$mean_mass))/sd(body_mass_genus$mean_mass)

# add body mass to physiological tolerance data frame
fulldata<-left_join(fulldata,body_mass_genus,by=c("family","genus"))


## Change in abundance over time

# read in bee count data for 2002-2019
bee_wide<-read.csv("SEVBeeData2002-2019_revised2023-07-19.csv")

# read in species list
specieslist<-read.csv("SEVBeeSpeciesList2002-2019_revised2023-07-19.csv")

# subset to exclude data from 2016 and 2017 due to incomplete sampling in these years
bee_wide<-subset(bee_wide,year!=2016 & year!=2017)

# convert the data to long format
bee_melt<-pivot_longer(data=bee_wide,cols=11:351,names_to="code",values_to="count")

# sum the bee counts for each year x transect x bee species combination
bee_melt_year <- bee_melt %>% group_by(year,ecosystem,transect,code) %>% summarise(count=sum(count))

# remove bee species only observed in 2016 and/or 2017
bee_melt_year <- subset(bee_melt_year, code!="APTRILUN" & code!="HALASDI57")

# read in list of species from thermal and desiccation tolerance trials
list_tol<-read.csv("speciessummaryctmaxdestol_2023-08-29.csv")

# filter data frame to just contain species for which we measured thermal and/or desiccation tolerance
bee_melt_year <- filter(bee_melt_year,code %in% list_tol$code)

# add family, genus, and species data
bee_melt_year<-left_join(bee_melt_year,specieslist,by="code")

# calculate genus sums
genus_count <- bee_melt_year %>% group_by(genus,year) %>% summarise(total_abund=sum(count))

# make genus data wide
genus_wide <- genus_count %>% pivot_wider(names_from = genus, values_from = total_abund)

# create a vector of genus names
genus<-colnames(genus_wide[,2:14])

# create a new matrix to put output of the below loop in
Year_output<-matrix(nrow=1,ncol=4)


# Loop through each genus, regressing abundance against time

for (i in 1:length(genus)){
  
  Genus_id<-genus[i]
  
  Year_reg<-lm(formula(paste(Genus_id,"~year")),data=genus_wide)
  
  # put slope, s.e., and p-value into new matrix
  
  slope<-summary(Year_reg)$coefficients[2,1]
  
  se<-summary(Year_reg)$coefficients[2,2]
  
  p_value<-summary(Year_reg)$coefficients[2,4]
  
  output1<-cbind(Genus_id,slope,se,p_value)
  
  Year_output <-rbind(Year_output,output1)
  
}

# create data frame of output, and remove blank first row
Year_output<-as.data.frame(Year_output)
Year_output<-Year_output[-1,]

# rename genus column
names(Year_output)[1]<-"genus"


# merge physiological tolerance data and change in abundance over time data
fulldata <-left_join(fulldata, Year_output, by=c("genus"))


##### Complete final data formatting and read in the bee phylogeny #####

# remove genera w/ NAs in some columns
fulldata <-fulldata[complete.cases(fulldata[ , 22]),]

# treat slope and se as numeric
fulldata$slope <- as.numeric(fulldata$slope)
fulldata$se <- as.numeric(fulldata$se)

# read in and format the genus-level bee phylogeny
sev_genus_tree<-read.tree("Sev_genus_tree_2021-09-27.tre")
plot(sev_genus_tree)

# subset the phylogeny to just include focal genera
sev_genus_tree_2<-keep.tip(sev_genus_tree,fulldata$genus)
plot(sev_genus_tree_2)



##### Phylogenetic mixed effects models: magnitude of change in abundance over time as a function of physiological or body mass metrics #####

# thermal tolerance (CTMax)
p1 <- pglmm(slope ~ mean_end_temp_z + (1|genus__), data = fulldata, family = "gaussian", cov_ranef = list(genus = sev_genus_tree_2),REML = FALSE)
p1 
R2_lik(p1)

# desiccation tolerance (duration of trial)
p1 <- pglmm(slope ~ mean_duration_hr_z + (1|genus__), data = fulldata, family = "gaussian", cov_ranef = list(genus = sev_genus_tree_2),REML = FALSE)
p1 
R2_lik(p1)

# critical water content
p1 <- pglmm(slope ~ mean_cwc_z + (1|genus__), data = fulldata, family = "gaussian", cov_ranef = list(genus = sev_genus_tree_2),REML = FALSE)
p1 
R2_lik(p1)  

# VPD
p1 <- pglmm(slope ~ mean_end_vpd_z + (1|genus__), data = fulldata, family = "gaussian", cov_ranef = list(genus = sev_genus_tree_2),REML = FALSE)
p1 
R2_lik(p1)

# combined tolerance metric
p1 <- pglmm(slope ~ combined_tol_duration_hr + (1|genus__), data = fulldata, family = "gaussian", cov_ranef = list(genus = sev_genus_tree_2),REML = FALSE)
p1 
R2_lik(p1)

# body mass
p1 <- pglmm(slope ~ mean_mass_z + (1|genus__), data = fulldata, family = "gaussian", cov_ranef = list(genus = sev_genus_tree_2),REML = FALSE)
p1  
R2_lik(p1) 


##### Phylogenetic mixed effects models: linear CSF parameter estimate as a function of physiological or body mass metrics #####

# thermal tolerance (CTMax)
p1 <- pglmm(lin_param_mean ~ mean_end_temp_z + (1|genus__), data = fulldata, family = "gaussian", cov_ranef = list(genus = sev_genus_tree_2),REML = FALSE)
p1 
R2_lik(p1)

# desiccation tolerance (duration of trial)
p1 <- pglmm(lin_param_mean ~ mean_duration_hr_z + (1|genus__), data = fulldata, family = "gaussian", cov_ranef = list(genus = sev_genus_tree_2),REML = FALSE)
p1 
R2_lik(p1)

# critical water content
p1 <- pglmm(lin_param_mean ~ mean_cwc_z + (1|genus__), data = fulldata, family = "gaussian", cov_ranef = list(genus = sev_genus_tree_2),REML = FALSE)
p1 
R2_lik(p1)  

# VPD
p1 <- pglmm(lin_param_mean ~ mean_end_vpd_z + (1|genus__), data = fulldata, family = "gaussian", cov_ranef = list(genus = sev_genus_tree_2),REML = FALSE)
p1 
R2_lik(p1)

# combined tolerance metric
p1 <- pglmm(lin_param_mean ~ combined_tol_duration_hr + (1|genus__), data = fulldata, family = "gaussian", cov_ranef = list(genus = sev_genus_tree_2),REML = FALSE)
p1 
R2_lik(p1)

# body mass
p1 <- pglmm(lin_param_mean ~ mean_mass_z + (1|genus__), data = fulldata, family = "gaussian", cov_ranef = list(genus = sev_genus_tree_2),REML = FALSE)
p1  
R2_lik(p1) 


##### Plot key trends #####

# Linear CSF parameter as a function of combined tolerance metric
csf_plot<-ggplot(data=fulldata, aes(x=combined_tol_duration_hr, y=lin_param_mean,fill=genus)) +
  geom_smooth(method=lm, color="azure4",inherit.aes=FALSE, data=fulldata, aes(x=combined_tol_duration_hr, y=lin_param_mean), size=0.75, fill="azure3") +
  geom_errorbarh(aes(xmax = combined_tol_duration_hr + combined_tol_duration_hr_se, xmin = combined_tol_duration_hr - combined_tol_duration_hr_se))+
  geom_linerange(aes(ymin=lin_param_mean-lin_param_se, ymax=lin_param_mean+lin_param_se)) +
  geom_point(size=4,shape=21) +
  xlab("Combined thermal and \ndesiccation tolerance") +
  ylab("Linear \nparameter \nfrom climate \nsensitivity \nfunction")+
  theme_classic()+
  theme(axis.text.x=element_blank(), axis.title.x=element_blank()) +
  theme(axis.text.y = element_text(color="black",size=14)) +
  theme(axis.title=element_text(size=18)) +
  scale_color_manual("black") +
  scale_fill_brewer(palette = "Paired") +
  theme(legend.text = element_text(face="italic"))+
  labs(fill="Genus")
csf_plot

# Magnitude of change in abundance a function of combined tolerance metric
y_title="<br><br>Magnitude <br>of change in <br>abundance <br>(<i>β</i><sub>time</sub>)"
slope_plot <- ggplot(data=fulldata, aes(x=combined_tol_duration_hr, y=slope,fill=genus)) +
  geom_smooth(method=lm, color="azure4",inherit.aes=FALSE, data=fulldata, aes(x=combined_tol_duration_hr, y=slope), size=0.75, fill="azure3") +
  geom_errorbarh(aes(xmax = combined_tol_duration_hr + combined_tol_duration_hr_se, xmin = combined_tol_duration_hr - combined_tol_duration_hr_se))+
  geom_linerange(aes(ymin=slope-se, ymax=slope+se)) +
  geom_point(size=4,shape=21) +
  theme_classic()+
  labs(y=y_title) +
  theme(axis.title.y = element_markdown()) +
  xlab("Combined thermal and \ndesiccation tolerance") +
  theme(axis.text.x = element_text(color="black",size=14)) +
  theme(axis.text.y = element_text(color="black",size=14)) +
  theme(axis.title=element_text(size=18)) +
  scale_color_manual("black") +
  scale_fill_brewer(palette = "Paired") +
  theme(legend.text = element_text(face="italic"))+
  theme(legend.position = "none") +
  labs(fill="Genus")
slope_plot

# Combine and save the plots
p1<-csf_plot + slope_plot + plot_layout(ncol = 1, guides = "collect")
p1
#ggsave("combined_tol_forfig2.png", p1, width=7,height=9,units = c("in"),dpi = 600)


##### Test for phylogenetic signal in physiological tolerance metrics and linear CSF parameter

# create separate data frames for each focal metric
combined_tol<-data.frame(fulldata[,"combined_tol_duration_hr"])
end_temp<-data.frame(fulldata[,"mean_end_temp"])
duration_hr<-data.frame(fulldata[,"mean_duration_hr"])
vpd<-data.frame(fulldata[,"mean_end_vpd"])
lin_param<-data.frame(fulldata[,"lin_param_mean"])

# add row names to data frames
rownames(combined_tol) <- fulldata$genus
rownames(end_temp) <- fulldata$genus
rownames(duration_hr) <- fulldata$genus
rownames(vpd) <- fulldata$genus
rownames(lin_param) <- fulldata$genus

# add node labels to tree
sev_genus_tree_2$node.label <- 1:11

# for each metric, combine the phylogeny and metric into a "phylo4d" object, and test for phylogenetic signal

# combined tolerance metric
p4d <- phylo4d(x=sev_genus_tree_2, tip.data=combined_tol, node.data=NULL)
plot(p4d)
phyloSignal(p4d = p4d, method = "all")

# thermal tolerance (CTMax)
p4d <- phylo4d(x=sev_genus_tree_2, tip.data=end_temp, node.data=NULL)
plot(p4d)
phyloSignal(p4d = p4d, method = "all")

# desiccation tolerance
p4d <- phylo4d(x=sev_genus_tree_2, tip.data=duration_hr, node.data=NULL)
plot(p4d)
phyloSignal(p4d = p4d, method = "all")

# VPD
p4d <- phylo4d(x=sev_genus_tree_2, tip.data=vpd, node.data=NULL)
plot(p4d)
phyloSignal(p4d = p4d, method = "all")

# linear CSF parameter estimate
p4d <- phylo4d(x=sev_genus_tree_2, tip.data=lin_param, node.data=NULL)
plot(p4d)
phyloSignal(p4d = p4d, method = "all")



##### Phylogenetic mixed effects models: quadratic CSF parameter estimate as a function of physiological or body mass metrics #####

# exclude Agapostemon data (outlier data point)
fulldata_noaga<-subset(fulldata,genus!="Agapostemon")

# drop Agapostemon from the phylogeny
sev_genus_tree_2<-keep.tip(sev_genus_tree,fulldata_noaga$genus)
plot(sev_genus_tree_2)

# thermal tolerance (CTMax)
p1 <- pglmm(quad_param_mean ~ mean_end_temp_z + (1|genus__), data = fulldata_noaga, family = "gaussian", cov_ranef = list(genus = sev_genus_tree_2),REML = FALSE)
p1 
R2_lik(p1)

# desiccation tolerance (duration of trial)
p1 <- pglmm(quad_param_mean ~ mean_duration_hr_z + (1|genus__), data = fulldata_noaga, family = "gaussian", cov_ranef = list(genus = sev_genus_tree_2),REML = FALSE)
p1 
R2_lik(p1)


# VPD
p1 <- pglmm(quad_param_mean ~ mean_end_vpd_z + (1|genus__), data = fulldata_noaga, family = "gaussian", cov_ranef = list(genus = sev_genus_tree_2),REML = FALSE)
p1 
R2_lik(p1)

# combined tolerance metric
p1 <- pglmm(quad_param_mean ~ combined_tol_duration_hr + (1|genus__), data = fulldata_noaga, family = "gaussian", cov_ranef = list(genus = sev_genus_tree_2),REML = FALSE)
p1 
R2_lik(p1)

# body mass
p1 <- pglmm(quad_param_mean ~ mean_mass_z + (1|genus__), data = fulldata_noaga, family = "gaussian", cov_ranef = list(genus = sev_genus_tree_2),REML = FALSE)
p1  
R2_lik(p1) 



##### Phylogenetic mixed effects models: how does physiological tolerance vary with body mass? #####

# subset the phylogeny to contain all focal genera
sev_genus_tree_2<-keep.tip(sev_genus_tree,fulldata$genus)
plot(sev_genus_tree_2)

# thermal tolerance
p1 <- pglmm(mean_end_temp_z ~ mean_mass_z + (1|genus__), data = fulldata, family = "gaussian", cov_ranef = list(genus = sev_genus_tree_2),REML = FALSE)
p1
R2_lik(p1)

# vpd
p1 <- pglmm(mean_end_vpd_z ~ mean_mass_z + (1|genus__), data = fulldata, family = "gaussian", cov_ranef = list(genus = sev_genus_tree_2),REML = FALSE)
p1
R2_lik(p1)

# des tol
p1 <- pglmm(mean_duration_hr_z ~ mean_mass_z + (1|genus__), data = fulldata, family = "gaussian", cov_ranef = list(genus = sev_genus_tree_2),REML = FALSE)
p1
R2_lik(p1)

# combined tolerance (ctmax + mean duration hr) 
p1 <- pglmm(combined_tol_duration_hr ~ mean_mass_z + (1|genus__), data = fulldata, family = "gaussian", cov_ranef = list(genus = sev_genus_tree_2),REML = FALSE)
p1
R2_lik(p1)


