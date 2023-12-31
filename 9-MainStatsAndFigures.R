################################################################################### 
# Main statistics and figures

# Heat and desiccation tolerances predict bee abundance under climate change
# Melanie R. Kazenel, Karen W. Wright, Terry Griswold, Kenneth D. Whitney, and Jennifer A. Rudgers

# Date: 2023-12-01
# Corresponding author's email: melanie.kazenel@gmail.com
################################################################################### 

# Load required packages
library(ggplot2)
library(dplyr)
library(lme4)
library(vegan)
library(reshape2)
library(car)
library(piecewiseSEM)
library(nlme)
library(MuMIn)
library(lubridate)
library(Rfast)
library(simr)
library(tidyr)
library(plantecophys)
library(ape)
library(ggtree)
library(phyr)
library(pez)
library(rr2)
library(phylosignal)
library(phylobase)
library(ggstance)
library(phytools)
library(codyn)
library(visreg)
library(cowplot)
library(textshape)
library(patchwork)
library(ggbreak)


##### CLIMATE SENSITIVITY FUNCTIONS: SUMMARIZING RESULTS #####

# Read in model output from climate sensitivity functions
CSFs <- read.csv("bee_CSFs_2023-12-01.csv")

### For how many bee populations did models run? ###
length(CSFs$code)

### What percentage of bee populations (species x ecosystem combinations) were sensitive to aridity? ###
CSFs_better_than_null<-subset(CSFs,dAICc_null!=0)
length(CSFs_better_than_null$code) # number of populations
length(CSFs_better_than_null$code)/length(CSFs$code) # percentage

### How many species did these climate-sensitive populations represent? ###
summary_species_betterthannull<-CSFs_better_than_null %>% group_by(code) %>% summarise(count=n())
length(summary_species_betterthannull$code)

### Amongst the populations that tracked aridity, how many had linear, quadratic, and cubic relationships? ###

### LINEAR

# Create a separate CSF results data frame for each model variant, containing data for the populations for which that model variant was superior, and format the data
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

# Bind together the separate data frames for the different model variants, and add a "model type" column indicating that the linear relationship was best
paramest_linear <- bind_rows(monsoon6mo,monsoon6mo_AR1,monsoon6mo_AR2,monsoon6mo_lag,monsoon6mo_lag_AR1,monsoon6mo_lag_AR2,spring6mo,spring6mo_AR1,spring6mo_AR2,spring6mo_lag,spring6mo_lag_AR1,spring6mo_lag_AR2)

paramest_linear$model_type <- "m1"

# Positive the data so that positive values indicate an increase in abundance under increasing aridity
paramest_linear$ParamEstimate<-paramest_linear$ParamEstimate*(-1)


### QUADRATIC

# Create a separate CSF results data frame for each model variant, containing data for the populations for which that model variant was superior, and format the data
monsoon6mo <- subset(CSFs, dAICc_m2_monsoon6mo == 0)
monsoon6mo <- monsoon6mo[,c("code","ecosystem","dAICc_m2_monsoon6mo","ParamEst_quad_m2_monsoon6mo","P_lin_m2_monsoon6mo","P_quad_m2_monsoon6mo")]
monsoon6mo$best_model <- "m2_monsoon6mo"
names(monsoon6mo)[3] <- "dAICc"
names(monsoon6mo)[4] <- "ParamEstimate"
names(monsoon6mo)[5] <- "p_lin"
names(monsoon6mo)[6] <- "p_quad"

monsoon6mo_AR1 <- subset(CSFs, dAICc_m2_monsoon6mo_AR1 == 0)
monsoon6mo_AR1 <- monsoon6mo_AR1[,c("code","ecosystem","dAICc_m2_monsoon6mo_AR1","ParamEst_quad_m2_monsoon6mo_AR1","P_lin_m2_monsoon6mo_AR1","P_quad_m2_monsoon6mo_AR1")]
monsoon6mo_AR1$best_model <- "m2_monsoon6mo_AR1"
names(monsoon6mo_AR1)[3] <- "dAICc"
names(monsoon6mo_AR1)[4] <- "ParamEstimate"
names(monsoon6mo_AR1)[5] <- "p_lin"
names(monsoon6mo_AR1)[6] <- "p_quad"

monsoon6mo_AR2 <- subset(CSFs, dAICc_m2_monsoon6mo_AR2 == 0)
monsoon6mo_AR2 <- monsoon6mo_AR2[,c("code","ecosystem","dAICc_m2_monsoon6mo_AR2","ParamEst_quad_m2_monsoon6mo_AR2","P_lin_m2_monsoon6mo_AR2","P_quad_m2_monsoon6mo_AR2")]
monsoon6mo_AR2$best_model <- "m2_monsoon6mo_AR2"
names(monsoon6mo_AR2)[3] <- "dAICc"
names(monsoon6mo_AR2)[4] <- "ParamEstimate"
names(monsoon6mo_AR2)[5] <- "p_lin"
names(monsoon6mo_AR2)[6] <- "p_quad"

monsoon6mo_lag <- subset(CSFs, dAICc_m2_monsoon6mo_lag == 0)
monsoon6mo_lag <- monsoon6mo_lag[,c("code","ecosystem","dAICc_m2_monsoon6mo_lag","ParamEst_quad_m2_monsoon6mo_lag","P_lin_m2_monsoon6mo_lag","P_quad_m2_monsoon6mo_lag")]
monsoon6mo_lag$best_model <- "m2_monsoon6mo_lag"
names(monsoon6mo_lag)[3] <- "dAICc"
names(monsoon6mo_lag)[4] <- "ParamEstimate"
names(monsoon6mo_lag)[5] <- "p_lin"
names(monsoon6mo_lag)[6] <- "p_quad"

monsoon6mo_lag_AR1 <- subset(CSFs, dAICc_m2_monsoon6mo_lag_AR1 == 0)
monsoon6mo_lag_AR1 <- monsoon6mo_lag_AR1[,c("code","ecosystem","dAICc_m2_monsoon6mo_lag_AR1","ParamEst_quad_m2_monsoon6mo_lag_AR1","P_lin_m2_monsoon6mo_lag_AR1","P_quad_m2_monsoon6mo_lag_AR1")]
monsoon6mo_lag_AR1$best_model <- "m2_monsoon6mo_lag_AR1"
names(monsoon6mo_lag_AR1)[3] <- "dAICc"
names(monsoon6mo_lag_AR1)[4] <- "ParamEstimate"
names(monsoon6mo_lag_AR1)[5] <- "p_lin"
names(monsoon6mo_lag_AR1)[6] <- "p_quad"

monsoon6mo_lag_AR2 <- subset(CSFs, dAICc_m2_monsoon6mo_lag_AR2 == 0)
monsoon6mo_lag_AR2 <- monsoon6mo_lag_AR2[,c("code","ecosystem","dAICc_m2_monsoon6mo_lag_AR2","ParamEst_quad_m2_monsoon6mo_lag_AR2","P_lin_m2_monsoon6mo_lag_AR2","P_quad_m2_monsoon6mo_lag_AR2")]
monsoon6mo_lag_AR2$best_model <- "m2_monsoon6mo_lag_AR2"
names(monsoon6mo_lag_AR2)[3] <- "dAICc"
names(monsoon6mo_lag_AR2)[4] <- "ParamEstimate"
names(monsoon6mo_lag_AR2)[5] <- "p_lin"
names(monsoon6mo_lag_AR2)[6] <- "p_quad"

spring6mo <- subset(CSFs, dAICc_m2_spring6mo == 0)
spring6mo <- spring6mo[,c("code","ecosystem","dAICc_m2_spring6mo","ParamEst_quad_m2_spring6mo","P_lin_m2_spring6mo","P_quad_m2_spring6mo")]
spring6mo$best_model <- "m2_spring6mo"
names(spring6mo)[3] <- "dAICc"
names(spring6mo)[4] <- "ParamEstimate"
names(spring6mo)[5] <- "p_lin"
names(spring6mo)[6] <- "p_quad"

spring6mo_AR1 <- subset(CSFs, dAICc_m2_spring6mo_AR1 == 0)
spring6mo_AR1 <- spring6mo_AR1[,c("code","ecosystem","dAICc_m2_spring6mo_AR1","ParamEst_quad_m2_spring6mo_AR1","P_lin_m2_spring6mo_AR1","P_quad_m2_spring6mo_AR1")]
spring6mo_AR1$best_model <- "m2_spring6mo_AR1"
names(spring6mo_AR1)[3] <- "dAICc"
names(spring6mo_AR1)[4] <- "ParamEstimate"
names(spring6mo_AR1)[5] <- "p_lin"
names(spring6mo_AR1)[6] <- "p_quad"

spring6mo_AR2 <- subset(CSFs, dAICc_m2_spring6mo_AR2 == 0)
spring6mo_AR2 <- spring6mo_AR2[,c("code","ecosystem","dAICc_m2_spring6mo_AR2","ParamEst_quad_m2_spring6mo_AR2", "P_lin_m2_spring6mo_AR2","P_quad_m2_spring6mo_AR2")]
spring6mo_AR2$best_model <- "m2_spring6mo_AR2"
names(spring6mo_AR2)[3] <- "dAICc"
names(spring6mo_AR2)[4] <- "ParamEstimate"
names(spring6mo_AR2)[5] <- "p_lin"
names(spring6mo_AR2)[6] <- "p_quad"

spring6mo_lag <- subset(CSFs, dAICc_m2_spring6mo_lag == 0)
spring6mo_lag <- spring6mo_lag[,c("code","ecosystem","dAICc_m2_spring6mo_lag","ParamEst_quad_m2_spring6mo_lag","P_lin_m2_spring6mo_lag","P_quad_m2_spring6mo_lag")]
spring6mo_lag$best_model <- "m2_spring6mo_lag"
names(spring6mo_lag)[3] <- "dAICc"
names(spring6mo_lag)[4] <- "ParamEstimate"
names(spring6mo_lag)[5] <- "p_lin"
names(spring6mo_lag)[6] <- "p_quad"

spring6mo_lag_AR1 <- subset(CSFs, dAICc_m2_spring6mo_lag_AR1 == 0)
spring6mo_lag_AR1 <- spring6mo_lag_AR1[,c("code","ecosystem","dAICc_m2_spring6mo_lag_AR1","ParamEst_quad_m2_spring6mo_lag_AR1","P_lin_m2_spring6mo_lag_AR1","P_quad_m2_spring6mo_lag_AR1")]
spring6mo_lag_AR1$best_model <- "m2_spring6mo_lag_AR1"
names(spring6mo_lag_AR1)[3] <- "dAICc"
names(spring6mo_lag_AR1)[4] <- "ParamEstimate"
names(spring6mo_lag_AR1)[5] <- "p_lin"
names(spring6mo_lag_AR1)[6] <- "p_quad"

spring6mo_lag_AR2 <- subset(CSFs, dAICc_m2_spring6mo_lag_AR2 == 0)
spring6mo_lag_AR2 <- spring6mo_lag_AR2[,c("code","ecosystem","dAICc_m2_spring6mo_lag_AR2","ParamEst_quad_m2_spring6mo_lag_AR2","P_lin_m2_spring6mo_lag_AR2","P_quad_m2_spring6mo_lag_AR2")]
spring6mo_lag_AR2$best_model <- "m2_spring6mo_lag_AR2"
names(spring6mo_lag_AR2)[3] <- "dAICc"
names(spring6mo_lag_AR2)[4] <- "ParamEstimate"
names(spring6mo_lag_AR2)[5] <- "p_lin"
names(spring6mo_lag_AR2)[6] <- "p_quad"

# Bind together the separate data frames for the different model variants, and add a "model type" column indicating that the quadratic relationship was best
paramest_quadratic <- bind_rows(monsoon6mo,monsoon6mo_AR1,monsoon6mo_AR2,monsoon6mo_lag,monsoon6mo_lag_AR1,monsoon6mo_lag_AR2,spring6mo,spring6mo_AR1,spring6mo_AR2,spring6mo_lag,spring6mo_lag_AR1,spring6mo_lag_AR2)

paramest_quadratic$model_type <- "m2"

# Split into separate datasets for cases in which the linear term of the quadratic model was significant vs. not
paramest_lin_sig_quadratic<-subset(paramest_quadratic,p_lin<=0.05)
paramest_lin_sig_quadratic[,5:6]<-NULL
paramest_lin_sig_quadratic$model_type <- "m2_linsig"

paramest_lin_nonsig_quadratic<-subset(paramest_quadratic,p_lin>0.05)
paramest_lin_nonsig_quadratic[,5:6]<-NULL
paramest_lin_nonsig_quadratic$model_type <- "m2_linnonsig"


### CUBIC

# Create a separate CSF results data frame for each model variant, containing data for the populations for which that model variant was superior, and format the data
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

# Bind together the separate data frames for the different model variants, and add a "model type" column indicating that the cubic relationship was best
paramest_cubic <- bind_rows(monsoon6mo,monsoon6mo_AR1,monsoon6mo_AR2,monsoon6mo_lag,monsoon6mo_lag_AR1,monsoon6mo_lag_AR2,spring6mo,spring6mo_AR1,spring6mo_AR2,spring6mo_lag,spring6mo_lag_AR1,spring6mo_lag_AR2)

paramest_cubic$model_type <- "m3"

# Positive the data so that positive values indicate an increase in abundance under increasing aridity
paramest_cubic$ParamEstimate<-paramest_cubic$ParamEstimate*(-1)


# Calculate the proportion of populations with linear, quadratic, and cubic relationships with aridity, respectively 
# linear
length(paramest_linear$code) / length(CSFs_better_than_null$code)
# quadratic, linear parameter nonsignificant
length(paramest_lin_nonsig_quadratic$code) / length(CSFs_better_than_null$code)
# quadratic, linear parameter significant
length(paramest_lin_sig_quadratic$code) / length(CSFs_better_than_null$code)
# cubic
length(paramest_cubic$code) / length(CSFs_better_than_null$code)


# Combine data frames across model types 
paramest_all <- bind_rows(paramest_linear,paramest_lin_sig_quadratic,paramest_lin_nonsig_quadratic,paramest_cubic)

# Add column related to the direction of the parameter estimate
for (i in 1:length(paramest_all[,1])){
  
  if (paramest_all$ParamEstimate[i] > 0) {
    paramest_all$Direction[i] <- "Positive" } 
  
  else if (paramest_all$ParamEstimate[i] < 0) {
    paramest_all$Direction[i] <- "Negative"}
}

# Add column indicating both model type and direction 
paramest_all$model_direction <- paste(paramest_all$model_type, "_", paramest_all$Direction)


### Graph the results (number of populations for which a given model shape was superior) ###

# Create count column
paramest_all$Count <- 1

# Get counts for each model type to be included in graph
param_summary <- paramest_all %>% group_by(best_model, model_type, model_direction) %>% summarize(total=sum(Count))

# Recode data
param_summary$model_direction<-as.factor(param_summary$model_direction)
levels(param_summary$model_direction)
param_summary$model_direction<-factor(param_summary$model_direction, levels = c("m1 _ Positive", "m1 _ Negative", "m2_linnonsig _ Positive", "m2_linnonsig _ Negative", "m2_linsig _ Positive", "m2_linsig _ Negative", "m3 _ Positive", "m3 _ Negative"))
levels(param_summary$model_direction)

# Graph the trends
p <- ggplot(aes(x=model_direction,y=total,fill=model_type), data=param_summary) +
  geom_col(width = 0.75) +
  xlab("Best model") +
  ylab("Count") +
  theme_classic()+ 
  theme(axis.text.y = element_text(color="black",size=12)) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  theme(axis.title=element_text(size=16)) +
  theme(legend.position = "none") +
  scale_fill_manual(values=c("#ff8243","#004458","#df73ff","#4c008b"))
p

# Save the graph
#ggsave("CSFdirectionmagnitude_2023-12-01.pdf", p ,width=5.5,height=3.25,units = c("in"),dpi = 600)


### For what proportion of populations did present vs. prior year's aridity best predict abundance? ###
# Count the number of populations for which each model variant was superior
summary<-paramest_all %>% group_by(best_model) %>% summarise(count=n())
# Extract rows for prior year models, and sum the counts
summary_lag<-summary[c(4:6,10:12,16:18,22:24,28:30,34:36),]
sum(summary_lag$count)
# Calculate proportion of populations for which prior year's aridity best predicted abundance
sum(summary_lag$count)/length(CSFs_better_than_null$code)
# Calculate proportion of populations for which present year's aridity best predicted abundance
1-(sum(summary_lag$count)/length(CSFs_better_than_null$code))


### For what proportion of populations did spring vs. monsoon season aridity best predict abundance? ###
# Add a column for season to the paramest_all data frame
paramest_all$season<-substr(paramest_all$best_model,4,6)
# Count the number of populations for which aridity from each season was most predictive
summary<-paramest_all %>% group_by(season) %>% summarise(count=n())
# Monsoon: 
summary$count[summary$season=="mon"]/length(CSFs_better_than_null$code)
# Spring:
summary$count[summary$season=="spr"]/length(CSFs_better_than_null$code)


### Differences among ecosystems in abundance relationships with aridity ###

# Get counts for each species of the number of different model types/directions there were across populations, and count the total number of populations for which models ran
shape_dir<-paramest_all %>% group_by(code, model_direction) %>% summarise(sum=sum(Count))
shape_dir2<-shape_dir %>% group_by(code) %>% summarise(diff_shape_dir=n(),total_pops=sum(sum))

# For how many species did models run in >1 ecosystem?
shape_dir3<-subset(shape_dir2,total_pops>1)
length(shape_dir3$code)

# For what proportion of these species did sensitivity to aridity (as indicated by model shape or direction) differ among ecosystems?
shape_dir4<-subset(shape_dir3,diff_shape_dir>1)
length(shape_dir4$code)
length(shape_dir4$code)/length(shape_dir3$code)

# Calculate the proportion of species with positive and negative relationships with aridity in each ecosystem
summary2<-paramest_all %>% group_by(ecosystem,Direction) %>% summarise(sum=sum(Count))
summary3<-paramest_all %>% group_by(ecosystem) %>% summarise(eco_total=sum(Count))
summary2<-left_join(summary2,summary3,by="ecosystem")
summary2$prop<-summary2$sum/summary2$eco_total





##### PREDICTED FUTURE CHANGE IN INDIVIDUAL SPECIES ABUNDANCES #####
## Data manipulation #####

## Format predicted future bee abundance data ##

# read in data from each GCM, and add a column to each data frame indicating the GCM
canesm<-read.csv("predicted_abundances_CanESM2-2023-12-01.csv")
canesm$gcm<-"CanESM2"

access<-read.csv("predicted_abundances_ACCESS1-0-2023-12-01.csv")
access$gcm<-"ACCESS1-0"

ccsm<-read.csv("predicted_abundances_CCSM4-2023-12-01.csv")
ccsm$gcm<-"CCSM4"

cnrm<-read.csv("predicted_abundances_CNRM-CM5-2023-12-01.csv")
cnrm$gcm<-"CNRM"

csiro<-read.csv("predicted_abundances_CSIRO-Mk3-6-0-2023-12-01.csv")
csiro$gcm<-"CSIRO"

inm<-read.csv("predicted_abundances_INM-CM4-2023-12-01.csv")
inm$gcm<-"INM"

# combine the data from the different GCMs
abund <- bind_rows(canesm,access,ccsm,cnrm,csiro,inm)

# remove NA and Inf values
abund<-subset(abund,predicted_max_abundance_per_transect>-Inf & predicted_max_abundance_per_transect<Inf)
abund$spei_value<-as.numeric(abund$spei_value)
abund<-na.omit(abund)

# remove outlier bee abundance values by excluding records that are more than five times the maximum number of bees of a single species recorded for a given transect x year combination 
# max. number recorded = 822 individuals; 822*5 = 4110
abund_subset<-subset(abund,predicted_max_abundance_per_transect<=4110)

# calculate mean and se for abundance across GCMs for each bee species x ecosystem x transect x scenario x year combination
abund_updated<-abund_subset %>% group_by(code,ecosystem,transect,scenario,year) %>% 
  summarise(mean_predicted_max_abundance_per_transect=mean(predicted_max_abundance_per_transect),
            se_predicted_max_abundance_per_transect=sd(predicted_max_abundance_per_transect)/sqrt(n()))

# replace NA values with zeros in se column
abund_updated$se_predicted_max_abundance_per_transect[is.na(abund_updated$se_predicted_max_abundance_per_transect)] <- 0

# create list of species for which we have predicted future abundance data (i.e., species for which CSFs ran)
species_future<-unique(abund_updated$code)


## Format historic bee abundance data ##

# read in data
abund_hist_wide<-read.csv("bee_wide_year_2002-2019_no2016or2017_maxabund_2023-12-01.csv")

# add "scenario" column
abund_hist_wide <- abund_hist_wide %>%  mutate(scenario="historic", .after = monsoon6SPEI_prioryear)

# pivot to long form
abund_historic_melt<-pivot_longer(data=abund_hist_wide, cols = 11:349, names_to = "code", values_to = "abund")

# subset to just include species in predicted future abundance dataset (i.e., species for which CSFs ran)
abund_historic_2_melt <- abund_historic_melt %>% filter(code %in% species_future)


## Combine historic and predicted future data frames ##

# select columns and format data frames in preparation for combining them
abund_historic_final<-abund_historic_2_melt[,c(2,3,10,1,11,12)]
abund_historic_final$abund_se<-0
abund_future<-abund_updated[,c(2:5,1,6,7)]
names(abund_future)[6:7]<-c("abund","abund_se")

# combine data frames
abund_combined<-bind_rows(abund_future,abund_historic_final)

# add station column
abund_combined$station<-abund_combined$ecosystem

# rename levels of the "station" variable
abund_combined$station<-as.factor(abund_combined$station)
levels(abund_combined$station)[levels(abund_combined$station)=="B"] <- "Blue Grama"
levels(abund_combined$station)[levels(abund_combined$station)=="C"] <- "Five Points"
levels(abund_combined$station)[levels(abund_combined$station)=="G"] <- "Five Points"

# sort the data frame
abund_combined_2<-abund_combined[order(abund_combined$code, abund_combined$ecosystem,abund_combined$transect,abund_combined$scenario,abund_combined$year),]


## Loop: Which species are predicted to increase, decrease, or remain stable in abundance over time? #####

# create a vector of species codes
species<-levels(as.factor(abund_combined_2$code))

# create a vector of scenario names
scenarios <- c("rcp2.6","rcp4.5","rcp8.5")

# create a new matrix to hold the output of the loop below
year_output<-matrix(nrow=1,ncol=5)

# LOOP through each species and scenario, regressing abundance against time
for (i in 1:length(species)){
  for (j in 1:length(scenarios)){
    
    # create species and scenario ID objects
    species_id<-species[i]
    scenario_id<-scenarios[j]
    
    # subset the data
    model_data<-subset(abund_combined_2,code==species_id)
    model_data<-subset(model_data,scenario==scenario_id | scenario=="historic")
    
    # regression abundance ~ year
    year_reg<-lm(abund~year,data=model_data)
    
    # get statistical output from the regression
    slope_year<-summary(year_reg)$coefficients[2,1]
    
    se_year<-summary(year_reg)$coefficients[2,2]
    
    p_value<-summary(year_reg)$coefficients[2,4]
    
    # bind ID and statistical output values
    output1<-cbind(species_id,scenario_id,slope_year,se_year,p_value
    )
    
    # bind species output to main output data frame
    year_output <-rbind(year_output,output1)
    
  }
  
}

# create and format data frame containing output
year_output<-as.data.frame(year_output)
names(year_output)[1:2]<-c("code","scenario")
year_output<-year_output[-1,]

# read in species list, and add species information to data frame
specieslist<-read.csv("SEVBeeSpeciesList2002-2019_revised2023-07-19.csv")
specieslist2<-left_join(specieslist[,1:5],year_output,by="code")

# write .csv file of data frame
# write.csv(specieslist2,"slopes_predicted_future_change_MeanAcrossGCMs_2023-12-01.csv",row.names = FALSE)



## Summary statistics: Which species predicted to increase, decrease, or remain stable in abundance over time? #####

# Read in slope data (calculated above)
slopes <- read.csv("slopes_predicted_future_change_MeanAcrossGCMs_2023-12-01.csv")

# RCP 2.6
# create separate data frames for species predicted to increase, decrease, or not change in abundance
decrease_2 <- subset(slopes, scenario=="rcp2.6" & slope_year<0 & p_value<=0.05)
increase_2<- subset(slopes, scenario=="rcp2.6" & slope_year>0 & p_value<=0.05)
no_change_2 <- subset(slopes, scenario=="rcp2.6" & p_value>0.05)

# create single data frame with counts of species predicted to increase, decrease, or not change in abundance
df<-data.frame(rbind(c("decrease",nrow(decrease_2)),c("increase",nrow(increase_2)),c("no_change",nrow(no_change_2))))
names(df)<-c("direction","count")
df$scenario<-"rcp2.6"


# RCP 4.5
# create separate data frames for species predicted to increase, decrease, or not change in abundance
decrease_4 <- subset(slopes, scenario=="rcp4.5" & slope_year<0 & p_value<=0.05)
increase_4<- subset(slopes, scenario=="rcp4.5" & slope_year>0 & p_value<=0.05)
no_change_4 <- subset(slopes, scenario=="rcp4.5" & p_value>0.05)

# create single data frame with counts of species predicted to increase, decrease, or not change in abundance
df2<-data.frame(rbind(c("decrease",nrow(decrease_4)),c("increase",nrow(increase_4)),c("no_change",nrow(no_change_4))))
names(df2)<-c("direction","count")
df2$scenario<-"rcp4.5"
df<-rbind(df,df2)

# RCP 8.5
# create separate data frames for species predicted to increase, decrease, or not change in abundance
decrease_8 <- subset(slopes, scenario=="rcp8.5" & slope_year<0 & p_value<=0.05)
increase_8<- subset(slopes, scenario=="rcp8.5" & slope_year>0 & p_value<=0.05)
no_change_8 <- subset(slopes, scenario=="rcp8.5" & p_value>0.05)

# create single data frame with counts of species predicted to increase, decrease, or not change in abundance
df2<-data.frame(rbind(c("decrease",nrow(decrease_8)),c("increase",nrow(increase_8)),c("no_change",nrow(no_change_8))))
names(df2)<-c("direction","count")
df2$scenario<-"rcp8.5"
df<-rbind(df,df2)

# Calculate proportion of species predicted to increase, decrease, or not change in each scenario
df$proportion=as.numeric(df$count)/243

# Format for table in manuscript
df$count<-NULL
df_wide<-pivot_wider(df,names_from="direction",values_from="proportion")
df_wide<-df_wide %>% mutate(GCM="Mean Across GCMs", .before="scenario") 
#write.csv(df_wide,"table_predicted_abundance_change_MeanAcrossGCMs_2023-12-01.csv",row.names=FALSE)


## Map climate change winners vs. losers onto the bee phylogeny #####

# read in and plot the genus-level tree
sev_genus_tree<-read.tree("Sev_genus_tree_2021-09-27.tre")
plot(sev_genus_tree)

# read in species list, and remove species only observed in 2016 and/or 2017
specieslist<-read.csv("SEVBeeSpeciesList2002-2019_revised2023-07-19.csv")
specieslist2 <- subset(specieslist, code!="APTRILUN" & code!="HALASDI57")

# rename one genus in the dataset for consistency with the tree
specieslist2$genus[specieslist2$genus=="Lithurgopsis"] <- "Lithurgus"

# make a data frame to hold species codes for subsequent steps
species_codes<-specieslist2

# add genus_species column
species_codes$genus_species<-paste(species_codes$genus,species_codes$species,sep="_")

# turn the phylogeny into a species level phylogeny by adding fake species epithets specified with "_zzz" to each genus
sev_genus_tree$tip.label = paste0(sev_genus_tree$tip.label,'_zzz')
# add species to the tree as polytomies
sev_species_tree<-congeneric.merge(species_codes$genus_species, tree = sev_genus_tree, split = "_")
plot(sev_species_tree)

# drop tips not included in dataset
sev_species_tree2<-keep.tip(phy=sev_species_tree,tip=species_codes$genus_species)
plot(sev_species_tree2)

# recode genus labels on tree
sev_species_tree2$genus<-species_codes$genus

# plot a circular tree and turn it into a cladogram via setting branch lengths to "none"
circ_spp_tree2 <- ggtree(sev_species_tree2, layout = "circular", branch.length="none") + 
  geom_tiplab2(size = 2, offset =1.5, aes(angle = angle))
plot(circ_spp_tree2)

# plot the tree with genus rather than species tip labels
circ_spp_tree <- ggtree(sev_species_tree2, layout = "circular", branch.length="none") + 
  theme(legend.position = "none") +
  geom_cladelab(node=353, label="Agapostemon", offset=1.5, angle="auto",fontface=3) +
  geom_cladelab(node=363, label="Andrena", offset=1.5, angle="auto",fontface=3) +
  geom_cladelab(node=429, label="Anthidium", offset=1.5, angle="auto",fontface=3) +
  geom_cladelab(node=378, label="Anthophora", offset=1.5, angle="auto",fontface=3) +
  geom_cladelab(node=405, label="Anthophorula", offset=1.5, angle="auto",fontface=3) +
  geom_cladelab(node=102, label="Apis", offset=1.5, angle="auto",fontface=3) +
  geom_cladelab(node=424, label="Ashmeadiella", offset=1.5, angle="auto",fontface=3) +
  geom_cladelab(node=425, label="Atoposmia", offset=1.5, angle="auto",fontface=3) +
  geom_cladelab(node=206, label="Augochlorella", offset=1.5, angle="auto",fontface=3) +
  geom_cladelab(node=409, label="Bombus", offset=1.5, angle="auto",fontface=3) +
  geom_cladelab(node=373, label="Calliopsis", offset=1.5, angle="auto",fontface=3) +
  geom_cladelab(node=186, label="Caupolicana", offset=1.5, angle="auto",fontface=3) +
  geom_cladelab(node=410, label="Centris", offset=1.5, angle="auto",fontface=3) +
  geom_cladelab(node=110, label="Ceratina", offset=1.5, angle="auto",fontface=3) +
  geom_cladelab(node=419, label="Coelioxys", offset=1.5, angle="auto",fontface=3) +
  geom_cladelab(node=347, label="Colletes", offset=1.5, angle="auto",fontface=3) +
  geom_cladelab(node=207, label="Conanthalictus", offset=1.5, angle="auto",fontface=3) +
  geom_cladelab(node=403, label="Diadasia", offset=1.5, angle="auto",fontface=3) +
  geom_cladelab(node=431, label="Dianthidium", offset=1.5, angle="auto",fontface=3) +
  geom_cladelab(node=414, label="Dioxys", offset=1.5, angle="auto",fontface=3) +
  geom_cladelab(node=361, label="Dufourea", offset=1.5, angle="auto",fontface=3) +
  geom_cladelab(node=387, label="Epeolus", offset=1.5, angle="auto",fontface=3) +
  geom_cladelab(node=122, label="Ericrocis", offset=1.5, angle="auto",fontface=3) +
  geom_cladelab(node=398, label="Eucera", offset=1.5, angle="auto",fontface=3) +
  geom_cladelab(node=132, label="Exomalopsis", offset=1.5, angle="auto",fontface=3) +
  geom_cladelab(node=377, label="Habropoda", offset=1.5, angle="auto",fontface=3) +
  geom_cladelab(node=357, label="Halictus", offset=1.5, angle="auto",fontface=3) +
  geom_cladelab(node=433, label="Hesperapis", offset=1.5, angle="auto",fontface=3) +
  geom_cladelab(node=388, label="Holcopasites", offset=1.5, angle="auto",fontface=3) +
  geom_cladelab(node=426, label="Hoplitis", offset=1.5, angle="auto",fontface=3) +
  geom_cladelab(node=346, label="Hylaeus", offset=1.5, angle="auto",fontface=3) +
  geom_cladelab(node=358, label="Lasioglossum", offset=1.5, angle="auto",fontface=3) +
  geom_cladelab(node=298, label="Lithurgopsis", offset=1.5, angle="auto",fontface=3) +
  geom_cladelab(node=368, label="Macrotera", offset=1.5, angle="auto",fontface=3) +
  geom_cladelab(node=137, label="Martinapis", offset=1.5, angle="auto",fontface=3) +
  geom_cladelab(node=418, label="Megachile", offset=1.5, angle="auto",fontface=3) +
  geom_cladelab(node=393, label="Melecta", offset=1.5, angle="auto",fontface=3) +
  geom_cladelab(node=402, label="Melissodes", offset=1.5, angle="auto",fontface=3) +
  geom_cladelab(node=163, label="Neolarra", offset=1.5, angle="auto",fontface=3) +
  geom_cladelab(node=164, label="Nomada", offset=1.5, angle="auto",fontface=3) +
  geom_cladelab(node=350, label="Nomia", offset=1.5, angle="auto",fontface=3) +
  geom_cladelab(node=422, label="Osmia", offset=1.5, angle="auto",fontface=3) +
  geom_cladelab(node=369, label="Perdita", offset=1.5, angle="auto",fontface=3) +
  geom_cladelab(node=371, label="Protandrena", offset=1.5, angle="auto",fontface=3) +
  geom_cladelab(node=82, label="Protoxaea", offset=1.5, angle="auto",fontface=3) +
  geom_cladelab(node=372, label="Pseudopanurgus", offset=1.5, angle="auto",fontface=3) +
  geom_cladelab(node=355, label="Sphecodes", offset=1.5, angle="auto",fontface=3) +
  geom_cladelab(node=263, label="Sphecodosoma", offset=1.5, angle="auto",fontface=3) +
  geom_cladelab(node=432, label="Stelis", offset=1.5, angle="auto",fontface=3) +
  geom_cladelab(node=401, label="Svastra", offset=1.5, angle="auto",fontface=3) +
  geom_cladelab(node=170, label="Townsendiella", offset=1.5, angle="auto",fontface=3) +
  geom_cladelab(node=337, label="Trachusa", offset=1.5, angle="auto",fontface=3) +
  geom_cladelab(node=386, label="Triepeolus", offset=1.5, angle="auto",fontface=3) +
  geom_cladelab(node=182, label="Xenoglossa", offset=1.5, angle="auto",fontface=3) +
  geom_cladelab(node=183, label="Xeromelecta", offset=1.5, angle="auto",fontface=3) +
  geom_cladelab(node=184, label="Xylocopa", offset=1.5, angle="auto",fontface=3) +
  geom_cladelab(node=185, label="Zacosmia", offset=1.5, angle="auto",fontface=3)

plot(circ_spp_tree)

# Add a heatmap to the tree indicating whether each species is predicted to increase, decrease, or not chnage in abundance over time

# read in predicted future abundance change data
slopes <- read.csv("slopes_predicted_future_change_MeanAcrossGCMs_2023-12-01.csv")

# format the data for mapping onto the phylogeny
decrease_4 <- subset(slopes, scenario=="rcp4.5" & slope_year<0 & p_value<=0.05)
decrease_4$direction_num<-"Decrease"
increase_4<- subset(slopes, scenario=="rcp4.5" & slope_year>0 & p_value<=0.05)
increase_4$direction_num<-"Increase"
no_change_4 <- subset(slopes, scenario=="rcp4.5" & p_value>0.05)
no_change_4$direction_num<-"No change"

slopes2<-rbind(decrease_4,increase_4,no_change_4)

slopes2$genus_species<-paste(slopes2$genus,slopes2$species,sep = "_")

slopes_for_tree_4.5 <- data.frame(slopes2[,"direction_num"])
rownames(slopes_for_tree_4.5) <- slopes2$genus_species

# plot the tree with genus tip labels and colors indicating predicted change in abundance
circ_spp_tree_heatmap <-gheatmap(circ_spp_tree, slopes_for_tree_4.5, width=0.05, color = "black", colnames = FALSE) + 
  theme(legend.position = c(0.38,0.595)) +
  scale_fill_manual(na.value="white",name = "Predicted change in \nabundance over time",values=c("orangered3","blue","gray82")) +
  theme(legend.text=element_text(size=14), legend.title=element_text(size=15))
circ_spp_tree_heatmap

# save the plot. 
# save_plot("bee_phylogeny_winnerslosers_MeanAcrossGCMs_2023-12-01_withlegend.pdf", circ_spp_tree_heatmap, base_width = 13, base_height = 13)

# plot the tree with species tip labels and colors indicating predicted change in abundance
circ_spp_tree_heatmap2 <-gheatmap(circ_spp_tree2, slopes_for_tree_4.5, width=0.05, color = "black", colnames = FALSE) + 
  theme(legend.position = c(0.38,0.595)) +
  scale_fill_manual(na.value="white",name = "Predicted change in \nabundance over time",values=c("orangered3","blue","gray82"))+
  theme(legend.text=element_text(size=14), legend.title=element_text(size=15))

# rename one species
circ_spp_tree_heatmap2$data$label[circ_spp_tree_heatmap2$data$label=="Lithurgus_apicalis"] <- "Lithurgopsis_apicalis"

# examine the tree
circ_spp_tree_heatmap2

# save the plot
# save_plot("bee_phylogeny_winnerslosers_MeanAcrossGCMs_speciesnames_2023-12-01_withlegend.jpg", circ_spp_tree_heatmap2, base_width = 13, base_height = 13,dpi=600)


## Do life history traits predict "winners" vs. "losers" under future climate change? #####

# Read in required data
# life history trait data:
lhtraits<-read.csv("lhtraits_2023-08-29_forpub.csv")
# predicted future change in abundance data:
predslopes<-read.csv("slopes_predicted_future_change_MeanAcrossGCMs_2023-12-01.csv")

# Subset predicted change in abundance data to just include points from the RCP 4.5 climate scenario, and create separate data frames for species predicted to increase, decrease, and not change in abundance over time, adding a column indicating direction of predicted change
decrease_4 <- subset(predslopes, scenario=="rcp4.5" & slope_year<0 & p_value<=0.05)
decrease_4$direction_num<-"Decrease"
increase_4<- subset(predslopes, scenario=="rcp4.5" & slope_year>0 & p_value<=0.05)
increase_4$direction_num<-"Increase"
no_change_4 <- subset(predslopes, scenario=="rcp4.5" & p_value>0.05)
no_change_4$direction_num<-"No change"

# Bind the separate data frames back together
predslopes2<-rbind(decrease_4,increase_4,no_change_4)

# Join the trait data to the predicted change in abundance data
focaldata<-left_join(predslopes2,lhtraits,by=c("family","genus","subgenus","species"))

# Does diet breadth predict winner vs. loser status? ###
# Remove species for which we lack trait data
diet<-filter(focaldata,diet_breadth!="unknown")
unique(diet$diet_breadth)
# Format data for analysis
diet$diet_breadth[diet$diet_breadth == "polylectic?"] <- "polylectic"
# Chi-squared test
chisq.test(y=diet$direction_num,x=diet$diet_breadth, simulate.p.value = TRUE)

# Does sociality predict winner vs. loser status? ###
# Remove species for which we lack trait data
soc<-filter(focaldata,sociality!="unknown")
unique(soc$sociality)
# Chi-squared test
chisq.test(y=soc$direction_num,x=soc$sociality, simulate.p.value = TRUE)

# Does overwintering stage predict winner vs. loser status? ###
# Remove species for which we lack trait data
over<-filter(focaldata,overwintering_stage=="adult" | overwintering_stage=="prepupae")
unique(over$overwintering_stage)
# Chi-squared test
chisq.test(y=over$direction_num,x=over$overwintering_stage, simulate.p.value = TRUE)


## For what proportion of species might future increases in climate variability buffer against declines? #####

# read in predicted future abundance change data
slopes <- read.csv("slopes_predicted_future_change_MeanAcrossGCMs_2023-12-01.csv")

# select data from RCP 4.5, and add a column indicating whether each species is predicted to increase, decrease, or not change in abundance over time
decrease_4 <- subset(slopes, scenario=="rcp4.5" & slope_year<0 & p_value<=0.05)
decrease_4$direction_num<-"Decrease"
increase_4<- subset(slopes, scenario=="rcp4.5" & slope_year>0 & p_value<=0.05)
increase_4$direction_num<-"Increase"
no_change_4 <- subset(slopes, scenario=="rcp4.5" & p_value>0.05)
no_change_4$direction_num<-"No change"
slopes2<-rbind(decrease_4,increase_4,no_change_4)

# read in climate sensitivity function results 
CSFs <- read.csv("bee_CSFs_2023-12-01.csv")

# create new dataset for calculations
CSFsquad <- CSFs

# add a column containing the AICc value for the best quadratic model for each population
CSFsquad$minAICcquad <- apply(CSFsquad[,c("dAICc_m2_monsoon6mo", "dAICc_m2_monsoon6mo_lag", "dAICc_m2_spring6mo", "dAICc_m2_spring6mo_lag","dAICc_m2_spring6mo_AR1", "dAICc_m2_spring6mo_AR2","dAICc_m2_monsoon6mo_AR1","dAICc_m2_monsoon6mo_AR2","dAICc_m2_spring6mo_lag_AR1", "dAICc_m2_spring6mo_lag_AR2","dAICc_m2_monsoon6mo_lag_AR1","dAICc_m2_monsoon6mo_lag_AR2")], 1, FUN=min)

# get all cases in which a given quadratic model was best, and create a data frame with the relevant parameters from that model, repeating this process for all quadratic model variants
monsoon6mo <- subset(CSFsquad, dAICc_m2_monsoon6mo == minAICcquad)
monsoon6mo <- monsoon6mo[,c("code","ecosystem","dAICc_m2_monsoon6mo","minAICcquad","ParamEst_quad_m2_monsoon6mo")]
monsoon6mo$best_model <- "monsoon6mo"
names(monsoon6mo)[3] <- "dAICc"
names(monsoon6mo)[4] <- "minAICc"
names(monsoon6mo)[5] <- "ParamEst_quad"

monsoon6mo_AR1 <- subset(CSFsquad, dAICc_m2_monsoon6mo_AR1 == minAICcquad)
monsoon6mo_AR1 <- monsoon6mo_AR1[,c("code","ecosystem","dAICc_m2_monsoon6mo_AR1","minAICcquad","ParamEst_quad_m2_monsoon6mo_AR1")]
monsoon6mo_AR1$best_model <- "monsoon6mo_AR1"
names(monsoon6mo_AR1)[3] <- "dAICc"
names(monsoon6mo_AR1)[4] <- "minAICc"
names(monsoon6mo_AR1)[5] <- "ParamEst_quad"

monsoon6mo_AR2 <- subset(CSFsquad, dAICc_m2_monsoon6mo_AR2 == minAICcquad)
monsoon6mo_AR2 <- monsoon6mo_AR2[,c("code","ecosystem","dAICc_m2_monsoon6mo_AR2","minAICcquad","ParamEst_quad_m2_monsoon6mo_AR2")]
monsoon6mo_AR2$best_model <- "monsoon6mo_AR2"
names(monsoon6mo_AR2)[3] <- "dAICc"
names(monsoon6mo_AR2)[4] <- "minAICc"
names(monsoon6mo_AR2)[5] <- "ParamEst_quad"

monsoon6mo_lag <- subset(CSFsquad, dAICc_m2_monsoon6mo_lag == minAICcquad)
monsoon6mo_lag <- monsoon6mo_lag[,c("code","ecosystem","dAICc_m2_monsoon6mo_lag","minAICcquad","ParamEst_quad_m2_monsoon6mo_lag")]
monsoon6mo_lag$best_model <- "monsoon6mo_lag"
names(monsoon6mo_lag)[3] <- "dAICc"
names(monsoon6mo_lag)[4] <- "minAICc"
names(monsoon6mo_lag)[5] <- "ParamEst_quad"

monsoon6mo_lag_AR1 <- subset(CSFsquad, dAICc_m2_monsoon6mo_lag_AR1 == minAICcquad)
monsoon6mo_lag_AR1 <- monsoon6mo_lag_AR1[,c("code","ecosystem","dAICc_m2_monsoon6mo_lag_AR1","minAICcquad","ParamEst_quad_m2_monsoon6mo_lag_AR1")]
monsoon6mo_lag_AR1$best_model <- "monsoon6mo_lag_AR1"
names(monsoon6mo_lag_AR1)[3] <- "dAICc"
names(monsoon6mo_lag_AR1)[4] <- "minAICc"
names(monsoon6mo_lag_AR1)[5] <- "ParamEst_quad"

monsoon6mo_lag_AR2 <- subset(CSFsquad, dAICc_m2_monsoon6mo_lag_AR2 == minAICcquad)
monsoon6mo_lag_AR2 <- monsoon6mo_lag_AR2[,c("code","ecosystem","dAICc_m2_monsoon6mo_lag_AR2","minAICcquad","ParamEst_quad_m2_monsoon6mo_lag_AR2")]
monsoon6mo_lag_AR2$best_model <- "monsoon6mo_lag_AR2"
names(monsoon6mo_lag_AR2)[3] <- "dAICc"
names(monsoon6mo_lag_AR2)[4] <- "minAICc"
names(monsoon6mo_lag_AR2)[5] <- "ParamEst_quad"

spring6mo <- subset(CSFsquad, dAICc_m2_spring6mo == minAICcquad)
spring6mo <- spring6mo[,c("code","ecosystem","dAICc_m2_spring6mo","minAICcquad","ParamEst_quad_m2_spring6mo")]
spring6mo$best_model <- "spring6mo"
names(spring6mo)[3] <- "dAICc"
names(spring6mo)[4] <- "minAICc"
names(spring6mo)[5] <- "ParamEst_quad"

spring6mo_AR1 <- subset(CSFsquad, dAICc_m2_spring6mo_AR1 == minAICcquad)
spring6mo_AR1 <- spring6mo_AR1[,c("code","ecosystem","dAICc_m2_spring6mo_AR1","minAICcquad","ParamEst_quad_m2_spring6mo_AR1")]
spring6mo_AR1$best_model <- "spring6mo_AR1"
names(spring6mo_AR1)[3] <- "dAICc"
names(spring6mo_AR1)[4] <- "minAICc"
names(spring6mo_AR1)[5] <- "ParamEst_quad"

spring6mo_AR2 <- subset(CSFsquad, dAICc_m2_spring6mo_AR2 == minAICcquad)
spring6mo_AR2 <- spring6mo_AR2[,c("code","ecosystem","dAICc_m2_spring6mo_AR2","minAICcquad","ParamEst_quad_m2_spring6mo_AR2")]
spring6mo_AR2$best_model <- "spring6mo_AR2"
names(spring6mo_AR2)[3] <- "dAICc"
names(spring6mo_AR2)[4] <- "minAICc"
names(spring6mo_AR2)[5] <- "ParamEst_quad"

spring6mo_lag <- subset(CSFsquad, dAICc_m2_spring6mo_lag == minAICcquad)
spring6mo_lag <- spring6mo_lag[,c("code","ecosystem","dAICc_m2_spring6mo_lag","minAICcquad","ParamEst_quad_m2_spring6mo_lag")]
spring6mo_lag$best_model <- "spring6mo_lag"
names(spring6mo_lag)[3] <- "dAICc"
names(spring6mo_lag)[4] <- "minAICc"
names(spring6mo_lag)[5] <- "ParamEst_quad"

spring6mo_lag_AR1 <- subset(CSFsquad, dAICc_m2_spring6mo_lag_AR1 == minAICcquad)
spring6mo_lag_AR1 <- spring6mo_lag_AR1[,c("code","ecosystem","dAICc_m2_spring6mo_lag_AR1","minAICcquad","ParamEst_quad_m2_spring6mo_lag_AR1")]
spring6mo_lag_AR1$best_model <- "spring6mo_lag_AR1"
names(spring6mo_lag_AR1)[3] <- "dAICc"
names(spring6mo_lag_AR1)[4] <- "minAICc"
names(spring6mo_lag_AR1)[5] <- "ParamEst_quad"

spring6mo_lag_AR2 <- subset(CSFsquad, dAICc_m2_spring6mo_lag_AR2 == minAICcquad)
spring6mo_lag_AR2 <- spring6mo_lag_AR2[,c("code","ecosystem","dAICc_m2_spring6mo_lag_AR2","minAICcquad","ParamEst_quad_m2_spring6mo_lag_AR2")]
spring6mo_lag_AR2$best_model <- "spring6mo_lag_AR2"
names(spring6mo_lag_AR2)[3] <- "dAICc"
names(spring6mo_lag_AR2)[4] <- "minAICc"
names(spring6mo_lag_AR2)[5] <- "ParamEst_quad"

# bind the separate data frames together
paramest <- bind_rows(monsoon6mo,monsoon6mo_AR1,monsoon6mo_AR2,monsoon6mo_lag,monsoon6mo_lag_AR1,monsoon6mo_lag_AR2,spring6mo,spring6mo_AR1,spring6mo_AR2,spring6mo_lag,spring6mo_lag_AR1,spring6mo_lag_AR2)

# calculate the mean quadratic parameter estimate for each species (average across ecosystems)
summaryquad<-paramest %>% group_by(code) %>% summarise(mean_quad_param=mean(ParamEst_quad))

# join the change in abundance over time and quadratic parameter estimate datasets
slopes3<-left_join(slopes2,summaryquad,by="code")

# filter the dataset to just include species with positive quadratic parameter estimates that are predicted to not decrease in abundance over time
subset2<-filter(slopes3, direction_num!="Decrease" & mean_quad_param>0)


##### PHYLOGENETIC SIGNAL TESTS #####

## Test for phylogenetic signal in change in abundance over time  #####

# subset change in abundance data frame to just contain future climate scenario of interest
slopes_4.5<-subset(slopes,scenario=="rcp4.5")

# rename one genus in the dataset for consistency with the tree
slopes_4.5$genus[slopes_4.5$genus=="Lithurgopsis"] <- "Lithurgus"

# add genus_species column
slopes_4.5$genus_species<-paste(slopes_4.5$genus,slopes_4.5$species,sep="_")

# drop tips from the phylogeny that are not included in slopes dataset
sev_species_tree_3<-keep.tip(phy=sev_species_tree2,tip=slopes_4.5$genus_species)
plot(sev_species_tree_3)

# format change in abundance data as a single column with row names
slopes_for_physig_4.5<-data.frame(slopes_4.5[,"slope_year"])
rownames(slopes_for_physig_4.5) <- slopes_4.5$genus_species

# add node numbers to the tree
sev_species_tree_3$node.label <- 1:85

# combine the phylogeny and trait data into a "phylo4d" object, to be used with the phylosignal package
p4d <- phylo4d(x=sev_species_tree_3, tip.data=slopes_for_physig_4.5, node.data=NULL)

# visualize data being mapped onto the phylogeny
barplot.phylo4d(p4d, tree.type = "phylo")

# test for phylogenetic signal
phyloSignal(p4d = p4d, method = "all")


##### PREDICTED FUTURE CHANGE IN TOTAL ABUNDANCE ACROSS SPECIES #####
## How will total bee abundance change over time? #####

# calculate total abundance across bee species for each ecosystem x year combination in the historic and projected future data
abund_summary <- abund_combined_2 %>% group_by(ecosystem,scenario,year) %>% summarise(total_abund=sum(abund))

# remove very high outlier values
abund_summary<-subset(abund_summary,total_abund<10000)


# add climate data:

# read in climate data from each GCM, and add a column indicating GCM identity
clim_canesm<-read.csv("spei_historic_and_future_byscenario_CanESM2_2023-12-01.csv")
clim_canesm$gcm<-"CanESM2"

clim_access<-read.csv("spei_historic_and_future_byscenario_ACCESS1-0_2023-12-01.csv")
clim_access$gcm<-"ACCESS1-0"

clim_ccsm<-read.csv("spei_historic_and_future_byscenario_CCSM4_2023-12-01.csv")
clim_ccsm$gcm<-"CCSM4"

clim_cnrm<-read.csv("spei_historic_and_future_byscenario_CNRM-CM5_2023-12-01.csv")
clim_cnrm$gcm<-"CNRM"

clim_csiro<-read.csv("spei_historic_and_future_byscenario_CSIRO-Mk3-6-0_2023-12-01.csv")
clim_csiro$gcm<-"CSIRO"

clim_inm<-read.csv("spei_historic_and_future_byscenario_INM-CM4_2023-12-01.csv")
clim_inm$gcm<-"INM"

# bind climate data from the GCMs together
clim_all<-bind_rows(clim_canesm,
                    clim_access,
                    clim_ccsm,
                    clim_cnrm,
                    clim_csiro,
                    clim_inm)

# calculate mean and se monsoon aridity across GCMs
clim_avg <- clim_all %>% group_by(year,source,station) %>%
  summarise(mean_monsoon6SPEI=mean(monsoon6SPEI),
            se_monsoon6SPEI=sd(monsoon6SPEI)/sqrt(n()))

# replace NA values with zeros in se column
clim_avg$se_monsoon6SPEI[is.na(clim_avg$se_monsoon6SPEI)] <- 0

# add station column
abund_summary$station<-abund_summary$ecosystem

# rename levels of the "station" variable
abund_summary$station<-as.factor(abund_summary$station)
levels(abund_summary$station)[levels(abund_summary$station)=="B"] <- "Blue Grama"
levels(abund_summary$station)[levels(abund_summary$station)=="C"] <- "Five Points"
levels(abund_summary$station)[levels(abund_summary$station)=="G"] <- "Five Points"

# create data frame for historic data, and add climate data to it
colnames(abund_summary)
colnames(clim_avg)
historic<-subset(abund_summary,year<=2020)
clim_historic<-subset(clim_avg,source=="SEV-met")
historic2<-left_join(historic,clim_historic,by=c("year","station"))
historic2$source<-NULL

# create data frame for future bee data, and add climate data to it
future<-subset(abund_summary,year>2020)
clim_future<-subset(clim_avg,source!="SEV-met")
names(clim_future)[2]<-"scenario"
future2<-left_join(future,clim_future,by=c("year","scenario","station"))

# combine the data frames
colnames(historic2)
colnames(future2)
abund_summary<-bind_rows(historic2,
                           future2)

# sort the resulting data frame
abund_summary<-abund_summary[order(abund_summary$ecosystem,abund_summary$scenario,abund_summary$year),]


# Plot predicted change in total abundance over time for RCP 4.5, with points colored by SPEI

# format data
rcp4.5<-subset(abund_summary, scenario=="historic" | scenario=="rcp4.5")
rcp4.5$ecosystem<-as.factor(rcp4.5$ecosystem)
levels(rcp4.5$ecosystem) <- c("Plains grassland","Chihuahuan Desert shrubland","Chihuahuan Desert grassland")
rcp4.5$scenario<-as.factor(rcp4.5$scenario)
levels(rcp4.5$scenario)<-c("long-term historic","predicted future")

# plot the trend
plot_abundchange<-ggplot(data=rcp4.5, aes(x=year,y=total_abund)) + 
  geom_point(size=3,aes(fill=mean_monsoon6SPEI*(-1),shape=scenario))+ 
  scale_shape_manual(values=c(22,21)) +
  xlab("Year") + ylab("Total abundance") + 
  theme_classic() + theme(axis.text.x = element_text(size=16))+ theme(axis.text.y = element_text(size=16)) + theme(legend.text=element_text(size=13), legend.title=element_text(size=15))+
  theme(axis.title.x = element_text(size=20))+ theme(axis.title.y = element_text(size=20)) +
  geom_smooth(method="lm",color="black") + 
  scale_fill_gradient2(high = "red3", low = "royalblue",mid = "white",midpoint = 0) + labs(fill="Aridity index", shape="Data type") +
  annotate("text", x = 2075, y = 3000, label = expression(atop(italic(P)*" = 0.1127",italic(R^2)=="0.14")),size=4.5) +
  scale_y_continuous(breaks = c(0,500,1500,2500,3500)) +
  theme(plot.margin = unit(c(0,0,0,20), "pt"))
plot_abundchange


# Plot predicted change in abundance over time for RCP 2.6 and RCP 8.5, with points colored by SPEI

# format data
rcp2.6<-subset(abund_summary, scenario=="historic" | scenario=="rcp2.6")
rcp2.6$scenario2<-"rcp2.6"

rcp8.5<-subset(abund_summary, scenario=="historic" | scenario=="rcp8.5")
rcp8.5$scenario2<-"rcp8.5"

rcp2_8<-bind_rows(rcp2.6,rcp8.5)

rcp2_8$ecosystem<-as.factor(rcp2_8$ecosystem)
levels(rcp2_8$ecosystem) <- c("Plains grassland","Chihuahuan Desert shrubland","Chihuahuan Desert grassland")
rcp2_8$scenario<-as.factor(rcp2_8$scenario)
rcp2_8$scenario2<-as.factor(rcp2_8$scenario2)
levels(rcp2_8$scenario)<-c("long-term historic","predicted future","predicted future")
levels(rcp2_8$scenario2)<-c("RCP 2.6","RCP 8.5")

# plot trends
plot<-ggplot(data=rcp2_8, aes(x=year,y=total_abund)) + 
  geom_point(size=3,aes(fill=mean_monsoon6SPEI*(-1),shape=scenario))+ 
  scale_shape_manual(values=c(22,21)) +
  xlab("Year") + ylab("Total abundance") + 
  theme_bw() + theme(axis.text.x = element_text(size=16))+ theme(axis.text.y = element_text(size=16)) + theme(legend.text=element_text(size=12), legend.title=element_text(size=15))+
  theme(axis.title.x = element_text(size=20))+ theme(axis.title.y = element_text(size=20)) +
  geom_smooth(method="lm",color="black") + 
  scale_fill_gradient2(high = "red3", low = "royalblue",mid = "white",midpoint = 0) + labs(fill="Aridity index", shape="Data type") + 
  #ylim(0,10000) + 
  facet_wrap(~scenario2) + theme(panel.spacing.x = unit(10, "mm")) + theme(strip.text.x = element_text(size = 18)) + theme(axis.line = element_line(color='black'), plot.background = element_blank(),panel.grid.minor = element_blank(),panel.grid.major = element_blank())
plot

# save plot
# ggsave("change_total_abundance_rcp2.6_and_8.5_2023-12-01.jpg", plot,
#        width=12,height=6.25,units = c("in"),
#        dpi = 300)


# Linear regressions corresponding with trends above (how has total abundance changed over time under each climate scenario?)

# RCP 4.5
rcp4.5<-subset(abund_summary, scenario=="historic" | scenario=="rcp4.5")
anova(lm(total_abund~year*ecosystem, data=rcp4.5))
summary(lm(total_abund~year*ecosystem, data=rcp4.5))

# RCP 2.6
rcp2.6<-subset(abund_summary, scenario=="historic" | scenario=="rcp2.6")
anova(lm(total_abund~year*ecosystem, data=rcp2.6))
summary(lm(total_abund~year*ecosystem, data=rcp2.6))

# RCP 8.5
rcp8.5<-subset(abund_summary, scenario=="historic" | scenario=="rcp8.5")
anova(lm(total_abund~year*ecosystem, data=rcp8.5))
summary(lm(total_abund~year*ecosystem, data=rcp8.5))


## How does current rank abundance relate to magnitude of predicted change in abundance over time? #####

# calculate abundance ranks during 2002-2019
historic<-subset(abund_combined, scenario == "historic")
rank_historic <- historic %>% group_by(code) %>% summarise(total_abund=sum(abund))
rank_historic <- rank_historic[order(-rank_historic$total_abund),]
rank_historic$abundance_rank_historic<-1:243

# read in predicted abundance change data
slopes <- read.csv("slopes_predicted_future_change_MeanAcrossGCMs_2023-12-01.csv")

# subset change in abundance data to just include RCP 4.5, and join to historic rank abundance data frame
slopes<-subset(slopes,scenario=="rcp4.5")
rank_historic<-left_join(rank_historic,slopes,by="code")

# add genus_species column
rank_historic$genus_species<-paste(rank_historic$genus, rank_historic$species, sep=" ")

# subset data for graph
for_graph<-subset(rank_historic,abundance_rank_historic<=50) # just plotting abundance ranks of top 50 most abundant species

# add column indicating indicating whether abundance change over time was significant
for_graph$signif<-ifelse(for_graph$p_value<0.05,"significant","nonsignificant")

# get a list of genera represented by the 50 most abundant species
summary<-for_graph%>%group_by(genus)%>%summarise(count=n())

# create y-axis label
y_label<-expression(atop("Magnitude of change in","abundance ("*italic("\u03B2")[time]*")"))

# plot magnitude of change in abundance over time as a function of abundance rank
plot_rank<-ggplot(data=for_graph, aes(x=abundance_rank_historic, y=slope_year, fill=genus)) +
  geom_errorbar(aes(ymax = slope_year + se_year, ymin = slope_year - se_year, color=signif))+
  geom_hline(yintercept=0, linetype="dashed", color = "darkgray") +
  geom_point(size=3,shape=21,aes(color=signif)) +
  
  theme_classic() + theme(axis.text.x = element_text(size=16))+ theme(axis.text.y = element_text(size=16)) + theme(legend.text=element_text(size=12), legend.title=element_text(size=15))+
  theme(axis.title.x = element_text(size=20))+ theme(axis.title.y = element_text(size=20)) +
  
  xlab("Abundance rank") + ylab(y_label) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_color_manual(values=c("gray","black")) +
  scale_fill_manual(values=c("#A6CEE3",
                             "#1F78B4",
                             "#B2DF8A",
                             "#33A02C",
                             "#0000FF",
                             "#FB9A99",
                             "#E31A1C",
                             "#663300",
                             "#003333",
                             "#FF00FF",
                             "#D95F02",
                             "#FDBF6F",
                             "#FF7F00",
                             "#666666",
                             "#CAB2D6",
                             "#D9D9D9",
                             "#6A3D9A",
                             "#BC80BD",
                             "#80B1D3",
                             "#FFFF99")) +
  theme(legend.text = element_text(face="italic"))+
  labs(fill="Genus") +
  guides(color = "none") +
  scale_y_break(breaks = c(-0.035, -0.11), scale=5, space=0.15, ticklabels=c(-0.025,0,0.025,0.05)) +
  scale_y_break(breaks = c(0.07, 0.09), scale=0.6, ticklabels=c(0.10, 0.12), space=0.15)
plot_rank

# combine graphs into a multi-figure plot
p1<- plot_rank + plot_abundchange + plot_layout(ncol = 2) + plot_annotation(tag_levels = 'a') & 
  theme(plot.tag = element_text(size = 30))
p1

# save plot
# ggsave("rank_abund_and_change_over_time_2023-12-01.pdf", p1,
#        width=40,height=18,units = c("cm"),
#        dpi = 600, device = cairo_pdf)


##### CHANGE IN COMMUNITY-WEIGHTED MEAN BODY MASS WITH ARIDITY IN HISTORIC DATASET #####
## Format data #####

# select columns and format data frames in preparation for combining them
abund_historic_final<-abund_historic_2_melt[,c(2,3,10,1,11,12)]

# add station column
abund_historic_final$station<-abund_historic_final$ecosystem

# rename levels of the "station" variable
abund_historic_final$station<-as.factor(abund_historic_final$station)
levels(abund_historic_final$station)[levels(abund_historic_final$station)=="B"] <- "Blue Grama"
levels(abund_historic_final$station)[levels(abund_historic_final$station)=="C"] <- "Five Points"
levels(abund_historic_final$station)[levels(abund_historic_final$station)=="G"] <- "Five Points"

## Add climate data ##

# read in SPEI data
clim_canesm<-read.csv("spei_historic_and_future_byscenario_CanESM2_2023-12-01.csv")

# create data frame for historic data, and add climate data to it
colnames(clim_canesm)
clim_historic<-subset(clim_canesm,source=="SEV-met")
historic2<-left_join(abund_historic_final,clim_historic,by=c("year","station"))
historic2$source<-NULL

# read in mass data for each bee species
mass<-read.csv("SEVBeeBodyMassData_2023-08-29_forpub.csv")

# calculate average mass for each species
avg_mass<-mass %>% group_by(code) %>% summarise(mean_mass_mg=mean(mass_mg),se_mass_mg=sd(mass_mg)/sqrt(n()))

# Calculate community-weighted mean body mass:

# merge historic bee abundance data with mass data
abund_focal_merged <- left_join(historic2, avg_mass, by = c("code"))

# remove rows with NAs (species for which we lack body mass data)
abund_focal_merged<-subset(abund_focal_merged,!is.na(mean_mass_mg))

# calculate total bee abundance for each ecosystem x transect x year combination
summary <- abund_focal_merged %>% group_by(ecosystem, transect, year) %>% summarise(site_year_total_abund = sum(abund))

# join the two data frames from above
cwdata <- left_join(abund_focal_merged, summary, by = c("ecosystem","transect", "year"))

# add a column of mass-weighted abundance 
cwdata$prop_abund <- cwdata$abund/cwdata$site_year_total_abund
cwdata$cw_mass <- cwdata$prop_abund*cwdata$mean_mass_mg

# calculate community-weighted mean (CWM) body mass
cwdata_final <- cwdata %>% group_by(ecosystem, transect, year, monsoon6SPEI) %>% summarise(cwm_mass = sum(cw_mass, na.rm = TRUE))

# create inverse monsoon SPEI column
cwdata_final$monsoon6SPEI_positivized<-cwdata_final$monsoon6SPEI*(-1)


## Relationship between CWM mass and SPEI in historic dataset #####

# Across ecosystems:

# mixed-effects models: how does CWM body mass vary with monsoon aridity?
# linear model
m3<-lmer(cwm_mass~monsoon6SPEI_positivized*ecosystem+(1|year)+(1|transect),data=cwdata_final,na.action=na.omit)
# quadratic model
m4<-lmer(cwm_mass~monsoon6SPEI_positivized*ecosystem+I(monsoon6SPEI_positivized^2)+(1|year)+(1|transect),data=cwdata_final,na.action=na.omit)

# compare the two models based on AICc
AICc(m3,m4) 
# models have very similar AICc values, with linear model slightly superior --> choose linear model as superior for parsimony

# get statistical values for superior model
Anova(m3,type = 3)
rsquared(m3)

# plot the trend using the visreg package
op <- par(mfrow = c(1, 1), 
          pty = "s",mgp=c(3.6,1,0))
plot<-visreg(m3,"monsoon6SPEI_positivized",type="conditional",points.par=list(cex=1.2,col="black",pch=15),
       cex.axis=1.4, line.par=list(col="black"),
       xlab=list("Aridity index", cex=1.8),
       ylab=list("Community-weighted \nmean body mass (mg)", cex=1.8),xlim=c(-2,1.85)) 

