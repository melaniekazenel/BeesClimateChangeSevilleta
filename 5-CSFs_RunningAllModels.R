################################################################################### 
# Climate sensitivity functions for bees in three ecosystems at the Sevilleta National Wildlife Refuge
# Code to build models, extract results, and plot trends for each bee population (species x ecosystem combination)

# Heat and desiccation tolerances predict bee abundance under climate change
# Melanie R. Kazenel, Karen W. Wright, Terry Griswold, Kenneth D. Whitney, and Jennifer A. Rudgers

# Date: 2023-08-29
# Corresponding author's email: melanie.kazenel@gmail.com
################################################################################### 


# Load required packages
library(dplyr)
library(car)
library(MuMIn)
library(ggplot2)
library(nlme)
library(piecewiseSEM)



##### Format data #####

# Read in data frame of square-root transformed bee count data and SPEI data
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

# Create a matrix to hold statistical results
beeCSF_output<-matrix(nrow=1,ncol=546,byrow=TRUE,dimnames=list(c("row1"),c("code","ecosystem","P_shapiro_null","dAICc_m1_monsoon6mo",
"dAICc_m1_monsoon6mo_lag",
"dAICc_m1_spring6mo",
"dAICc_m1_spring6mo_lag",
"dAICc_m2_monsoon6mo",
"dAICc_m2_monsoon6mo_lag",
"dAICc_m2_spring6mo",
"dAICc_m2_spring6mo_lag",
"dAICc_m3_monsoon6mo",
"dAICc_m3_monsoon6mo_lag",
"dAICc_m3_spring6mo",
"dAICc_m3_spring6mo_lag",
"dAICc_null",
"dAICc_m1_spring6mo_AR1",
"dAICc_m1_spring6mo_AR2",
"dAICc_m2_spring6mo_AR1",
"dAICc_m2_spring6mo_AR2",
"dAICc_m3_spring6mo_AR1",
"dAICc_m3_spring6mo_AR2",
"dAICc_m1_monsoon6mo_AR1",
"dAICc_m1_monsoon6mo_AR2",
"dAICc_m2_monsoon6mo_AR1",
"dAICc_m2_monsoon6mo_AR2",
"dAICc_m3_monsoon6mo_AR1",
"dAICc_m3_monsoon6mo_AR2",
"dAICc_m1_spring6mo_lag_AR1",
"dAICc_m1_spring6mo_lag_AR2",
"dAICc_m2_spring6mo_lag_AR1",
"dAICc_m2_spring6mo_lag_AR2",
"dAICc_m3_spring6mo_lag_AR1",
"dAICc_m3_spring6mo_lag_AR2",
"dAICc_m1_monsoon6mo_lag_AR1",
"dAICc_m1_monsoon6mo_lag_AR2",
"dAICc_m2_monsoon6mo_lag_AR1",
"dAICc_m2_monsoon6mo_lag_AR2",
"dAICc_m3_monsoon6mo_lag_AR1",
"dAICc_m3_monsoon6mo_lag_AR2",
"numDF_lin_m1_spring6mo",
"denDF_lin_m1_spring6mo", 
"F_lin_m1_spring6mo",
"P_lin_m1_spring6mo",
"numDF_lin_m2_spring6mo",
"denDF_lin_m2_spring6mo",
"F_lin_m2_spring6mo",
"P_lin_m2_spring6mo",
"numDF_lin_m3_spring6mo",
"denDF_lin_m3_spring6mo",
"F_lin_m3_spring6mo",
"P_lin_m3_spring6mo",
"numDF_lin_m1_monsoon6mo",
"denDF_lin_m1_monsoon6mo",
"F_lin_m1_monsoon6mo",
"P_lin_m1_monsoon6mo",
"numDF_lin_m2_monsoon6mo",
"denDF_lin_m2_monsoon6mo",
"F_lin_m2_monsoon6mo",
"P_lin_m2_monsoon6mo",
"numDF_lin_m3_monsoon6mo",
"denDF_lin_m3_monsoon6mo",
"F_lin_m3_monsoon6mo",
"P_lin_m3_monsoon6mo",
"numDF_lin_m1_spring6mo_lag",
"denDF_lin_m1_spring6mo_lag",
"F_lin_m1_spring6mo_lag",
"P_lin_m1_spring6mo_lag",
"numDF_lin_m2_spring6mo_lag",
"denDF_lin_m2_spring6mo_lag",
"F_lin_m2_spring6mo_lag",
"P_lin_m2_spring6mo_lag",
"numDF_lin_m3_spring6mo_lag",
"denDF_lin_m3_spring6mo_lag",
"F_lin_m3_spring6mo_lag",
"P_lin_m3_spring6mo_lag",
"numDF_lin_m1_monsoon6mo_lag",
"denDF_lin_m1_monsoon6mo_lag",
"F_lin_m1_monsoon6mo_lag",
"P_lin_m1_monsoon6mo_lag",
"numDF_lin_m2_monsoon6mo_lag",
"denDF_lin_m2_monsoon6mo_lag",
"F_lin_m2_monsoon6mo_lag",
"P_lin_m2_monsoon6mo_lag",
"numDF_lin_m3_monsoon6mo_lag",
"denDF_lin_m3_monsoon6mo_lag",
"F_lin_m3_monsoon6mo_lag",
"P_lin_m3_monsoon6mo_lag",
"ParamEst_lin_m1_monsoon6mo",
"ParamEst_lin_m1_monsoon6mo_lag",
"ParamEst_lin_m1_spring6mo",
"ParamEst_lin_m1_spring6mo_lag",
"ParamEst_lin_m2_monsoon6mo",
"ParamEst_lin_m2_monsoon6mo_lag",
"ParamEst_lin_m2_spring6mo",
"ParamEst_lin_m2_spring6mo_lag",
"ParamEst_lin_m3_monsoon6mo",
"ParamEst_lin_m3_monsoon6mo_lag",
"ParamEst_lin_m3_spring6mo",
"ParamEst_lin_m3_spring6mo_lag",
"SE_lin_m1_monsoon6mo",
"SE_lin_m1_monsoon6mo_lag",
"SE_lin_m1_spring6mo",
"SE_lin_m1_spring6mo_lag",
"SE_lin_m2_monsoon6mo",
"SE_lin_m2_monsoon6mo_lag",
"SE_lin_m2_spring6mo",
"SE_lin_m2_spring6mo_lag",
"SE_lin_m3_monsoon6mo",
"SE_lin_m3_monsoon6mo_lag",
"SE_lin_m3_spring6mo",
"SE_lin_m3_spring6mo_lag",
"numDF_lin_m1_spring6mo_AR1",
"denDF_lin_m1_spring6mo_AR1",
"F_lin_m1_spring6mo_AR1",
"P_lin_m1_spring6mo_AR1",
"ParamEst_lin_m1_spring6mo_AR1",
"SE_lin_m1_spring6mo_AR1",
"numDF_lin_m1_spring6mo_AR2",
"denDF_lin_m1_spring6mo_AR2",
"F_lin_m1_spring6mo_AR2",
"P_lin_m1_spring6mo_AR2",
"ParamEst_lin_m1_spring6mo_AR2",
"SE_lin_m1_spring6mo_AR2",
"numDF_lin_m2_spring6mo_AR1",
"denDF_lin_m2_spring6mo_AR1",
"F_lin_m2_spring6mo_AR1",
"P_lin_m2_spring6mo_AR1",
"ParamEst_lin_m2_spring6mo_AR1",
"SE_lin_m2_spring6mo_AR1",
"numDF_lin_m2_spring6mo_AR2",
"denDF_lin_m2_spring6mo_AR2",
"F_lin_m2_spring6mo_AR2",
"P_lin_m2_spring6mo_AR2",
"ParamEst_lin_m2_spring6mo_AR2",
"SE_lin_m2_spring6mo_AR2",
"numDF_lin_m3_spring6mo_AR1",
"denDF_lin_m3_spring6mo_AR1",
"F_lin_m3_spring6mo_AR1",
"P_lin_m3_spring6mo_AR1",
"ParamEst_lin_m3_spring6mo_AR1",
"SE_lin_m3_spring6mo_AR1",
"numDF_lin_m3_spring6mo_AR2",
"denDF_lin_m3_spring6mo_AR2",
"F_lin_m3_spring6mo_AR2",
"P_lin_m3_spring6mo_AR2",
"ParamEst_lin_m3_spring6mo_AR2",
"SE_lin_m3_spring6mo_AR2",
"numDF_lin_m1_monsoon6mo_AR1",
"denDF_lin_m1_monsoon6mo_AR1",
"F_lin_m1_monsoon6mo_AR1",
"P_lin_m1_monsoon6mo_AR1",
"ParamEst_lin_m1_monsoon6mo_AR1",
"SE_lin_m1_monsoon6mo_AR1",
"numDF_lin_m1_monsoon6mo_AR2",
"denDF_lin_m1_monsoon6mo_AR2",
"F_lin_m1_monsoon6mo_AR2",
"P_lin_m1_monsoon6mo_AR2",
"ParamEst_lin_m1_monsoon6mo_AR2",
"SE_lin_m1_monsoon6mo_AR2",
"numDF_lin_m2_monsoon6mo_AR1",
"denDF_lin_m2_monsoon6mo_AR1",
"F_lin_m2_monsoon6mo_AR1",
"P_lin_m2_monsoon6mo_AR1",
"ParamEst_lin_m2_monsoon6mo_AR1",
"SE_lin_m2_monsoon6mo_AR1",
"numDF_lin_m2_monsoon6mo_AR2",
"denDF_lin_m2_monsoon6mo_AR2",
"F_lin_m2_monsoon6mo_AR2",
"P_lin_m2_monsoon6mo_AR2",
"ParamEst_lin_m2_monsoon6mo_AR2",
"SE_lin_m2_monsoon6mo_AR2",
"numDF_lin_m3_monsoon6mo_AR1",
"denDF_lin_m3_monsoon6mo_AR1",
"F_lin_m3_monsoon6mo_AR1",
"P_lin_m3_monsoon6mo_AR1",
"ParamEst_lin_m3_monsoon6mo_AR1",
"SE_lin_m3_monsoon6mo_AR1",
"numDF_lin_m3_monsoon6mo_AR2",
"denDF_lin_m3_monsoon6mo_AR2",
"F_lin_m3_monsoon6mo_AR2",
"P_lin_m3_monsoon6mo_AR2",
"ParamEst_lin_m3_monsoon6mo_AR2",
"SE_lin_m3_monsoon6mo_AR2",
"numDF_lin_m1_spring6mo_lag_AR1",
"denDF_lin_m1_spring6mo_lag_AR1",
"F_lin_m1_spring6mo_lag_AR1",
"P_lin_m1_spring6mo_lag_AR1",
"ParamEst_lin_m1_spring6mo_lag_AR1",
"SE_lin_m1_spring6mo_lag_AR1",
"numDF_lin_m1_spring6mo_lag_AR2",
"denDF_lin_m1_spring6mo_lag_AR2",
"F_lin_m1_spring6mo_lag_AR2",
"P_lin_m1_spring6mo_lag_AR2",
"ParamEst_lin_m1_spring6mo_lag_AR2",
"SE_lin_m1_spring6mo_lag_AR2",
"numDF_lin_m2_spring6mo_lag_AR1",
"denDF_lin_m2_spring6mo_lag_AR1",
"F_lin_m2_spring6mo_lag_AR1",
"P_lin_m2_spring6mo_lag_AR1",
"ParamEst_lin_m2_spring6mo_lag_AR1",
"SE_lin_m2_spring6mo_lag_AR1",
"numDF_lin_m2_spring6mo_lag_AR2",
"denDF_lin_m2_spring6mo_lag_AR2",
"F_lin_m2_spring6mo_lag_AR2",
"P_lin_m2_spring6mo_lag_AR2",
"ParamEst_lin_m2_spring6mo_lag_AR2",
"SE_lin_m2_spring6mo_lag_AR2",
"numDF_lin_m3_spring6mo_lag_AR1",
"denDF_lin_m3_spring6mo_lag_AR1",
"F_lin_m3_spring6mo_lag_AR1",
"P_lin_m3_spring6mo_lag_AR1",
"ParamEst_lin_m3_spring6mo_lag_AR1",
"SE_lin_m3_spring6mo_lag_AR1",
"numDF_lin_m3_spring6mo_lag_AR2",
"denDF_lin_m3_spring6mo_lag_AR2",
"F_lin_m3_spring6mo_lag_AR2",
"P_lin_m3_spring6mo_lag_AR2",
"ParamEst_lin_m3_spring6mo_lag_AR2",
"SE_lin_m3_spring6mo_lag_AR2",
"numDF_lin_m1_monsoon6mo_lag_AR1",
"denDF_lin_m1_monsoon6mo_lag_AR1",
"F_lin_m1_monsoon6mo_lag_AR1",
"P_lin_m1_monsoon6mo_lag_AR1",
"ParamEst_lin_m1_monsoon6mo_lag_AR1",
"SE_lin_m1_monsoon6mo_lag_AR1",
"numDF_lin_m1_monsoon6mo_lag_AR2",
"denDF_lin_m1_monsoon6mo_lag_AR2",
"F_lin_m1_monsoon6mo_lag_AR2",
"P_lin_m1_monsoon6mo_lag_AR2",
"ParamEst_lin_m1_monsoon6mo_lag_AR2",
"SE_lin_m1_monsoon6mo_lag_AR2",
"numDF_lin_m2_monsoon6mo_lag_AR1",
"denDF_lin_m2_monsoon6mo_lag_AR1",
"F_lin_m2_monsoon6mo_lag_AR1",
"P_lin_m2_monsoon6mo_lag_AR1",
"ParamEst_lin_m2_monsoon6mo_lag_AR1",
"SE_lin_m2_monsoon6mo_lag_AR1",
"numDF_lin_m2_monsoon6mo_lag_AR2",
"denDF_lin_m2_monsoon6mo_lag_AR2",
"F_lin_m2_monsoon6mo_lag_AR2",
"P_lin_m2_monsoon6mo_lag_AR2",
"ParamEst_lin_m2_monsoon6mo_lag_AR2",
"SE_lin_m2_monsoon6mo_lag_AR2",
"numDF_lin_m3_monsoon6mo_lag_AR1",
"denDF_lin_m3_monsoon6mo_lag_AR1",
"F_lin_m3_monsoon6mo_lag_AR1",
"P_lin_m3_monsoon6mo_lag_AR1",
"ParamEst_lin_m3_monsoon6mo_lag_AR1",
"SE_lin_m3_monsoon6mo_lag_AR1",
"numDF_lin_m3_monsoon6mo_lag_AR2",
"denDF_lin_m3_monsoon6mo_lag_AR2",
"F_lin_m3_monsoon6mo_lag_AR2",
"P_lin_m3_monsoon6mo_lag_AR2",
"ParamEst_lin_m3_monsoon6mo_lag_AR2",
"SE_lin_m3_monsoon6mo_lag_AR2",
"Rsquared_marginal_m_null",
"Rsquared_marginal_m1_spring6mo",
"Rsquared_marginal_m2_spring6mo",
"Rsquared_marginal_m3_spring6mo",
"Rsquared_marginal_m1_spring6mo_AR1",
"Rsquared_marginal_m1_spring6mo_AR2",
"Rsquared_marginal_m2_spring6mo_AR1",
"Rsquared_marginal_m2_spring6mo_AR2",
"Rsquared_marginal_m3_spring6mo_AR1",
"Rsquared_marginal_m3_spring6mo_AR2",
"Rsquared_marginal_m1_monsoon6mo",
"Rsquared_marginal_m2_monsoon6mo",
"Rsquared_marginal_m3_monsoon6mo",
"Rsquared_marginal_m1_monsoon6mo_AR1",
"Rsquared_marginal_m1_monsoon6mo_AR2",
"Rsquared_marginal_m2_monsoon6mo_AR1",
"Rsquared_marginal_m2_monsoon6mo_AR2",
"Rsquared_marginal_m3_monsoon6mo_AR1",
"Rsquared_marginal_m3_monsoon6mo_AR2",
"Rsquared_marginal_m1_spring6mo_lag",
"Rsquared_marginal_m2_spring6mo_lag",
"Rsquared_marginal_m3_spring6mo_lag",
"Rsquared_marginal_m1_spring6mo_lag_AR1",
"Rsquared_marginal_m1_spring6mo_lag_AR2",
"Rsquared_marginal_m2_spring6mo_lag_AR1",
"Rsquared_marginal_m2_spring6mo_lag_AR2",
"Rsquared_marginal_m3_spring6mo_lag_AR1",
"Rsquared_marginal_m3_spring6mo_lag_AR2",
"Rsquared_marginal_m1_monsoon6mo_lag",
"Rsquared_marginal_m2_monsoon6mo_lag",
"Rsquared_marginal_m3_monsoon6mo_lag",
"Rsquared_marginal_m1_monsoon6mo_lag_AR1",
"Rsquared_marginal_m1_monsoon6mo_lag_AR2",
"Rsquared_marginal_m2_monsoon6mo_lag_AR1",
"Rsquared_marginal_m2_monsoon6mo_lag_AR2",
"Rsquared_marginal_m3_monsoon6mo_lag_AR1",
"Rsquared_marginal_m3_monsoon6mo_lag_AR2",
"Rsquared_conditional_m_null",
"Rsquared_conditional_m1_spring6mo",
"Rsquared_conditional_m2_spring6mo",
"Rsquared_conditional_m3_spring6mo",
"Rsquared_conditional_m1_spring6mo_AR1",
"Rsquared_conditional_m1_spring6mo_AR2",
"Rsquared_conditional_m2_spring6mo_AR1",
"Rsquared_conditional_m2_spring6mo_AR2",
"Rsquared_conditional_m3_spring6mo_AR1",
"Rsquared_conditional_m3_spring6mo_AR2",
"Rsquared_conditional_m1_monsoon6mo",
"Rsquared_conditional_m2_monsoon6mo",
"Rsquared_conditional_m3_monsoon6mo",
"Rsquared_conditional_m1_monsoon6mo_AR1",
"Rsquared_conditional_m1_monsoon6mo_AR2",
"Rsquared_conditional_m2_monsoon6mo_AR1",
"Rsquared_conditional_m2_monsoon6mo_AR2",
"Rsquared_conditional_m3_monsoon6mo_AR1",
"Rsquared_conditional_m3_monsoon6mo_AR2",
"Rsquared_conditional_m1_spring6mo_lag",
"Rsquared_conditional_m2_spring6mo_lag",
"Rsquared_conditional_m3_spring6mo_lag",
"Rsquared_conditional_m1_spring6mo_lag_AR1",
"Rsquared_conditional_m1_spring6mo_lag_AR2",
"Rsquared_conditional_m2_spring6mo_lag_AR1",
"Rsquared_conditional_m2_spring6mo_lag_AR2",
"Rsquared_conditional_m3_spring6mo_lag_AR1",
"Rsquared_conditional_m3_spring6mo_lag_AR2",
"Rsquared_conditional_m1_monsoon6mo_lag",
"Rsquared_conditional_m2_monsoon6mo_lag",
"Rsquared_conditional_m3_monsoon6mo_lag",
"Rsquared_conditional_m1_monsoon6mo_lag_AR1",
"Rsquared_conditional_m1_monsoon6mo_lag_AR2",
"Rsquared_conditional_m2_monsoon6mo_lag_AR1",
"Rsquared_conditional_m2_monsoon6mo_lag_AR2",
"Rsquared_conditional_m3_monsoon6mo_lag_AR1",
"Rsquared_conditional_m3_monsoon6mo_lag_AR2",
"numDF_quad_m2_spring6mo",
"numDF_quad_m2_spring6mo_AR1",
"numDF_quad_m2_spring6mo_AR2",
"numDF_quad_m2_monsoon6mo",
"numDF_quad_m2_monsoon6mo_AR1",
"numDF_quad_m2_monsoon6mo_AR2",
"numDF_quad_m2_spring6mo_lag",
"numDF_quad_m2_spring6mo_lag_AR1",
"numDF_quad_m2_spring6mo_lag_AR2",
"numDF_quad_m2_monsoon6mo_lag",
"numDF_quad_m2_monsoon6mo_lag_AR1",
"numDF_quad_m2_monsoon6mo_lag_AR2",
"numDF_quad_m3_spring6mo",
"numDF_quad_m3_spring6mo_AR1",
"numDF_quad_m3_spring6mo_AR2",
"numDF_quad_m3_monsoon6mo",
"numDF_quad_m3_monsoon6mo_AR1",
"numDF_quad_m3_monsoon6mo_AR2",
"numDF_quad_m3_spring6mo_lag",
"numDF_quad_m3_spring6mo_lag_AR1",
"numDF_quad_m3_spring6mo_lag_AR2",
"numDF_quad_m3_monsoon6mo_lag",
"numDF_quad_m3_monsoon6mo_lag_AR1",
"numDF_quad_m3_monsoon6mo_lag_AR2",
"denDF_quad_m2_spring6mo",
"denDF_quad_m2_spring6mo_AR1",
"denDF_quad_m2_spring6mo_AR2",
"denDF_quad_m2_monsoon6mo",
"denDF_quad_m2_monsoon6mo_AR1",
"denDF_quad_m2_monsoon6mo_AR2",
"denDF_quad_m2_spring6mo_lag",
"denDF_quad_m2_spring6mo_lag_AR1",
"denDF_quad_m2_spring6mo_lag_AR2",
"denDF_quad_m2_monsoon6mo_lag",
"denDF_quad_m2_monsoon6mo_lag_AR1",
"denDF_quad_m2_monsoon6mo_lag_AR2",
"denDF_quad_m3_spring6mo",
"denDF_quad_m3_spring6mo_AR1",
"denDF_quad_m3_spring6mo_AR2",
"denDF_quad_m3_monsoon6mo",
"denDF_quad_m3_monsoon6mo_AR1",
"denDF_quad_m3_monsoon6mo_AR2",
"denDF_quad_m3_spring6mo_lag",
"denDF_quad_m3_spring6mo_lag_AR1",
"denDF_quad_m3_spring6mo_lag_AR2",
"denDF_quad_m3_monsoon6mo_lag",
"denDF_quad_m3_monsoon6mo_lag_AR1",
"denDF_quad_m3_monsoon6mo_lag_AR2",
"F_quad_m2_spring6mo",
"F_quad_m2_spring6mo_AR1",
"F_quad_m2_spring6mo_AR2",
"F_quad_m2_monsoon6mo",
"F_quad_m2_monsoon6mo_AR1",
"F_quad_m2_monsoon6mo_AR2",
"F_quad_m2_spring6mo_lag",
"F_quad_m2_spring6mo_lag_AR1",
"F_quad_m2_spring6mo_lag_AR2",
"F_quad_m2_monsoon6mo_lag",
"F_quad_m2_monsoon6mo_lag_AR1",
"F_quad_m2_monsoon6mo_lag_AR2",
"F_quad_m3_spring6mo",
"F_quad_m3_spring6mo_AR1",
"F_quad_m3_spring6mo_AR2",
"F_quad_m3_monsoon6mo",
"F_quad_m3_monsoon6mo_AR1",
"F_quad_m3_monsoon6mo_AR2",
"F_quad_m3_spring6mo_lag",
"F_quad_m3_spring6mo_lag_AR1",
"F_quad_m3_spring6mo_lag_AR2",
"F_quad_m3_monsoon6mo_lag",
"F_quad_m3_monsoon6mo_lag_AR1",
"F_quad_m3_monsoon6mo_lag_AR2",
"P_quad_m2_spring6mo",
"P_quad_m2_spring6mo_AR1",
"P_quad_m2_spring6mo_AR2",
"P_quad_m2_monsoon6mo",
"P_quad_m2_monsoon6mo_AR1",
"P_quad_m2_monsoon6mo_AR2",
"P_quad_m2_spring6mo_lag",
"P_quad_m2_spring6mo_lag_AR1",
"P_quad_m2_spring6mo_lag_AR2",
"P_quad_m2_monsoon6mo_lag",
"P_quad_m2_monsoon6mo_lag_AR1",
"P_quad_m2_monsoon6mo_lag_AR2",
"P_quad_m3_spring6mo",
"P_quad_m3_spring6mo_AR1",
"P_quad_m3_spring6mo_AR2",
"P_quad_m3_monsoon6mo",
"P_quad_m3_monsoon6mo_AR1",
"P_quad_m3_monsoon6mo_AR2",
"P_quad_m3_spring6mo_lag",
"P_quad_m3_spring6mo_lag_AR1",
"P_quad_m3_spring6mo_lag_AR2",
"P_quad_m3_monsoon6mo_lag",
"P_quad_m3_monsoon6mo_lag_AR1",
"P_quad_m3_monsoon6mo_lag_AR2",
"ParamEst_quad_m2_spring6mo",
"ParamEst_quad_m2_spring6mo_AR1",
"ParamEst_quad_m2_spring6mo_AR2",
"ParamEst_quad_m2_monsoon6mo",
"ParamEst_quad_m2_monsoon6mo_AR1",
"ParamEst_quad_m2_monsoon6mo_AR2",
"ParamEst_quad_m2_spring6mo_lag",
"ParamEst_quad_m2_spring6mo_lag_AR1",
"ParamEst_quad_m2_spring6mo_lag_AR2",
"ParamEst_quad_m2_monsoon6mo_lag",
"ParamEst_quad_m2_monsoon6mo_lag_AR1",
"ParamEst_quad_m2_monsoon6mo_lag_AR2",
"ParamEst_quad_m3_spring6mo",
"ParamEst_quad_m3_spring6mo_AR1",
"ParamEst_quad_m3_spring6mo_AR2",
"ParamEst_quad_m3_monsoon6mo",
"ParamEst_quad_m3_monsoon6mo_AR1",
"ParamEst_quad_m3_monsoon6mo_AR2",
"ParamEst_quad_m3_spring6mo_lag",
"ParamEst_quad_m3_spring6mo_lag_AR1",
"ParamEst_quad_m3_spring6mo_lag_AR2",
"ParamEst_quad_m3_monsoon6mo_lag",
"ParamEst_quad_m3_monsoon6mo_lag_AR1",
"ParamEst_quad_m3_monsoon6mo_lag_AR2",
"SE_quad_m2_spring6mo",
"SE_quad_m2_spring6mo_AR1",
"SE_quad_m2_spring6mo_AR2",
"SE_quad_m2_monsoon6mo",
"SE_quad_m2_monsoon6mo_AR1",
"SE_quad_m2_monsoon6mo_AR2",
"SE_quad_m2_spring6mo_lag",
"SE_quad_m2_spring6mo_lag_AR1",
"SE_quad_m2_spring6mo_lag_AR2",
"SE_quad_m2_monsoon6mo_lag",
"SE_quad_m2_monsoon6mo_lag_AR1",
"SE_quad_m2_monsoon6mo_lag_AR2",
"SE_quad_m3_spring6mo",
"SE_quad_m3_spring6mo_AR1",
"SE_quad_m3_spring6mo_AR2",
"SE_quad_m3_monsoon6mo",
"SE_quad_m3_monsoon6mo_AR1",
"SE_quad_m3_monsoon6mo_AR2",
"SE_quad_m3_spring6mo_lag",
"SE_quad_m3_spring6mo_lag_AR1",
"SE_quad_m3_spring6mo_lag_AR2",
"SE_quad_m3_monsoon6mo_lag",
"SE_quad_m3_monsoon6mo_lag_AR1",
"SE_quad_m3_monsoon6mo_lag_AR2",
"numDF_cub_m3_spring6mo",
"numDF_cub_m3_spring6mo_AR1",
"numDF_cub_m3_spring6mo_AR2",
"numDF_cub_m3_monsoon6mo",
"numDF_cub_m3_monsoon6mo_AR1",
"numDF_cub_m3_monsoon6mo_AR2",
"numDF_cub_m3_spring6mo_lag",
"numDF_cub_m3_spring6mo_lag_AR1",
"numDF_cub_m3_spring6mo_lag_AR2",
"numDF_cub_m3_monsoon6mo_lag",
"numDF_cub_m3_monsoon6mo_lag_AR1",
"numDF_cub_m3_monsoon6mo_lag_AR2",
"denDF_cub_m3_spring6mo",
"denDF_cub_m3_spring6mo_AR1",
"denDF_cub_m3_spring6mo_AR2",
"denDF_cub_m3_monsoon6mo",
"denDF_cub_m3_monsoon6mo_AR1",
"denDF_cub_m3_monsoon6mo_AR2",
"denDF_cub_m3_spring6mo_lag",
"denDF_cub_m3_spring6mo_lag_AR1",
"denDF_cub_m3_spring6mo_lag_AR2",
"denDF_cub_m3_monsoon6mo_lag",
"denDF_cub_m3_monsoon6mo_lag_AR1",
"denDF_cub_m3_monsoon6mo_lag_AR2",
"F_cub_m3_spring6mo",
"F_cub_m3_spring6mo_AR1",
"F_cub_m3_spring6mo_AR2",
"F_cub_m3_monsoon6mo",
"F_cub_m3_monsoon6mo_AR1",
"F_cub_m3_monsoon6mo_AR2",
"F_cub_m3_spring6mo_lag",
"F_cub_m3_spring6mo_lag_AR1",
"F_cub_m3_spring6mo_lag_AR2",
"F_cub_m3_monsoon6mo_lag",
"F_cub_m3_monsoon6mo_lag_AR1",
"F_cub_m3_monsoon6mo_lag_AR2",
"P_cub_m3_spring6mo",
"P_cub_m3_spring6mo_AR1",
"P_cub_m3_spring6mo_AR2",
"P_cub_m3_monsoon6mo",
"P_cub_m3_monsoon6mo_AR1",
"P_cub_m3_monsoon6mo_AR2",
"P_cub_m3_spring6mo_lag",
"P_cub_m3_spring6mo_lag_AR1",
"P_cub_m3_spring6mo_lag_AR2",
"P_cub_m3_monsoon6mo_lag",
"P_cub_m3_monsoon6mo_lag_AR1",
"P_cub_m3_monsoon6mo_lag_AR2",
"ParamEst_cub_m3_spring6mo",
"ParamEst_cub_m3_spring6mo_AR1",
"ParamEst_cub_m3_spring6mo_AR2",
"ParamEst_cub_m3_monsoon6mo",
"ParamEst_cub_m3_monsoon6mo_AR1",
"ParamEst_cub_m3_monsoon6mo_AR2",
"ParamEst_cub_m3_spring6mo_lag",
"ParamEst_cub_m3_spring6mo_lag_AR1",
"ParamEst_cub_m3_spring6mo_lag_AR2",
"ParamEst_cub_m3_monsoon6mo_lag",
"ParamEst_cub_m3_monsoon6mo_lag_AR1",
"ParamEst_cub_m3_monsoon6mo_lag_AR2",
"SE_cub_m3_spring6mo",
"SE_cub_m3_spring6mo_AR1",
"SE_cub_m3_spring6mo_AR2",
"SE_cub_m3_monsoon6mo",
"SE_cub_m3_monsoon6mo_AR1",
"SE_cub_m3_monsoon6mo_AR2",
"SE_cub_m3_spring6mo_lag",
"SE_cub_m3_spring6mo_lag_AR1",
"SE_cub_m3_spring6mo_lag_AR2",
"SE_cub_m3_monsoon6mo_lag",
"SE_cub_m3_monsoon6mo_lag_AR1",
"SE_cub_m3_monsoon6mo_lag_AR2")))

##### CSFs: Chihuahuan Desert Shrubland #####

# Create a new data frame of the original data
creo_original <- creo

# Create a data frame of just climate data
creo_climate <- creo[,5:8] 

# Create a data frame of just the bee abundance matrix (descriptor variables removed)
speciesMatrix <- creo[,9:229]

# Create a vector of species codes
speciesCodes <- colnames(speciesMatrix)

# Create a vector with a number corresponding to each species
number <-1:length(speciesCodes)

### Loop through each column of speciesMatrix, running CSFs and putting the model output in the beeCSF_output matrix ###

for (i in 1:length(speciesMatrix[1,])) {
  
  # save the species code for column i
  speciesCode <- speciesCodes[i]
  
  # create an object with the name of the ecosystem type
  ecosystem <-"C"
  
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
  
  
  # normality test for null model
  shapiro_null<-shapiro.test(resid(m_null))
  P_shapiro_null<-shapiro_null$p.value
  print(P_shapiro_null)
  
  
  # calculate delta AICc for each model
  
  # delta AICc for AR models
  dAICc_m1_spring6mo_AR1<-(AICc(m1_spring6mo_AR1))-min(AICc(m_null),AICc(m1_spring6mo),AICc(m2_spring6mo),AICc(m3_spring6mo), AICc(m1_monsoon6mo), AICc(m2_monsoon6mo),AICc(m3_monsoon6mo),AICc(m1_spring6mo_lag),AICc(m2_spring6mo_lag),AICc(m3_spring6mo_lag),AICc(m1_monsoon6mo_lag),AICc(m2_monsoon6mo_lag),AICc(m3_monsoon6mo_lag),AICc(m1_spring6mo_AR1),AICc(m1_spring6mo_AR2),AICc(m2_spring6mo_AR1),AICc(m2_spring6mo_AR2),AICc(m3_spring6mo_AR1),AICc(m3_spring6mo_AR2),AICc(m1_monsoon6mo_AR1),AICc(m1_monsoon6mo_AR2),AICc(m2_monsoon6mo_AR1),AICc(m2_monsoon6mo_AR2),AICc(m3_monsoon6mo_AR1),AICc(m3_monsoon6mo_AR2),AICc(m1_spring6mo_lag_AR1),AICc(m1_spring6mo_lag_AR2),AICc(m2_spring6mo_lag_AR1),AICc(m2_spring6mo_lag_AR2),AICc(m3_spring6mo_lag_AR1),AICc(m3_spring6mo_lag_AR2),AICc(m1_monsoon6mo_lag_AR1),AICc(m1_monsoon6mo_lag_AR2),AICc(m2_monsoon6mo_lag_AR1),AICc(m2_monsoon6mo_lag_AR2),AICc(m3_monsoon6mo_lag_AR1),AICc(m3_monsoon6mo_lag_AR2))
  
  dAICc_m1_spring6mo_AR2<-(AICc(m1_spring6mo_AR2))-min(AICc(m_null),AICc(m1_spring6mo),AICc(m2_spring6mo),AICc(m3_spring6mo), AICc(m1_monsoon6mo), AICc(m2_monsoon6mo),AICc(m3_monsoon6mo),AICc(m1_spring6mo_lag),AICc(m2_spring6mo_lag),AICc(m3_spring6mo_lag),AICc(m1_monsoon6mo_lag),AICc(m2_monsoon6mo_lag),AICc(m3_monsoon6mo_lag),AICc(m1_spring6mo_AR1),AICc(m1_spring6mo_AR2),AICc(m2_spring6mo_AR1),AICc(m2_spring6mo_AR2),AICc(m3_spring6mo_AR1),AICc(m3_spring6mo_AR2),AICc(m1_monsoon6mo_AR1),AICc(m1_monsoon6mo_AR2),AICc(m2_monsoon6mo_AR1),AICc(m2_monsoon6mo_AR2),AICc(m3_monsoon6mo_AR1),AICc(m3_monsoon6mo_AR2),AICc(m1_spring6mo_lag_AR1),AICc(m1_spring6mo_lag_AR2),AICc(m2_spring6mo_lag_AR1),AICc(m2_spring6mo_lag_AR2),AICc(m3_spring6mo_lag_AR1),AICc(m3_spring6mo_lag_AR2),AICc(m1_monsoon6mo_lag_AR1),AICc(m1_monsoon6mo_lag_AR2),AICc(m2_monsoon6mo_lag_AR1),AICc(m2_monsoon6mo_lag_AR2),AICc(m3_monsoon6mo_lag_AR1),AICc(m3_monsoon6mo_lag_AR2))
  
  dAICc_m2_spring6mo_AR1<-(AICc(m2_spring6mo_AR1))-min(AICc(m_null),AICc(m1_spring6mo),AICc(m2_spring6mo),AICc(m3_spring6mo), AICc(m1_monsoon6mo), AICc(m2_monsoon6mo),AICc(m3_monsoon6mo),AICc(m1_spring6mo_lag),AICc(m2_spring6mo_lag),AICc(m3_spring6mo_lag),AICc(m1_monsoon6mo_lag),AICc(m2_monsoon6mo_lag),AICc(m3_monsoon6mo_lag),AICc(m1_spring6mo_AR1),AICc(m1_spring6mo_AR2),AICc(m2_spring6mo_AR1),AICc(m2_spring6mo_AR2),AICc(m3_spring6mo_AR1),AICc(m3_spring6mo_AR2),AICc(m1_monsoon6mo_AR1),AICc(m1_monsoon6mo_AR2),AICc(m2_monsoon6mo_AR1),AICc(m2_monsoon6mo_AR2),AICc(m3_monsoon6mo_AR1),AICc(m3_monsoon6mo_AR2),AICc(m1_spring6mo_lag_AR1),AICc(m1_spring6mo_lag_AR2),AICc(m2_spring6mo_lag_AR1),AICc(m2_spring6mo_lag_AR2),AICc(m3_spring6mo_lag_AR1),AICc(m3_spring6mo_lag_AR2),AICc(m1_monsoon6mo_lag_AR1),AICc(m1_monsoon6mo_lag_AR2),AICc(m2_monsoon6mo_lag_AR1),AICc(m2_monsoon6mo_lag_AR2),AICc(m3_monsoon6mo_lag_AR1),AICc(m3_monsoon6mo_lag_AR2))
  
  dAICc_m2_spring6mo_AR2<-(AICc(m2_spring6mo_AR2))-min(AICc(m_null),AICc(m1_spring6mo),AICc(m2_spring6mo),AICc(m3_spring6mo), AICc(m1_monsoon6mo), AICc(m2_monsoon6mo),AICc(m3_monsoon6mo),AICc(m1_spring6mo_lag),AICc(m2_spring6mo_lag),AICc(m3_spring6mo_lag),AICc(m1_monsoon6mo_lag),AICc(m2_monsoon6mo_lag),AICc(m3_monsoon6mo_lag),AICc(m1_spring6mo_AR1),AICc(m1_spring6mo_AR2),AICc(m2_spring6mo_AR1),AICc(m2_spring6mo_AR2),AICc(m3_spring6mo_AR1),AICc(m3_spring6mo_AR2),AICc(m1_monsoon6mo_AR1),AICc(m1_monsoon6mo_AR2),AICc(m2_monsoon6mo_AR1),AICc(m2_monsoon6mo_AR2),AICc(m3_monsoon6mo_AR1),AICc(m3_monsoon6mo_AR2),AICc(m1_spring6mo_lag_AR1),AICc(m1_spring6mo_lag_AR2),AICc(m2_spring6mo_lag_AR1),AICc(m2_spring6mo_lag_AR2),AICc(m3_spring6mo_lag_AR1),AICc(m3_spring6mo_lag_AR2),AICc(m1_monsoon6mo_lag_AR1),AICc(m1_monsoon6mo_lag_AR2),AICc(m2_monsoon6mo_lag_AR1),AICc(m2_monsoon6mo_lag_AR2),AICc(m3_monsoon6mo_lag_AR1),AICc(m3_monsoon6mo_lag_AR2))
  
  dAICc_m3_spring6mo_AR1<-(AICc(m3_spring6mo_AR1))-min(AICc(m_null),AICc(m1_spring6mo),AICc(m2_spring6mo),AICc(m3_spring6mo), AICc(m1_monsoon6mo), AICc(m2_monsoon6mo),AICc(m3_monsoon6mo),AICc(m1_spring6mo_lag),AICc(m2_spring6mo_lag),AICc(m3_spring6mo_lag),AICc(m1_monsoon6mo_lag),AICc(m2_monsoon6mo_lag),AICc(m3_monsoon6mo_lag),AICc(m1_spring6mo_AR1),AICc(m1_spring6mo_AR2),AICc(m2_spring6mo_AR1),AICc(m2_spring6mo_AR2),AICc(m3_spring6mo_AR1),AICc(m3_spring6mo_AR2),AICc(m1_monsoon6mo_AR1),AICc(m1_monsoon6mo_AR2),AICc(m2_monsoon6mo_AR1),AICc(m2_monsoon6mo_AR2),AICc(m3_monsoon6mo_AR1),AICc(m3_monsoon6mo_AR2),AICc(m1_spring6mo_lag_AR1),AICc(m1_spring6mo_lag_AR2),AICc(m2_spring6mo_lag_AR1),AICc(m2_spring6mo_lag_AR2),AICc(m3_spring6mo_lag_AR1),AICc(m3_spring6mo_lag_AR2),AICc(m1_monsoon6mo_lag_AR1),AICc(m1_monsoon6mo_lag_AR2),AICc(m2_monsoon6mo_lag_AR1),AICc(m2_monsoon6mo_lag_AR2),AICc(m3_monsoon6mo_lag_AR1),AICc(m3_monsoon6mo_lag_AR2))
  
  dAICc_m3_spring6mo_AR2<-(AICc(m3_spring6mo_AR2))-min(AICc(m_null),AICc(m1_spring6mo),AICc(m2_spring6mo),AICc(m3_spring6mo), AICc(m1_monsoon6mo), AICc(m2_monsoon6mo),AICc(m3_monsoon6mo),AICc(m1_spring6mo_lag),AICc(m2_spring6mo_lag),AICc(m3_spring6mo_lag),AICc(m1_monsoon6mo_lag),AICc(m2_monsoon6mo_lag),AICc(m3_monsoon6mo_lag),AICc(m1_spring6mo_AR1),AICc(m1_spring6mo_AR2),AICc(m2_spring6mo_AR1),AICc(m2_spring6mo_AR2),AICc(m3_spring6mo_AR1),AICc(m3_spring6mo_AR2),AICc(m1_monsoon6mo_AR1),AICc(m1_monsoon6mo_AR2),AICc(m2_monsoon6mo_AR1),AICc(m2_monsoon6mo_AR2),AICc(m3_monsoon6mo_AR1),AICc(m3_monsoon6mo_AR2),AICc(m1_spring6mo_lag_AR1),AICc(m1_spring6mo_lag_AR2),AICc(m2_spring6mo_lag_AR1),AICc(m2_spring6mo_lag_AR2),AICc(m3_spring6mo_lag_AR1),AICc(m3_spring6mo_lag_AR2),AICc(m1_monsoon6mo_lag_AR1),AICc(m1_monsoon6mo_lag_AR2),AICc(m2_monsoon6mo_lag_AR1),AICc(m2_monsoon6mo_lag_AR2),AICc(m3_monsoon6mo_lag_AR1),AICc(m3_monsoon6mo_lag_AR2))
  
  dAICc_m1_monsoon6mo_AR1<-(AICc(m1_monsoon6mo_AR1))-min(AICc(m_null),AICc(m1_spring6mo),AICc(m2_spring6mo),AICc(m3_spring6mo), AICc(m1_monsoon6mo), AICc(m2_monsoon6mo),AICc(m3_monsoon6mo),AICc(m1_spring6mo_lag),AICc(m2_spring6mo_lag),AICc(m3_spring6mo_lag),AICc(m1_monsoon6mo_lag),AICc(m2_monsoon6mo_lag),AICc(m3_monsoon6mo_lag),AICc(m1_spring6mo_AR1),AICc(m1_spring6mo_AR2),AICc(m2_spring6mo_AR1),AICc(m2_spring6mo_AR2),AICc(m3_spring6mo_AR1),AICc(m3_spring6mo_AR2),AICc(m1_monsoon6mo_AR1),AICc(m1_monsoon6mo_AR2),AICc(m2_monsoon6mo_AR1),AICc(m2_monsoon6mo_AR2),AICc(m3_monsoon6mo_AR1),AICc(m3_monsoon6mo_AR2),AICc(m1_spring6mo_lag_AR1),AICc(m1_spring6mo_lag_AR2),AICc(m2_spring6mo_lag_AR1),AICc(m2_spring6mo_lag_AR2),AICc(m3_spring6mo_lag_AR1),AICc(m3_spring6mo_lag_AR2),AICc(m1_monsoon6mo_lag_AR1),AICc(m1_monsoon6mo_lag_AR2),AICc(m2_monsoon6mo_lag_AR1),AICc(m2_monsoon6mo_lag_AR2),AICc(m3_monsoon6mo_lag_AR1),AICc(m3_monsoon6mo_lag_AR2))
  
  dAICc_m1_monsoon6mo_AR2<-(AICc(m1_monsoon6mo_AR2))-min(AICc(m_null),AICc(m1_spring6mo),AICc(m2_spring6mo),AICc(m3_spring6mo), AICc(m1_monsoon6mo), AICc(m2_monsoon6mo),AICc(m3_monsoon6mo),AICc(m1_spring6mo_lag),AICc(m2_spring6mo_lag),AICc(m3_spring6mo_lag),AICc(m1_monsoon6mo_lag),AICc(m2_monsoon6mo_lag),AICc(m3_monsoon6mo_lag),AICc(m1_spring6mo_AR1),AICc(m1_spring6mo_AR2),AICc(m2_spring6mo_AR1),AICc(m2_spring6mo_AR2),AICc(m3_spring6mo_AR1),AICc(m3_spring6mo_AR2),AICc(m1_monsoon6mo_AR1),AICc(m1_monsoon6mo_AR2),AICc(m2_monsoon6mo_AR1),AICc(m2_monsoon6mo_AR2),AICc(m3_monsoon6mo_AR1),AICc(m3_monsoon6mo_AR2),AICc(m1_spring6mo_lag_AR1),AICc(m1_spring6mo_lag_AR2),AICc(m2_spring6mo_lag_AR1),AICc(m2_spring6mo_lag_AR2),AICc(m3_spring6mo_lag_AR1),AICc(m3_spring6mo_lag_AR2),AICc(m1_monsoon6mo_lag_AR1),AICc(m1_monsoon6mo_lag_AR2),AICc(m2_monsoon6mo_lag_AR1),AICc(m2_monsoon6mo_lag_AR2),AICc(m3_monsoon6mo_lag_AR1),AICc(m3_monsoon6mo_lag_AR2))
  
  dAICc_m2_monsoon6mo_AR1<-(AICc(m2_monsoon6mo_AR1))-min(AICc(m_null),AICc(m1_spring6mo),AICc(m2_spring6mo),AICc(m3_spring6mo), AICc(m1_monsoon6mo), AICc(m2_monsoon6mo),AICc(m3_monsoon6mo),AICc(m1_spring6mo_lag),AICc(m2_spring6mo_lag),AICc(m3_spring6mo_lag),AICc(m1_monsoon6mo_lag),AICc(m2_monsoon6mo_lag),AICc(m3_monsoon6mo_lag),AICc(m1_spring6mo_AR1),AICc(m1_spring6mo_AR2),AICc(m2_spring6mo_AR1),AICc(m2_spring6mo_AR2),AICc(m3_spring6mo_AR1),AICc(m3_spring6mo_AR2),AICc(m1_monsoon6mo_AR1),AICc(m1_monsoon6mo_AR2),AICc(m2_monsoon6mo_AR1),AICc(m2_monsoon6mo_AR2),AICc(m3_monsoon6mo_AR1),AICc(m3_monsoon6mo_AR2),AICc(m1_spring6mo_lag_AR1),AICc(m1_spring6mo_lag_AR2),AICc(m2_spring6mo_lag_AR1),AICc(m2_spring6mo_lag_AR2),AICc(m3_spring6mo_lag_AR1),AICc(m3_spring6mo_lag_AR2),AICc(m1_monsoon6mo_lag_AR1),AICc(m1_monsoon6mo_lag_AR2),AICc(m2_monsoon6mo_lag_AR1),AICc(m2_monsoon6mo_lag_AR2),AICc(m3_monsoon6mo_lag_AR1),AICc(m3_monsoon6mo_lag_AR2))
  
  dAICc_m2_monsoon6mo_AR2<-(AICc(m2_monsoon6mo_AR2))-min(AICc(m_null),AICc(m1_spring6mo),AICc(m2_spring6mo),AICc(m3_spring6mo), AICc(m1_monsoon6mo), AICc(m2_monsoon6mo),AICc(m3_monsoon6mo),AICc(m1_spring6mo_lag),AICc(m2_spring6mo_lag),AICc(m3_spring6mo_lag),AICc(m1_monsoon6mo_lag),AICc(m2_monsoon6mo_lag),AICc(m3_monsoon6mo_lag),AICc(m1_spring6mo_AR1),AICc(m1_spring6mo_AR2),AICc(m2_spring6mo_AR1),AICc(m2_spring6mo_AR2),AICc(m3_spring6mo_AR1),AICc(m3_spring6mo_AR2),AICc(m1_monsoon6mo_AR1),AICc(m1_monsoon6mo_AR2),AICc(m2_monsoon6mo_AR1),AICc(m2_monsoon6mo_AR2),AICc(m3_monsoon6mo_AR1),AICc(m3_monsoon6mo_AR2),AICc(m1_spring6mo_lag_AR1),AICc(m1_spring6mo_lag_AR2),AICc(m2_spring6mo_lag_AR1),AICc(m2_spring6mo_lag_AR2),AICc(m3_spring6mo_lag_AR1),AICc(m3_spring6mo_lag_AR2),AICc(m1_monsoon6mo_lag_AR1),AICc(m1_monsoon6mo_lag_AR2),AICc(m2_monsoon6mo_lag_AR1),AICc(m2_monsoon6mo_lag_AR2),AICc(m3_monsoon6mo_lag_AR1),AICc(m3_monsoon6mo_lag_AR2))
  
  dAICc_m3_monsoon6mo_AR1<-(AICc(m3_monsoon6mo_AR1))-min(AICc(m_null),AICc(m1_spring6mo),AICc(m2_spring6mo),AICc(m3_spring6mo), AICc(m1_monsoon6mo), AICc(m2_monsoon6mo),AICc(m3_monsoon6mo),AICc(m1_spring6mo_lag),AICc(m2_spring6mo_lag),AICc(m3_spring6mo_lag),AICc(m1_monsoon6mo_lag),AICc(m2_monsoon6mo_lag),AICc(m3_monsoon6mo_lag),AICc(m1_spring6mo_AR1),AICc(m1_spring6mo_AR2),AICc(m2_spring6mo_AR1),AICc(m2_spring6mo_AR2),AICc(m3_spring6mo_AR1),AICc(m3_spring6mo_AR2),AICc(m1_monsoon6mo_AR1),AICc(m1_monsoon6mo_AR2),AICc(m2_monsoon6mo_AR1),AICc(m2_monsoon6mo_AR2),AICc(m3_monsoon6mo_AR1),AICc(m3_monsoon6mo_AR2),AICc(m1_spring6mo_lag_AR1),AICc(m1_spring6mo_lag_AR2),AICc(m2_spring6mo_lag_AR1),AICc(m2_spring6mo_lag_AR2),AICc(m3_spring6mo_lag_AR1),AICc(m3_spring6mo_lag_AR2),AICc(m1_monsoon6mo_lag_AR1),AICc(m1_monsoon6mo_lag_AR2),AICc(m2_monsoon6mo_lag_AR1),AICc(m2_monsoon6mo_lag_AR2),AICc(m3_monsoon6mo_lag_AR1),AICc(m3_monsoon6mo_lag_AR2))
  
  dAICc_m3_monsoon6mo_AR2<-(AICc(m3_monsoon6mo_AR2))-min(AICc(m_null),AICc(m1_spring6mo),AICc(m2_spring6mo),AICc(m3_spring6mo), AICc(m1_monsoon6mo), AICc(m2_monsoon6mo),AICc(m3_monsoon6mo),AICc(m1_spring6mo_lag),AICc(m2_spring6mo_lag),AICc(m3_spring6mo_lag),AICc(m1_monsoon6mo_lag),AICc(m2_monsoon6mo_lag),AICc(m3_monsoon6mo_lag),AICc(m1_spring6mo_AR1),AICc(m1_spring6mo_AR2),AICc(m2_spring6mo_AR1),AICc(m2_spring6mo_AR2),AICc(m3_spring6mo_AR1),AICc(m3_spring6mo_AR2),AICc(m1_monsoon6mo_AR1),AICc(m1_monsoon6mo_AR2),AICc(m2_monsoon6mo_AR1),AICc(m2_monsoon6mo_AR2),AICc(m3_monsoon6mo_AR1),AICc(m3_monsoon6mo_AR2),AICc(m1_spring6mo_lag_AR1),AICc(m1_spring6mo_lag_AR2),AICc(m2_spring6mo_lag_AR1),AICc(m2_spring6mo_lag_AR2),AICc(m3_spring6mo_lag_AR1),AICc(m3_spring6mo_lag_AR2),AICc(m1_monsoon6mo_lag_AR1),AICc(m1_monsoon6mo_lag_AR2),AICc(m2_monsoon6mo_lag_AR1),AICc(m2_monsoon6mo_lag_AR2),AICc(m3_monsoon6mo_lag_AR1),AICc(m3_monsoon6mo_lag_AR2))
  
  dAICc_m1_spring6mo_lag_AR1<-(AICc(m1_spring6mo_lag_AR1))-min(AICc(m_null),AICc(m1_spring6mo),AICc(m2_spring6mo),AICc(m3_spring6mo), AICc(m1_monsoon6mo), AICc(m2_monsoon6mo),AICc(m3_monsoon6mo),AICc(m1_spring6mo_lag),AICc(m2_spring6mo_lag),AICc(m3_spring6mo_lag),AICc(m1_monsoon6mo_lag),AICc(m2_monsoon6mo_lag),AICc(m3_monsoon6mo_lag),AICc(m1_spring6mo_AR1),AICc(m1_spring6mo_AR2),AICc(m2_spring6mo_AR1),AICc(m2_spring6mo_AR2),AICc(m3_spring6mo_AR1),AICc(m3_spring6mo_AR2),AICc(m1_monsoon6mo_AR1),AICc(m1_monsoon6mo_AR2),AICc(m2_monsoon6mo_AR1),AICc(m2_monsoon6mo_AR2),AICc(m3_monsoon6mo_AR1),AICc(m3_monsoon6mo_AR2),AICc(m1_spring6mo_lag_AR1),AICc(m1_spring6mo_lag_AR2),AICc(m2_spring6mo_lag_AR1),AICc(m2_spring6mo_lag_AR2),AICc(m3_spring6mo_lag_AR1),AICc(m3_spring6mo_lag_AR2),AICc(m1_monsoon6mo_lag_AR1),AICc(m1_monsoon6mo_lag_AR2),AICc(m2_monsoon6mo_lag_AR1),AICc(m2_monsoon6mo_lag_AR2),AICc(m3_monsoon6mo_lag_AR1),AICc(m3_monsoon6mo_lag_AR2))
  
  dAICc_m1_spring6mo_lag_AR2<-(AICc(m1_spring6mo_lag_AR2))-min(AICc(m_null),AICc(m1_spring6mo),AICc(m2_spring6mo),AICc(m3_spring6mo), AICc(m1_monsoon6mo), AICc(m2_monsoon6mo),AICc(m3_monsoon6mo),AICc(m1_spring6mo_lag),AICc(m2_spring6mo_lag),AICc(m3_spring6mo_lag),AICc(m1_monsoon6mo_lag),AICc(m2_monsoon6mo_lag),AICc(m3_monsoon6mo_lag),AICc(m1_spring6mo_AR1),AICc(m1_spring6mo_AR2),AICc(m2_spring6mo_AR1),AICc(m2_spring6mo_AR2),AICc(m3_spring6mo_AR1),AICc(m3_spring6mo_AR2),AICc(m1_monsoon6mo_AR1),AICc(m1_monsoon6mo_AR2),AICc(m2_monsoon6mo_AR1),AICc(m2_monsoon6mo_AR2),AICc(m3_monsoon6mo_AR1),AICc(m3_monsoon6mo_AR2),AICc(m1_spring6mo_lag_AR1),AICc(m1_spring6mo_lag_AR2),AICc(m2_spring6mo_lag_AR1),AICc(m2_spring6mo_lag_AR2),AICc(m3_spring6mo_lag_AR1),AICc(m3_spring6mo_lag_AR2),AICc(m1_monsoon6mo_lag_AR1),AICc(m1_monsoon6mo_lag_AR2),AICc(m2_monsoon6mo_lag_AR1),AICc(m2_monsoon6mo_lag_AR2),AICc(m3_monsoon6mo_lag_AR1),AICc(m3_monsoon6mo_lag_AR2))
  
  dAICc_m2_spring6mo_lag_AR1<-(AICc(m2_spring6mo_lag_AR1))-min(AICc(m_null),AICc(m1_spring6mo),AICc(m2_spring6mo),AICc(m3_spring6mo), AICc(m1_monsoon6mo), AICc(m2_monsoon6mo),AICc(m3_monsoon6mo),AICc(m1_spring6mo_lag),AICc(m2_spring6mo_lag),AICc(m3_spring6mo_lag),AICc(m1_monsoon6mo_lag),AICc(m2_monsoon6mo_lag),AICc(m3_monsoon6mo_lag),AICc(m1_spring6mo_AR1),AICc(m1_spring6mo_AR2),AICc(m2_spring6mo_AR1),AICc(m2_spring6mo_AR2),AICc(m3_spring6mo_AR1),AICc(m3_spring6mo_AR2),AICc(m1_monsoon6mo_AR1),AICc(m1_monsoon6mo_AR2),AICc(m2_monsoon6mo_AR1),AICc(m2_monsoon6mo_AR2),AICc(m3_monsoon6mo_AR1),AICc(m3_monsoon6mo_AR2),AICc(m1_spring6mo_lag_AR1),AICc(m1_spring6mo_lag_AR2),AICc(m2_spring6mo_lag_AR1),AICc(m2_spring6mo_lag_AR2),AICc(m3_spring6mo_lag_AR1),AICc(m3_spring6mo_lag_AR2),AICc(m1_monsoon6mo_lag_AR1),AICc(m1_monsoon6mo_lag_AR2),AICc(m2_monsoon6mo_lag_AR1),AICc(m2_monsoon6mo_lag_AR2),AICc(m3_monsoon6mo_lag_AR1),AICc(m3_monsoon6mo_lag_AR2))
  
  dAICc_m2_spring6mo_lag_AR2<-(AICc(m2_spring6mo_lag_AR2))-min(AICc(m_null),AICc(m1_spring6mo),AICc(m2_spring6mo),AICc(m3_spring6mo), AICc(m1_monsoon6mo), AICc(m2_monsoon6mo),AICc(m3_monsoon6mo),AICc(m1_spring6mo_lag),AICc(m2_spring6mo_lag),AICc(m3_spring6mo_lag),AICc(m1_monsoon6mo_lag),AICc(m2_monsoon6mo_lag),AICc(m3_monsoon6mo_lag),AICc(m1_spring6mo_AR1),AICc(m1_spring6mo_AR2),AICc(m2_spring6mo_AR1),AICc(m2_spring6mo_AR2),AICc(m3_spring6mo_AR1),AICc(m3_spring6mo_AR2),AICc(m1_monsoon6mo_AR1),AICc(m1_monsoon6mo_AR2),AICc(m2_monsoon6mo_AR1),AICc(m2_monsoon6mo_AR2),AICc(m3_monsoon6mo_AR1),AICc(m3_monsoon6mo_AR2),AICc(m1_spring6mo_lag_AR1),AICc(m1_spring6mo_lag_AR2),AICc(m2_spring6mo_lag_AR1),AICc(m2_spring6mo_lag_AR2),AICc(m3_spring6mo_lag_AR1),AICc(m3_spring6mo_lag_AR2),AICc(m1_monsoon6mo_lag_AR1),AICc(m1_monsoon6mo_lag_AR2),AICc(m2_monsoon6mo_lag_AR1),AICc(m2_monsoon6mo_lag_AR2),AICc(m3_monsoon6mo_lag_AR1),AICc(m3_monsoon6mo_lag_AR2))
  
  dAICc_m3_spring6mo_lag_AR1<-(AICc(m3_spring6mo_lag_AR1))-min(AICc(m_null),AICc(m1_spring6mo),AICc(m2_spring6mo),AICc(m3_spring6mo), AICc(m1_monsoon6mo), AICc(m2_monsoon6mo),AICc(m3_monsoon6mo),AICc(m1_spring6mo_lag),AICc(m2_spring6mo_lag),AICc(m3_spring6mo_lag),AICc(m1_monsoon6mo_lag),AICc(m2_monsoon6mo_lag),AICc(m3_monsoon6mo_lag),AICc(m1_spring6mo_AR1),AICc(m1_spring6mo_AR2),AICc(m2_spring6mo_AR1),AICc(m2_spring6mo_AR2),AICc(m3_spring6mo_AR1),AICc(m3_spring6mo_AR2),AICc(m1_monsoon6mo_AR1),AICc(m1_monsoon6mo_AR2),AICc(m2_monsoon6mo_AR1),AICc(m2_monsoon6mo_AR2),AICc(m3_monsoon6mo_AR1),AICc(m3_monsoon6mo_AR2),AICc(m1_spring6mo_lag_AR1),AICc(m1_spring6mo_lag_AR2),AICc(m2_spring6mo_lag_AR1),AICc(m2_spring6mo_lag_AR2),AICc(m3_spring6mo_lag_AR1),AICc(m3_spring6mo_lag_AR2),AICc(m1_monsoon6mo_lag_AR1),AICc(m1_monsoon6mo_lag_AR2),AICc(m2_monsoon6mo_lag_AR1),AICc(m2_monsoon6mo_lag_AR2),AICc(m3_monsoon6mo_lag_AR1),AICc(m3_monsoon6mo_lag_AR2))
  
  dAICc_m3_spring6mo_lag_AR2<-(AICc(m3_spring6mo_lag_AR2))-min(AICc(m_null),AICc(m1_spring6mo),AICc(m2_spring6mo),AICc(m3_spring6mo), AICc(m1_monsoon6mo), AICc(m2_monsoon6mo),AICc(m3_monsoon6mo),AICc(m1_spring6mo_lag),AICc(m2_spring6mo_lag),AICc(m3_spring6mo_lag),AICc(m1_monsoon6mo_lag),AICc(m2_monsoon6mo_lag),AICc(m3_monsoon6mo_lag),AICc(m1_spring6mo_AR1),AICc(m1_spring6mo_AR2),AICc(m2_spring6mo_AR1),AICc(m2_spring6mo_AR2),AICc(m3_spring6mo_AR1),AICc(m3_spring6mo_AR2),AICc(m1_monsoon6mo_AR1),AICc(m1_monsoon6mo_AR2),AICc(m2_monsoon6mo_AR1),AICc(m2_monsoon6mo_AR2),AICc(m3_monsoon6mo_AR1),AICc(m3_monsoon6mo_AR2),AICc(m1_spring6mo_lag_AR1),AICc(m1_spring6mo_lag_AR2),AICc(m2_spring6mo_lag_AR1),AICc(m2_spring6mo_lag_AR2),AICc(m3_spring6mo_lag_AR1),AICc(m3_spring6mo_lag_AR2),AICc(m1_monsoon6mo_lag_AR1),AICc(m1_monsoon6mo_lag_AR2),AICc(m2_monsoon6mo_lag_AR1),AICc(m2_monsoon6mo_lag_AR2),AICc(m3_monsoon6mo_lag_AR1),AICc(m3_monsoon6mo_lag_AR2))
  
  dAICc_m1_monsoon6mo_lag_AR1<-(AICc(m1_monsoon6mo_lag_AR1))-min(AICc(m_null),AICc(m1_spring6mo),AICc(m2_spring6mo),AICc(m3_spring6mo), AICc(m1_monsoon6mo), AICc(m2_monsoon6mo),AICc(m3_monsoon6mo),AICc(m1_spring6mo_lag),AICc(m2_spring6mo_lag),AICc(m3_spring6mo_lag),AICc(m1_monsoon6mo_lag),AICc(m2_monsoon6mo_lag),AICc(m3_monsoon6mo_lag),AICc(m1_spring6mo_AR1),AICc(m1_spring6mo_AR2),AICc(m2_spring6mo_AR1),AICc(m2_spring6mo_AR2),AICc(m3_spring6mo_AR1),AICc(m3_spring6mo_AR2),AICc(m1_monsoon6mo_AR1),AICc(m1_monsoon6mo_AR2),AICc(m2_monsoon6mo_AR1),AICc(m2_monsoon6mo_AR2),AICc(m3_monsoon6mo_AR1),AICc(m3_monsoon6mo_AR2),AICc(m1_spring6mo_lag_AR1),AICc(m1_spring6mo_lag_AR2),AICc(m2_spring6mo_lag_AR1),AICc(m2_spring6mo_lag_AR2),AICc(m3_spring6mo_lag_AR1),AICc(m3_spring6mo_lag_AR2),AICc(m1_monsoon6mo_lag_AR1),AICc(m1_monsoon6mo_lag_AR2),AICc(m2_monsoon6mo_lag_AR1),AICc(m2_monsoon6mo_lag_AR2),AICc(m3_monsoon6mo_lag_AR1),AICc(m3_monsoon6mo_lag_AR2))
  
  dAICc_m1_monsoon6mo_lag_AR2<-(AICc(m1_monsoon6mo_lag_AR2))-min(AICc(m_null),AICc(m1_spring6mo),AICc(m2_spring6mo),AICc(m3_spring6mo), AICc(m1_monsoon6mo), AICc(m2_monsoon6mo),AICc(m3_monsoon6mo),AICc(m1_spring6mo_lag),AICc(m2_spring6mo_lag),AICc(m3_spring6mo_lag),AICc(m1_monsoon6mo_lag),AICc(m2_monsoon6mo_lag),AICc(m3_monsoon6mo_lag),AICc(m1_spring6mo_AR1),AICc(m1_spring6mo_AR2),AICc(m2_spring6mo_AR1),AICc(m2_spring6mo_AR2),AICc(m3_spring6mo_AR1),AICc(m3_spring6mo_AR2),AICc(m1_monsoon6mo_AR1),AICc(m1_monsoon6mo_AR2),AICc(m2_monsoon6mo_AR1),AICc(m2_monsoon6mo_AR2),AICc(m3_monsoon6mo_AR1),AICc(m3_monsoon6mo_AR2),AICc(m1_spring6mo_lag_AR1),AICc(m1_spring6mo_lag_AR2),AICc(m2_spring6mo_lag_AR1),AICc(m2_spring6mo_lag_AR2),AICc(m3_spring6mo_lag_AR1),AICc(m3_spring6mo_lag_AR2),AICc(m1_monsoon6mo_lag_AR1),AICc(m1_monsoon6mo_lag_AR2),AICc(m2_monsoon6mo_lag_AR1),AICc(m2_monsoon6mo_lag_AR2),AICc(m3_monsoon6mo_lag_AR1),AICc(m3_monsoon6mo_lag_AR2))
  
  dAICc_m2_monsoon6mo_lag_AR1<-(AICc(m2_monsoon6mo_lag_AR1))-min(AICc(m_null),AICc(m1_spring6mo),AICc(m2_spring6mo),AICc(m3_spring6mo), AICc(m1_monsoon6mo), AICc(m2_monsoon6mo),AICc(m3_monsoon6mo),AICc(m1_spring6mo_lag),AICc(m2_spring6mo_lag),AICc(m3_spring6mo_lag),AICc(m1_monsoon6mo_lag),AICc(m2_monsoon6mo_lag),AICc(m3_monsoon6mo_lag),AICc(m1_spring6mo_AR1),AICc(m1_spring6mo_AR2),AICc(m2_spring6mo_AR1),AICc(m2_spring6mo_AR2),AICc(m3_spring6mo_AR1),AICc(m3_spring6mo_AR2),AICc(m1_monsoon6mo_AR1),AICc(m1_monsoon6mo_AR2),AICc(m2_monsoon6mo_AR1),AICc(m2_monsoon6mo_AR2),AICc(m3_monsoon6mo_AR1),AICc(m3_monsoon6mo_AR2),AICc(m1_spring6mo_lag_AR1),AICc(m1_spring6mo_lag_AR2),AICc(m2_spring6mo_lag_AR1),AICc(m2_spring6mo_lag_AR2),AICc(m3_spring6mo_lag_AR1),AICc(m3_spring6mo_lag_AR2),AICc(m1_monsoon6mo_lag_AR1),AICc(m1_monsoon6mo_lag_AR2),AICc(m2_monsoon6mo_lag_AR1),AICc(m2_monsoon6mo_lag_AR2),AICc(m3_monsoon6mo_lag_AR1),AICc(m3_monsoon6mo_lag_AR2))
  
  dAICc_m2_monsoon6mo_lag_AR2<-(AICc(m2_monsoon6mo_lag_AR2))-min(AICc(m_null),AICc(m1_spring6mo),AICc(m2_spring6mo),AICc(m3_spring6mo), AICc(m1_monsoon6mo), AICc(m2_monsoon6mo),AICc(m3_monsoon6mo),AICc(m1_spring6mo_lag),AICc(m2_spring6mo_lag),AICc(m3_spring6mo_lag),AICc(m1_monsoon6mo_lag),AICc(m2_monsoon6mo_lag),AICc(m3_monsoon6mo_lag),AICc(m1_spring6mo_AR1),AICc(m1_spring6mo_AR2),AICc(m2_spring6mo_AR1),AICc(m2_spring6mo_AR2),AICc(m3_spring6mo_AR1),AICc(m3_spring6mo_AR2),AICc(m1_monsoon6mo_AR1),AICc(m1_monsoon6mo_AR2),AICc(m2_monsoon6mo_AR1),AICc(m2_monsoon6mo_AR2),AICc(m3_monsoon6mo_AR1),AICc(m3_monsoon6mo_AR2),AICc(m1_spring6mo_lag_AR1),AICc(m1_spring6mo_lag_AR2),AICc(m2_spring6mo_lag_AR1),AICc(m2_spring6mo_lag_AR2),AICc(m3_spring6mo_lag_AR1),AICc(m3_spring6mo_lag_AR2),AICc(m1_monsoon6mo_lag_AR1),AICc(m1_monsoon6mo_lag_AR2),AICc(m2_monsoon6mo_lag_AR1),AICc(m2_monsoon6mo_lag_AR2),AICc(m3_monsoon6mo_lag_AR1),AICc(m3_monsoon6mo_lag_AR2))
  
  dAICc_m3_monsoon6mo_lag_AR1<-(AICc(m3_monsoon6mo_lag_AR1))-min(AICc(m_null),AICc(m1_spring6mo),AICc(m2_spring6mo),AICc(m3_spring6mo), AICc(m1_monsoon6mo), AICc(m2_monsoon6mo),AICc(m3_monsoon6mo),AICc(m1_spring6mo_lag),AICc(m2_spring6mo_lag),AICc(m3_spring6mo_lag),AICc(m1_monsoon6mo_lag),AICc(m2_monsoon6mo_lag),AICc(m3_monsoon6mo_lag),AICc(m1_spring6mo_AR1),AICc(m1_spring6mo_AR2),AICc(m2_spring6mo_AR1),AICc(m2_spring6mo_AR2),AICc(m3_spring6mo_AR1),AICc(m3_spring6mo_AR2),AICc(m1_monsoon6mo_AR1),AICc(m1_monsoon6mo_AR2),AICc(m2_monsoon6mo_AR1),AICc(m2_monsoon6mo_AR2),AICc(m3_monsoon6mo_AR1),AICc(m3_monsoon6mo_AR2),AICc(m1_spring6mo_lag_AR1),AICc(m1_spring6mo_lag_AR2),AICc(m2_spring6mo_lag_AR1),AICc(m2_spring6mo_lag_AR2),AICc(m3_spring6mo_lag_AR1),AICc(m3_spring6mo_lag_AR2),AICc(m1_monsoon6mo_lag_AR1),AICc(m1_monsoon6mo_lag_AR2),AICc(m2_monsoon6mo_lag_AR1),AICc(m2_monsoon6mo_lag_AR2),AICc(m3_monsoon6mo_lag_AR1),AICc(m3_monsoon6mo_lag_AR2))
  
  dAICc_m3_monsoon6mo_lag_AR2<-(AICc(m3_monsoon6mo_lag_AR2))-min(AICc(m_null),AICc(m1_spring6mo),AICc(m2_spring6mo),AICc(m3_spring6mo), AICc(m1_monsoon6mo), AICc(m2_monsoon6mo),AICc(m3_monsoon6mo),AICc(m1_spring6mo_lag),AICc(m2_spring6mo_lag),AICc(m3_spring6mo_lag),AICc(m1_monsoon6mo_lag),AICc(m2_monsoon6mo_lag),AICc(m3_monsoon6mo_lag),AICc(m1_spring6mo_AR1),AICc(m1_spring6mo_AR2),AICc(m2_spring6mo_AR1),AICc(m2_spring6mo_AR2),AICc(m3_spring6mo_AR1),AICc(m3_spring6mo_AR2),AICc(m1_monsoon6mo_AR1),AICc(m1_monsoon6mo_AR2),AICc(m2_monsoon6mo_AR1),AICc(m2_monsoon6mo_AR2),AICc(m3_monsoon6mo_AR1),AICc(m3_monsoon6mo_AR2),AICc(m1_spring6mo_lag_AR1),AICc(m1_spring6mo_lag_AR2),AICc(m2_spring6mo_lag_AR1),AICc(m2_spring6mo_lag_AR2),AICc(m3_spring6mo_lag_AR1),AICc(m3_spring6mo_lag_AR2),AICc(m1_monsoon6mo_lag_AR1),AICc(m1_monsoon6mo_lag_AR2),AICc(m2_monsoon6mo_lag_AR1),AICc(m2_monsoon6mo_lag_AR2),AICc(m3_monsoon6mo_lag_AR1),AICc(m3_monsoon6mo_lag_AR2))
  
  # delta AICc for other models
  dAICc_null<-(AICc(m_null))-min(AICc(m_null),AICc(m1_spring6mo),AICc(m2_spring6mo),AICc(m3_spring6mo), AICc(m1_monsoon6mo), AICc(m2_monsoon6mo),AICc(m3_monsoon6mo),AICc(m1_spring6mo_lag),AICc(m2_spring6mo_lag),AICc(m3_spring6mo_lag),AICc(m1_monsoon6mo_lag),AICc(m2_monsoon6mo_lag),AICc(m3_monsoon6mo_lag),AICc(m1_spring6mo_AR1),AICc(m1_spring6mo_AR2),AICc(m2_spring6mo_AR1),AICc(m2_spring6mo_AR2),AICc(m3_spring6mo_AR1),AICc(m3_spring6mo_AR2),AICc(m1_monsoon6mo_AR1),AICc(m1_monsoon6mo_AR2),AICc(m2_monsoon6mo_AR1),AICc(m2_monsoon6mo_AR2),AICc(m3_monsoon6mo_AR1),AICc(m3_monsoon6mo_AR2),AICc(m1_spring6mo_lag_AR1),AICc(m1_spring6mo_lag_AR2),AICc(m2_spring6mo_lag_AR1),AICc(m2_spring6mo_lag_AR2),AICc(m3_spring6mo_lag_AR1),AICc(m3_spring6mo_lag_AR2),AICc(m1_monsoon6mo_lag_AR1),AICc(m1_monsoon6mo_lag_AR2),AICc(m2_monsoon6mo_lag_AR1),AICc(m2_monsoon6mo_lag_AR2),AICc(m3_monsoon6mo_lag_AR1),AICc(m3_monsoon6mo_lag_AR2))
  
  dAICc_m1_spring6mo<-(AICc(m1_spring6mo))-min(AICc(m_null),AICc(m1_spring6mo),AICc(m2_spring6mo),AICc(m3_spring6mo), AICc(m1_monsoon6mo), AICc(m2_monsoon6mo),AICc(m3_monsoon6mo),AICc(m1_spring6mo_lag),AICc(m2_spring6mo_lag),AICc(m3_spring6mo_lag),AICc(m1_monsoon6mo_lag),AICc(m2_monsoon6mo_lag),AICc(m3_monsoon6mo_lag),AICc(m1_spring6mo_AR1),AICc(m1_spring6mo_AR2),AICc(m2_spring6mo_AR1),AICc(m2_spring6mo_AR2),AICc(m3_spring6mo_AR1),AICc(m3_spring6mo_AR2),AICc(m1_monsoon6mo_AR1),AICc(m1_monsoon6mo_AR2),AICc(m2_monsoon6mo_AR1),AICc(m2_monsoon6mo_AR2),AICc(m3_monsoon6mo_AR1),AICc(m3_monsoon6mo_AR2),AICc(m1_spring6mo_lag_AR1),AICc(m1_spring6mo_lag_AR2),AICc(m2_spring6mo_lag_AR1),AICc(m2_spring6mo_lag_AR2),AICc(m3_spring6mo_lag_AR1),AICc(m3_spring6mo_lag_AR2),AICc(m1_monsoon6mo_lag_AR1),AICc(m1_monsoon6mo_lag_AR2),AICc(m2_monsoon6mo_lag_AR1),AICc(m2_monsoon6mo_lag_AR2),AICc(m3_monsoon6mo_lag_AR1),AICc(m3_monsoon6mo_lag_AR2))
  
  dAICc_m2_spring6mo<-(AICc(m2_spring6mo))-min(AICc(m_null),AICc(m1_spring6mo),AICc(m2_spring6mo),AICc(m3_spring6mo), AICc(m1_monsoon6mo), AICc(m2_monsoon6mo),AICc(m3_monsoon6mo),AICc(m1_spring6mo_lag),AICc(m2_spring6mo_lag),AICc(m3_spring6mo_lag),AICc(m1_monsoon6mo_lag),AICc(m2_monsoon6mo_lag),AICc(m3_monsoon6mo_lag),AICc(m1_spring6mo_AR1),AICc(m1_spring6mo_AR2),AICc(m2_spring6mo_AR1),AICc(m2_spring6mo_AR2),AICc(m3_spring6mo_AR1),AICc(m3_spring6mo_AR2),AICc(m1_monsoon6mo_AR1),AICc(m1_monsoon6mo_AR2),AICc(m2_monsoon6mo_AR1),AICc(m2_monsoon6mo_AR2),AICc(m3_monsoon6mo_AR1),AICc(m3_monsoon6mo_AR2),AICc(m1_spring6mo_lag_AR1),AICc(m1_spring6mo_lag_AR2),AICc(m2_spring6mo_lag_AR1),AICc(m2_spring6mo_lag_AR2),AICc(m3_spring6mo_lag_AR1),AICc(m3_spring6mo_lag_AR2),AICc(m1_monsoon6mo_lag_AR1),AICc(m1_monsoon6mo_lag_AR2),AICc(m2_monsoon6mo_lag_AR1),AICc(m2_monsoon6mo_lag_AR2),AICc(m3_monsoon6mo_lag_AR1),AICc(m3_monsoon6mo_lag_AR2))
  
  dAICc_m3_spring6mo<-(AICc(m3_spring6mo))-min(AICc(m_null),AICc(m1_spring6mo),AICc(m2_spring6mo),AICc(m3_spring6mo), AICc(m1_monsoon6mo), AICc(m2_monsoon6mo),AICc(m3_monsoon6mo),AICc(m1_spring6mo_lag),AICc(m2_spring6mo_lag),AICc(m3_spring6mo_lag),AICc(m1_monsoon6mo_lag),AICc(m2_monsoon6mo_lag),AICc(m3_monsoon6mo_lag),AICc(m1_spring6mo_AR1),AICc(m1_spring6mo_AR2),AICc(m2_spring6mo_AR1),AICc(m2_spring6mo_AR2),AICc(m3_spring6mo_AR1),AICc(m3_spring6mo_AR2),AICc(m1_monsoon6mo_AR1),AICc(m1_monsoon6mo_AR2),AICc(m2_monsoon6mo_AR1),AICc(m2_monsoon6mo_AR2),AICc(m3_monsoon6mo_AR1),AICc(m3_monsoon6mo_AR2),AICc(m1_spring6mo_lag_AR1),AICc(m1_spring6mo_lag_AR2),AICc(m2_spring6mo_lag_AR1),AICc(m2_spring6mo_lag_AR2),AICc(m3_spring6mo_lag_AR1),AICc(m3_spring6mo_lag_AR2),AICc(m1_monsoon6mo_lag_AR1),AICc(m1_monsoon6mo_lag_AR2),AICc(m2_monsoon6mo_lag_AR1),AICc(m2_monsoon6mo_lag_AR2),AICc(m3_monsoon6mo_lag_AR1),AICc(m3_monsoon6mo_lag_AR2))
  
  dAICc_m1_monsoon6mo<-(AICc(m1_monsoon6mo))-min(AICc(m_null),AICc(m1_spring6mo),AICc(m2_spring6mo),AICc(m3_spring6mo), AICc(m1_monsoon6mo), AICc(m2_monsoon6mo),AICc(m3_monsoon6mo),AICc(m1_spring6mo_lag),AICc(m2_spring6mo_lag),AICc(m3_spring6mo_lag),AICc(m1_monsoon6mo_lag),AICc(m2_monsoon6mo_lag),AICc(m3_monsoon6mo_lag),AICc(m1_spring6mo_AR1),AICc(m1_spring6mo_AR2),AICc(m2_spring6mo_AR1),AICc(m2_spring6mo_AR2),AICc(m3_spring6mo_AR1),AICc(m3_spring6mo_AR2),AICc(m1_monsoon6mo_AR1),AICc(m1_monsoon6mo_AR2),AICc(m2_monsoon6mo_AR1),AICc(m2_monsoon6mo_AR2),AICc(m3_monsoon6mo_AR1),AICc(m3_monsoon6mo_AR2),AICc(m1_spring6mo_lag_AR1),AICc(m1_spring6mo_lag_AR2),AICc(m2_spring6mo_lag_AR1),AICc(m2_spring6mo_lag_AR2),AICc(m3_spring6mo_lag_AR1),AICc(m3_spring6mo_lag_AR2),AICc(m1_monsoon6mo_lag_AR1),AICc(m1_monsoon6mo_lag_AR2),AICc(m2_monsoon6mo_lag_AR1),AICc(m2_monsoon6mo_lag_AR2),AICc(m3_monsoon6mo_lag_AR1),AICc(m3_monsoon6mo_lag_AR2))
  
  dAICc_m2_monsoon6mo<-(AICc(m2_monsoon6mo))-min(AICc(m_null),AICc(m1_spring6mo),AICc(m2_spring6mo),AICc(m3_spring6mo), AICc(m1_monsoon6mo), AICc(m2_monsoon6mo),AICc(m3_monsoon6mo),AICc(m1_spring6mo_lag),AICc(m2_spring6mo_lag),AICc(m3_spring6mo_lag),AICc(m1_monsoon6mo_lag),AICc(m2_monsoon6mo_lag),AICc(m3_monsoon6mo_lag),AICc(m1_spring6mo_AR1),AICc(m1_spring6mo_AR2),AICc(m2_spring6mo_AR1),AICc(m2_spring6mo_AR2),AICc(m3_spring6mo_AR1),AICc(m3_spring6mo_AR2),AICc(m1_monsoon6mo_AR1),AICc(m1_monsoon6mo_AR2),AICc(m2_monsoon6mo_AR1),AICc(m2_monsoon6mo_AR2),AICc(m3_monsoon6mo_AR1),AICc(m3_monsoon6mo_AR2),AICc(m1_spring6mo_lag_AR1),AICc(m1_spring6mo_lag_AR2),AICc(m2_spring6mo_lag_AR1),AICc(m2_spring6mo_lag_AR2),AICc(m3_spring6mo_lag_AR1),AICc(m3_spring6mo_lag_AR2),AICc(m1_monsoon6mo_lag_AR1),AICc(m1_monsoon6mo_lag_AR2),AICc(m2_monsoon6mo_lag_AR1),AICc(m2_monsoon6mo_lag_AR2),AICc(m3_monsoon6mo_lag_AR1),AICc(m3_monsoon6mo_lag_AR2))
  
  dAICc_m3_monsoon6mo<-(AICc(m3_monsoon6mo))-min(AICc(m_null),AICc(m1_spring6mo),AICc(m2_spring6mo),AICc(m3_spring6mo), AICc(m1_monsoon6mo), AICc(m2_monsoon6mo),AICc(m3_monsoon6mo),AICc(m1_spring6mo_lag),AICc(m2_spring6mo_lag),AICc(m3_spring6mo_lag),AICc(m1_monsoon6mo_lag),AICc(m2_monsoon6mo_lag),AICc(m3_monsoon6mo_lag),AICc(m1_spring6mo_AR1),AICc(m1_spring6mo_AR2),AICc(m2_spring6mo_AR1),AICc(m2_spring6mo_AR2),AICc(m3_spring6mo_AR1),AICc(m3_spring6mo_AR2),AICc(m1_monsoon6mo_AR1),AICc(m1_monsoon6mo_AR2),AICc(m2_monsoon6mo_AR1),AICc(m2_monsoon6mo_AR2),AICc(m3_monsoon6mo_AR1),AICc(m3_monsoon6mo_AR2),AICc(m1_spring6mo_lag_AR1),AICc(m1_spring6mo_lag_AR2),AICc(m2_spring6mo_lag_AR1),AICc(m2_spring6mo_lag_AR2),AICc(m3_spring6mo_lag_AR1),AICc(m3_spring6mo_lag_AR2),AICc(m1_monsoon6mo_lag_AR1),AICc(m1_monsoon6mo_lag_AR2),AICc(m2_monsoon6mo_lag_AR1),AICc(m2_monsoon6mo_lag_AR2),AICc(m3_monsoon6mo_lag_AR1),AICc(m3_monsoon6mo_lag_AR2))
  
  dAICc_m1_spring6mo_lag<-(AICc(m1_spring6mo_lag))-min(AICc(m_null),AICc(m1_spring6mo),AICc(m2_spring6mo),AICc(m3_spring6mo), AICc(m1_monsoon6mo), AICc(m2_monsoon6mo),AICc(m3_monsoon6mo),AICc(m1_spring6mo_lag),AICc(m2_spring6mo_lag),AICc(m3_spring6mo_lag),AICc(m1_monsoon6mo_lag),AICc(m2_monsoon6mo_lag),AICc(m3_monsoon6mo_lag),AICc(m1_spring6mo_AR1),AICc(m1_spring6mo_AR2),AICc(m2_spring6mo_AR1),AICc(m2_spring6mo_AR2),AICc(m3_spring6mo_AR1),AICc(m3_spring6mo_AR2),AICc(m1_monsoon6mo_AR1),AICc(m1_monsoon6mo_AR2),AICc(m2_monsoon6mo_AR1),AICc(m2_monsoon6mo_AR2),AICc(m3_monsoon6mo_AR1),AICc(m3_monsoon6mo_AR2),AICc(m1_spring6mo_lag_AR1),AICc(m1_spring6mo_lag_AR2),AICc(m2_spring6mo_lag_AR1),AICc(m2_spring6mo_lag_AR2),AICc(m3_spring6mo_lag_AR1),AICc(m3_spring6mo_lag_AR2),AICc(m1_monsoon6mo_lag_AR1),AICc(m1_monsoon6mo_lag_AR2),AICc(m2_monsoon6mo_lag_AR1),AICc(m2_monsoon6mo_lag_AR2),AICc(m3_monsoon6mo_lag_AR1),AICc(m3_monsoon6mo_lag_AR2))
  
  dAICc_m2_spring6mo_lag<-(AICc(m2_spring6mo_lag))-min(AICc(m_null),AICc(m1_spring6mo),AICc(m2_spring6mo),AICc(m3_spring6mo), AICc(m1_monsoon6mo), AICc(m2_monsoon6mo),AICc(m3_monsoon6mo),AICc(m1_spring6mo_lag),AICc(m2_spring6mo_lag),AICc(m3_spring6mo_lag),AICc(m1_monsoon6mo_lag),AICc(m2_monsoon6mo_lag),AICc(m3_monsoon6mo_lag),AICc(m1_spring6mo_AR1),AICc(m1_spring6mo_AR2),AICc(m2_spring6mo_AR1),AICc(m2_spring6mo_AR2),AICc(m3_spring6mo_AR1),AICc(m3_spring6mo_AR2),AICc(m1_monsoon6mo_AR1),AICc(m1_monsoon6mo_AR2),AICc(m2_monsoon6mo_AR1),AICc(m2_monsoon6mo_AR2),AICc(m3_monsoon6mo_AR1),AICc(m3_monsoon6mo_AR2),AICc(m1_spring6mo_lag_AR1),AICc(m1_spring6mo_lag_AR2),AICc(m2_spring6mo_lag_AR1),AICc(m2_spring6mo_lag_AR2),AICc(m3_spring6mo_lag_AR1),AICc(m3_spring6mo_lag_AR2),AICc(m1_monsoon6mo_lag_AR1),AICc(m1_monsoon6mo_lag_AR2),AICc(m2_monsoon6mo_lag_AR1),AICc(m2_monsoon6mo_lag_AR2),AICc(m3_monsoon6mo_lag_AR1),AICc(m3_monsoon6mo_lag_AR2))
  
  dAICc_m3_spring6mo_lag<-(AICc(m3_spring6mo_lag))-min(AICc(m_null),AICc(m1_spring6mo),AICc(m2_spring6mo),AICc(m3_spring6mo), AICc(m1_monsoon6mo), AICc(m2_monsoon6mo),AICc(m3_monsoon6mo),AICc(m1_spring6mo_lag),AICc(m2_spring6mo_lag),AICc(m3_spring6mo_lag),AICc(m1_monsoon6mo_lag),AICc(m2_monsoon6mo_lag),AICc(m3_monsoon6mo_lag),AICc(m1_spring6mo_AR1),AICc(m1_spring6mo_AR2),AICc(m2_spring6mo_AR1),AICc(m2_spring6mo_AR2),AICc(m3_spring6mo_AR1),AICc(m3_spring6mo_AR2),AICc(m1_monsoon6mo_AR1),AICc(m1_monsoon6mo_AR2),AICc(m2_monsoon6mo_AR1),AICc(m2_monsoon6mo_AR2),AICc(m3_monsoon6mo_AR1),AICc(m3_monsoon6mo_AR2),AICc(m1_spring6mo_lag_AR1),AICc(m1_spring6mo_lag_AR2),AICc(m2_spring6mo_lag_AR1),AICc(m2_spring6mo_lag_AR2),AICc(m3_spring6mo_lag_AR1),AICc(m3_spring6mo_lag_AR2),AICc(m1_monsoon6mo_lag_AR1),AICc(m1_monsoon6mo_lag_AR2),AICc(m2_monsoon6mo_lag_AR1),AICc(m2_monsoon6mo_lag_AR2),AICc(m3_monsoon6mo_lag_AR1),AICc(m3_monsoon6mo_lag_AR2))
  
  dAICc_m1_monsoon6mo_lag<-(AICc(m1_monsoon6mo_lag))-min(AICc(m_null),AICc(m1_spring6mo),AICc(m2_spring6mo),AICc(m3_spring6mo), AICc(m1_monsoon6mo), AICc(m2_monsoon6mo),AICc(m3_monsoon6mo),AICc(m1_spring6mo_lag),AICc(m2_spring6mo_lag),AICc(m3_spring6mo_lag),AICc(m1_monsoon6mo_lag),AICc(m2_monsoon6mo_lag),AICc(m3_monsoon6mo_lag),AICc(m1_spring6mo_AR1),AICc(m1_spring6mo_AR2),AICc(m2_spring6mo_AR1),AICc(m2_spring6mo_AR2),AICc(m3_spring6mo_AR1),AICc(m3_spring6mo_AR2),AICc(m1_monsoon6mo_AR1),AICc(m1_monsoon6mo_AR2),AICc(m2_monsoon6mo_AR1),AICc(m2_monsoon6mo_AR2),AICc(m3_monsoon6mo_AR1),AICc(m3_monsoon6mo_AR2),AICc(m1_spring6mo_lag_AR1),AICc(m1_spring6mo_lag_AR2),AICc(m2_spring6mo_lag_AR1),AICc(m2_spring6mo_lag_AR2),AICc(m3_spring6mo_lag_AR1),AICc(m3_spring6mo_lag_AR2),AICc(m1_monsoon6mo_lag_AR1),AICc(m1_monsoon6mo_lag_AR2),AICc(m2_monsoon6mo_lag_AR1),AICc(m2_monsoon6mo_lag_AR2),AICc(m3_monsoon6mo_lag_AR1),AICc(m3_monsoon6mo_lag_AR2))
  
  dAICc_m2_monsoon6mo_lag<-(AICc(m2_monsoon6mo_lag))-min(AICc(m_null),AICc(m1_spring6mo),AICc(m2_spring6mo),AICc(m3_spring6mo), AICc(m1_monsoon6mo), AICc(m2_monsoon6mo),AICc(m3_monsoon6mo),AICc(m1_spring6mo_lag),AICc(m2_spring6mo_lag),AICc(m3_spring6mo_lag),AICc(m1_monsoon6mo_lag),AICc(m2_monsoon6mo_lag),AICc(m3_monsoon6mo_lag),AICc(m1_spring6mo_AR1),AICc(m1_spring6mo_AR2),AICc(m2_spring6mo_AR1),AICc(m2_spring6mo_AR2),AICc(m3_spring6mo_AR1),AICc(m3_spring6mo_AR2),AICc(m1_monsoon6mo_AR1),AICc(m1_monsoon6mo_AR2),AICc(m2_monsoon6mo_AR1),AICc(m2_monsoon6mo_AR2),AICc(m3_monsoon6mo_AR1),AICc(m3_monsoon6mo_AR2),AICc(m1_spring6mo_lag_AR1),AICc(m1_spring6mo_lag_AR2),AICc(m2_spring6mo_lag_AR1),AICc(m2_spring6mo_lag_AR2),AICc(m3_spring6mo_lag_AR1),AICc(m3_spring6mo_lag_AR2),AICc(m1_monsoon6mo_lag_AR1),AICc(m1_monsoon6mo_lag_AR2),AICc(m2_monsoon6mo_lag_AR1),AICc(m2_monsoon6mo_lag_AR2),AICc(m3_monsoon6mo_lag_AR1),AICc(m3_monsoon6mo_lag_AR2))
  
  dAICc_m3_monsoon6mo_lag<-(AICc(m3_monsoon6mo_lag))-min(AICc(m_null),AICc(m1_spring6mo),AICc(m2_spring6mo),AICc(m3_spring6mo), AICc(m1_monsoon6mo), AICc(m2_monsoon6mo),AICc(m3_monsoon6mo),AICc(m1_spring6mo_lag),AICc(m2_spring6mo_lag),AICc(m3_spring6mo_lag),AICc(m1_monsoon6mo_lag),AICc(m2_monsoon6mo_lag),AICc(m3_monsoon6mo_lag),AICc(m1_spring6mo_AR1),AICc(m1_spring6mo_AR2),AICc(m2_spring6mo_AR1),AICc(m2_spring6mo_AR2),AICc(m3_spring6mo_AR1),AICc(m3_spring6mo_AR2),AICc(m1_monsoon6mo_AR1),AICc(m1_monsoon6mo_AR2),AICc(m2_monsoon6mo_AR1),AICc(m2_monsoon6mo_AR2),AICc(m3_monsoon6mo_AR1),AICc(m3_monsoon6mo_AR2),AICc(m1_spring6mo_lag_AR1),AICc(m1_spring6mo_lag_AR2),AICc(m2_spring6mo_lag_AR1),AICc(m2_spring6mo_lag_AR2),AICc(m3_spring6mo_lag_AR1),AICc(m3_spring6mo_lag_AR2),AICc(m1_monsoon6mo_lag_AR1),AICc(m1_monsoon6mo_lag_AR2),AICc(m2_monsoon6mo_lag_AR1),AICc(m2_monsoon6mo_lag_AR2),AICc(m3_monsoon6mo_lag_AR1),AICc(m3_monsoon6mo_lag_AR2))
  
  # graph the best CSF for the species (maximum abundance as a function of SPEI), if the best model is not the null, and save the graph
  
  if (dAICc_m1_spring6mo == 0 | dAICc_m1_spring6mo_AR1 == 0 | dAICc_m1_spring6mo_AR2 == 0) {
    
    speciesData <- creo_original[speciesCode]
    speciesData <-cbind(speciesData,creo_climate)
    names(speciesData)[1] <- "speciesCode"
    
    p <- ggplot(speciesData,aes(x=spring6SPEI,y=speciesCode^2))+
      geom_point(size=3,color="darkgreen")+
      geom_smooth(formula=y~poly(x,1),method="lm",color="black")+
      xlab("Spring SPEI\n")+
      ylab("Maximum abundance") +       ggtitle(speciesCode) +
      ggtitle(speciesCode) +
      theme(legend.position = "none") +
      theme_classic() +
      theme(axis.text.x = element_text(color="black",size=14)) +
      theme(axis.text.y = element_text(color="black",size=14)) +
      theme(axis.title=element_text(size=18))+
      theme(plot.title = element_text(size=20))
    print(p)
    
    plot_name <- paste("C",speciesCodes[i], "plot",number[i],".pdf", sep="_")
    
    ggsave(paste(plot_name), p,
           width=4,height=4,units = c("in"),
           dpi = 600)
    
    speciesData$speciesCode <-NULL
  }
  
  if (dAICc_m2_spring6mo == 0 | dAICc_m2_spring6mo_AR1 == 0 | dAICc_m2_spring6mo_AR2 == 0) {
    speciesData <- creo_original[speciesCode]
    speciesData <-cbind(speciesData,creo_climate)
    names(speciesData)[1] <- "speciesCode"
    
    p <- ggplot(speciesData,aes(x=spring6SPEI,y=speciesCode^2))+
      geom_point(size=3,color="darkgreen")+
      geom_smooth(formula=y~poly(x,2),method="lm",color="black")+
      xlab("Spring SPEI\n")+
      ylab("Maximum abundance") +       ggtitle(speciesCode) +
      theme(legend.position = "none") +
      theme_classic() +
      theme(axis.text.x = element_text(color="black",size=14)) +
      theme(axis.text.y = element_text(color="black",size=14)) +
      theme(axis.title=element_text(size=18))+
      theme(plot.title = element_text(size=20))
    print(p)
    
    plot_name <- paste("C",speciesCodes[i], "plot",number[i],".pdf", sep="_")
    
    ggsave(paste(plot_name), p,
           width=4,height=4,units = c("in"),
           dpi = 600)
    
    speciesData$speciesCode <-NULL
  }
  
  if (dAICc_m3_spring6mo == 0 | dAICc_m3_spring6mo_AR1 == 0 | dAICc_m3_spring6mo_AR2 == 0) {
    speciesData <- creo_original[speciesCode]
    speciesData <-cbind(speciesData,creo_climate)
    names(speciesData)[1] <- "speciesCode"
    
    p <- ggplot(speciesData,aes(x=spring6SPEI,y=speciesCode^2))+
      geom_point(size=3,color="darkgreen")+
      geom_smooth(formula=y~poly(x,3),method="lm",color="black")+
      xlab("Spring SPEI\n")+
      ylab("Maximum abundance") +       ggtitle(speciesCode) +
      theme(legend.position = "none") +
      theme_classic() +
      theme(axis.text.x = element_text(color="black",size=14)) +
      theme(axis.text.y = element_text(color="black",size=14)) +
      theme(axis.title=element_text(size=18))+
      theme(plot.title = element_text(size=20))
    print(p)
    
    plot_name <- paste("C",speciesCodes[i], "plot",number[i],".pdf", sep="_")
    
    ggsave(paste(plot_name), p,
           width=4,height=4,units = c("in"),
           dpi = 600)    
    speciesData$speciesCode <-NULL
  }
  
  if (dAICc_m1_monsoon6mo == 0 | dAICc_m1_monsoon6mo_AR1 == 0 | dAICc_m1_monsoon6mo_AR2 == 0) {
    speciesData <- creo_original[speciesCode]
    speciesData <-cbind(speciesData,creo_climate)
    names(speciesData)[1] <- "speciesCode"
    
    p <- ggplot(speciesData,aes(x=monsoon6SPEI,y=speciesCode^2))+
      geom_point(size=3,color="darkgreen")+
      geom_smooth(formula=y~poly(x,1),method="lm",color="black")+
      xlab("Monsoon SPEI\n")+
      ylab("Maximum abundance") +       ggtitle(speciesCode) +
      theme(legend.position = "none") +
      theme_classic() +
      theme(axis.text.x = element_text(color="black",size=14)) +
      theme(axis.text.y = element_text(color="black",size=14)) +
      theme(axis.title=element_text(size=18))+
      theme(plot.title = element_text(size=20))
    print(p)
    
    plot_name <- paste("C",speciesCodes[i], "plot",number[i],".pdf", sep="_")
    
    ggsave(paste(plot_name), p,
           width=4,height=4,units = c("in"),
           dpi = 600)
    
    speciesData$speciesCode <-NULL
  }
  
  if (dAICc_m2_monsoon6mo == 0 | dAICc_m2_monsoon6mo_AR1 == 0 | dAICc_m2_monsoon6mo_AR2 == 0) {
    speciesData <- creo_original[speciesCode]
    speciesData <-cbind(speciesData,creo_climate)
    names(speciesData)[1] <- "speciesCode"
    
    p <- ggplot(speciesData,aes(x=monsoon6SPEI,y=speciesCode^2))+
      geom_point(size=3,color="darkgreen")+
      geom_smooth(formula=y~poly(x,2),method="lm",color="black")+
      xlab("Monsoon SPEI\n")+
      ylab("Maximum abundance") +       ggtitle(speciesCode) +
      theme(legend.position = "none") +
      theme_classic() +
      theme(axis.text.x = element_text(color="black",size=14)) +
      theme(axis.text.y = element_text(color="black",size=14)) +
      theme(axis.title=element_text(size=18))+
      theme(plot.title = element_text(size=20))
    print(p)
    
    plot_name <- paste("C",speciesCodes[i], "plot",number[i],".pdf", sep="_")
    
    ggsave(paste(plot_name), p,
           width=4,height=4,units = c("in"),
           dpi = 600)    
    speciesData$speciesCode <-NULL
  }
  
  if (dAICc_m3_monsoon6mo == 0 | dAICc_m3_monsoon6mo_AR1 == 0 | dAICc_m3_monsoon6mo_AR2 == 0) {
    speciesData <- creo_original[speciesCode]
    speciesData <-cbind(speciesData,creo_climate)
    names(speciesData)[1] <- "speciesCode"
    
    p <- ggplot(speciesData,aes(x=monsoon6SPEI,y=speciesCode^2))+
      geom_point(size=3,color="darkgreen")+
      geom_smooth(formula=y~poly(x,3),method="lm",color="black")+
      xlab("Monsoon SPEI\n")+
      ylab("Maximum abundance") +       ggtitle(speciesCode) +
      theme(legend.position = "none") +
      theme_classic() +
      theme(axis.text.x = element_text(color="black",size=14)) +
      theme(axis.text.y = element_text(color="black",size=14)) +
      theme(axis.title=element_text(size=18))+
      theme(plot.title = element_text(size=20))
    print(p)
    
    plot_name <- paste("C",speciesCodes[i], "plot",number[i],".pdf", sep="_")
    
    ggsave(paste(plot_name), p,
           width=4,height=4,units = c("in"),
           dpi = 600)
    
    speciesData$speciesCode <-NULL
  }
  
  if (dAICc_m1_spring6mo_lag == 0 | dAICc_m1_spring6mo_lag_AR1 == 0 | dAICc_m1_spring6mo_lag_AR2 == 0) {
    speciesData <- creo_original[speciesCode]
    speciesData <-cbind(speciesData,creo_climate)
    names(speciesData)[1] <- "speciesCode"
    
    p <- ggplot(speciesData,aes(x=spring6SPEI_prioryear,y=speciesCode^2))+
      geom_point(size=3,color="darkgreen")+
      geom_smooth(formula=y~poly(x,1),method="lm",color="black")+
      xlab("Spring SPEI \n(previous year)")+
      ylab("Maximum abundance") +       ggtitle(speciesCode) +
      theme(legend.position = "none") +
      theme_classic() +
      theme(axis.text.x = element_text(color="black",size=14)) +
      theme(axis.text.y = element_text(color="black",size=14)) +
      theme(axis.title=element_text(size=18))+
      theme(plot.title = element_text(size=20))
    print(p)
    
    plot_name <- paste("C",speciesCodes[i], "plot",number[i],".pdf", sep="_")
    
    ggsave(paste(plot_name), p,
           width=4,height=4,units = c("in"),
           dpi = 600)
    
    speciesData$speciesCode <-NULL
  }
  
  if (dAICc_m2_spring6mo_lag == 0 | dAICc_m2_spring6mo_lag_AR1 == 0 | dAICc_m2_spring6mo_lag_AR2 == 0) {
    speciesData <- creo_original[speciesCode]
    speciesData <-cbind(speciesData,creo_climate)
    names(speciesData)[1] <- "speciesCode"
    
    p <- ggplot(speciesData,aes(x=spring6SPEI_prioryear,y=speciesCode^2))+
      geom_point(size=3,color="darkgreen")+
      geom_smooth(formula=y~poly(x,2),method="lm",color="black")+
      xlab("Spring SPEI \n(previous year)")+
      ylab("Maximum abundance") +       ggtitle(speciesCode) +
      theme(legend.position = "none") +
      theme_classic() +
      theme(axis.text.x = element_text(color="black",size=14)) +
      theme(axis.text.y = element_text(color="black",size=14)) +
      theme(axis.title=element_text(size=18))+
      theme(plot.title = element_text(size=20))
    print(p)
    
    plot_name <- paste("C",speciesCodes[i], "plot",number[i],".pdf", sep="_")
    
    ggsave(paste(plot_name), p,
           width=4,height=4,units = c("in"),
           dpi = 600)
    
    speciesData$speciesCode <-NULL
  }
  
  if (dAICc_m3_spring6mo_lag == 0 | dAICc_m3_spring6mo_lag_AR1 == 0 | dAICc_m3_spring6mo_lag_AR2 == 0) {
    speciesData <- creo_original[speciesCode]
    speciesData <-cbind(speciesData,creo_climate)
    names(speciesData)[1] <- "speciesCode"
    
    p <- ggplot(speciesData,aes(x=spring6SPEI_prioryear,y=speciesCode^2))+
      geom_point(size=3,color="darkgreen")+
      geom_smooth(formula=y~poly(x,3),method="lm",color="black")+
      xlab("Spring SPEI \n(previous year)")+
      ylab("Maximum abundance") +       ggtitle(speciesCode) +
      theme(legend.position = "none") +
      theme_classic() +
      theme(axis.text.x = element_text(color="black",size=14)) +
      theme(axis.text.y = element_text(color="black",size=14)) +
      theme(axis.title=element_text(size=18))+
      theme(plot.title = element_text(size=20))
    print(p)
    
    plot_name <- paste("C",speciesCodes[i], "plot",number[i],".pdf", sep="_")
    
    ggsave(paste(plot_name), p,
           width=4,height=4,units = c("in"),
           dpi = 600)
    
    speciesData$speciesCode <-NULL
  }
  
  if (dAICc_m1_monsoon6mo_lag == 0 | dAICc_m1_monsoon6mo_lag_AR1 == 0 | dAICc_m1_monsoon6mo_lag_AR2 == 0) {
    speciesData <- creo_original[speciesCode]
    speciesData <-cbind(speciesData,creo_climate)
    names(speciesData)[1] <- "speciesCode"
    
    p <- ggplot(speciesData,aes(x=monsoon6SPEI_prioryear,y=speciesCode^2))+
      geom_point(size=3,color="darkgreen")+
      geom_smooth(formula=y~poly(x,1),method="lm",color="black")+
      xlab("Monsoon SPEI \n(previous year)")+
      ylab("Maximum abundance") +       ggtitle(speciesCode) +
      theme(legend.position = "none") +
      theme_classic() +
      theme(axis.text.x = element_text(color="black",size=14)) +
      theme(axis.text.y = element_text(color="black",size=14)) +
      theme(axis.title=element_text(size=18))+
      theme(plot.title = element_text(size=20))
    print(p)
    
    plot_name <- paste("C",speciesCodes[i], "plot",number[i],".pdf", sep="_")
    
    ggsave(paste(plot_name), p,
           width=4,height=4,units = c("in"),
           dpi = 600)
    
    speciesData$speciesCode <-NULL
  }
  
  if (dAICc_m2_monsoon6mo_lag == 0 | dAICc_m2_monsoon6mo_lag_AR1 == 0 | dAICc_m2_monsoon6mo_lag_AR2 == 0) {
    speciesData <- creo_original[speciesCode]
    speciesData <-cbind(speciesData,creo_climate)
    names(speciesData)[1] <- "speciesCode"
    
    p <- ggplot(speciesData,aes(x=monsoon6SPEI_prioryear,y=speciesCode^2))+
      geom_point(size=3,color="darkgreen")+
      geom_smooth(formula=y~poly(x,2),method="lm",color="black")+
      xlab("Monsoon SPEI \n(previous year)")+
      ylab("Maximum abundance") +       ggtitle(speciesCode) +
      theme(legend.position = "none") +
      theme_classic() +
      theme(axis.text.x = element_text(color="black",size=14)) +
      theme(axis.text.y = element_text(color="black",size=14)) +
      theme(axis.title=element_text(size=18))+
      theme(plot.title = element_text(size=20))
    print(p)
    
    plot_name <- paste("C",speciesCodes[i], "plot",number[i],".pdf", sep="_")
    
    ggsave(paste(plot_name), p,
           width=4,height=4,units = c("in"),
           dpi = 600)
    
    speciesData$speciesCode <-NULL
  }
  
  if (dAICc_m3_monsoon6mo_lag == 0 | dAICc_m3_monsoon6mo_lag_AR1 == 0 | dAICc_m3_monsoon6mo_lag_AR2 == 0) {
    speciesData <- creo_original[speciesCode]
    speciesData <-cbind(speciesData,creo_climate)
    names(speciesData)[1] <- "speciesCode"
    
    p<-ggplot(speciesData,aes(x=monsoon6SPEI_prioryear,y=speciesCode^2))+
      geom_point(size=3,color="darkgreen")+
      geom_smooth(formula=y~poly(x,3),method="lm",color="black")+
      xlab("Monsoon SPEI \n(previous year)")+
      ylab("Maximum abundance") +       ggtitle(speciesCode) +
      theme(legend.position = "none") +
      theme_classic() +
      theme(axis.text.x = element_text(color="black",size=14)) +
      theme(axis.text.y = element_text(color="black",size=14)) +
      theme(axis.title=element_text(size=18))+
      theme(plot.title = element_text(size=20))
    print(p)
    
    plot_name <- paste("C",speciesCodes[i], "plot",number[i],".pdf", sep="_")
    
    ggsave(paste(plot_name), p,
           width=4,height=4,units = c("in"),
           dpi = 600)
    
    speciesData$speciesCode <-NULL
  }
  
  
  # create objects containing statistical values related to each model
  
  numDF_lin_m1_spring6mo<-anova(m1_spring6mo,type="marginal")[2,1]
  denDF_lin_m1_spring6mo<-anova(m1_spring6mo,type="marginal")[2,2]
  F_lin_m1_spring6mo<-anova(m1_spring6mo,type="marginal")[2,3]
  P_lin_m1_spring6mo<-anova(m1_spring6mo,type="marginal")[2,4]
  ParamEst_lin_m1_spring6mo<-coef(summary(m1_spring6mo))[2,1]
  SE_lin_m1_spring6mo<-coef(summary(m1_spring6mo))[2,2]
  
  numDF_lin_m2_spring6mo<-anova(m2_spring6mo,type="marginal")[2,1]
  denDF_lin_m2_spring6mo<-anova(m2_spring6mo,type="marginal")[2,2]
  F_lin_m2_spring6mo<-anova(m2_spring6mo,type="marginal")[2,3]
  P_lin_m2_spring6mo<-anova(m2_spring6mo,type="marginal")[2,4]
  ParamEst_lin_m2_spring6mo<-coef(summary(m2_spring6mo))[2,1]
  SE_lin_m2_spring6mo<-coef(summary(m2_spring6mo))[2,2]
  
  numDF_lin_m3_spring6mo<-anova(m3_spring6mo,type="marginal")[2,1]
  denDF_lin_m3_spring6mo<-anova(m3_spring6mo,type="marginal")[2,2]
  F_lin_m3_spring6mo<-anova(m3_spring6mo,type="marginal")[2,3]
  P_lin_m3_spring6mo<-anova(m3_spring6mo,type="marginal")[2,4]
  ParamEst_lin_m3_spring6mo<-coef(summary(m3_spring6mo))[2,1]
  SE_lin_m3_spring6mo<-coef(summary(m3_spring6mo))[2,2]
  
  numDF_lin_m1_monsoon6mo<-anova(m1_monsoon6mo,type="marginal")[2,1]
  denDF_lin_m1_monsoon6mo<-anova(m1_monsoon6mo,type="marginal")[2,2]
  F_lin_m1_monsoon6mo<-anova(m1_monsoon6mo,type="marginal")[2,3]
  P_lin_m1_monsoon6mo<-anova(m1_monsoon6mo,type="marginal")[2,4]
  ParamEst_lin_m1_monsoon6mo<-coef(summary(m1_monsoon6mo))[2,1]
  SE_lin_m1_monsoon6mo<-coef(summary(m1_monsoon6mo))[2,2]
  
  numDF_lin_m2_monsoon6mo<-anova(m2_monsoon6mo,type="marginal")[2,1]
  denDF_lin_m2_monsoon6mo<-anova(m2_monsoon6mo,type="marginal")[2,2]
  F_lin_m2_monsoon6mo<-anova(m2_monsoon6mo,type="marginal")[2,3]
  P_lin_m2_monsoon6mo<-anova(m2_monsoon6mo,type="marginal")[2,4]
  ParamEst_lin_m2_monsoon6mo<-coef(summary(m2_monsoon6mo))[2,1]
  SE_lin_m2_monsoon6mo<-coef(summary(m2_monsoon6mo))[2,2]
  
  numDF_lin_m3_monsoon6mo<-anova(m3_monsoon6mo,type="marginal")[2,1]
  denDF_lin_m3_monsoon6mo<-anova(m3_monsoon6mo,type="marginal")[2,2]
  F_lin_m3_monsoon6mo<-anova(m3_monsoon6mo,type="marginal")[2,3]
  P_lin_m3_monsoon6mo<-anova(m3_monsoon6mo,type="marginal")[2,4]
  ParamEst_lin_m3_monsoon6mo<-coef(summary(m3_monsoon6mo))[2,1]
  SE_lin_m3_monsoon6mo<-coef(summary(m3_monsoon6mo))[2,2]
  
  numDF_lin_m1_spring6mo_lag<-anova(m1_spring6mo_lag,type="marginal")[2,1]
  denDF_lin_m1_spring6mo_lag<-anova(m1_spring6mo_lag,type="marginal")[2,2]
  F_lin_m1_spring6mo_lag<-anova(m1_spring6mo_lag,type="marginal")[2,3]
  P_lin_m1_spring6mo_lag<-anova(m1_spring6mo_lag,type="marginal")[2,4]
  ParamEst_lin_m1_spring6mo_lag<-coef(summary(m1_spring6mo_lag))[2,1]
  SE_lin_m1_spring6mo_lag<-coef(summary(m1_spring6mo_lag))[2,2]
  
  numDF_lin_m2_spring6mo_lag<-anova(m2_spring6mo_lag,type="marginal")[2,1]
  denDF_lin_m2_spring6mo_lag<-anova(m2_spring6mo_lag,type="marginal")[2,2]
  F_lin_m2_spring6mo_lag<-anova(m2_spring6mo_lag,type="marginal")[2,3]
  P_lin_m2_spring6mo_lag<-anova(m2_spring6mo_lag,type="marginal")[2,4]
  ParamEst_lin_m2_spring6mo_lag<-coef(summary(m2_spring6mo_lag))[2,1]
  SE_lin_m2_spring6mo_lag<-coef(summary(m2_spring6mo_lag))[2,2]
  
  numDF_lin_m3_spring6mo_lag<-anova(m3_spring6mo_lag,type="marginal")[2,1]
  denDF_lin_m3_spring6mo_lag<-anova(m3_spring6mo_lag,type="marginal")[2,2]
  F_lin_m3_spring6mo_lag<-anova(m3_spring6mo_lag,type="marginal")[2,3]
  P_lin_m3_spring6mo_lag<-anova(m3_spring6mo_lag,type="marginal")[2,4]
  ParamEst_lin_m3_spring6mo_lag<-coef(summary(m3_spring6mo_lag))[2,1]
  SE_lin_m3_spring6mo_lag<-coef(summary(m3_spring6mo_lag))[2,2]
  
  numDF_lin_m1_monsoon6mo_lag<-anova(m1_monsoon6mo_lag,type="marginal")[2,1]
  denDF_lin_m1_monsoon6mo_lag<-anova(m1_monsoon6mo_lag,type="marginal")[2,2]
  F_lin_m1_monsoon6mo_lag<-anova(m1_monsoon6mo_lag,type="marginal")[2,3]
  P_lin_m1_monsoon6mo_lag<-anova(m1_monsoon6mo_lag,type="marginal")[2,4]
  ParamEst_lin_m1_monsoon6mo_lag<-coef(summary(m1_monsoon6mo_lag))[2,1]
  SE_lin_m1_monsoon6mo_lag<-coef(summary(m1_monsoon6mo_lag))[2,2]
  
  numDF_lin_m2_monsoon6mo_lag<-anova(m2_monsoon6mo_lag,type="marginal")[2,1]
  denDF_lin_m2_monsoon6mo_lag<-anova(m2_monsoon6mo_lag,type="marginal")[2,2]
  F_lin_m2_monsoon6mo_lag<-anova(m2_monsoon6mo_lag,type="marginal")[2,3]
  P_lin_m2_monsoon6mo_lag<-anova(m2_monsoon6mo_lag,type="marginal")[2,4]
  ParamEst_lin_m2_monsoon6mo_lag<-coef(summary(m2_monsoon6mo_lag))[2,1]
  SE_lin_m2_monsoon6mo_lag<-coef(summary(m2_monsoon6mo_lag))[2,2]
  
  numDF_lin_m3_monsoon6mo_lag<-anova(m3_monsoon6mo_lag,type="marginal")[2,1]
  denDF_lin_m3_monsoon6mo_lag<-anova(m3_monsoon6mo_lag,type="marginal")[2,2]
  F_lin_m3_monsoon6mo_lag<-anova(m3_monsoon6mo_lag,type="marginal")[2,3]
  P_lin_m3_monsoon6mo_lag<-anova(m3_monsoon6mo_lag,type="marginal")[2,4]
  ParamEst_lin_m3_monsoon6mo_lag<-coef(summary(m3_monsoon6mo_lag))[2,1]
  SE_lin_m3_monsoon6mo_lag<-coef(summary(m3_monsoon6mo_lag))[2,2]
  
  numDF_lin_m1_spring6mo_AR1<-anova(m1_spring6mo_AR1,type="marginal")[2,1]
  denDF_lin_m1_spring6mo_AR1<-anova(m1_spring6mo_AR1,type="marginal")[2,2]
  F_lin_m1_spring6mo_AR1<-anova(m1_spring6mo_AR1,type="marginal")[2,3]
  P_lin_m1_spring6mo_AR1<-anova(m1_spring6mo_AR1,type="marginal")[2,4]
  ParamEst_lin_m1_spring6mo_AR1<-coef(summary(m1_spring6mo_AR1))[2,1]
  SE_lin_m1_spring6mo_AR1<-coef(summary(m1_spring6mo_AR1))[2,2]
  
  numDF_lin_m1_spring6mo_AR2<-anova(m1_spring6mo_AR2,type="marginal")[2,1]
  denDF_lin_m1_spring6mo_AR2<-anova(m1_spring6mo_AR2,type="marginal")[2,2]
  F_lin_m1_spring6mo_AR2<-anova(m1_spring6mo_AR2,type="marginal")[2,3]
  P_lin_m1_spring6mo_AR2<-anova(m1_spring6mo_AR2,type="marginal")[2,4]
  ParamEst_lin_m1_spring6mo_AR2<-coef(summary(m1_spring6mo_AR2))[2,1]
  SE_lin_m1_spring6mo_AR2<-coef(summary(m1_spring6mo_AR2))[2,2]
  
  numDF_lin_m2_spring6mo_AR1<-anova(m2_spring6mo_AR1,type="marginal")[2,1]
  denDF_lin_m2_spring6mo_AR1<-anova(m2_spring6mo_AR1,type="marginal")[2,2]
  F_lin_m2_spring6mo_AR1<-anova(m2_spring6mo_AR1,type="marginal")[2,3]
  P_lin_m2_spring6mo_AR1<-anova(m2_spring6mo_AR1,type="marginal")[2,4]
  ParamEst_lin_m2_spring6mo_AR1<-coef(summary(m2_spring6mo_AR1))[2,1]
  SE_lin_m2_spring6mo_AR1<-coef(summary(m2_spring6mo_AR1))[2,2]
  
  numDF_lin_m2_spring6mo_AR2<-anova(m2_spring6mo_AR2,type="marginal")[2,1]
  denDF_lin_m2_spring6mo_AR2<-anova(m2_spring6mo_AR2,type="marginal")[2,2]
  F_lin_m2_spring6mo_AR2<-anova(m2_spring6mo_AR2,type="marginal")[2,3]
  P_lin_m2_spring6mo_AR2<-anova(m2_spring6mo_AR2,type="marginal")[2,4]
  ParamEst_lin_m2_spring6mo_AR2<-coef(summary(m2_spring6mo_AR2))[2,1]
  SE_lin_m2_spring6mo_AR2<-coef(summary(m2_spring6mo_AR2))[2,2]
  
  numDF_lin_m3_spring6mo_AR1<-anova(m3_spring6mo_AR1,type="marginal")[2,1]
  denDF_lin_m3_spring6mo_AR1<-anova(m3_spring6mo_AR1,type="marginal")[2,2]
  F_lin_m3_spring6mo_AR1<-anova(m3_spring6mo_AR1,type="marginal")[2,3]
  P_lin_m3_spring6mo_AR1<-anova(m3_spring6mo_AR1,type="marginal")[2,4]
  ParamEst_lin_m3_spring6mo_AR1<-coef(summary(m3_spring6mo_AR1))[2,1]
  SE_lin_m3_spring6mo_AR1<-coef(summary(m3_spring6mo_AR1))[2,2]
  
  numDF_lin_m3_spring6mo_AR2<-anova(m3_spring6mo_AR2,type="marginal")[2,1]
  denDF_lin_m3_spring6mo_AR2<-anova(m3_spring6mo_AR2,type="marginal")[2,2]
  F_lin_m3_spring6mo_AR2<-anova(m3_spring6mo_AR2,type="marginal")[2,3]
  P_lin_m3_spring6mo_AR2<-anova(m3_spring6mo_AR2,type="marginal")[2,4]
  ParamEst_lin_m3_spring6mo_AR2<-coef(summary(m3_spring6mo_AR2))[2,1]
  SE_lin_m3_spring6mo_AR2<-coef(summary(m3_spring6mo_AR2))[2,2]
  
  numDF_lin_m1_monsoon6mo_AR1<-anova(m1_monsoon6mo_AR1,type="marginal")[2,1]
  denDF_lin_m1_monsoon6mo_AR1<-anova(m1_monsoon6mo_AR1,type="marginal")[2,2]
  F_lin_m1_monsoon6mo_AR1<-anova(m1_monsoon6mo_AR1,type="marginal")[2,3]
  P_lin_m1_monsoon6mo_AR1<-anova(m1_monsoon6mo_AR1,type="marginal")[2,4]
  ParamEst_lin_m1_monsoon6mo_AR1<-coef(summary(m1_monsoon6mo_AR1))[2,1]
  SE_lin_m1_monsoon6mo_AR1<-coef(summary(m1_monsoon6mo_AR1))[2,2]
  
  numDF_lin_m1_monsoon6mo_AR2<-anova(m1_monsoon6mo_AR2,type="marginal")[2,1]
  denDF_lin_m1_monsoon6mo_AR2<-anova(m1_monsoon6mo_AR2,type="marginal")[2,2]
  F_lin_m1_monsoon6mo_AR2<-anova(m1_monsoon6mo_AR2,type="marginal")[2,3]
  P_lin_m1_monsoon6mo_AR2<-anova(m1_monsoon6mo_AR2,type="marginal")[2,4]
  ParamEst_lin_m1_monsoon6mo_AR2<-coef(summary(m1_monsoon6mo_AR2))[2,1]
  SE_lin_m1_monsoon6mo_AR2<-coef(summary(m1_monsoon6mo_AR2))[2,2]
  
  numDF_lin_m2_monsoon6mo_AR1<-anova(m2_monsoon6mo_AR1,type="marginal")[2,1]
  denDF_lin_m2_monsoon6mo_AR1<-anova(m2_monsoon6mo_AR1,type="marginal")[2,2]
  F_lin_m2_monsoon6mo_AR1<-anova(m2_monsoon6mo_AR1,type="marginal")[2,3]
  P_lin_m2_monsoon6mo_AR1<-anova(m2_monsoon6mo_AR1,type="marginal")[2,4]
  ParamEst_lin_m2_monsoon6mo_AR1<-coef(summary(m2_monsoon6mo_AR1))[2,1]
  SE_lin_m2_monsoon6mo_AR1<-coef(summary(m2_monsoon6mo_AR1))[2,2]
  
  numDF_lin_m2_monsoon6mo_AR2<-anova(m2_monsoon6mo_AR2,type="marginal")[2,1]
  denDF_lin_m2_monsoon6mo_AR2<-anova(m2_monsoon6mo_AR2,type="marginal")[2,2]
  F_lin_m2_monsoon6mo_AR2<-anova(m2_monsoon6mo_AR2,type="marginal")[2,3]
  P_lin_m2_monsoon6mo_AR2<-anova(m2_monsoon6mo_AR2,type="marginal")[2,4]
  ParamEst_lin_m2_monsoon6mo_AR2<-coef(summary(m2_monsoon6mo_AR2))[2,1]
  SE_lin_m2_monsoon6mo_AR2<-coef(summary(m2_monsoon6mo_AR2))[2,2]
  
  numDF_lin_m3_monsoon6mo_AR1<-anova(m3_monsoon6mo_AR1,type="marginal")[2,1]
  denDF_lin_m3_monsoon6mo_AR1<-anova(m3_monsoon6mo_AR1,type="marginal")[2,2]
  F_lin_m3_monsoon6mo_AR1<-anova(m3_monsoon6mo_AR1,type="marginal")[2,3]
  P_lin_m3_monsoon6mo_AR1<-anova(m3_monsoon6mo_AR1,type="marginal")[2,4]
  ParamEst_lin_m3_monsoon6mo_AR1<-coef(summary(m3_monsoon6mo_AR1))[2,1]
  SE_lin_m3_monsoon6mo_AR1<-coef(summary(m3_monsoon6mo_AR1))[2,2]
  
  numDF_lin_m3_monsoon6mo_AR2<-anova(m3_monsoon6mo_AR2,type="marginal")[2,1]
  denDF_lin_m3_monsoon6mo_AR2<-anova(m3_monsoon6mo_AR2,type="marginal")[2,2]
  F_lin_m3_monsoon6mo_AR2<-anova(m3_monsoon6mo_AR2,type="marginal")[2,3]
  P_lin_m3_monsoon6mo_AR2<-anova(m3_monsoon6mo_AR2,type="marginal")[2,4]
  ParamEst_lin_m3_monsoon6mo_AR2<-coef(summary(m3_monsoon6mo_AR2))[2,1]
  SE_lin_m3_monsoon6mo_AR2<-coef(summary(m3_monsoon6mo_AR2))[2,2]
  
  numDF_lin_m1_spring6mo_lag_AR1<-anova(m1_spring6mo_lag_AR1,type="marginal")[2,1]
  denDF_lin_m1_spring6mo_lag_AR1<-anova(m1_spring6mo_lag_AR1,type="marginal")[2,2]
  F_lin_m1_spring6mo_lag_AR1<-anova(m1_spring6mo_lag_AR1,type="marginal")[2,3]
  P_lin_m1_spring6mo_lag_AR1<-anova(m1_spring6mo_lag_AR1,type="marginal")[2,4]
  ParamEst_lin_m1_spring6mo_lag_AR1<-coef(summary(m1_spring6mo_lag_AR1))[2,1]
  SE_lin_m1_spring6mo_lag_AR1<-coef(summary(m1_spring6mo_lag_AR1))[2,2]
  
  numDF_lin_m1_spring6mo_lag_AR2<-anova(m1_spring6mo_lag_AR2,type="marginal")[2,1]
  denDF_lin_m1_spring6mo_lag_AR2<-anova(m1_spring6mo_lag_AR2,type="marginal")[2,2]
  F_lin_m1_spring6mo_lag_AR2<-anova(m1_spring6mo_lag_AR2,type="marginal")[2,3]
  P_lin_m1_spring6mo_lag_AR2<-anova(m1_spring6mo_lag_AR2,type="marginal")[2,4]
  ParamEst_lin_m1_spring6mo_lag_AR2<-coef(summary(m1_spring6mo_lag_AR2))[2,1]
  SE_lin_m1_spring6mo_lag_AR2<-coef(summary(m1_spring6mo_lag_AR2))[2,2]
  
  numDF_lin_m2_spring6mo_lag_AR1<-anova(m2_spring6mo_lag_AR1,type="marginal")[2,1]
  denDF_lin_m2_spring6mo_lag_AR1<-anova(m2_spring6mo_lag_AR1,type="marginal")[2,2]
  F_lin_m2_spring6mo_lag_AR1<-anova(m2_spring6mo_lag_AR1,type="marginal")[2,3]
  P_lin_m2_spring6mo_lag_AR1<-anova(m2_spring6mo_lag_AR1,type="marginal")[2,4]
  ParamEst_lin_m2_spring6mo_lag_AR1<-coef(summary(m2_spring6mo_lag_AR1))[2,1]
  SE_lin_m2_spring6mo_lag_AR1<-coef(summary(m2_spring6mo_lag_AR1))[2,2]
  
  numDF_lin_m2_spring6mo_lag_AR2<-anova(m2_spring6mo_lag_AR2,type="marginal")[2,1]
  denDF_lin_m2_spring6mo_lag_AR2<-anova(m2_spring6mo_lag_AR2,type="marginal")[2,2]
  F_lin_m2_spring6mo_lag_AR2<-anova(m2_spring6mo_lag_AR2,type="marginal")[2,3]
  P_lin_m2_spring6mo_lag_AR2<-anova(m2_spring6mo_lag_AR2,type="marginal")[2,4]
  ParamEst_lin_m2_spring6mo_lag_AR2<-coef(summary(m2_spring6mo_lag_AR2))[2,1]
  SE_lin_m2_spring6mo_lag_AR2<-coef(summary(m2_spring6mo_lag_AR2))[2,2]
  
  numDF_lin_m3_spring6mo_lag_AR1<-anova(m3_spring6mo_lag_AR1,type="marginal")[2,1]
  denDF_lin_m3_spring6mo_lag_AR1<-anova(m3_spring6mo_lag_AR1,type="marginal")[2,2]
  F_lin_m3_spring6mo_lag_AR1<-anova(m3_spring6mo_lag_AR1,type="marginal")[2,3]
  P_lin_m3_spring6mo_lag_AR1<-anova(m3_spring6mo_lag_AR1,type="marginal")[2,4]
  ParamEst_lin_m3_spring6mo_lag_AR1<-coef(summary(m3_spring6mo_lag_AR1))[2,1]
  SE_lin_m3_spring6mo_lag_AR1<-coef(summary(m3_spring6mo_lag_AR1))[2,2]
  
  numDF_lin_m3_spring6mo_lag_AR2<-anova(m3_spring6mo_lag_AR2,type="marginal")[2,1]
  denDF_lin_m3_spring6mo_lag_AR2<-anova(m3_spring6mo_lag_AR2,type="marginal")[2,2]
  F_lin_m3_spring6mo_lag_AR2<-anova(m3_spring6mo_lag_AR2,type="marginal")[2,3]
  P_lin_m3_spring6mo_lag_AR2<-anova(m3_spring6mo_lag_AR2,type="marginal")[2,4]
  ParamEst_lin_m3_spring6mo_lag_AR2<-coef(summary(m3_spring6mo_lag_AR2))[2,1]
  SE_lin_m3_spring6mo_lag_AR2<-coef(summary(m3_spring6mo_lag_AR2))[2,2]
  
  numDF_lin_m1_monsoon6mo_lag_AR1<-anova(m1_monsoon6mo_lag_AR1,type="marginal")[2,1]
  denDF_lin_m1_monsoon6mo_lag_AR1<-anova(m1_monsoon6mo_lag_AR1,type="marginal")[2,2]
  F_lin_m1_monsoon6mo_lag_AR1<-anova(m1_monsoon6mo_lag_AR1,type="marginal")[2,3]
  P_lin_m1_monsoon6mo_lag_AR1<-anova(m1_monsoon6mo_lag_AR1,type="marginal")[2,4]
  ParamEst_lin_m1_monsoon6mo_lag_AR1<-coef(summary(m1_monsoon6mo_lag_AR1))[2,1]
  SE_lin_m1_monsoon6mo_lag_AR1<-coef(summary(m1_monsoon6mo_lag_AR1))[2,2]
  
  numDF_lin_m1_monsoon6mo_lag_AR2<-anova(m1_monsoon6mo_lag_AR2,type="marginal")[2,1]
  denDF_lin_m1_monsoon6mo_lag_AR2<-anova(m1_monsoon6mo_lag_AR2,type="marginal")[2,2]
  F_lin_m1_monsoon6mo_lag_AR2<-anova(m1_monsoon6mo_lag_AR2,type="marginal")[2,3]
  P_lin_m1_monsoon6mo_lag_AR2<-anova(m1_monsoon6mo_lag_AR2,type="marginal")[2,4]
  ParamEst_lin_m1_monsoon6mo_lag_AR2<-coef(summary(m1_monsoon6mo_lag_AR2))[2,1]
  SE_lin_m1_monsoon6mo_lag_AR2<-coef(summary(m1_monsoon6mo_lag_AR2))[2,2]
  
  numDF_lin_m2_monsoon6mo_lag_AR1<-anova(m2_monsoon6mo_lag_AR1,type="marginal")[2,1]
  denDF_lin_m2_monsoon6mo_lag_AR1<-anova(m2_monsoon6mo_lag_AR1,type="marginal")[2,2]
  F_lin_m2_monsoon6mo_lag_AR1<-anova(m2_monsoon6mo_lag_AR1,type="marginal")[2,3]
  P_lin_m2_monsoon6mo_lag_AR1<-anova(m2_monsoon6mo_lag_AR1,type="marginal")[2,4]
  ParamEst_lin_m2_monsoon6mo_lag_AR1<-coef(summary(m2_monsoon6mo_lag_AR1))[2,1]
  SE_lin_m2_monsoon6mo_lag_AR1<-coef(summary(m2_monsoon6mo_lag_AR1))[2,2]
  
  numDF_lin_m2_monsoon6mo_lag_AR2<-anova(m2_monsoon6mo_lag_AR2,type="marginal")[2,1]
  denDF_lin_m2_monsoon6mo_lag_AR2<-anova(m2_monsoon6mo_lag_AR2,type="marginal")[2,2]
  F_lin_m2_monsoon6mo_lag_AR2<-anova(m2_monsoon6mo_lag_AR2,type="marginal")[2,3]
  P_lin_m2_monsoon6mo_lag_AR2<-anova(m2_monsoon6mo_lag_AR2,type="marginal")[2,4]
  ParamEst_lin_m2_monsoon6mo_lag_AR2<-coef(summary(m2_monsoon6mo_lag_AR2))[2,1]
  SE_lin_m2_monsoon6mo_lag_AR2<-coef(summary(m2_monsoon6mo_lag_AR2))[2,2]
  
  numDF_lin_m3_monsoon6mo_lag_AR1<-anova(m3_monsoon6mo_lag_AR1,type="marginal")[2,1]
  denDF_lin_m3_monsoon6mo_lag_AR1<-anova(m3_monsoon6mo_lag_AR1,type="marginal")[2,2]
  F_lin_m3_monsoon6mo_lag_AR1<-anova(m3_monsoon6mo_lag_AR1,type="marginal")[2,3]
  P_lin_m3_monsoon6mo_lag_AR1<-anova(m3_monsoon6mo_lag_AR1,type="marginal")[2,4]
  ParamEst_lin_m3_monsoon6mo_lag_AR1<-coef(summary(m3_monsoon6mo_lag_AR1))[2,1]
  SE_lin_m3_monsoon6mo_lag_AR1<-coef(summary(m3_monsoon6mo_lag_AR1))[2,2]
  
  numDF_lin_m3_monsoon6mo_lag_AR2<-anova(m3_monsoon6mo_lag_AR2,type="marginal")[2,1]
  denDF_lin_m3_monsoon6mo_lag_AR2<-anova(m3_monsoon6mo_lag_AR2,type="marginal")[2,2]
  F_lin_m3_monsoon6mo_lag_AR2<-anova(m3_monsoon6mo_lag_AR2,type="marginal")[2,3]
  P_lin_m3_monsoon6mo_lag_AR2<-anova(m3_monsoon6mo_lag_AR2,type="marginal")[2,4]
  ParamEst_lin_m3_monsoon6mo_lag_AR2<-coef(summary(m3_monsoon6mo_lag_AR2))[2,1]
  SE_lin_m3_monsoon6mo_lag_AR2<-coef(summary(m3_monsoon6mo_lag_AR2))[2,2]
  
  Rsquared_marginal_m_null<-rsquared(m_null)[1,5]
  Rsquared_marginal_m1_spring6mo<-rsquared(m1_spring6mo)[1,5]
  Rsquared_marginal_m2_spring6mo<-rsquared(m2_spring6mo)[1,5]
  Rsquared_marginal_m3_spring6mo<-rsquared(m3_spring6mo)[1,5]
  Rsquared_marginal_m1_spring6mo_AR1<-rsquared(m1_spring6mo_AR1)[1,5]
  Rsquared_marginal_m1_spring6mo_AR2<-rsquared(m1_spring6mo_AR2)[1,5]
  Rsquared_marginal_m2_spring6mo_AR1<-rsquared(m2_spring6mo_AR1)[1,5]
  Rsquared_marginal_m2_spring6mo_AR2<-rsquared(m2_spring6mo_AR2)[1,5]
  Rsquared_marginal_m3_spring6mo_AR1<-rsquared(m3_spring6mo_AR1)[1,5]
  Rsquared_marginal_m3_spring6mo_AR2<-rsquared(m3_spring6mo_AR2)[1,5]
  Rsquared_marginal_m1_monsoon6mo<-rsquared(m1_monsoon6mo)[1,5]
  Rsquared_marginal_m2_monsoon6mo<-rsquared(m2_monsoon6mo)[1,5]
  Rsquared_marginal_m3_monsoon6mo<-rsquared(m3_monsoon6mo)[1,5]
  Rsquared_marginal_m1_monsoon6mo_AR1<-rsquared(m1_monsoon6mo_AR1)[1,5]
  Rsquared_marginal_m1_monsoon6mo_AR2<-rsquared(m1_monsoon6mo_AR2)[1,5]
  Rsquared_marginal_m2_monsoon6mo_AR1<-rsquared(m2_monsoon6mo_AR1)[1,5]
  Rsquared_marginal_m2_monsoon6mo_AR2<-rsquared(m2_monsoon6mo_AR2)[1,5]
  Rsquared_marginal_m3_monsoon6mo_AR1<-rsquared(m3_monsoon6mo_AR1)[1,5]
  Rsquared_marginal_m3_monsoon6mo_AR2<-rsquared(m3_monsoon6mo_AR2)[1,5]
  Rsquared_marginal_m1_spring6mo_lag<-rsquared(m1_spring6mo_lag)[1,5]
  Rsquared_marginal_m2_spring6mo_lag<-rsquared(m2_spring6mo_lag)[1,5]
  Rsquared_marginal_m3_spring6mo_lag<-rsquared(m3_spring6mo_lag)[1,5]
  Rsquared_marginal_m1_spring6mo_lag_AR1<-rsquared(m1_spring6mo_lag_AR1)[1,5]
  Rsquared_marginal_m1_spring6mo_lag_AR2<-rsquared(m1_spring6mo_lag_AR2)[1,5]
  Rsquared_marginal_m2_spring6mo_lag_AR1<-rsquared(m2_spring6mo_lag_AR1)[1,5]
  Rsquared_marginal_m2_spring6mo_lag_AR2<-rsquared(m2_spring6mo_lag_AR2)[1,5]
  Rsquared_marginal_m3_spring6mo_lag_AR1<-rsquared(m3_spring6mo_lag_AR1)[1,5]
  Rsquared_marginal_m3_spring6mo_lag_AR2<-rsquared(m3_spring6mo_lag_AR2)[1,5]
  Rsquared_marginal_m1_monsoon6mo_lag<-rsquared(m1_monsoon6mo_lag)[1,5]
  Rsquared_marginal_m2_monsoon6mo_lag<-rsquared(m2_monsoon6mo_lag)[1,5]
  Rsquared_marginal_m3_monsoon6mo_lag<-rsquared(m3_monsoon6mo_lag)[1,5]
  Rsquared_marginal_m1_monsoon6mo_lag_AR1<-rsquared(m1_monsoon6mo_lag_AR1)[1,5]
  Rsquared_marginal_m1_monsoon6mo_lag_AR2<-rsquared(m1_monsoon6mo_lag_AR2)[1,5]
  Rsquared_marginal_m2_monsoon6mo_lag_AR1<-rsquared(m2_monsoon6mo_lag_AR1)[1,5]
  Rsquared_marginal_m2_monsoon6mo_lag_AR2<-rsquared(m2_monsoon6mo_lag_AR2)[1,5]
  Rsquared_marginal_m3_monsoon6mo_lag_AR1<-rsquared(m3_monsoon6mo_lag_AR1)[1,5]
  Rsquared_marginal_m3_monsoon6mo_lag_AR2<-rsquared(m3_monsoon6mo_lag_AR2)[1,5]
  
  Rsquared_conditional_m_null<-rsquared(m_null)[1,6]
  Rsquared_conditional_m1_spring6mo<-rsquared(m1_spring6mo)[1,6]
  Rsquared_conditional_m2_spring6mo<-rsquared(m2_spring6mo)[1,6]
  Rsquared_conditional_m3_spring6mo<-rsquared(m3_spring6mo)[1,6]
  Rsquared_conditional_m1_spring6mo_AR1<-rsquared(m1_spring6mo_AR1)[1,6]
  Rsquared_conditional_m1_spring6mo_AR2<-rsquared(m1_spring6mo_AR2)[1,6]
  Rsquared_conditional_m2_spring6mo_AR1<-rsquared(m2_spring6mo_AR1)[1,6]
  Rsquared_conditional_m2_spring6mo_AR2<-rsquared(m2_spring6mo_AR2)[1,6]
  Rsquared_conditional_m3_spring6mo_AR1<-rsquared(m3_spring6mo_AR1)[1,6]
  Rsquared_conditional_m3_spring6mo_AR2<-rsquared(m3_spring6mo_AR2)[1,6]
  Rsquared_conditional_m1_monsoon6mo<-rsquared(m1_monsoon6mo)[1,6]
  Rsquared_conditional_m2_monsoon6mo<-rsquared(m2_monsoon6mo)[1,6]
  Rsquared_conditional_m3_monsoon6mo<-rsquared(m3_monsoon6mo)[1,6]
  Rsquared_conditional_m1_monsoon6mo_AR1<-rsquared(m1_monsoon6mo_AR1)[1,6]
  Rsquared_conditional_m1_monsoon6mo_AR2<-rsquared(m1_monsoon6mo_AR2)[1,6]
  Rsquared_conditional_m2_monsoon6mo_AR1<-rsquared(m2_monsoon6mo_AR1)[1,6]
  Rsquared_conditional_m2_monsoon6mo_AR2<-rsquared(m2_monsoon6mo_AR2)[1,6]
  Rsquared_conditional_m3_monsoon6mo_AR1<-rsquared(m3_monsoon6mo_AR1)[1,6]
  Rsquared_conditional_m3_monsoon6mo_AR2<-rsquared(m3_monsoon6mo_AR2)[1,6]
  Rsquared_conditional_m1_spring6mo_lag<-rsquared(m1_spring6mo_lag)[1,6]
  Rsquared_conditional_m2_spring6mo_lag<-rsquared(m2_spring6mo_lag)[1,6]
  Rsquared_conditional_m3_spring6mo_lag<-rsquared(m3_spring6mo_lag)[1,6]
  Rsquared_conditional_m1_spring6mo_lag_AR1<-rsquared(m1_spring6mo_lag_AR1)[1,6]
  Rsquared_conditional_m1_spring6mo_lag_AR2<-rsquared(m1_spring6mo_lag_AR2)[1,6]
  Rsquared_conditional_m2_spring6mo_lag_AR1<-rsquared(m2_spring6mo_lag_AR1)[1,6]
  Rsquared_conditional_m2_spring6mo_lag_AR2<-rsquared(m2_spring6mo_lag_AR2)[1,6]
  Rsquared_conditional_m3_spring6mo_lag_AR1<-rsquared(m3_spring6mo_lag_AR1)[1,6]
  Rsquared_conditional_m3_spring6mo_lag_AR2<-rsquared(m3_spring6mo_lag_AR2)[1,6]
  Rsquared_conditional_m1_monsoon6mo_lag<-rsquared(m1_monsoon6mo_lag)[1,6]
  Rsquared_conditional_m2_monsoon6mo_lag<-rsquared(m2_monsoon6mo_lag)[1,6]
  Rsquared_conditional_m3_monsoon6mo_lag<-rsquared(m3_monsoon6mo_lag)[1,6]
  Rsquared_conditional_m1_monsoon6mo_lag_AR1<-rsquared(m1_monsoon6mo_lag_AR1)[1,6]
  Rsquared_conditional_m1_monsoon6mo_lag_AR2<-rsquared(m1_monsoon6mo_lag_AR2)[1,6]
  Rsquared_conditional_m2_monsoon6mo_lag_AR1<-rsquared(m2_monsoon6mo_lag_AR1)[1,6]
  Rsquared_conditional_m2_monsoon6mo_lag_AR2<-rsquared(m2_monsoon6mo_lag_AR2)[1,6]
  Rsquared_conditional_m3_monsoon6mo_lag_AR1<-rsquared(m3_monsoon6mo_lag_AR1)[1,6]
  Rsquared_conditional_m3_monsoon6mo_lag_AR2<-rsquared(m3_monsoon6mo_lag_AR2)[1,6]
  
  numDF_quad_m2_spring6mo<-anova(m2_spring6mo,type="marginal")[3,1]
  numDF_quad_m2_spring6mo_AR1<-anova(m2_spring6mo_AR1,type="marginal")[3,1]
  numDF_quad_m2_spring6mo_AR2<-anova(m2_spring6mo_AR2,type="marginal")[3,1]
  numDF_quad_m2_monsoon6mo<-anova(m2_monsoon6mo,type="marginal")[3,1]
  numDF_quad_m2_monsoon6mo_AR1<-anova(m2_monsoon6mo_AR1,type="marginal")[3,1]
  numDF_quad_m2_monsoon6mo_AR2<-anova(m2_monsoon6mo_AR2,type="marginal")[3,1]
  numDF_quad_m2_spring6mo_lag<-anova(m2_spring6mo_lag,type="marginal")[3,1]
  numDF_quad_m2_spring6mo_lag_AR1<-anova(m2_spring6mo_lag_AR1,type="marginal")[3,1]
  numDF_quad_m2_spring6mo_lag_AR2<-anova(m2_spring6mo_lag_AR2,type="marginal")[3,1]
  numDF_quad_m2_monsoon6mo_lag<-anova(m2_monsoon6mo_lag,type="marginal")[3,1]
  numDF_quad_m2_monsoon6mo_lag_AR1<-anova(m2_monsoon6mo_lag_AR1,type="marginal")[3,1]
  numDF_quad_m2_monsoon6mo_lag_AR2<-anova(m2_monsoon6mo_lag_AR2,type="marginal")[3,1]
  
  numDF_quad_m3_spring6mo<-anova(m3_spring6mo,type="marginal")[3,1]
  numDF_quad_m3_spring6mo_AR1<-anova(m3_spring6mo_AR1,type="marginal")[3,1]
  numDF_quad_m3_spring6mo_AR2<-anova(m3_spring6mo_AR2,type="marginal")[3,1]
  numDF_quad_m3_monsoon6mo<-anova(m3_monsoon6mo,type="marginal")[3,1]
  numDF_quad_m3_monsoon6mo_AR1<-anova(m3_monsoon6mo_AR1,type="marginal")[3,1]
  numDF_quad_m3_monsoon6mo_AR2<-anova(m3_monsoon6mo_AR2,type="marginal")[3,1]
  numDF_quad_m3_spring6mo_lag<-anova(m3_spring6mo_lag,type="marginal")[3,1]
  numDF_quad_m3_spring6mo_lag_AR1<-anova(m3_spring6mo_lag_AR1,type="marginal")[3,1]
  numDF_quad_m3_spring6mo_lag_AR2<-anova(m3_spring6mo_lag_AR2,type="marginal")[3,1]
  numDF_quad_m3_monsoon6mo_lag<-anova(m3_monsoon6mo_lag,type="marginal")[3,1]
  numDF_quad_m3_monsoon6mo_lag_AR1<-anova(m3_monsoon6mo_lag_AR1,type="marginal")[3,1]
  numDF_quad_m3_monsoon6mo_lag_AR2<-anova(m3_monsoon6mo_lag_AR2,type="marginal")[3,1]
  
  denDF_quad_m2_spring6mo<-anova(m2_spring6mo,type="marginal")[3,2]
  denDF_quad_m2_spring6mo_AR1<-anova(m2_spring6mo_AR1,type="marginal")[3,2]
  denDF_quad_m2_spring6mo_AR2<-anova(m2_spring6mo_AR2,type="marginal")[3,2]
  denDF_quad_m2_monsoon6mo<-anova(m2_monsoon6mo,type="marginal")[3,2]
  denDF_quad_m2_monsoon6mo_AR1<-anova(m2_monsoon6mo_AR1,type="marginal")[3,2]
  denDF_quad_m2_monsoon6mo_AR2<-anova(m2_monsoon6mo_AR2,type="marginal")[3,2]
  denDF_quad_m2_spring6mo_lag<-anova(m2_spring6mo_lag,type="marginal")[3,2]
  denDF_quad_m2_spring6mo_lag_AR1<-anova(m2_spring6mo_lag_AR1,type="marginal")[3,2]
  denDF_quad_m2_spring6mo_lag_AR2<-anova(m2_spring6mo_lag_AR2,type="marginal")[3,2]
  denDF_quad_m2_monsoon6mo_lag<-anova(m2_monsoon6mo_lag,type="marginal")[3,2]
  denDF_quad_m2_monsoon6mo_lag_AR1<-anova(m2_monsoon6mo_lag_AR1,type="marginal")[3,2]
  denDF_quad_m2_monsoon6mo_lag_AR2<-anova(m2_monsoon6mo_lag_AR2,type="marginal")[3,2]
  
  denDF_quad_m3_spring6mo<-anova(m3_spring6mo,type="marginal")[3,2]
  denDF_quad_m3_spring6mo_AR1<-anova(m3_spring6mo_AR1,type="marginal")[3,2]
  denDF_quad_m3_spring6mo_AR2<-anova(m3_spring6mo_AR2,type="marginal")[3,2]
  denDF_quad_m3_monsoon6mo<-anova(m3_monsoon6mo,type="marginal")[3,2]
  denDF_quad_m3_monsoon6mo_AR1<-anova(m3_monsoon6mo_AR1,type="marginal")[3,2]
  denDF_quad_m3_monsoon6mo_AR2<-anova(m3_monsoon6mo_AR2,type="marginal")[3,2]
  denDF_quad_m3_spring6mo_lag<-anova(m3_spring6mo_lag,type="marginal")[3,2]
  denDF_quad_m3_spring6mo_lag_AR1<-anova(m3_spring6mo_lag_AR1,type="marginal")[3,2]
  denDF_quad_m3_spring6mo_lag_AR2<-anova(m3_spring6mo_lag_AR2,type="marginal")[3,2]
  denDF_quad_m3_monsoon6mo_lag<-anova(m3_monsoon6mo_lag,type="marginal")[3,2]
  denDF_quad_m3_monsoon6mo_lag_AR1<-anova(m3_monsoon6mo_lag_AR1,type="marginal")[3,2]
  denDF_quad_m3_monsoon6mo_lag_AR2<-anova(m3_monsoon6mo_lag_AR2,type="marginal")[3,2]
  
  F_quad_m2_spring6mo<-anova(m2_spring6mo,type="marginal")[3,3]
  F_quad_m2_spring6mo_AR1<-anova(m2_spring6mo_AR1,type="marginal")[3,3]
  F_quad_m2_spring6mo_AR2<-anova(m2_spring6mo_AR2,type="marginal")[3,3]
  F_quad_m2_monsoon6mo<-anova(m2_monsoon6mo,type="marginal")[3,3]
  F_quad_m2_monsoon6mo_AR1<-anova(m2_monsoon6mo_AR1,type="marginal")[3,3]
  F_quad_m2_monsoon6mo_AR2<-anova(m2_monsoon6mo_AR2,type="marginal")[3,3]
  F_quad_m2_spring6mo_lag<-anova(m2_spring6mo_lag,type="marginal")[3,3]
  F_quad_m2_spring6mo_lag_AR1<-anova(m2_spring6mo_lag_AR1,type="marginal")[3,3]
  F_quad_m2_spring6mo_lag_AR2<-anova(m2_spring6mo_lag_AR2,type="marginal")[3,3]
  F_quad_m2_monsoon6mo_lag<-anova(m2_monsoon6mo_lag,type="marginal")[3,3]
  F_quad_m2_monsoon6mo_lag_AR1<-anova(m2_monsoon6mo_lag_AR1,type="marginal")[3,3]
  F_quad_m2_monsoon6mo_lag_AR2<-anova(m2_monsoon6mo_lag_AR2,type="marginal")[3,3]
  
  F_quad_m3_spring6mo<-anova(m3_spring6mo,type="marginal")[3,3]
  F_quad_m3_spring6mo_AR1<-anova(m3_spring6mo_AR1,type="marginal")[3,3]
  F_quad_m3_spring6mo_AR2<-anova(m3_spring6mo_AR2,type="marginal")[3,3]
  F_quad_m3_monsoon6mo<-anova(m3_monsoon6mo,type="marginal")[3,3]
  F_quad_m3_monsoon6mo_AR1<-anova(m3_monsoon6mo_AR1,type="marginal")[3,3]
  F_quad_m3_monsoon6mo_AR2<-anova(m3_monsoon6mo_AR2,type="marginal")[3,3]
  F_quad_m3_spring6mo_lag<-anova(m3_spring6mo_lag,type="marginal")[3,3]
  F_quad_m3_spring6mo_lag_AR1<-anova(m3_spring6mo_lag_AR1,type="marginal")[3,3]
  F_quad_m3_spring6mo_lag_AR2<-anova(m3_spring6mo_lag_AR2,type="marginal")[3,3]
  F_quad_m3_monsoon6mo_lag<-anova(m3_monsoon6mo_lag,type="marginal")[3,3]
  F_quad_m3_monsoon6mo_lag_AR1<-anova(m3_monsoon6mo_lag_AR1,type="marginal")[3,3]
  F_quad_m3_monsoon6mo_lag_AR2<-anova(m3_monsoon6mo_lag_AR2,type="marginal")[3,3]
  
  P_quad_m2_spring6mo<-anova(m2_spring6mo,type="marginal")[3,4]
  P_quad_m2_spring6mo_AR1<-anova(m2_spring6mo_AR1,type="marginal")[3,4]
  P_quad_m2_spring6mo_AR2<-anova(m2_spring6mo_AR2,type="marginal")[3,4]
  P_quad_m2_monsoon6mo<-anova(m2_monsoon6mo,type="marginal")[3,4]
  P_quad_m2_monsoon6mo_AR1<-anova(m2_monsoon6mo_AR1,type="marginal")[3,4]
  P_quad_m2_monsoon6mo_AR2<-anova(m2_monsoon6mo_AR2,type="marginal")[3,4]
  P_quad_m2_spring6mo_lag<-anova(m2_spring6mo_lag,type="marginal")[3,4]
  P_quad_m2_spring6mo_lag_AR1<-anova(m2_spring6mo_lag_AR1,type="marginal")[3,4]
  P_quad_m2_spring6mo_lag_AR2<-anova(m2_spring6mo_lag_AR2,type="marginal")[3,4]
  P_quad_m2_monsoon6mo_lag<-anova(m2_monsoon6mo_lag,type="marginal")[3,4]
  P_quad_m2_monsoon6mo_lag_AR1<-anova(m2_monsoon6mo_lag_AR1,type="marginal")[3,4]
  P_quad_m2_monsoon6mo_lag_AR2<-anova(m2_monsoon6mo_lag_AR2,type="marginal")[3,4]
  
  P_quad_m3_spring6mo<-anova(m3_spring6mo,type="marginal")[3,4]
  P_quad_m3_spring6mo_AR1<-anova(m3_spring6mo_AR1,type="marginal")[3,4]
  P_quad_m3_spring6mo_AR2<-anova(m3_spring6mo_AR2,type="marginal")[3,4]
  P_quad_m3_monsoon6mo<-anova(m3_monsoon6mo,type="marginal")[3,4]
  P_quad_m3_monsoon6mo_AR1<-anova(m3_monsoon6mo_AR1,type="marginal")[3,4]
  P_quad_m3_monsoon6mo_AR2<-anova(m3_monsoon6mo_AR2,type="marginal")[3,4]
  P_quad_m3_spring6mo_lag<-anova(m3_spring6mo_lag,type="marginal")[3,4]
  P_quad_m3_spring6mo_lag_AR1<-anova(m3_spring6mo_lag_AR1,type="marginal")[3,4]
  P_quad_m3_spring6mo_lag_AR2<-anova(m3_spring6mo_lag_AR2,type="marginal")[3,4]
  P_quad_m3_monsoon6mo_lag<-anova(m3_monsoon6mo_lag,type="marginal")[3,4]
  P_quad_m3_monsoon6mo_lag_AR1<-anova(m3_monsoon6mo_lag_AR1,type="marginal")[3,4]
  P_quad_m3_monsoon6mo_lag_AR2<-anova(m3_monsoon6mo_lag_AR2,type="marginal")[3,4]
  
  ParamEst_quad_m2_spring6mo<-coef(summary(m2_spring6mo))[3,1]
  ParamEst_quad_m2_spring6mo_AR1<-coef(summary(m2_spring6mo_AR1))[3,1]
  ParamEst_quad_m2_spring6mo_AR2<-coef(summary(m2_spring6mo_AR2))[3,1]
  ParamEst_quad_m2_monsoon6mo<-coef(summary(m2_monsoon6mo))[3,1]
  ParamEst_quad_m2_monsoon6mo_AR1<-coef(summary(m2_monsoon6mo_AR1))[3,1]
  ParamEst_quad_m2_monsoon6mo_AR2<-coef(summary(m2_monsoon6mo_AR2))[3,1]
  ParamEst_quad_m2_spring6mo_lag<-coef(summary(m2_spring6mo_lag))[3,1]
  ParamEst_quad_m2_spring6mo_lag_AR1<-coef(summary(m2_spring6mo_lag_AR1))[3,1]
  ParamEst_quad_m2_spring6mo_lag_AR2<-coef(summary(m2_spring6mo_lag_AR2))[3,1]
  ParamEst_quad_m2_monsoon6mo_lag<-coef(summary(m2_monsoon6mo_lag))[3,1]
  ParamEst_quad_m2_monsoon6mo_lag_AR1<-coef(summary(m2_monsoon6mo_lag_AR1))[3,1]
  ParamEst_quad_m2_monsoon6mo_lag_AR2<-coef(summary(m2_monsoon6mo_lag_AR2))[3,1]
  
  ParamEst_quad_m3_spring6mo<-coef(summary(m3_spring6mo))[3,1]
  ParamEst_quad_m3_spring6mo_AR1<-coef(summary(m3_spring6mo_AR1))[3,1]
  ParamEst_quad_m3_spring6mo_AR2<-coef(summary(m3_spring6mo_AR2))[3,1]
  ParamEst_quad_m3_monsoon6mo<-coef(summary(m3_monsoon6mo))[3,1]
  ParamEst_quad_m3_monsoon6mo_AR1<-coef(summary(m3_monsoon6mo_AR1))[3,1]
  ParamEst_quad_m3_monsoon6mo_AR2<-coef(summary(m3_monsoon6mo_AR2))[3,1]
  ParamEst_quad_m3_spring6mo_lag<-coef(summary(m3_spring6mo_lag))[3,1]
  ParamEst_quad_m3_spring6mo_lag_AR1<-coef(summary(m3_spring6mo_lag_AR1))[3,1]
  ParamEst_quad_m3_spring6mo_lag_AR2<-coef(summary(m3_spring6mo_lag_AR2))[3,1]
  ParamEst_quad_m3_monsoon6mo_lag<-coef(summary(m3_monsoon6mo_lag))[3,1]
  ParamEst_quad_m3_monsoon6mo_lag_AR1<-coef(summary(m3_monsoon6mo_lag_AR1))[3,1]
  ParamEst_quad_m3_monsoon6mo_lag_AR2<-coef(summary(m3_monsoon6mo_lag_AR2))[3,1]
  
  SE_quad_m2_spring6mo<-coef(summary(m2_spring6mo))[3,2]
  SE_quad_m2_spring6mo_AR1<-coef(summary(m2_spring6mo_AR1))[3,2]
  SE_quad_m2_spring6mo_AR2<-coef(summary(m2_spring6mo_AR2))[3,2]
  SE_quad_m2_monsoon6mo<-coef(summary(m2_monsoon6mo))[3,2]
  SE_quad_m2_monsoon6mo_AR1<-coef(summary(m2_monsoon6mo_AR1))[3,2]
  SE_quad_m2_monsoon6mo_AR2<-coef(summary(m2_monsoon6mo_AR2))[3,2]
  SE_quad_m2_spring6mo_lag<-coef(summary(m2_spring6mo_lag))[3,2]
  SE_quad_m2_spring6mo_lag_AR1<-coef(summary(m2_spring6mo_lag_AR1))[3,2]
  SE_quad_m2_spring6mo_lag_AR2<-coef(summary(m2_spring6mo_lag_AR2))[3,2]
  SE_quad_m2_monsoon6mo_lag<-coef(summary(m2_monsoon6mo_lag))[3,2]
  SE_quad_m2_monsoon6mo_lag_AR1<-coef(summary(m2_monsoon6mo_lag_AR1))[3,2]
  SE_quad_m2_monsoon6mo_lag_AR2<-coef(summary(m2_monsoon6mo_lag_AR2))[3,2]
  
  SE_quad_m3_spring6mo<-coef(summary(m3_spring6mo))[3,2]
  SE_quad_m3_spring6mo_AR1<-coef(summary(m3_spring6mo_AR1))[3,2]
  SE_quad_m3_spring6mo_AR2<-coef(summary(m3_spring6mo_AR2))[3,2]
  SE_quad_m3_monsoon6mo<-coef(summary(m3_monsoon6mo))[3,2]
  SE_quad_m3_monsoon6mo_AR1<-coef(summary(m3_monsoon6mo_AR1))[3,2]
  SE_quad_m3_monsoon6mo_AR2<-coef(summary(m3_monsoon6mo_AR2))[3,2]
  SE_quad_m3_spring6mo_lag<-coef(summary(m3_spring6mo_lag))[3,2]
  SE_quad_m3_spring6mo_lag_AR1<-coef(summary(m3_spring6mo_lag_AR1))[3,2]
  SE_quad_m3_spring6mo_lag_AR2<-coef(summary(m3_spring6mo_lag_AR2))[3,2]
  SE_quad_m3_monsoon6mo_lag<-coef(summary(m3_monsoon6mo_lag))[3,2]
  SE_quad_m3_monsoon6mo_lag_AR1<-coef(summary(m3_monsoon6mo_lag_AR1))[3,2]
  SE_quad_m3_monsoon6mo_lag_AR2<-coef(summary(m3_monsoon6mo_lag_AR2))[3,2]
  
  numDF_cub_m3_spring6mo<-anova(m3_spring6mo,type="marginal")[4,1]
  numDF_cub_m3_spring6mo_AR1<-anova(m3_spring6mo_AR1,type="marginal")[4,1]
  numDF_cub_m3_spring6mo_AR2<-anova(m3_spring6mo_AR2,type="marginal")[4,1]
  numDF_cub_m3_monsoon6mo<-anova(m3_monsoon6mo,type="marginal")[4,1]
  numDF_cub_m3_monsoon6mo_AR1<-anova(m3_monsoon6mo_AR1,type="marginal")[4,1]
  numDF_cub_m3_monsoon6mo_AR2<-anova(m3_monsoon6mo_AR2,type="marginal")[4,1]
  numDF_cub_m3_spring6mo_lag<-anova(m3_spring6mo_lag,type="marginal")[4,1]
  numDF_cub_m3_spring6mo_lag_AR1<-anova(m3_spring6mo_lag_AR1,type="marginal")[4,1]
  numDF_cub_m3_spring6mo_lag_AR2<-anova(m3_spring6mo_lag_AR2,type="marginal")[4,1]
  numDF_cub_m3_monsoon6mo_lag<-anova(m3_monsoon6mo_lag,type="marginal")[4,1]
  numDF_cub_m3_monsoon6mo_lag_AR1<-anova(m3_monsoon6mo_lag_AR1,type="marginal")[4,1]
  numDF_cub_m3_monsoon6mo_lag_AR2<-anova(m3_monsoon6mo_lag_AR2,type="marginal")[4,1]
  
  denDF_cub_m3_spring6mo<-anova(m3_spring6mo,type="marginal")[4,2]
  denDF_cub_m3_spring6mo_AR1<-anova(m3_spring6mo_AR1,type="marginal")[4,2]
  denDF_cub_m3_spring6mo_AR2<-anova(m3_spring6mo_AR2,type="marginal")[4,2]
  denDF_cub_m3_monsoon6mo<-anova(m3_monsoon6mo,type="marginal")[4,2]
  denDF_cub_m3_monsoon6mo_AR1<-anova(m3_monsoon6mo_AR1,type="marginal")[4,2]
  denDF_cub_m3_monsoon6mo_AR2<-anova(m3_monsoon6mo_AR2,type="marginal")[4,2]
  denDF_cub_m3_spring6mo_lag<-anova(m3_spring6mo_lag,type="marginal")[4,2]
  denDF_cub_m3_spring6mo_lag_AR1<-anova(m3_spring6mo_lag_AR1,type="marginal")[4,2]
  denDF_cub_m3_spring6mo_lag_AR2<-anova(m3_spring6mo_lag_AR2,type="marginal")[4,2]
  denDF_cub_m3_monsoon6mo_lag<-anova(m3_monsoon6mo_lag,type="marginal")[4,2]
  denDF_cub_m3_monsoon6mo_lag_AR1<-anova(m3_monsoon6mo_lag_AR1,type="marginal")[4,2]
  denDF_cub_m3_monsoon6mo_lag_AR2<-anova(m3_monsoon6mo_lag_AR2,type="marginal")[4,2]
  
  F_cub_m3_spring6mo<-anova(m3_spring6mo,type="marginal")[4,3]
  F_cub_m3_spring6mo_AR1<-anova(m3_spring6mo_AR1,type="marginal")[4,3]
  F_cub_m3_spring6mo_AR2<-anova(m3_spring6mo_AR2,type="marginal")[4,3]
  F_cub_m3_monsoon6mo<-anova(m3_monsoon6mo,type="marginal")[4,3]
  F_cub_m3_monsoon6mo_AR1<-anova(m3_monsoon6mo_AR1,type="marginal")[4,3]
  F_cub_m3_monsoon6mo_AR2<-anova(m3_monsoon6mo_AR2,type="marginal")[4,3]
  F_cub_m3_spring6mo_lag<-anova(m3_spring6mo_lag,type="marginal")[4,3]
  F_cub_m3_spring6mo_lag_AR1<-anova(m3_spring6mo_lag_AR1,type="marginal")[4,3]
  F_cub_m3_spring6mo_lag_AR2<-anova(m3_spring6mo_lag_AR2,type="marginal")[4,3]
  F_cub_m3_monsoon6mo_lag<-anova(m3_monsoon6mo_lag,type="marginal")[4,3]
  F_cub_m3_monsoon6mo_lag_AR1<-anova(m3_monsoon6mo_lag_AR1,type="marginal")[4,3]
  F_cub_m3_monsoon6mo_lag_AR2<-anova(m3_monsoon6mo_lag_AR2,type="marginal")[4,3]
  
  P_cub_m3_spring6mo<-anova(m3_spring6mo,type="marginal")[4,4]
  P_cub_m3_spring6mo_AR1<-anova(m3_spring6mo_AR1,type="marginal")[4,4]
  P_cub_m3_spring6mo_AR2<-anova(m3_spring6mo_AR2,type="marginal")[4,4]
  P_cub_m3_monsoon6mo<-anova(m3_monsoon6mo,type="marginal")[4,4]
  P_cub_m3_monsoon6mo_AR1<-anova(m3_monsoon6mo_AR1,type="marginal")[4,4]
  P_cub_m3_monsoon6mo_AR2<-anova(m3_monsoon6mo_AR2,type="marginal")[4,4]
  P_cub_m3_spring6mo_lag<-anova(m3_spring6mo_lag,type="marginal")[4,4]
  P_cub_m3_spring6mo_lag_AR1<-anova(m3_spring6mo_lag_AR1,type="marginal")[4,4]
  P_cub_m3_spring6mo_lag_AR2<-anova(m3_spring6mo_lag_AR2,type="marginal")[4,4]
  P_cub_m3_monsoon6mo_lag<-anova(m3_monsoon6mo_lag,type="marginal")[4,4]
  P_cub_m3_monsoon6mo_lag_AR1<-anova(m3_monsoon6mo_lag_AR1,type="marginal")[4,4]
  P_cub_m3_monsoon6mo_lag_AR2<-anova(m3_monsoon6mo_lag_AR2,type="marginal")[4,4]
  
  ParamEst_cub_m3_spring6mo<-coef(summary(m3_spring6mo))[4,1]
  ParamEst_cub_m3_spring6mo_AR1<-coef(summary(m3_spring6mo_AR1))[4,1]
  ParamEst_cub_m3_spring6mo_AR2<-coef(summary(m3_spring6mo_AR2))[4,1]
  ParamEst_cub_m3_monsoon6mo<-coef(summary(m3_monsoon6mo))[4,1]
  ParamEst_cub_m3_monsoon6mo_AR1<-coef(summary(m3_monsoon6mo_AR1))[4,1]
  ParamEst_cub_m3_monsoon6mo_AR2<-coef(summary(m3_monsoon6mo_AR2))[4,1]
  ParamEst_cub_m3_spring6mo_lag<-coef(summary(m3_spring6mo_lag))[4,1]
  ParamEst_cub_m3_spring6mo_lag_AR1<-coef(summary(m3_spring6mo_lag_AR1))[4,1]
  ParamEst_cub_m3_spring6mo_lag_AR2<-coef(summary(m3_spring6mo_lag_AR2))[4,1]
  ParamEst_cub_m3_monsoon6mo_lag<-coef(summary(m3_monsoon6mo_lag))[4,1]
  ParamEst_cub_m3_monsoon6mo_lag_AR1<-coef(summary(m3_monsoon6mo_lag_AR1))[4,1]
  ParamEst_cub_m3_monsoon6mo_lag_AR2<-coef(summary(m3_monsoon6mo_lag_AR2))[4,1]
  
  SE_cub_m3_spring6mo<-coef(summary(m3_spring6mo))[4,2]
  SE_cub_m3_spring6mo_AR1<-coef(summary(m3_spring6mo_AR1))[4,2]
  SE_cub_m3_spring6mo_AR2<-coef(summary(m3_spring6mo_AR2))[4,2]
  SE_cub_m3_monsoon6mo<-coef(summary(m3_monsoon6mo))[4,2]
  SE_cub_m3_monsoon6mo_AR1<-coef(summary(m3_monsoon6mo_AR1))[4,2]
  SE_cub_m3_monsoon6mo_AR2<-coef(summary(m3_monsoon6mo_AR2))[4,2]
  SE_cub_m3_spring6mo_lag<-coef(summary(m3_spring6mo_lag))[4,2]
  SE_cub_m3_spring6mo_lag_AR1<-coef(summary(m3_spring6mo_lag_AR1))[4,2]
  SE_cub_m3_spring6mo_lag_AR2<-coef(summary(m3_spring6mo_lag_AR2))[4,2]
  SE_cub_m3_monsoon6mo_lag<-coef(summary(m3_monsoon6mo_lag))[4,2]
  SE_cub_m3_monsoon6mo_lag_AR1<-coef(summary(m3_monsoon6mo_lag_AR1))[4,2]
  SE_cub_m3_monsoon6mo_lag_AR2<-coef(summary(m3_monsoon6mo_lag_AR2))[4,2]
  
  # bind all target output values together
  output_id<-cbind(speciesCode, 
                   ecosystem, 
                   P_shapiro_null,
                   dAICc_m1_monsoon6mo,
                   dAICc_m1_monsoon6mo_lag,
                   dAICc_m1_spring6mo,
                   dAICc_m1_spring6mo_lag,
                   dAICc_m2_monsoon6mo,
                   dAICc_m2_monsoon6mo_lag,
                   dAICc_m2_spring6mo,
                   dAICc_m2_spring6mo_lag,
                   dAICc_m3_monsoon6mo,
                   dAICc_m3_monsoon6mo_lag,
                   dAICc_m3_spring6mo,
                   dAICc_m3_spring6mo_lag,
                   dAICc_null,
                   dAICc_m1_spring6mo_AR1,
                   dAICc_m1_spring6mo_AR2,
                   dAICc_m2_spring6mo_AR1,
                   dAICc_m2_spring6mo_AR2,
                   dAICc_m3_spring6mo_AR1,
                   dAICc_m3_spring6mo_AR2,
                   dAICc_m1_monsoon6mo_AR1,
                   dAICc_m1_monsoon6mo_AR2,
                   dAICc_m2_monsoon6mo_AR1,
                   dAICc_m2_monsoon6mo_AR2,
                   dAICc_m3_monsoon6mo_AR1,
                   dAICc_m3_monsoon6mo_AR2,
                   dAICc_m1_spring6mo_lag_AR1,
                   dAICc_m1_spring6mo_lag_AR2,
                   dAICc_m2_spring6mo_lag_AR1,
                   dAICc_m2_spring6mo_lag_AR2,
                   dAICc_m3_spring6mo_lag_AR1,
                   dAICc_m3_spring6mo_lag_AR2,
                   dAICc_m1_monsoon6mo_lag_AR1,
                   dAICc_m1_monsoon6mo_lag_AR2,
                   dAICc_m2_monsoon6mo_lag_AR1,
                   dAICc_m2_monsoon6mo_lag_AR2,
                   dAICc_m3_monsoon6mo_lag_AR1,
                   dAICc_m3_monsoon6mo_lag_AR2,
                   numDF_lin_m1_spring6mo,
                   denDF_lin_m1_spring6mo,
                   F_lin_m1_spring6mo,
                   P_lin_m1_spring6mo,
                   numDF_lin_m2_spring6mo,
                   denDF_lin_m2_spring6mo,
                   F_lin_m2_spring6mo,
                   P_lin_m2_spring6mo,
                   numDF_lin_m3_spring6mo,
                   denDF_lin_m3_spring6mo,
                   F_lin_m3_spring6mo,
                   P_lin_m3_spring6mo,
                   numDF_lin_m1_monsoon6mo,
                   denDF_lin_m1_monsoon6mo,
                   F_lin_m1_monsoon6mo,
                   P_lin_m1_monsoon6mo,
                   numDF_lin_m2_monsoon6mo,
                   denDF_lin_m2_monsoon6mo,
                   F_lin_m2_monsoon6mo,
                   P_lin_m2_monsoon6mo,
                   numDF_lin_m3_monsoon6mo,
                   denDF_lin_m3_monsoon6mo,
                   F_lin_m3_monsoon6mo,
                   P_lin_m3_monsoon6mo,
                   numDF_lin_m1_spring6mo_lag,
                   denDF_lin_m1_spring6mo_lag,
                   F_lin_m1_spring6mo_lag,
                   P_lin_m1_spring6mo_lag,
                   numDF_lin_m2_spring6mo_lag,
                   denDF_lin_m2_spring6mo_lag,
                   F_lin_m2_spring6mo_lag,
                   P_lin_m2_spring6mo_lag,
                   numDF_lin_m3_spring6mo_lag,
                   denDF_lin_m3_spring6mo_lag,
                   F_lin_m3_spring6mo_lag,
                   P_lin_m3_spring6mo_lag,
                   numDF_lin_m1_monsoon6mo_lag,
                   denDF_lin_m1_monsoon6mo_lag,
                   F_lin_m1_monsoon6mo_lag,
                   P_lin_m1_monsoon6mo_lag,
                   numDF_lin_m2_monsoon6mo_lag,
                   denDF_lin_m2_monsoon6mo_lag,
                   F_lin_m2_monsoon6mo_lag,
                   P_lin_m2_monsoon6mo_lag,
                   numDF_lin_m3_monsoon6mo_lag,
                   denDF_lin_m3_monsoon6mo_lag,
                   F_lin_m3_monsoon6mo_lag,
                   P_lin_m3_monsoon6mo_lag,
                   ParamEst_lin_m1_monsoon6mo,
                   ParamEst_lin_m1_monsoon6mo_lag,
                   ParamEst_lin_m1_spring6mo,
                   ParamEst_lin_m1_spring6mo_lag,
                   ParamEst_lin_m2_monsoon6mo,
                   ParamEst_lin_m2_monsoon6mo_lag,
                   ParamEst_lin_m2_spring6mo,
                   ParamEst_lin_m2_spring6mo_lag,
                   ParamEst_lin_m3_monsoon6mo,
                   ParamEst_lin_m3_monsoon6mo_lag,
                   ParamEst_lin_m3_spring6mo,
                   ParamEst_lin_m3_spring6mo_lag,
                   SE_lin_m1_monsoon6mo,
                   SE_lin_m1_monsoon6mo_lag,
                   SE_lin_m1_spring6mo,
                   SE_lin_m1_spring6mo_lag,
                   SE_lin_m2_monsoon6mo,
                   SE_lin_m2_monsoon6mo_lag,
                   SE_lin_m2_spring6mo,
                   SE_lin_m2_spring6mo_lag,
                   SE_lin_m3_monsoon6mo,
                   SE_lin_m3_monsoon6mo_lag,
                   SE_lin_m3_spring6mo,
                   SE_lin_m3_spring6mo_lag,
                   numDF_lin_m1_spring6mo_AR1,
                   denDF_lin_m1_spring6mo_AR1,
                   F_lin_m1_spring6mo_AR1,
                   P_lin_m1_spring6mo_AR1,
                   ParamEst_lin_m1_spring6mo_AR1,
                   SE_lin_m1_spring6mo_AR1,
                   numDF_lin_m1_spring6mo_AR2,
                   denDF_lin_m1_spring6mo_AR2,
                   F_lin_m1_spring6mo_AR2,
                   P_lin_m1_spring6mo_AR2,
                   ParamEst_lin_m1_spring6mo_AR2,
                   SE_lin_m1_spring6mo_AR2,
                   numDF_lin_m2_spring6mo_AR1,
                   denDF_lin_m2_spring6mo_AR1,
                   F_lin_m2_spring6mo_AR1,
                   P_lin_m2_spring6mo_AR1,
                   ParamEst_lin_m2_spring6mo_AR1,
                   SE_lin_m2_spring6mo_AR1,
                   numDF_lin_m2_spring6mo_AR2,
                   denDF_lin_m2_spring6mo_AR2,
                   F_lin_m2_spring6mo_AR2,
                   P_lin_m2_spring6mo_AR2,
                   ParamEst_lin_m2_spring6mo_AR2,
                   SE_lin_m2_spring6mo_AR2,
                   numDF_lin_m3_spring6mo_AR1,
                   denDF_lin_m3_spring6mo_AR1,
                   F_lin_m3_spring6mo_AR1,
                   P_lin_m3_spring6mo_AR1,
                   ParamEst_lin_m3_spring6mo_AR1,
                   SE_lin_m3_spring6mo_AR1,
                   numDF_lin_m3_spring6mo_AR2,
                   denDF_lin_m3_spring6mo_AR2,
                   F_lin_m3_spring6mo_AR2,
                   P_lin_m3_spring6mo_AR2,
                   ParamEst_lin_m3_spring6mo_AR2,
                   SE_lin_m3_spring6mo_AR2,
                   numDF_lin_m1_monsoon6mo_AR1,
                   denDF_lin_m1_monsoon6mo_AR1,
                   F_lin_m1_monsoon6mo_AR1,
                   P_lin_m1_monsoon6mo_AR1,
                   ParamEst_lin_m1_monsoon6mo_AR1,
                   SE_lin_m1_monsoon6mo_AR1,
                   numDF_lin_m1_monsoon6mo_AR2,
                   denDF_lin_m1_monsoon6mo_AR2,
                   F_lin_m1_monsoon6mo_AR2,
                   P_lin_m1_monsoon6mo_AR2,
                   ParamEst_lin_m1_monsoon6mo_AR2,
                   SE_lin_m1_monsoon6mo_AR2,
                   numDF_lin_m2_monsoon6mo_AR1,
                   denDF_lin_m2_monsoon6mo_AR1,
                   F_lin_m2_monsoon6mo_AR1,
                   P_lin_m2_monsoon6mo_AR1,
                   ParamEst_lin_m2_monsoon6mo_AR1,
                   SE_lin_m2_monsoon6mo_AR1,
                   numDF_lin_m2_monsoon6mo_AR2,
                   denDF_lin_m2_monsoon6mo_AR2,
                   F_lin_m2_monsoon6mo_AR2,
                   P_lin_m2_monsoon6mo_AR2,
                   ParamEst_lin_m2_monsoon6mo_AR2,
                   SE_lin_m2_monsoon6mo_AR2,
                   numDF_lin_m3_monsoon6mo_AR1,
                   denDF_lin_m3_monsoon6mo_AR1,
                   F_lin_m3_monsoon6mo_AR1,
                   P_lin_m3_monsoon6mo_AR1,
                   ParamEst_lin_m3_monsoon6mo_AR1,
                   SE_lin_m3_monsoon6mo_AR1,
                   numDF_lin_m3_monsoon6mo_AR2,
                   denDF_lin_m3_monsoon6mo_AR2,
                   F_lin_m3_monsoon6mo_AR2,
                   P_lin_m3_monsoon6mo_AR2,
                   ParamEst_lin_m3_monsoon6mo_AR2,
                   SE_lin_m3_monsoon6mo_AR2,
                   numDF_lin_m1_spring6mo_lag_AR1,
                   denDF_lin_m1_spring6mo_lag_AR1,
                   F_lin_m1_spring6mo_lag_AR1,
                   P_lin_m1_spring6mo_lag_AR1,
                   ParamEst_lin_m1_spring6mo_lag_AR1,
                   SE_lin_m1_spring6mo_lag_AR1,
                   numDF_lin_m1_spring6mo_lag_AR2,
                   denDF_lin_m1_spring6mo_lag_AR2,
                   F_lin_m1_spring6mo_lag_AR2,
                   P_lin_m1_spring6mo_lag_AR2,
                   ParamEst_lin_m1_spring6mo_lag_AR2,
                   SE_lin_m1_spring6mo_lag_AR2,
                   numDF_lin_m2_spring6mo_lag_AR1,
                   denDF_lin_m2_spring6mo_lag_AR1,
                   F_lin_m2_spring6mo_lag_AR1,
                   P_lin_m2_spring6mo_lag_AR1,
                   ParamEst_lin_m2_spring6mo_lag_AR1,
                   SE_lin_m2_spring6mo_lag_AR1,
                   numDF_lin_m2_spring6mo_lag_AR2,
                   denDF_lin_m2_spring6mo_lag_AR2,
                   F_lin_m2_spring6mo_lag_AR2,
                   P_lin_m2_spring6mo_lag_AR2,
                   ParamEst_lin_m2_spring6mo_lag_AR2,
                   SE_lin_m2_spring6mo_lag_AR2,
                   numDF_lin_m3_spring6mo_lag_AR1,
                   denDF_lin_m3_spring6mo_lag_AR1,
                   F_lin_m3_spring6mo_lag_AR1,
                   P_lin_m3_spring6mo_lag_AR1,
                   ParamEst_lin_m3_spring6mo_lag_AR1,
                   SE_lin_m3_spring6mo_lag_AR1,
                   numDF_lin_m3_spring6mo_lag_AR2,
                   denDF_lin_m3_spring6mo_lag_AR2,
                   F_lin_m3_spring6mo_lag_AR2,
                   P_lin_m3_spring6mo_lag_AR2,
                   ParamEst_lin_m3_spring6mo_lag_AR2,
                   SE_lin_m3_spring6mo_lag_AR2,
                   numDF_lin_m1_monsoon6mo_lag_AR1,
                   denDF_lin_m1_monsoon6mo_lag_AR1,
                   F_lin_m1_monsoon6mo_lag_AR1,
                   P_lin_m1_monsoon6mo_lag_AR1,
                   ParamEst_lin_m1_monsoon6mo_lag_AR1,
                   SE_lin_m1_monsoon6mo_lag_AR1,
                   numDF_lin_m1_monsoon6mo_lag_AR2,
                   denDF_lin_m1_monsoon6mo_lag_AR2,
                   F_lin_m1_monsoon6mo_lag_AR2,
                   P_lin_m1_monsoon6mo_lag_AR2,
                   ParamEst_lin_m1_monsoon6mo_lag_AR2,
                   SE_lin_m1_monsoon6mo_lag_AR2,
                   numDF_lin_m2_monsoon6mo_lag_AR1,
                   denDF_lin_m2_monsoon6mo_lag_AR1,
                   F_lin_m2_monsoon6mo_lag_AR1,
                   P_lin_m2_monsoon6mo_lag_AR1,
                   ParamEst_lin_m2_monsoon6mo_lag_AR1,
                   SE_lin_m2_monsoon6mo_lag_AR1,
                   numDF_lin_m2_monsoon6mo_lag_AR2,
                   denDF_lin_m2_monsoon6mo_lag_AR2,
                   F_lin_m2_monsoon6mo_lag_AR2,
                   P_lin_m2_monsoon6mo_lag_AR2,
                   ParamEst_lin_m2_monsoon6mo_lag_AR2,
                   SE_lin_m2_monsoon6mo_lag_AR2,
                   numDF_lin_m3_monsoon6mo_lag_AR1,
                   denDF_lin_m3_monsoon6mo_lag_AR1,
                   F_lin_m3_monsoon6mo_lag_AR1,
                   P_lin_m3_monsoon6mo_lag_AR1,
                   ParamEst_lin_m3_monsoon6mo_lag_AR1,
                   SE_lin_m3_monsoon6mo_lag_AR1,
                   numDF_lin_m3_monsoon6mo_lag_AR2,
                   denDF_lin_m3_monsoon6mo_lag_AR2,
                   F_lin_m3_monsoon6mo_lag_AR2,
                   P_lin_m3_monsoon6mo_lag_AR2,
                   ParamEst_lin_m3_monsoon6mo_lag_AR2,
                   SE_lin_m3_monsoon6mo_lag_AR2,
                   Rsquared_marginal_m_null,
                   Rsquared_marginal_m1_spring6mo,
                   Rsquared_marginal_m2_spring6mo,
                   Rsquared_marginal_m3_spring6mo,
                   Rsquared_marginal_m1_spring6mo_AR1,
                   Rsquared_marginal_m1_spring6mo_AR2,
                   Rsquared_marginal_m2_spring6mo_AR1,
                   Rsquared_marginal_m2_spring6mo_AR2,
                   Rsquared_marginal_m3_spring6mo_AR1,
                   Rsquared_marginal_m3_spring6mo_AR2,
                   Rsquared_marginal_m1_monsoon6mo,
                   Rsquared_marginal_m2_monsoon6mo,
                   Rsquared_marginal_m3_monsoon6mo,
                   Rsquared_marginal_m1_monsoon6mo_AR1,
                   Rsquared_marginal_m1_monsoon6mo_AR2,
                   Rsquared_marginal_m2_monsoon6mo_AR1,
                   Rsquared_marginal_m2_monsoon6mo_AR2,
                   Rsquared_marginal_m3_monsoon6mo_AR1,
                   Rsquared_marginal_m3_monsoon6mo_AR2,
                   Rsquared_marginal_m1_spring6mo_lag,
                   Rsquared_marginal_m2_spring6mo_lag,
                   Rsquared_marginal_m3_spring6mo_lag,
                   Rsquared_marginal_m1_spring6mo_lag_AR1,
                   Rsquared_marginal_m1_spring6mo_lag_AR2,
                   Rsquared_marginal_m2_spring6mo_lag_AR1,
                   Rsquared_marginal_m2_spring6mo_lag_AR2,
                   Rsquared_marginal_m3_spring6mo_lag_AR1,
                   Rsquared_marginal_m3_spring6mo_lag_AR2,
                   Rsquared_marginal_m1_monsoon6mo_lag,
                   Rsquared_marginal_m2_monsoon6mo_lag,
                   Rsquared_marginal_m3_monsoon6mo_lag,
                   Rsquared_marginal_m1_monsoon6mo_lag_AR1,
                   Rsquared_marginal_m1_monsoon6mo_lag_AR2,
                   Rsquared_marginal_m2_monsoon6mo_lag_AR1,
                   Rsquared_marginal_m2_monsoon6mo_lag_AR2,
                   Rsquared_marginal_m3_monsoon6mo_lag_AR1,
                   Rsquared_marginal_m3_monsoon6mo_lag_AR2,
                   Rsquared_conditional_m_null,
                   Rsquared_conditional_m1_spring6mo,
                   Rsquared_conditional_m2_spring6mo,
                   Rsquared_conditional_m3_spring6mo,
                   Rsquared_conditional_m1_spring6mo_AR1,
                   Rsquared_conditional_m1_spring6mo_AR2,
                   Rsquared_conditional_m2_spring6mo_AR1,
                   Rsquared_conditional_m2_spring6mo_AR2,
                   Rsquared_conditional_m3_spring6mo_AR1,
                   Rsquared_conditional_m3_spring6mo_AR2,
                   Rsquared_conditional_m1_monsoon6mo,
                   Rsquared_conditional_m2_monsoon6mo,
                   Rsquared_conditional_m3_monsoon6mo,
                   Rsquared_conditional_m1_monsoon6mo_AR1,
                   Rsquared_conditional_m1_monsoon6mo_AR2,
                   Rsquared_conditional_m2_monsoon6mo_AR1,
                   Rsquared_conditional_m2_monsoon6mo_AR2,
                   Rsquared_conditional_m3_monsoon6mo_AR1,
                   Rsquared_conditional_m3_monsoon6mo_AR2,
                   Rsquared_conditional_m1_spring6mo_lag,
                   Rsquared_conditional_m2_spring6mo_lag,
                   Rsquared_conditional_m3_spring6mo_lag,
                   Rsquared_conditional_m1_spring6mo_lag_AR1,
                   Rsquared_conditional_m1_spring6mo_lag_AR2,
                   Rsquared_conditional_m2_spring6mo_lag_AR1,
                   Rsquared_conditional_m2_spring6mo_lag_AR2,
                   Rsquared_conditional_m3_spring6mo_lag_AR1,
                   Rsquared_conditional_m3_spring6mo_lag_AR2,
                   Rsquared_conditional_m1_monsoon6mo_lag,
                   Rsquared_conditional_m2_monsoon6mo_lag,
                   Rsquared_conditional_m3_monsoon6mo_lag,
                   Rsquared_conditional_m1_monsoon6mo_lag_AR1,
                   Rsquared_conditional_m1_monsoon6mo_lag_AR2,
                   Rsquared_conditional_m2_monsoon6mo_lag_AR1,
                   Rsquared_conditional_m2_monsoon6mo_lag_AR2,
                   Rsquared_conditional_m3_monsoon6mo_lag_AR1,
                   Rsquared_conditional_m3_monsoon6mo_lag_AR2,
                   numDF_quad_m2_spring6mo,
                   numDF_quad_m2_spring6mo_AR1,
                   numDF_quad_m2_spring6mo_AR2,
                   numDF_quad_m2_monsoon6mo,
                   numDF_quad_m2_monsoon6mo_AR1,
                   numDF_quad_m2_monsoon6mo_AR2,
                   numDF_quad_m2_spring6mo_lag,
                   numDF_quad_m2_spring6mo_lag_AR1,
                   numDF_quad_m2_spring6mo_lag_AR2,
                   numDF_quad_m2_monsoon6mo_lag,
                   numDF_quad_m2_monsoon6mo_lag_AR1,
                   numDF_quad_m2_monsoon6mo_lag_AR2,
                   numDF_quad_m3_spring6mo,
                   numDF_quad_m3_spring6mo_AR1,
                   numDF_quad_m3_spring6mo_AR2,
                   numDF_quad_m3_monsoon6mo,
                   numDF_quad_m3_monsoon6mo_AR1,
                   numDF_quad_m3_monsoon6mo_AR2,
                   numDF_quad_m3_spring6mo_lag,
                   numDF_quad_m3_spring6mo_lag_AR1,
                   numDF_quad_m3_spring6mo_lag_AR2,
                   numDF_quad_m3_monsoon6mo_lag,
                   numDF_quad_m3_monsoon6mo_lag_AR1,
                   numDF_quad_m3_monsoon6mo_lag_AR2,
                   denDF_quad_m2_spring6mo,
                   denDF_quad_m2_spring6mo_AR1,
                   denDF_quad_m2_spring6mo_AR2,
                   denDF_quad_m2_monsoon6mo,
                   denDF_quad_m2_monsoon6mo_AR1,
                   denDF_quad_m2_monsoon6mo_AR2,
                   denDF_quad_m2_spring6mo_lag,
                   denDF_quad_m2_spring6mo_lag_AR1,
                   denDF_quad_m2_spring6mo_lag_AR2,
                   denDF_quad_m2_monsoon6mo_lag,
                   denDF_quad_m2_monsoon6mo_lag_AR1,
                   denDF_quad_m2_monsoon6mo_lag_AR2,
                   denDF_quad_m3_spring6mo,
                   denDF_quad_m3_spring6mo_AR1,
                   denDF_quad_m3_spring6mo_AR2,
                   denDF_quad_m3_monsoon6mo,
                   denDF_quad_m3_monsoon6mo_AR1,
                   denDF_quad_m3_monsoon6mo_AR2,
                   denDF_quad_m3_spring6mo_lag,
                   denDF_quad_m3_spring6mo_lag_AR1,
                   denDF_quad_m3_spring6mo_lag_AR2,
                   denDF_quad_m3_monsoon6mo_lag,
                   denDF_quad_m3_monsoon6mo_lag_AR1,
                   denDF_quad_m3_monsoon6mo_lag_AR2,
                   F_quad_m2_spring6mo,
                   F_quad_m2_spring6mo_AR1,
                   F_quad_m2_spring6mo_AR2,
                   F_quad_m2_monsoon6mo,
                   F_quad_m2_monsoon6mo_AR1,
                   F_quad_m2_monsoon6mo_AR2,
                   F_quad_m2_spring6mo_lag,
                   F_quad_m2_spring6mo_lag_AR1,
                   F_quad_m2_spring6mo_lag_AR2,
                   F_quad_m2_monsoon6mo_lag,
                   F_quad_m2_monsoon6mo_lag_AR1,
                   F_quad_m2_monsoon6mo_lag_AR2,
                   F_quad_m3_spring6mo,
                   F_quad_m3_spring6mo_AR1,
                   F_quad_m3_spring6mo_AR2,
                   F_quad_m3_monsoon6mo,
                   F_quad_m3_monsoon6mo_AR1,
                   F_quad_m3_monsoon6mo_AR2,
                   F_quad_m3_spring6mo_lag,
                   F_quad_m3_spring6mo_lag_AR1,
                   F_quad_m3_spring6mo_lag_AR2,
                   F_quad_m3_monsoon6mo_lag,
                   F_quad_m3_monsoon6mo_lag_AR1,
                   F_quad_m3_monsoon6mo_lag_AR2,
                   P_quad_m2_spring6mo,
                   P_quad_m2_spring6mo_AR1,
                   P_quad_m2_spring6mo_AR2,
                   P_quad_m2_monsoon6mo,
                   P_quad_m2_monsoon6mo_AR1,
                   P_quad_m2_monsoon6mo_AR2,
                   P_quad_m2_spring6mo_lag,
                   P_quad_m2_spring6mo_lag_AR1,
                   P_quad_m2_spring6mo_lag_AR2,
                   P_quad_m2_monsoon6mo_lag,
                   P_quad_m2_monsoon6mo_lag_AR1,
                   P_quad_m2_monsoon6mo_lag_AR2,
                   P_quad_m3_spring6mo,
                   P_quad_m3_spring6mo_AR1,
                   P_quad_m3_spring6mo_AR2,
                   P_quad_m3_monsoon6mo,
                   P_quad_m3_monsoon6mo_AR1,
                   P_quad_m3_monsoon6mo_AR2,
                   P_quad_m3_spring6mo_lag,
                   P_quad_m3_spring6mo_lag_AR1,
                   P_quad_m3_spring6mo_lag_AR2,
                   P_quad_m3_monsoon6mo_lag,
                   P_quad_m3_monsoon6mo_lag_AR1,
                   P_quad_m3_monsoon6mo_lag_AR2,
                   ParamEst_quad_m2_spring6mo,
                   ParamEst_quad_m2_spring6mo_AR1,
                   ParamEst_quad_m2_spring6mo_AR2,
                   ParamEst_quad_m2_monsoon6mo,
                   ParamEst_quad_m2_monsoon6mo_AR1,
                   ParamEst_quad_m2_monsoon6mo_AR2,
                   ParamEst_quad_m2_spring6mo_lag,
                   ParamEst_quad_m2_spring6mo_lag_AR1,
                   ParamEst_quad_m2_spring6mo_lag_AR2,
                   ParamEst_quad_m2_monsoon6mo_lag,
                   ParamEst_quad_m2_monsoon6mo_lag_AR1,
                   ParamEst_quad_m2_monsoon6mo_lag_AR2,
                   ParamEst_quad_m3_spring6mo,
                   ParamEst_quad_m3_spring6mo_AR1,
                   ParamEst_quad_m3_spring6mo_AR2,
                   ParamEst_quad_m3_monsoon6mo,
                   ParamEst_quad_m3_monsoon6mo_AR1,
                   ParamEst_quad_m3_monsoon6mo_AR2,
                   ParamEst_quad_m3_spring6mo_lag,
                   ParamEst_quad_m3_spring6mo_lag_AR1,
                   ParamEst_quad_m3_spring6mo_lag_AR2,
                   ParamEst_quad_m3_monsoon6mo_lag,
                   ParamEst_quad_m3_monsoon6mo_lag_AR1,
                   ParamEst_quad_m3_monsoon6mo_lag_AR2,
                   SE_quad_m2_spring6mo,
                   SE_quad_m2_spring6mo_AR1,
                   SE_quad_m2_spring6mo_AR2,
                   SE_quad_m2_monsoon6mo,
                   SE_quad_m2_monsoon6mo_AR1,
                   SE_quad_m2_monsoon6mo_AR2,
                   SE_quad_m2_spring6mo_lag,
                   SE_quad_m2_spring6mo_lag_AR1,
                   SE_quad_m2_spring6mo_lag_AR2,
                   SE_quad_m2_monsoon6mo_lag,
                   SE_quad_m2_monsoon6mo_lag_AR1,
                   SE_quad_m2_monsoon6mo_lag_AR2,
                   SE_quad_m3_spring6mo,
                   SE_quad_m3_spring6mo_AR1,
                   SE_quad_m3_spring6mo_AR2,
                   SE_quad_m3_monsoon6mo,
                   SE_quad_m3_monsoon6mo_AR1,
                   SE_quad_m3_monsoon6mo_AR2,
                   SE_quad_m3_spring6mo_lag,
                   SE_quad_m3_spring6mo_lag_AR1,
                   SE_quad_m3_spring6mo_lag_AR2,
                   SE_quad_m3_monsoon6mo_lag,
                   SE_quad_m3_monsoon6mo_lag_AR1,
                   SE_quad_m3_monsoon6mo_lag_AR2,
                   numDF_cub_m3_spring6mo,
                   numDF_cub_m3_spring6mo_AR1,
                   numDF_cub_m3_spring6mo_AR2,
                   numDF_cub_m3_monsoon6mo,
                   numDF_cub_m3_monsoon6mo_AR1,
                   numDF_cub_m3_monsoon6mo_AR2,
                   numDF_cub_m3_spring6mo_lag,
                   numDF_cub_m3_spring6mo_lag_AR1,
                   numDF_cub_m3_spring6mo_lag_AR2,
                   numDF_cub_m3_monsoon6mo_lag,
                   numDF_cub_m3_monsoon6mo_lag_AR1,
                   numDF_cub_m3_monsoon6mo_lag_AR2,
                   denDF_cub_m3_spring6mo,
                   denDF_cub_m3_spring6mo_AR1,
                   denDF_cub_m3_spring6mo_AR2,
                   denDF_cub_m3_monsoon6mo,
                   denDF_cub_m3_monsoon6mo_AR1,
                   denDF_cub_m3_monsoon6mo_AR2,
                   denDF_cub_m3_spring6mo_lag,
                   denDF_cub_m3_spring6mo_lag_AR1,
                   denDF_cub_m3_spring6mo_lag_AR2,
                   denDF_cub_m3_monsoon6mo_lag,
                   denDF_cub_m3_monsoon6mo_lag_AR1,
                   denDF_cub_m3_monsoon6mo_lag_AR2,
                   F_cub_m3_spring6mo,
                   F_cub_m3_spring6mo_AR1,
                   F_cub_m3_spring6mo_AR2,
                   F_cub_m3_monsoon6mo,
                   F_cub_m3_monsoon6mo_AR1,
                   F_cub_m3_monsoon6mo_AR2,
                   F_cub_m3_spring6mo_lag,
                   F_cub_m3_spring6mo_lag_AR1,
                   F_cub_m3_spring6mo_lag_AR2,
                   F_cub_m3_monsoon6mo_lag,
                   F_cub_m3_monsoon6mo_lag_AR1,
                   F_cub_m3_monsoon6mo_lag_AR2,
                   P_cub_m3_spring6mo,
                   P_cub_m3_spring6mo_AR1,
                   P_cub_m3_spring6mo_AR2,
                   P_cub_m3_monsoon6mo,
                   P_cub_m3_monsoon6mo_AR1,
                   P_cub_m3_monsoon6mo_AR2,
                   P_cub_m3_spring6mo_lag,
                   P_cub_m3_spring6mo_lag_AR1,
                   P_cub_m3_spring6mo_lag_AR2,
                   P_cub_m3_monsoon6mo_lag,
                   P_cub_m3_monsoon6mo_lag_AR1,
                   P_cub_m3_monsoon6mo_lag_AR2,
                   ParamEst_cub_m3_spring6mo,
                   ParamEst_cub_m3_spring6mo_AR1,
                   ParamEst_cub_m3_spring6mo_AR2,
                   ParamEst_cub_m3_monsoon6mo,
                   ParamEst_cub_m3_monsoon6mo_AR1,
                   ParamEst_cub_m3_monsoon6mo_AR2,
                   ParamEst_cub_m3_spring6mo_lag,
                   ParamEst_cub_m3_spring6mo_lag_AR1,
                   ParamEst_cub_m3_spring6mo_lag_AR2,
                   ParamEst_cub_m3_monsoon6mo_lag,
                   ParamEst_cub_m3_monsoon6mo_lag_AR1,
                   ParamEst_cub_m3_monsoon6mo_lag_AR2,
                   SE_cub_m3_spring6mo,
                   SE_cub_m3_spring6mo_AR1,
                   SE_cub_m3_spring6mo_AR2,
                   SE_cub_m3_monsoon6mo,
                   SE_cub_m3_monsoon6mo_AR1,
                   SE_cub_m3_monsoon6mo_AR2,
                   SE_cub_m3_spring6mo_lag,
                   SE_cub_m3_spring6mo_lag_AR1,
                   SE_cub_m3_spring6mo_lag_AR2,
                   SE_cub_m3_monsoon6mo_lag,
                   SE_cub_m3_monsoon6mo_lag_AR1,
                   SE_cub_m3_monsoon6mo_lag_AR2)
  
  # append results for the bee species to the output dataset
  beeCSF_output<-rbind(beeCSF_output,output_id)
  
}


##### CSFs: Chihuahuan Desert Grassland #####

# Create a new data frame of the original data
black_original <- black

# Create a data frame of just climate data
black_climate <- black[,5:8]

# Create a data frame of just the bee abundance matrix (descriptor variables removed)
speciesMatrix <- black[,9:224]

# Create a vector of species codes
speciesCodes <- colnames(speciesMatrix)

# Create a vector with a number for each species
number <-1:length(speciesCodes)


for (i in 1:length(speciesMatrix[1,])) {
  
  # save the species code for column i
  speciesCode <- speciesCodes[i]
  
  # create an object with the name of the ecosystem type
  ecosystem <-"G"
  
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
  
  
  # normality test for null model
  shapiro_null<-shapiro.test(resid(m_null))
  P_shapiro_null<-shapiro_null$p.value
  print(P_shapiro_null)
  
  
  # calculate delta AICc for each model
  
  # delta AICc for AR models
  dAICc_m1_spring6mo_AR1<-(AICc(m1_spring6mo_AR1))-min(AICc(m_null),AICc(m1_spring6mo),AICc(m2_spring6mo),AICc(m3_spring6mo), AICc(m1_monsoon6mo), AICc(m2_monsoon6mo),AICc(m3_monsoon6mo),AICc(m1_spring6mo_lag),AICc(m2_spring6mo_lag),AICc(m3_spring6mo_lag),AICc(m1_monsoon6mo_lag),AICc(m2_monsoon6mo_lag),AICc(m3_monsoon6mo_lag),AICc(m1_spring6mo_AR1),AICc(m1_spring6mo_AR2),AICc(m2_spring6mo_AR1),AICc(m2_spring6mo_AR2),AICc(m3_spring6mo_AR1),AICc(m3_spring6mo_AR2),AICc(m1_monsoon6mo_AR1),AICc(m1_monsoon6mo_AR2),AICc(m2_monsoon6mo_AR1),AICc(m2_monsoon6mo_AR2),AICc(m3_monsoon6mo_AR1),AICc(m3_monsoon6mo_AR2),AICc(m1_spring6mo_lag_AR1),AICc(m1_spring6mo_lag_AR2),AICc(m2_spring6mo_lag_AR1),AICc(m2_spring6mo_lag_AR2),AICc(m3_spring6mo_lag_AR1),AICc(m3_spring6mo_lag_AR2),AICc(m1_monsoon6mo_lag_AR1),AICc(m1_monsoon6mo_lag_AR2),AICc(m2_monsoon6mo_lag_AR1),AICc(m2_monsoon6mo_lag_AR2),AICc(m3_monsoon6mo_lag_AR1),AICc(m3_monsoon6mo_lag_AR2))
  
  dAICc_m1_spring6mo_AR2<-(AICc(m1_spring6mo_AR2))-min(AICc(m_null),AICc(m1_spring6mo),AICc(m2_spring6mo),AICc(m3_spring6mo), AICc(m1_monsoon6mo), AICc(m2_monsoon6mo),AICc(m3_monsoon6mo),AICc(m1_spring6mo_lag),AICc(m2_spring6mo_lag),AICc(m3_spring6mo_lag),AICc(m1_monsoon6mo_lag),AICc(m2_monsoon6mo_lag),AICc(m3_monsoon6mo_lag),AICc(m1_spring6mo_AR1),AICc(m1_spring6mo_AR2),AICc(m2_spring6mo_AR1),AICc(m2_spring6mo_AR2),AICc(m3_spring6mo_AR1),AICc(m3_spring6mo_AR2),AICc(m1_monsoon6mo_AR1),AICc(m1_monsoon6mo_AR2),AICc(m2_monsoon6mo_AR1),AICc(m2_monsoon6mo_AR2),AICc(m3_monsoon6mo_AR1),AICc(m3_monsoon6mo_AR2),AICc(m1_spring6mo_lag_AR1),AICc(m1_spring6mo_lag_AR2),AICc(m2_spring6mo_lag_AR1),AICc(m2_spring6mo_lag_AR2),AICc(m3_spring6mo_lag_AR1),AICc(m3_spring6mo_lag_AR2),AICc(m1_monsoon6mo_lag_AR1),AICc(m1_monsoon6mo_lag_AR2),AICc(m2_monsoon6mo_lag_AR1),AICc(m2_monsoon6mo_lag_AR2),AICc(m3_monsoon6mo_lag_AR1),AICc(m3_monsoon6mo_lag_AR2))
  
  dAICc_m2_spring6mo_AR1<-(AICc(m2_spring6mo_AR1))-min(AICc(m_null),AICc(m1_spring6mo),AICc(m2_spring6mo),AICc(m3_spring6mo), AICc(m1_monsoon6mo), AICc(m2_monsoon6mo),AICc(m3_monsoon6mo),AICc(m1_spring6mo_lag),AICc(m2_spring6mo_lag),AICc(m3_spring6mo_lag),AICc(m1_monsoon6mo_lag),AICc(m2_monsoon6mo_lag),AICc(m3_monsoon6mo_lag),AICc(m1_spring6mo_AR1),AICc(m1_spring6mo_AR2),AICc(m2_spring6mo_AR1),AICc(m2_spring6mo_AR2),AICc(m3_spring6mo_AR1),AICc(m3_spring6mo_AR2),AICc(m1_monsoon6mo_AR1),AICc(m1_monsoon6mo_AR2),AICc(m2_monsoon6mo_AR1),AICc(m2_monsoon6mo_AR2),AICc(m3_monsoon6mo_AR1),AICc(m3_monsoon6mo_AR2),AICc(m1_spring6mo_lag_AR1),AICc(m1_spring6mo_lag_AR2),AICc(m2_spring6mo_lag_AR1),AICc(m2_spring6mo_lag_AR2),AICc(m3_spring6mo_lag_AR1),AICc(m3_spring6mo_lag_AR2),AICc(m1_monsoon6mo_lag_AR1),AICc(m1_monsoon6mo_lag_AR2),AICc(m2_monsoon6mo_lag_AR1),AICc(m2_monsoon6mo_lag_AR2),AICc(m3_monsoon6mo_lag_AR1),AICc(m3_monsoon6mo_lag_AR2))
  
  dAICc_m2_spring6mo_AR2<-(AICc(m2_spring6mo_AR2))-min(AICc(m_null),AICc(m1_spring6mo),AICc(m2_spring6mo),AICc(m3_spring6mo), AICc(m1_monsoon6mo), AICc(m2_monsoon6mo),AICc(m3_monsoon6mo),AICc(m1_spring6mo_lag),AICc(m2_spring6mo_lag),AICc(m3_spring6mo_lag),AICc(m1_monsoon6mo_lag),AICc(m2_monsoon6mo_lag),AICc(m3_monsoon6mo_lag),AICc(m1_spring6mo_AR1),AICc(m1_spring6mo_AR2),AICc(m2_spring6mo_AR1),AICc(m2_spring6mo_AR2),AICc(m3_spring6mo_AR1),AICc(m3_spring6mo_AR2),AICc(m1_monsoon6mo_AR1),AICc(m1_monsoon6mo_AR2),AICc(m2_monsoon6mo_AR1),AICc(m2_monsoon6mo_AR2),AICc(m3_monsoon6mo_AR1),AICc(m3_monsoon6mo_AR2),AICc(m1_spring6mo_lag_AR1),AICc(m1_spring6mo_lag_AR2),AICc(m2_spring6mo_lag_AR1),AICc(m2_spring6mo_lag_AR2),AICc(m3_spring6mo_lag_AR1),AICc(m3_spring6mo_lag_AR2),AICc(m1_monsoon6mo_lag_AR1),AICc(m1_monsoon6mo_lag_AR2),AICc(m2_monsoon6mo_lag_AR1),AICc(m2_monsoon6mo_lag_AR2),AICc(m3_monsoon6mo_lag_AR1),AICc(m3_monsoon6mo_lag_AR2))
  
  dAICc_m3_spring6mo_AR1<-(AICc(m3_spring6mo_AR1))-min(AICc(m_null),AICc(m1_spring6mo),AICc(m2_spring6mo),AICc(m3_spring6mo), AICc(m1_monsoon6mo), AICc(m2_monsoon6mo),AICc(m3_monsoon6mo),AICc(m1_spring6mo_lag),AICc(m2_spring6mo_lag),AICc(m3_spring6mo_lag),AICc(m1_monsoon6mo_lag),AICc(m2_monsoon6mo_lag),AICc(m3_monsoon6mo_lag),AICc(m1_spring6mo_AR1),AICc(m1_spring6mo_AR2),AICc(m2_spring6mo_AR1),AICc(m2_spring6mo_AR2),AICc(m3_spring6mo_AR1),AICc(m3_spring6mo_AR2),AICc(m1_monsoon6mo_AR1),AICc(m1_monsoon6mo_AR2),AICc(m2_monsoon6mo_AR1),AICc(m2_monsoon6mo_AR2),AICc(m3_monsoon6mo_AR1),AICc(m3_monsoon6mo_AR2),AICc(m1_spring6mo_lag_AR1),AICc(m1_spring6mo_lag_AR2),AICc(m2_spring6mo_lag_AR1),AICc(m2_spring6mo_lag_AR2),AICc(m3_spring6mo_lag_AR1),AICc(m3_spring6mo_lag_AR2),AICc(m1_monsoon6mo_lag_AR1),AICc(m1_monsoon6mo_lag_AR2),AICc(m2_monsoon6mo_lag_AR1),AICc(m2_monsoon6mo_lag_AR2),AICc(m3_monsoon6mo_lag_AR1),AICc(m3_monsoon6mo_lag_AR2))
  
  dAICc_m3_spring6mo_AR2<-(AICc(m3_spring6mo_AR2))-min(AICc(m_null),AICc(m1_spring6mo),AICc(m2_spring6mo),AICc(m3_spring6mo), AICc(m1_monsoon6mo), AICc(m2_monsoon6mo),AICc(m3_monsoon6mo),AICc(m1_spring6mo_lag),AICc(m2_spring6mo_lag),AICc(m3_spring6mo_lag),AICc(m1_monsoon6mo_lag),AICc(m2_monsoon6mo_lag),AICc(m3_monsoon6mo_lag),AICc(m1_spring6mo_AR1),AICc(m1_spring6mo_AR2),AICc(m2_spring6mo_AR1),AICc(m2_spring6mo_AR2),AICc(m3_spring6mo_AR1),AICc(m3_spring6mo_AR2),AICc(m1_monsoon6mo_AR1),AICc(m1_monsoon6mo_AR2),AICc(m2_monsoon6mo_AR1),AICc(m2_monsoon6mo_AR2),AICc(m3_monsoon6mo_AR1),AICc(m3_monsoon6mo_AR2),AICc(m1_spring6mo_lag_AR1),AICc(m1_spring6mo_lag_AR2),AICc(m2_spring6mo_lag_AR1),AICc(m2_spring6mo_lag_AR2),AICc(m3_spring6mo_lag_AR1),AICc(m3_spring6mo_lag_AR2),AICc(m1_monsoon6mo_lag_AR1),AICc(m1_monsoon6mo_lag_AR2),AICc(m2_monsoon6mo_lag_AR1),AICc(m2_monsoon6mo_lag_AR2),AICc(m3_monsoon6mo_lag_AR1),AICc(m3_monsoon6mo_lag_AR2))
  
  dAICc_m1_monsoon6mo_AR1<-(AICc(m1_monsoon6mo_AR1))-min(AICc(m_null),AICc(m1_spring6mo),AICc(m2_spring6mo),AICc(m3_spring6mo), AICc(m1_monsoon6mo), AICc(m2_monsoon6mo),AICc(m3_monsoon6mo),AICc(m1_spring6mo_lag),AICc(m2_spring6mo_lag),AICc(m3_spring6mo_lag),AICc(m1_monsoon6mo_lag),AICc(m2_monsoon6mo_lag),AICc(m3_monsoon6mo_lag),AICc(m1_spring6mo_AR1),AICc(m1_spring6mo_AR2),AICc(m2_spring6mo_AR1),AICc(m2_spring6mo_AR2),AICc(m3_spring6mo_AR1),AICc(m3_spring6mo_AR2),AICc(m1_monsoon6mo_AR1),AICc(m1_monsoon6mo_AR2),AICc(m2_monsoon6mo_AR1),AICc(m2_monsoon6mo_AR2),AICc(m3_monsoon6mo_AR1),AICc(m3_monsoon6mo_AR2),AICc(m1_spring6mo_lag_AR1),AICc(m1_spring6mo_lag_AR2),AICc(m2_spring6mo_lag_AR1),AICc(m2_spring6mo_lag_AR2),AICc(m3_spring6mo_lag_AR1),AICc(m3_spring6mo_lag_AR2),AICc(m1_monsoon6mo_lag_AR1),AICc(m1_monsoon6mo_lag_AR2),AICc(m2_monsoon6mo_lag_AR1),AICc(m2_monsoon6mo_lag_AR2),AICc(m3_monsoon6mo_lag_AR1),AICc(m3_monsoon6mo_lag_AR2))
  
  dAICc_m1_monsoon6mo_AR2<-(AICc(m1_monsoon6mo_AR2))-min(AICc(m_null),AICc(m1_spring6mo),AICc(m2_spring6mo),AICc(m3_spring6mo), AICc(m1_monsoon6mo), AICc(m2_monsoon6mo),AICc(m3_monsoon6mo),AICc(m1_spring6mo_lag),AICc(m2_spring6mo_lag),AICc(m3_spring6mo_lag),AICc(m1_monsoon6mo_lag),AICc(m2_monsoon6mo_lag),AICc(m3_monsoon6mo_lag),AICc(m1_spring6mo_AR1),AICc(m1_spring6mo_AR2),AICc(m2_spring6mo_AR1),AICc(m2_spring6mo_AR2),AICc(m3_spring6mo_AR1),AICc(m3_spring6mo_AR2),AICc(m1_monsoon6mo_AR1),AICc(m1_monsoon6mo_AR2),AICc(m2_monsoon6mo_AR1),AICc(m2_monsoon6mo_AR2),AICc(m3_monsoon6mo_AR1),AICc(m3_monsoon6mo_AR2),AICc(m1_spring6mo_lag_AR1),AICc(m1_spring6mo_lag_AR2),AICc(m2_spring6mo_lag_AR1),AICc(m2_spring6mo_lag_AR2),AICc(m3_spring6mo_lag_AR1),AICc(m3_spring6mo_lag_AR2),AICc(m1_monsoon6mo_lag_AR1),AICc(m1_monsoon6mo_lag_AR2),AICc(m2_monsoon6mo_lag_AR1),AICc(m2_monsoon6mo_lag_AR2),AICc(m3_monsoon6mo_lag_AR1),AICc(m3_monsoon6mo_lag_AR2))
  
  dAICc_m2_monsoon6mo_AR1<-(AICc(m2_monsoon6mo_AR1))-min(AICc(m_null),AICc(m1_spring6mo),AICc(m2_spring6mo),AICc(m3_spring6mo), AICc(m1_monsoon6mo), AICc(m2_monsoon6mo),AICc(m3_monsoon6mo),AICc(m1_spring6mo_lag),AICc(m2_spring6mo_lag),AICc(m3_spring6mo_lag),AICc(m1_monsoon6mo_lag),AICc(m2_monsoon6mo_lag),AICc(m3_monsoon6mo_lag),AICc(m1_spring6mo_AR1),AICc(m1_spring6mo_AR2),AICc(m2_spring6mo_AR1),AICc(m2_spring6mo_AR2),AICc(m3_spring6mo_AR1),AICc(m3_spring6mo_AR2),AICc(m1_monsoon6mo_AR1),AICc(m1_monsoon6mo_AR2),AICc(m2_monsoon6mo_AR1),AICc(m2_monsoon6mo_AR2),AICc(m3_monsoon6mo_AR1),AICc(m3_monsoon6mo_AR2),AICc(m1_spring6mo_lag_AR1),AICc(m1_spring6mo_lag_AR2),AICc(m2_spring6mo_lag_AR1),AICc(m2_spring6mo_lag_AR2),AICc(m3_spring6mo_lag_AR1),AICc(m3_spring6mo_lag_AR2),AICc(m1_monsoon6mo_lag_AR1),AICc(m1_monsoon6mo_lag_AR2),AICc(m2_monsoon6mo_lag_AR1),AICc(m2_monsoon6mo_lag_AR2),AICc(m3_monsoon6mo_lag_AR1),AICc(m3_monsoon6mo_lag_AR2))
  
  dAICc_m2_monsoon6mo_AR2<-(AICc(m2_monsoon6mo_AR2))-min(AICc(m_null),AICc(m1_spring6mo),AICc(m2_spring6mo),AICc(m3_spring6mo), AICc(m1_monsoon6mo), AICc(m2_monsoon6mo),AICc(m3_monsoon6mo),AICc(m1_spring6mo_lag),AICc(m2_spring6mo_lag),AICc(m3_spring6mo_lag),AICc(m1_monsoon6mo_lag),AICc(m2_monsoon6mo_lag),AICc(m3_monsoon6mo_lag),AICc(m1_spring6mo_AR1),AICc(m1_spring6mo_AR2),AICc(m2_spring6mo_AR1),AICc(m2_spring6mo_AR2),AICc(m3_spring6mo_AR1),AICc(m3_spring6mo_AR2),AICc(m1_monsoon6mo_AR1),AICc(m1_monsoon6mo_AR2),AICc(m2_monsoon6mo_AR1),AICc(m2_monsoon6mo_AR2),AICc(m3_monsoon6mo_AR1),AICc(m3_monsoon6mo_AR2),AICc(m1_spring6mo_lag_AR1),AICc(m1_spring6mo_lag_AR2),AICc(m2_spring6mo_lag_AR1),AICc(m2_spring6mo_lag_AR2),AICc(m3_spring6mo_lag_AR1),AICc(m3_spring6mo_lag_AR2),AICc(m1_monsoon6mo_lag_AR1),AICc(m1_monsoon6mo_lag_AR2),AICc(m2_monsoon6mo_lag_AR1),AICc(m2_monsoon6mo_lag_AR2),AICc(m3_monsoon6mo_lag_AR1),AICc(m3_monsoon6mo_lag_AR2))
  
  dAICc_m3_monsoon6mo_AR1<-(AICc(m3_monsoon6mo_AR1))-min(AICc(m_null),AICc(m1_spring6mo),AICc(m2_spring6mo),AICc(m3_spring6mo), AICc(m1_monsoon6mo), AICc(m2_monsoon6mo),AICc(m3_monsoon6mo),AICc(m1_spring6mo_lag),AICc(m2_spring6mo_lag),AICc(m3_spring6mo_lag),AICc(m1_monsoon6mo_lag),AICc(m2_monsoon6mo_lag),AICc(m3_monsoon6mo_lag),AICc(m1_spring6mo_AR1),AICc(m1_spring6mo_AR2),AICc(m2_spring6mo_AR1),AICc(m2_spring6mo_AR2),AICc(m3_spring6mo_AR1),AICc(m3_spring6mo_AR2),AICc(m1_monsoon6mo_AR1),AICc(m1_monsoon6mo_AR2),AICc(m2_monsoon6mo_AR1),AICc(m2_monsoon6mo_AR2),AICc(m3_monsoon6mo_AR1),AICc(m3_monsoon6mo_AR2),AICc(m1_spring6mo_lag_AR1),AICc(m1_spring6mo_lag_AR2),AICc(m2_spring6mo_lag_AR1),AICc(m2_spring6mo_lag_AR2),AICc(m3_spring6mo_lag_AR1),AICc(m3_spring6mo_lag_AR2),AICc(m1_monsoon6mo_lag_AR1),AICc(m1_monsoon6mo_lag_AR2),AICc(m2_monsoon6mo_lag_AR1),AICc(m2_monsoon6mo_lag_AR2),AICc(m3_monsoon6mo_lag_AR1),AICc(m3_monsoon6mo_lag_AR2))
  
  dAICc_m3_monsoon6mo_AR2<-(AICc(m3_monsoon6mo_AR2))-min(AICc(m_null),AICc(m1_spring6mo),AICc(m2_spring6mo),AICc(m3_spring6mo), AICc(m1_monsoon6mo), AICc(m2_monsoon6mo),AICc(m3_monsoon6mo),AICc(m1_spring6mo_lag),AICc(m2_spring6mo_lag),AICc(m3_spring6mo_lag),AICc(m1_monsoon6mo_lag),AICc(m2_monsoon6mo_lag),AICc(m3_monsoon6mo_lag),AICc(m1_spring6mo_AR1),AICc(m1_spring6mo_AR2),AICc(m2_spring6mo_AR1),AICc(m2_spring6mo_AR2),AICc(m3_spring6mo_AR1),AICc(m3_spring6mo_AR2),AICc(m1_monsoon6mo_AR1),AICc(m1_monsoon6mo_AR2),AICc(m2_monsoon6mo_AR1),AICc(m2_monsoon6mo_AR2),AICc(m3_monsoon6mo_AR1),AICc(m3_monsoon6mo_AR2),AICc(m1_spring6mo_lag_AR1),AICc(m1_spring6mo_lag_AR2),AICc(m2_spring6mo_lag_AR1),AICc(m2_spring6mo_lag_AR2),AICc(m3_spring6mo_lag_AR1),AICc(m3_spring6mo_lag_AR2),AICc(m1_monsoon6mo_lag_AR1),AICc(m1_monsoon6mo_lag_AR2),AICc(m2_monsoon6mo_lag_AR1),AICc(m2_monsoon6mo_lag_AR2),AICc(m3_monsoon6mo_lag_AR1),AICc(m3_monsoon6mo_lag_AR2))
  
  dAICc_m1_spring6mo_lag_AR1<-(AICc(m1_spring6mo_lag_AR1))-min(AICc(m_null),AICc(m1_spring6mo),AICc(m2_spring6mo),AICc(m3_spring6mo), AICc(m1_monsoon6mo), AICc(m2_monsoon6mo),AICc(m3_monsoon6mo),AICc(m1_spring6mo_lag),AICc(m2_spring6mo_lag),AICc(m3_spring6mo_lag),AICc(m1_monsoon6mo_lag),AICc(m2_monsoon6mo_lag),AICc(m3_monsoon6mo_lag),AICc(m1_spring6mo_AR1),AICc(m1_spring6mo_AR2),AICc(m2_spring6mo_AR1),AICc(m2_spring6mo_AR2),AICc(m3_spring6mo_AR1),AICc(m3_spring6mo_AR2),AICc(m1_monsoon6mo_AR1),AICc(m1_monsoon6mo_AR2),AICc(m2_monsoon6mo_AR1),AICc(m2_monsoon6mo_AR2),AICc(m3_monsoon6mo_AR1),AICc(m3_monsoon6mo_AR2),AICc(m1_spring6mo_lag_AR1),AICc(m1_spring6mo_lag_AR2),AICc(m2_spring6mo_lag_AR1),AICc(m2_spring6mo_lag_AR2),AICc(m3_spring6mo_lag_AR1),AICc(m3_spring6mo_lag_AR2),AICc(m1_monsoon6mo_lag_AR1),AICc(m1_monsoon6mo_lag_AR2),AICc(m2_monsoon6mo_lag_AR1),AICc(m2_monsoon6mo_lag_AR2),AICc(m3_monsoon6mo_lag_AR1),AICc(m3_monsoon6mo_lag_AR2))
  
  dAICc_m1_spring6mo_lag_AR2<-(AICc(m1_spring6mo_lag_AR2))-min(AICc(m_null),AICc(m1_spring6mo),AICc(m2_spring6mo),AICc(m3_spring6mo), AICc(m1_monsoon6mo), AICc(m2_monsoon6mo),AICc(m3_monsoon6mo),AICc(m1_spring6mo_lag),AICc(m2_spring6mo_lag),AICc(m3_spring6mo_lag),AICc(m1_monsoon6mo_lag),AICc(m2_monsoon6mo_lag),AICc(m3_monsoon6mo_lag),AICc(m1_spring6mo_AR1),AICc(m1_spring6mo_AR2),AICc(m2_spring6mo_AR1),AICc(m2_spring6mo_AR2),AICc(m3_spring6mo_AR1),AICc(m3_spring6mo_AR2),AICc(m1_monsoon6mo_AR1),AICc(m1_monsoon6mo_AR2),AICc(m2_monsoon6mo_AR1),AICc(m2_monsoon6mo_AR2),AICc(m3_monsoon6mo_AR1),AICc(m3_monsoon6mo_AR2),AICc(m1_spring6mo_lag_AR1),AICc(m1_spring6mo_lag_AR2),AICc(m2_spring6mo_lag_AR1),AICc(m2_spring6mo_lag_AR2),AICc(m3_spring6mo_lag_AR1),AICc(m3_spring6mo_lag_AR2),AICc(m1_monsoon6mo_lag_AR1),AICc(m1_monsoon6mo_lag_AR2),AICc(m2_monsoon6mo_lag_AR1),AICc(m2_monsoon6mo_lag_AR2),AICc(m3_monsoon6mo_lag_AR1),AICc(m3_monsoon6mo_lag_AR2))
  
  dAICc_m2_spring6mo_lag_AR1<-(AICc(m2_spring6mo_lag_AR1))-min(AICc(m_null),AICc(m1_spring6mo),AICc(m2_spring6mo),AICc(m3_spring6mo), AICc(m1_monsoon6mo), AICc(m2_monsoon6mo),AICc(m3_monsoon6mo),AICc(m1_spring6mo_lag),AICc(m2_spring6mo_lag),AICc(m3_spring6mo_lag),AICc(m1_monsoon6mo_lag),AICc(m2_monsoon6mo_lag),AICc(m3_monsoon6mo_lag),AICc(m1_spring6mo_AR1),AICc(m1_spring6mo_AR2),AICc(m2_spring6mo_AR1),AICc(m2_spring6mo_AR2),AICc(m3_spring6mo_AR1),AICc(m3_spring6mo_AR2),AICc(m1_monsoon6mo_AR1),AICc(m1_monsoon6mo_AR2),AICc(m2_monsoon6mo_AR1),AICc(m2_monsoon6mo_AR2),AICc(m3_monsoon6mo_AR1),AICc(m3_monsoon6mo_AR2),AICc(m1_spring6mo_lag_AR1),AICc(m1_spring6mo_lag_AR2),AICc(m2_spring6mo_lag_AR1),AICc(m2_spring6mo_lag_AR2),AICc(m3_spring6mo_lag_AR1),AICc(m3_spring6mo_lag_AR2),AICc(m1_monsoon6mo_lag_AR1),AICc(m1_monsoon6mo_lag_AR2),AICc(m2_monsoon6mo_lag_AR1),AICc(m2_monsoon6mo_lag_AR2),AICc(m3_monsoon6mo_lag_AR1),AICc(m3_monsoon6mo_lag_AR2))
  
  dAICc_m2_spring6mo_lag_AR2<-(AICc(m2_spring6mo_lag_AR2))-min(AICc(m_null),AICc(m1_spring6mo),AICc(m2_spring6mo),AICc(m3_spring6mo), AICc(m1_monsoon6mo), AICc(m2_monsoon6mo),AICc(m3_monsoon6mo),AICc(m1_spring6mo_lag),AICc(m2_spring6mo_lag),AICc(m3_spring6mo_lag),AICc(m1_monsoon6mo_lag),AICc(m2_monsoon6mo_lag),AICc(m3_monsoon6mo_lag),AICc(m1_spring6mo_AR1),AICc(m1_spring6mo_AR2),AICc(m2_spring6mo_AR1),AICc(m2_spring6mo_AR2),AICc(m3_spring6mo_AR1),AICc(m3_spring6mo_AR2),AICc(m1_monsoon6mo_AR1),AICc(m1_monsoon6mo_AR2),AICc(m2_monsoon6mo_AR1),AICc(m2_monsoon6mo_AR2),AICc(m3_monsoon6mo_AR1),AICc(m3_monsoon6mo_AR2),AICc(m1_spring6mo_lag_AR1),AICc(m1_spring6mo_lag_AR2),AICc(m2_spring6mo_lag_AR1),AICc(m2_spring6mo_lag_AR2),AICc(m3_spring6mo_lag_AR1),AICc(m3_spring6mo_lag_AR2),AICc(m1_monsoon6mo_lag_AR1),AICc(m1_monsoon6mo_lag_AR2),AICc(m2_monsoon6mo_lag_AR1),AICc(m2_monsoon6mo_lag_AR2),AICc(m3_monsoon6mo_lag_AR1),AICc(m3_monsoon6mo_lag_AR2))
  
  dAICc_m3_spring6mo_lag_AR1<-(AICc(m3_spring6mo_lag_AR1))-min(AICc(m_null),AICc(m1_spring6mo),AICc(m2_spring6mo),AICc(m3_spring6mo), AICc(m1_monsoon6mo), AICc(m2_monsoon6mo),AICc(m3_monsoon6mo),AICc(m1_spring6mo_lag),AICc(m2_spring6mo_lag),AICc(m3_spring6mo_lag),AICc(m1_monsoon6mo_lag),AICc(m2_monsoon6mo_lag),AICc(m3_monsoon6mo_lag),AICc(m1_spring6mo_AR1),AICc(m1_spring6mo_AR2),AICc(m2_spring6mo_AR1),AICc(m2_spring6mo_AR2),AICc(m3_spring6mo_AR1),AICc(m3_spring6mo_AR2),AICc(m1_monsoon6mo_AR1),AICc(m1_monsoon6mo_AR2),AICc(m2_monsoon6mo_AR1),AICc(m2_monsoon6mo_AR2),AICc(m3_monsoon6mo_AR1),AICc(m3_monsoon6mo_AR2),AICc(m1_spring6mo_lag_AR1),AICc(m1_spring6mo_lag_AR2),AICc(m2_spring6mo_lag_AR1),AICc(m2_spring6mo_lag_AR2),AICc(m3_spring6mo_lag_AR1),AICc(m3_spring6mo_lag_AR2),AICc(m1_monsoon6mo_lag_AR1),AICc(m1_monsoon6mo_lag_AR2),AICc(m2_monsoon6mo_lag_AR1),AICc(m2_monsoon6mo_lag_AR2),AICc(m3_monsoon6mo_lag_AR1),AICc(m3_monsoon6mo_lag_AR2))
  
  dAICc_m3_spring6mo_lag_AR2<-(AICc(m3_spring6mo_lag_AR2))-min(AICc(m_null),AICc(m1_spring6mo),AICc(m2_spring6mo),AICc(m3_spring6mo), AICc(m1_monsoon6mo), AICc(m2_monsoon6mo),AICc(m3_monsoon6mo),AICc(m1_spring6mo_lag),AICc(m2_spring6mo_lag),AICc(m3_spring6mo_lag),AICc(m1_monsoon6mo_lag),AICc(m2_monsoon6mo_lag),AICc(m3_monsoon6mo_lag),AICc(m1_spring6mo_AR1),AICc(m1_spring6mo_AR2),AICc(m2_spring6mo_AR1),AICc(m2_spring6mo_AR2),AICc(m3_spring6mo_AR1),AICc(m3_spring6mo_AR2),AICc(m1_monsoon6mo_AR1),AICc(m1_monsoon6mo_AR2),AICc(m2_monsoon6mo_AR1),AICc(m2_monsoon6mo_AR2),AICc(m3_monsoon6mo_AR1),AICc(m3_monsoon6mo_AR2),AICc(m1_spring6mo_lag_AR1),AICc(m1_spring6mo_lag_AR2),AICc(m2_spring6mo_lag_AR1),AICc(m2_spring6mo_lag_AR2),AICc(m3_spring6mo_lag_AR1),AICc(m3_spring6mo_lag_AR2),AICc(m1_monsoon6mo_lag_AR1),AICc(m1_monsoon6mo_lag_AR2),AICc(m2_monsoon6mo_lag_AR1),AICc(m2_monsoon6mo_lag_AR2),AICc(m3_monsoon6mo_lag_AR1),AICc(m3_monsoon6mo_lag_AR2))
  
  dAICc_m1_monsoon6mo_lag_AR1<-(AICc(m1_monsoon6mo_lag_AR1))-min(AICc(m_null),AICc(m1_spring6mo),AICc(m2_spring6mo),AICc(m3_spring6mo), AICc(m1_monsoon6mo), AICc(m2_monsoon6mo),AICc(m3_monsoon6mo),AICc(m1_spring6mo_lag),AICc(m2_spring6mo_lag),AICc(m3_spring6mo_lag),AICc(m1_monsoon6mo_lag),AICc(m2_monsoon6mo_lag),AICc(m3_monsoon6mo_lag),AICc(m1_spring6mo_AR1),AICc(m1_spring6mo_AR2),AICc(m2_spring6mo_AR1),AICc(m2_spring6mo_AR2),AICc(m3_spring6mo_AR1),AICc(m3_spring6mo_AR2),AICc(m1_monsoon6mo_AR1),AICc(m1_monsoon6mo_AR2),AICc(m2_monsoon6mo_AR1),AICc(m2_monsoon6mo_AR2),AICc(m3_monsoon6mo_AR1),AICc(m3_monsoon6mo_AR2),AICc(m1_spring6mo_lag_AR1),AICc(m1_spring6mo_lag_AR2),AICc(m2_spring6mo_lag_AR1),AICc(m2_spring6mo_lag_AR2),AICc(m3_spring6mo_lag_AR1),AICc(m3_spring6mo_lag_AR2),AICc(m1_monsoon6mo_lag_AR1),AICc(m1_monsoon6mo_lag_AR2),AICc(m2_monsoon6mo_lag_AR1),AICc(m2_monsoon6mo_lag_AR2),AICc(m3_monsoon6mo_lag_AR1),AICc(m3_monsoon6mo_lag_AR2))
  
  dAICc_m1_monsoon6mo_lag_AR2<-(AICc(m1_monsoon6mo_lag_AR2))-min(AICc(m_null),AICc(m1_spring6mo),AICc(m2_spring6mo),AICc(m3_spring6mo), AICc(m1_monsoon6mo), AICc(m2_monsoon6mo),AICc(m3_monsoon6mo),AICc(m1_spring6mo_lag),AICc(m2_spring6mo_lag),AICc(m3_spring6mo_lag),AICc(m1_monsoon6mo_lag),AICc(m2_monsoon6mo_lag),AICc(m3_monsoon6mo_lag),AICc(m1_spring6mo_AR1),AICc(m1_spring6mo_AR2),AICc(m2_spring6mo_AR1),AICc(m2_spring6mo_AR2),AICc(m3_spring6mo_AR1),AICc(m3_spring6mo_AR2),AICc(m1_monsoon6mo_AR1),AICc(m1_monsoon6mo_AR2),AICc(m2_monsoon6mo_AR1),AICc(m2_monsoon6mo_AR2),AICc(m3_monsoon6mo_AR1),AICc(m3_monsoon6mo_AR2),AICc(m1_spring6mo_lag_AR1),AICc(m1_spring6mo_lag_AR2),AICc(m2_spring6mo_lag_AR1),AICc(m2_spring6mo_lag_AR2),AICc(m3_spring6mo_lag_AR1),AICc(m3_spring6mo_lag_AR2),AICc(m1_monsoon6mo_lag_AR1),AICc(m1_monsoon6mo_lag_AR2),AICc(m2_monsoon6mo_lag_AR1),AICc(m2_monsoon6mo_lag_AR2),AICc(m3_monsoon6mo_lag_AR1),AICc(m3_monsoon6mo_lag_AR2))
  
  dAICc_m2_monsoon6mo_lag_AR1<-(AICc(m2_monsoon6mo_lag_AR1))-min(AICc(m_null),AICc(m1_spring6mo),AICc(m2_spring6mo),AICc(m3_spring6mo), AICc(m1_monsoon6mo), AICc(m2_monsoon6mo),AICc(m3_monsoon6mo),AICc(m1_spring6mo_lag),AICc(m2_spring6mo_lag),AICc(m3_spring6mo_lag),AICc(m1_monsoon6mo_lag),AICc(m2_monsoon6mo_lag),AICc(m3_monsoon6mo_lag),AICc(m1_spring6mo_AR1),AICc(m1_spring6mo_AR2),AICc(m2_spring6mo_AR1),AICc(m2_spring6mo_AR2),AICc(m3_spring6mo_AR1),AICc(m3_spring6mo_AR2),AICc(m1_monsoon6mo_AR1),AICc(m1_monsoon6mo_AR2),AICc(m2_monsoon6mo_AR1),AICc(m2_monsoon6mo_AR2),AICc(m3_monsoon6mo_AR1),AICc(m3_monsoon6mo_AR2),AICc(m1_spring6mo_lag_AR1),AICc(m1_spring6mo_lag_AR2),AICc(m2_spring6mo_lag_AR1),AICc(m2_spring6mo_lag_AR2),AICc(m3_spring6mo_lag_AR1),AICc(m3_spring6mo_lag_AR2),AICc(m1_monsoon6mo_lag_AR1),AICc(m1_monsoon6mo_lag_AR2),AICc(m2_monsoon6mo_lag_AR1),AICc(m2_monsoon6mo_lag_AR2),AICc(m3_monsoon6mo_lag_AR1),AICc(m3_monsoon6mo_lag_AR2))
  
  dAICc_m2_monsoon6mo_lag_AR2<-(AICc(m2_monsoon6mo_lag_AR2))-min(AICc(m_null),AICc(m1_spring6mo),AICc(m2_spring6mo),AICc(m3_spring6mo), AICc(m1_monsoon6mo), AICc(m2_monsoon6mo),AICc(m3_monsoon6mo),AICc(m1_spring6mo_lag),AICc(m2_spring6mo_lag),AICc(m3_spring6mo_lag),AICc(m1_monsoon6mo_lag),AICc(m2_monsoon6mo_lag),AICc(m3_monsoon6mo_lag),AICc(m1_spring6mo_AR1),AICc(m1_spring6mo_AR2),AICc(m2_spring6mo_AR1),AICc(m2_spring6mo_AR2),AICc(m3_spring6mo_AR1),AICc(m3_spring6mo_AR2),AICc(m1_monsoon6mo_AR1),AICc(m1_monsoon6mo_AR2),AICc(m2_monsoon6mo_AR1),AICc(m2_monsoon6mo_AR2),AICc(m3_monsoon6mo_AR1),AICc(m3_monsoon6mo_AR2),AICc(m1_spring6mo_lag_AR1),AICc(m1_spring6mo_lag_AR2),AICc(m2_spring6mo_lag_AR1),AICc(m2_spring6mo_lag_AR2),AICc(m3_spring6mo_lag_AR1),AICc(m3_spring6mo_lag_AR2),AICc(m1_monsoon6mo_lag_AR1),AICc(m1_monsoon6mo_lag_AR2),AICc(m2_monsoon6mo_lag_AR1),AICc(m2_monsoon6mo_lag_AR2),AICc(m3_monsoon6mo_lag_AR1),AICc(m3_monsoon6mo_lag_AR2))
  
  dAICc_m3_monsoon6mo_lag_AR1<-(AICc(m3_monsoon6mo_lag_AR1))-min(AICc(m_null),AICc(m1_spring6mo),AICc(m2_spring6mo),AICc(m3_spring6mo), AICc(m1_monsoon6mo), AICc(m2_monsoon6mo),AICc(m3_monsoon6mo),AICc(m1_spring6mo_lag),AICc(m2_spring6mo_lag),AICc(m3_spring6mo_lag),AICc(m1_monsoon6mo_lag),AICc(m2_monsoon6mo_lag),AICc(m3_monsoon6mo_lag),AICc(m1_spring6mo_AR1),AICc(m1_spring6mo_AR2),AICc(m2_spring6mo_AR1),AICc(m2_spring6mo_AR2),AICc(m3_spring6mo_AR1),AICc(m3_spring6mo_AR2),AICc(m1_monsoon6mo_AR1),AICc(m1_monsoon6mo_AR2),AICc(m2_monsoon6mo_AR1),AICc(m2_monsoon6mo_AR2),AICc(m3_monsoon6mo_AR1),AICc(m3_monsoon6mo_AR2),AICc(m1_spring6mo_lag_AR1),AICc(m1_spring6mo_lag_AR2),AICc(m2_spring6mo_lag_AR1),AICc(m2_spring6mo_lag_AR2),AICc(m3_spring6mo_lag_AR1),AICc(m3_spring6mo_lag_AR2),AICc(m1_monsoon6mo_lag_AR1),AICc(m1_monsoon6mo_lag_AR2),AICc(m2_monsoon6mo_lag_AR1),AICc(m2_monsoon6mo_lag_AR2),AICc(m3_monsoon6mo_lag_AR1),AICc(m3_monsoon6mo_lag_AR2))
  
  dAICc_m3_monsoon6mo_lag_AR2<-(AICc(m3_monsoon6mo_lag_AR2))-min(AICc(m_null),AICc(m1_spring6mo),AICc(m2_spring6mo),AICc(m3_spring6mo), AICc(m1_monsoon6mo), AICc(m2_monsoon6mo),AICc(m3_monsoon6mo),AICc(m1_spring6mo_lag),AICc(m2_spring6mo_lag),AICc(m3_spring6mo_lag),AICc(m1_monsoon6mo_lag),AICc(m2_monsoon6mo_lag),AICc(m3_monsoon6mo_lag),AICc(m1_spring6mo_AR1),AICc(m1_spring6mo_AR2),AICc(m2_spring6mo_AR1),AICc(m2_spring6mo_AR2),AICc(m3_spring6mo_AR1),AICc(m3_spring6mo_AR2),AICc(m1_monsoon6mo_AR1),AICc(m1_monsoon6mo_AR2),AICc(m2_monsoon6mo_AR1),AICc(m2_monsoon6mo_AR2),AICc(m3_monsoon6mo_AR1),AICc(m3_monsoon6mo_AR2),AICc(m1_spring6mo_lag_AR1),AICc(m1_spring6mo_lag_AR2),AICc(m2_spring6mo_lag_AR1),AICc(m2_spring6mo_lag_AR2),AICc(m3_spring6mo_lag_AR1),AICc(m3_spring6mo_lag_AR2),AICc(m1_monsoon6mo_lag_AR1),AICc(m1_monsoon6mo_lag_AR2),AICc(m2_monsoon6mo_lag_AR1),AICc(m2_monsoon6mo_lag_AR2),AICc(m3_monsoon6mo_lag_AR1),AICc(m3_monsoon6mo_lag_AR2))
  
  # delta AICc for other models
  dAICc_null<-(AICc(m_null))-min(AICc(m_null),AICc(m1_spring6mo),AICc(m2_spring6mo),AICc(m3_spring6mo), AICc(m1_monsoon6mo), AICc(m2_monsoon6mo),AICc(m3_monsoon6mo),AICc(m1_spring6mo_lag),AICc(m2_spring6mo_lag),AICc(m3_spring6mo_lag),AICc(m1_monsoon6mo_lag),AICc(m2_monsoon6mo_lag),AICc(m3_monsoon6mo_lag),AICc(m1_spring6mo_AR1),AICc(m1_spring6mo_AR2),AICc(m2_spring6mo_AR1),AICc(m2_spring6mo_AR2),AICc(m3_spring6mo_AR1),AICc(m3_spring6mo_AR2),AICc(m1_monsoon6mo_AR1),AICc(m1_monsoon6mo_AR2),AICc(m2_monsoon6mo_AR1),AICc(m2_monsoon6mo_AR2),AICc(m3_monsoon6mo_AR1),AICc(m3_monsoon6mo_AR2),AICc(m1_spring6mo_lag_AR1),AICc(m1_spring6mo_lag_AR2),AICc(m2_spring6mo_lag_AR1),AICc(m2_spring6mo_lag_AR2),AICc(m3_spring6mo_lag_AR1),AICc(m3_spring6mo_lag_AR2),AICc(m1_monsoon6mo_lag_AR1),AICc(m1_monsoon6mo_lag_AR2),AICc(m2_monsoon6mo_lag_AR1),AICc(m2_monsoon6mo_lag_AR2),AICc(m3_monsoon6mo_lag_AR1),AICc(m3_monsoon6mo_lag_AR2))
  
  dAICc_m1_spring6mo<-(AICc(m1_spring6mo))-min(AICc(m_null),AICc(m1_spring6mo),AICc(m2_spring6mo),AICc(m3_spring6mo), AICc(m1_monsoon6mo), AICc(m2_monsoon6mo),AICc(m3_monsoon6mo),AICc(m1_spring6mo_lag),AICc(m2_spring6mo_lag),AICc(m3_spring6mo_lag),AICc(m1_monsoon6mo_lag),AICc(m2_monsoon6mo_lag),AICc(m3_monsoon6mo_lag),AICc(m1_spring6mo_AR1),AICc(m1_spring6mo_AR2),AICc(m2_spring6mo_AR1),AICc(m2_spring6mo_AR2),AICc(m3_spring6mo_AR1),AICc(m3_spring6mo_AR2),AICc(m1_monsoon6mo_AR1),AICc(m1_monsoon6mo_AR2),AICc(m2_monsoon6mo_AR1),AICc(m2_monsoon6mo_AR2),AICc(m3_monsoon6mo_AR1),AICc(m3_monsoon6mo_AR2),AICc(m1_spring6mo_lag_AR1),AICc(m1_spring6mo_lag_AR2),AICc(m2_spring6mo_lag_AR1),AICc(m2_spring6mo_lag_AR2),AICc(m3_spring6mo_lag_AR1),AICc(m3_spring6mo_lag_AR2),AICc(m1_monsoon6mo_lag_AR1),AICc(m1_monsoon6mo_lag_AR2),AICc(m2_monsoon6mo_lag_AR1),AICc(m2_monsoon6mo_lag_AR2),AICc(m3_monsoon6mo_lag_AR1),AICc(m3_monsoon6mo_lag_AR2))
  
  dAICc_m2_spring6mo<-(AICc(m2_spring6mo))-min(AICc(m_null),AICc(m1_spring6mo),AICc(m2_spring6mo),AICc(m3_spring6mo), AICc(m1_monsoon6mo), AICc(m2_monsoon6mo),AICc(m3_monsoon6mo),AICc(m1_spring6mo_lag),AICc(m2_spring6mo_lag),AICc(m3_spring6mo_lag),AICc(m1_monsoon6mo_lag),AICc(m2_monsoon6mo_lag),AICc(m3_monsoon6mo_lag),AICc(m1_spring6mo_AR1),AICc(m1_spring6mo_AR2),AICc(m2_spring6mo_AR1),AICc(m2_spring6mo_AR2),AICc(m3_spring6mo_AR1),AICc(m3_spring6mo_AR2),AICc(m1_monsoon6mo_AR1),AICc(m1_monsoon6mo_AR2),AICc(m2_monsoon6mo_AR1),AICc(m2_monsoon6mo_AR2),AICc(m3_monsoon6mo_AR1),AICc(m3_monsoon6mo_AR2),AICc(m1_spring6mo_lag_AR1),AICc(m1_spring6mo_lag_AR2),AICc(m2_spring6mo_lag_AR1),AICc(m2_spring6mo_lag_AR2),AICc(m3_spring6mo_lag_AR1),AICc(m3_spring6mo_lag_AR2),AICc(m1_monsoon6mo_lag_AR1),AICc(m1_monsoon6mo_lag_AR2),AICc(m2_monsoon6mo_lag_AR1),AICc(m2_monsoon6mo_lag_AR2),AICc(m3_monsoon6mo_lag_AR1),AICc(m3_monsoon6mo_lag_AR2))
  
  dAICc_m3_spring6mo<-(AICc(m3_spring6mo))-min(AICc(m_null),AICc(m1_spring6mo),AICc(m2_spring6mo),AICc(m3_spring6mo), AICc(m1_monsoon6mo), AICc(m2_monsoon6mo),AICc(m3_monsoon6mo),AICc(m1_spring6mo_lag),AICc(m2_spring6mo_lag),AICc(m3_spring6mo_lag),AICc(m1_monsoon6mo_lag),AICc(m2_monsoon6mo_lag),AICc(m3_monsoon6mo_lag),AICc(m1_spring6mo_AR1),AICc(m1_spring6mo_AR2),AICc(m2_spring6mo_AR1),AICc(m2_spring6mo_AR2),AICc(m3_spring6mo_AR1),AICc(m3_spring6mo_AR2),AICc(m1_monsoon6mo_AR1),AICc(m1_monsoon6mo_AR2),AICc(m2_monsoon6mo_AR1),AICc(m2_monsoon6mo_AR2),AICc(m3_monsoon6mo_AR1),AICc(m3_monsoon6mo_AR2),AICc(m1_spring6mo_lag_AR1),AICc(m1_spring6mo_lag_AR2),AICc(m2_spring6mo_lag_AR1),AICc(m2_spring6mo_lag_AR2),AICc(m3_spring6mo_lag_AR1),AICc(m3_spring6mo_lag_AR2),AICc(m1_monsoon6mo_lag_AR1),AICc(m1_monsoon6mo_lag_AR2),AICc(m2_monsoon6mo_lag_AR1),AICc(m2_monsoon6mo_lag_AR2),AICc(m3_monsoon6mo_lag_AR1),AICc(m3_monsoon6mo_lag_AR2))
  
  dAICc_m1_monsoon6mo<-(AICc(m1_monsoon6mo))-min(AICc(m_null),AICc(m1_spring6mo),AICc(m2_spring6mo),AICc(m3_spring6mo), AICc(m1_monsoon6mo), AICc(m2_monsoon6mo),AICc(m3_monsoon6mo),AICc(m1_spring6mo_lag),AICc(m2_spring6mo_lag),AICc(m3_spring6mo_lag),AICc(m1_monsoon6mo_lag),AICc(m2_monsoon6mo_lag),AICc(m3_monsoon6mo_lag),AICc(m1_spring6mo_AR1),AICc(m1_spring6mo_AR2),AICc(m2_spring6mo_AR1),AICc(m2_spring6mo_AR2),AICc(m3_spring6mo_AR1),AICc(m3_spring6mo_AR2),AICc(m1_monsoon6mo_AR1),AICc(m1_monsoon6mo_AR2),AICc(m2_monsoon6mo_AR1),AICc(m2_monsoon6mo_AR2),AICc(m3_monsoon6mo_AR1),AICc(m3_monsoon6mo_AR2),AICc(m1_spring6mo_lag_AR1),AICc(m1_spring6mo_lag_AR2),AICc(m2_spring6mo_lag_AR1),AICc(m2_spring6mo_lag_AR2),AICc(m3_spring6mo_lag_AR1),AICc(m3_spring6mo_lag_AR2),AICc(m1_monsoon6mo_lag_AR1),AICc(m1_monsoon6mo_lag_AR2),AICc(m2_monsoon6mo_lag_AR1),AICc(m2_monsoon6mo_lag_AR2),AICc(m3_monsoon6mo_lag_AR1),AICc(m3_monsoon6mo_lag_AR2))
  
  dAICc_m2_monsoon6mo<-(AICc(m2_monsoon6mo))-min(AICc(m_null),AICc(m1_spring6mo),AICc(m2_spring6mo),AICc(m3_spring6mo), AICc(m1_monsoon6mo), AICc(m2_monsoon6mo),AICc(m3_monsoon6mo),AICc(m1_spring6mo_lag),AICc(m2_spring6mo_lag),AICc(m3_spring6mo_lag),AICc(m1_monsoon6mo_lag),AICc(m2_monsoon6mo_lag),AICc(m3_monsoon6mo_lag),AICc(m1_spring6mo_AR1),AICc(m1_spring6mo_AR2),AICc(m2_spring6mo_AR1),AICc(m2_spring6mo_AR2),AICc(m3_spring6mo_AR1),AICc(m3_spring6mo_AR2),AICc(m1_monsoon6mo_AR1),AICc(m1_monsoon6mo_AR2),AICc(m2_monsoon6mo_AR1),AICc(m2_monsoon6mo_AR2),AICc(m3_monsoon6mo_AR1),AICc(m3_monsoon6mo_AR2),AICc(m1_spring6mo_lag_AR1),AICc(m1_spring6mo_lag_AR2),AICc(m2_spring6mo_lag_AR1),AICc(m2_spring6mo_lag_AR2),AICc(m3_spring6mo_lag_AR1),AICc(m3_spring6mo_lag_AR2),AICc(m1_monsoon6mo_lag_AR1),AICc(m1_monsoon6mo_lag_AR2),AICc(m2_monsoon6mo_lag_AR1),AICc(m2_monsoon6mo_lag_AR2),AICc(m3_monsoon6mo_lag_AR1),AICc(m3_monsoon6mo_lag_AR2))
  
  dAICc_m3_monsoon6mo<-(AICc(m3_monsoon6mo))-min(AICc(m_null),AICc(m1_spring6mo),AICc(m2_spring6mo),AICc(m3_spring6mo), AICc(m1_monsoon6mo), AICc(m2_monsoon6mo),AICc(m3_monsoon6mo),AICc(m1_spring6mo_lag),AICc(m2_spring6mo_lag),AICc(m3_spring6mo_lag),AICc(m1_monsoon6mo_lag),AICc(m2_monsoon6mo_lag),AICc(m3_monsoon6mo_lag),AICc(m1_spring6mo_AR1),AICc(m1_spring6mo_AR2),AICc(m2_spring6mo_AR1),AICc(m2_spring6mo_AR2),AICc(m3_spring6mo_AR1),AICc(m3_spring6mo_AR2),AICc(m1_monsoon6mo_AR1),AICc(m1_monsoon6mo_AR2),AICc(m2_monsoon6mo_AR1),AICc(m2_monsoon6mo_AR2),AICc(m3_monsoon6mo_AR1),AICc(m3_monsoon6mo_AR2),AICc(m1_spring6mo_lag_AR1),AICc(m1_spring6mo_lag_AR2),AICc(m2_spring6mo_lag_AR1),AICc(m2_spring6mo_lag_AR2),AICc(m3_spring6mo_lag_AR1),AICc(m3_spring6mo_lag_AR2),AICc(m1_monsoon6mo_lag_AR1),AICc(m1_monsoon6mo_lag_AR2),AICc(m2_monsoon6mo_lag_AR1),AICc(m2_monsoon6mo_lag_AR2),AICc(m3_monsoon6mo_lag_AR1),AICc(m3_monsoon6mo_lag_AR2))
  
  dAICc_m1_spring6mo_lag<-(AICc(m1_spring6mo_lag))-min(AICc(m_null),AICc(m1_spring6mo),AICc(m2_spring6mo),AICc(m3_spring6mo), AICc(m1_monsoon6mo), AICc(m2_monsoon6mo),AICc(m3_monsoon6mo),AICc(m1_spring6mo_lag),AICc(m2_spring6mo_lag),AICc(m3_spring6mo_lag),AICc(m1_monsoon6mo_lag),AICc(m2_monsoon6mo_lag),AICc(m3_monsoon6mo_lag),AICc(m1_spring6mo_AR1),AICc(m1_spring6mo_AR2),AICc(m2_spring6mo_AR1),AICc(m2_spring6mo_AR2),AICc(m3_spring6mo_AR1),AICc(m3_spring6mo_AR2),AICc(m1_monsoon6mo_AR1),AICc(m1_monsoon6mo_AR2),AICc(m2_monsoon6mo_AR1),AICc(m2_monsoon6mo_AR2),AICc(m3_monsoon6mo_AR1),AICc(m3_monsoon6mo_AR2),AICc(m1_spring6mo_lag_AR1),AICc(m1_spring6mo_lag_AR2),AICc(m2_spring6mo_lag_AR1),AICc(m2_spring6mo_lag_AR2),AICc(m3_spring6mo_lag_AR1),AICc(m3_spring6mo_lag_AR2),AICc(m1_monsoon6mo_lag_AR1),AICc(m1_monsoon6mo_lag_AR2),AICc(m2_monsoon6mo_lag_AR1),AICc(m2_monsoon6mo_lag_AR2),AICc(m3_monsoon6mo_lag_AR1),AICc(m3_monsoon6mo_lag_AR2))
  
  dAICc_m2_spring6mo_lag<-(AICc(m2_spring6mo_lag))-min(AICc(m_null),AICc(m1_spring6mo),AICc(m2_spring6mo),AICc(m3_spring6mo), AICc(m1_monsoon6mo), AICc(m2_monsoon6mo),AICc(m3_monsoon6mo),AICc(m1_spring6mo_lag),AICc(m2_spring6mo_lag),AICc(m3_spring6mo_lag),AICc(m1_monsoon6mo_lag),AICc(m2_monsoon6mo_lag),AICc(m3_monsoon6mo_lag),AICc(m1_spring6mo_AR1),AICc(m1_spring6mo_AR2),AICc(m2_spring6mo_AR1),AICc(m2_spring6mo_AR2),AICc(m3_spring6mo_AR1),AICc(m3_spring6mo_AR2),AICc(m1_monsoon6mo_AR1),AICc(m1_monsoon6mo_AR2),AICc(m2_monsoon6mo_AR1),AICc(m2_monsoon6mo_AR2),AICc(m3_monsoon6mo_AR1),AICc(m3_monsoon6mo_AR2),AICc(m1_spring6mo_lag_AR1),AICc(m1_spring6mo_lag_AR2),AICc(m2_spring6mo_lag_AR1),AICc(m2_spring6mo_lag_AR2),AICc(m3_spring6mo_lag_AR1),AICc(m3_spring6mo_lag_AR2),AICc(m1_monsoon6mo_lag_AR1),AICc(m1_monsoon6mo_lag_AR2),AICc(m2_monsoon6mo_lag_AR1),AICc(m2_monsoon6mo_lag_AR2),AICc(m3_monsoon6mo_lag_AR1),AICc(m3_monsoon6mo_lag_AR2))
  
  dAICc_m3_spring6mo_lag<-(AICc(m3_spring6mo_lag))-min(AICc(m_null),AICc(m1_spring6mo),AICc(m2_spring6mo),AICc(m3_spring6mo), AICc(m1_monsoon6mo), AICc(m2_monsoon6mo),AICc(m3_monsoon6mo),AICc(m1_spring6mo_lag),AICc(m2_spring6mo_lag),AICc(m3_spring6mo_lag),AICc(m1_monsoon6mo_lag),AICc(m2_monsoon6mo_lag),AICc(m3_monsoon6mo_lag),AICc(m1_spring6mo_AR1),AICc(m1_spring6mo_AR2),AICc(m2_spring6mo_AR1),AICc(m2_spring6mo_AR2),AICc(m3_spring6mo_AR1),AICc(m3_spring6mo_AR2),AICc(m1_monsoon6mo_AR1),AICc(m1_monsoon6mo_AR2),AICc(m2_monsoon6mo_AR1),AICc(m2_monsoon6mo_AR2),AICc(m3_monsoon6mo_AR1),AICc(m3_monsoon6mo_AR2),AICc(m1_spring6mo_lag_AR1),AICc(m1_spring6mo_lag_AR2),AICc(m2_spring6mo_lag_AR1),AICc(m2_spring6mo_lag_AR2),AICc(m3_spring6mo_lag_AR1),AICc(m3_spring6mo_lag_AR2),AICc(m1_monsoon6mo_lag_AR1),AICc(m1_monsoon6mo_lag_AR2),AICc(m2_monsoon6mo_lag_AR1),AICc(m2_monsoon6mo_lag_AR2),AICc(m3_monsoon6mo_lag_AR1),AICc(m3_monsoon6mo_lag_AR2))
  
  dAICc_m1_monsoon6mo_lag<-(AICc(m1_monsoon6mo_lag))-min(AICc(m_null),AICc(m1_spring6mo),AICc(m2_spring6mo),AICc(m3_spring6mo), AICc(m1_monsoon6mo), AICc(m2_monsoon6mo),AICc(m3_monsoon6mo),AICc(m1_spring6mo_lag),AICc(m2_spring6mo_lag),AICc(m3_spring6mo_lag),AICc(m1_monsoon6mo_lag),AICc(m2_monsoon6mo_lag),AICc(m3_monsoon6mo_lag),AICc(m1_spring6mo_AR1),AICc(m1_spring6mo_AR2),AICc(m2_spring6mo_AR1),AICc(m2_spring6mo_AR2),AICc(m3_spring6mo_AR1),AICc(m3_spring6mo_AR2),AICc(m1_monsoon6mo_AR1),AICc(m1_monsoon6mo_AR2),AICc(m2_monsoon6mo_AR1),AICc(m2_monsoon6mo_AR2),AICc(m3_monsoon6mo_AR1),AICc(m3_monsoon6mo_AR2),AICc(m1_spring6mo_lag_AR1),AICc(m1_spring6mo_lag_AR2),AICc(m2_spring6mo_lag_AR1),AICc(m2_spring6mo_lag_AR2),AICc(m3_spring6mo_lag_AR1),AICc(m3_spring6mo_lag_AR2),AICc(m1_monsoon6mo_lag_AR1),AICc(m1_monsoon6mo_lag_AR2),AICc(m2_monsoon6mo_lag_AR1),AICc(m2_monsoon6mo_lag_AR2),AICc(m3_monsoon6mo_lag_AR1),AICc(m3_monsoon6mo_lag_AR2))
  
  dAICc_m2_monsoon6mo_lag<-(AICc(m2_monsoon6mo_lag))-min(AICc(m_null),AICc(m1_spring6mo),AICc(m2_spring6mo),AICc(m3_spring6mo), AICc(m1_monsoon6mo), AICc(m2_monsoon6mo),AICc(m3_monsoon6mo),AICc(m1_spring6mo_lag),AICc(m2_spring6mo_lag),AICc(m3_spring6mo_lag),AICc(m1_monsoon6mo_lag),AICc(m2_monsoon6mo_lag),AICc(m3_monsoon6mo_lag),AICc(m1_spring6mo_AR1),AICc(m1_spring6mo_AR2),AICc(m2_spring6mo_AR1),AICc(m2_spring6mo_AR2),AICc(m3_spring6mo_AR1),AICc(m3_spring6mo_AR2),AICc(m1_monsoon6mo_AR1),AICc(m1_monsoon6mo_AR2),AICc(m2_monsoon6mo_AR1),AICc(m2_monsoon6mo_AR2),AICc(m3_monsoon6mo_AR1),AICc(m3_monsoon6mo_AR2),AICc(m1_spring6mo_lag_AR1),AICc(m1_spring6mo_lag_AR2),AICc(m2_spring6mo_lag_AR1),AICc(m2_spring6mo_lag_AR2),AICc(m3_spring6mo_lag_AR1),AICc(m3_spring6mo_lag_AR2),AICc(m1_monsoon6mo_lag_AR1),AICc(m1_monsoon6mo_lag_AR2),AICc(m2_monsoon6mo_lag_AR1),AICc(m2_monsoon6mo_lag_AR2),AICc(m3_monsoon6mo_lag_AR1),AICc(m3_monsoon6mo_lag_AR2))
  
  dAICc_m3_monsoon6mo_lag<-(AICc(m3_monsoon6mo_lag))-min(AICc(m_null),AICc(m1_spring6mo),AICc(m2_spring6mo),AICc(m3_spring6mo), AICc(m1_monsoon6mo), AICc(m2_monsoon6mo),AICc(m3_monsoon6mo),AICc(m1_spring6mo_lag),AICc(m2_spring6mo_lag),AICc(m3_spring6mo_lag),AICc(m1_monsoon6mo_lag),AICc(m2_monsoon6mo_lag),AICc(m3_monsoon6mo_lag),AICc(m1_spring6mo_AR1),AICc(m1_spring6mo_AR2),AICc(m2_spring6mo_AR1),AICc(m2_spring6mo_AR2),AICc(m3_spring6mo_AR1),AICc(m3_spring6mo_AR2),AICc(m1_monsoon6mo_AR1),AICc(m1_monsoon6mo_AR2),AICc(m2_monsoon6mo_AR1),AICc(m2_monsoon6mo_AR2),AICc(m3_monsoon6mo_AR1),AICc(m3_monsoon6mo_AR2),AICc(m1_spring6mo_lag_AR1),AICc(m1_spring6mo_lag_AR2),AICc(m2_spring6mo_lag_AR1),AICc(m2_spring6mo_lag_AR2),AICc(m3_spring6mo_lag_AR1),AICc(m3_spring6mo_lag_AR2),AICc(m1_monsoon6mo_lag_AR1),AICc(m1_monsoon6mo_lag_AR2),AICc(m2_monsoon6mo_lag_AR1),AICc(m2_monsoon6mo_lag_AR2),AICc(m3_monsoon6mo_lag_AR1),AICc(m3_monsoon6mo_lag_AR2))
  
  # graph the best CSF for the species (maximum abundance as a function of SPEI), if the best model is not the null, and save the graph
  
  if (dAICc_m1_spring6mo == 0 | dAICc_m1_spring6mo_AR1 == 0 | dAICc_m1_spring6mo_AR2 == 0) {
    
    speciesData <- black_original[speciesCode]
    speciesData <-cbind(speciesData,black_climate)
    names(speciesData)[1] <- "speciesCode"
    
    p <- ggplot(speciesData,aes(x=spring6SPEI,y=speciesCode^2))+
      geom_point(size=3,color="black")+
      geom_smooth(formula=y~poly(x,1),method="lm",color="black")+
      xlab("Spring SPEI\n")+
      ylab("Maximum abundance") +       ggtitle(speciesCode) +
      ggtitle(speciesCode) +
      theme(legend.position = "none") +
      theme_classic() +
      theme(axis.text.x = element_text(color="black",size=14)) +
      theme(axis.text.y = element_text(color="black",size=14)) +
      theme(axis.title=element_text(size=18))+
      theme(plot.title = element_text(size=20))
    print(p)
    
    plot_name <- paste("G",speciesCodes[i], "plot",number[i],".pdf", sep="_")
    
    ggsave(paste(plot_name), p,
           width=4,height=4,units = c("in"),
           dpi = 600)
    
    speciesData$speciesCode <-NULL
  }
  
  if (dAICc_m2_spring6mo == 0 | dAICc_m2_spring6mo_AR1 == 0 | dAICc_m2_spring6mo_AR2 == 0) {
    speciesData <- black_original[speciesCode]
    speciesData <-cbind(speciesData,black_climate)
    names(speciesData)[1] <- "speciesCode"
    
    p <- ggplot(speciesData,aes(x=spring6SPEI,y=speciesCode^2))+
      geom_point(size=3,color="black")+
      geom_smooth(formula=y~poly(x,2),method="lm",color="black")+
      xlab("Spring SPEI\n")+
      ylab("Maximum abundance") +       ggtitle(speciesCode) +
      theme(legend.position = "none") +
      theme_classic() +
      theme(axis.text.x = element_text(color="black",size=14)) +
      theme(axis.text.y = element_text(color="black",size=14)) +
      theme(axis.title=element_text(size=18))+
      theme(plot.title = element_text(size=20))
    print(p)
    
    plot_name <- paste("G",speciesCodes[i], "plot",number[i],".pdf", sep="_")
    
    ggsave(paste(plot_name), p,
           width=4,height=4,units = c("in"),
           dpi = 600)
    
    speciesData$speciesCode <-NULL
  }
  
  if (dAICc_m3_spring6mo == 0 | dAICc_m3_spring6mo_AR1 == 0 | dAICc_m3_spring6mo_AR2 == 0) {
    speciesData <- black_original[speciesCode]
    speciesData <-cbind(speciesData,black_climate)
    names(speciesData)[1] <- "speciesCode"
    
    p <- ggplot(speciesData,aes(x=spring6SPEI,y=speciesCode^2))+
      geom_point(size=3,color="black")+
      geom_smooth(formula=y~poly(x,3),method="lm",color="black")+
      xlab("Spring SPEI\n")+
      ylab("Maximum abundance") +       ggtitle(speciesCode) +
      theme(legend.position = "none") +
      theme_classic() +
      theme(axis.text.x = element_text(color="black",size=14)) +
      theme(axis.text.y = element_text(color="black",size=14)) +
      theme(axis.title=element_text(size=18))+
      theme(plot.title = element_text(size=20))
    print(p)
    
    plot_name <- paste("G",speciesCodes[i], "plot",number[i],".pdf", sep="_")
    
    ggsave(paste(plot_name), p,
           width=4,height=4,units = c("in"),
           dpi = 600)    
    speciesData$speciesCode <-NULL
  }
  
  if (dAICc_m1_monsoon6mo == 0 | dAICc_m1_monsoon6mo_AR1 == 0 | dAICc_m1_monsoon6mo_AR2 == 0) {
    speciesData <- black_original[speciesCode]
    speciesData <-cbind(speciesData,black_climate)
    names(speciesData)[1] <- "speciesCode"
    
    p <- ggplot(speciesData,aes(x=monsoon6SPEI,y=speciesCode^2))+
      geom_point(size=3,color="black")+
      geom_smooth(formula=y~poly(x,1),method="lm",color="black")+
      xlab("Monsoon SPEI\n")+
      ylab("Maximum abundance") +       ggtitle(speciesCode) +
      theme(legend.position = "none") +
      theme_classic() +
      theme(axis.text.x = element_text(color="black",size=14)) +
      theme(axis.text.y = element_text(color="black",size=14)) +
      theme(axis.title=element_text(size=18))+
      theme(plot.title = element_text(size=20))
    print(p)
    
    plot_name <- paste("G",speciesCodes[i], "plot",number[i],".pdf", sep="_")
    
    ggsave(paste(plot_name), p,
           width=4,height=4,units = c("in"),
           dpi = 600)
    
    speciesData$speciesCode <-NULL
  }
  
  if (dAICc_m2_monsoon6mo == 0 | dAICc_m2_monsoon6mo_AR1 == 0 | dAICc_m2_monsoon6mo_AR2 == 0) {
    speciesData <- black_original[speciesCode]
    speciesData <-cbind(speciesData,black_climate)
    names(speciesData)[1] <- "speciesCode"
    
    p <- ggplot(speciesData,aes(x=monsoon6SPEI,y=speciesCode^2))+
      geom_point(size=3,color="black")+
      geom_smooth(formula=y~poly(x,2),method="lm",color="black")+
      xlab("Monsoon SPEI\n")+
      ylab("Maximum abundance") +       ggtitle(speciesCode) +
      theme(legend.position = "none") +
      theme_classic() +
      theme(axis.text.x = element_text(color="black",size=14)) +
      theme(axis.text.y = element_text(color="black",size=14)) +
      theme(axis.title=element_text(size=18))+
      theme(plot.title = element_text(size=20))
    print(p)
    
    plot_name <- paste("G",speciesCodes[i], "plot",number[i],".pdf", sep="_")
    
    ggsave(paste(plot_name), p,
           width=4,height=4,units = c("in"),
           dpi = 600)    
    speciesData$speciesCode <-NULL
  }
  
  if (dAICc_m3_monsoon6mo == 0 | dAICc_m3_monsoon6mo_AR1 == 0 | dAICc_m3_monsoon6mo_AR2 == 0) {
    speciesData <- black_original[speciesCode]
    speciesData <-cbind(speciesData,black_climate)
    names(speciesData)[1] <- "speciesCode"
    
    p <- ggplot(speciesData,aes(x=monsoon6SPEI,y=speciesCode^2))+
      geom_point(size=3,color="black")+
      geom_smooth(formula=y~poly(x,3),method="lm",color="black")+
      xlab("Monsoon SPEI\n")+
      ylab("Maximum abundance") +       ggtitle(speciesCode) +
      theme(legend.position = "none") +
      theme_classic() +
      theme(axis.text.x = element_text(color="black",size=14)) +
      theme(axis.text.y = element_text(color="black",size=14)) +
      theme(axis.title=element_text(size=18))+
      theme(plot.title = element_text(size=20))
    print(p)
    
    plot_name <- paste("G",speciesCodes[i], "plot",number[i],".pdf", sep="_")
    
    ggsave(paste(plot_name), p,
           width=4,height=4,units = c("in"),
           dpi = 600)
    
    speciesData$speciesCode <-NULL
  }
  
  if (dAICc_m1_spring6mo_lag == 0 | dAICc_m1_spring6mo_lag_AR1 == 0 | dAICc_m1_spring6mo_lag_AR2 == 0) {
    speciesData <- black_original[speciesCode]
    speciesData <-cbind(speciesData,black_climate)
    names(speciesData)[1] <- "speciesCode"
    
    p <- ggplot(speciesData,aes(x=spring6SPEI_prioryear,y=speciesCode^2))+
      geom_point(size=3,color="black")+
      geom_smooth(formula=y~poly(x,1),method="lm",color="black")+
      xlab("Spring SPEI \n(previous year)")+
      ylab("Maximum abundance") +       ggtitle(speciesCode) +
      theme(legend.position = "none") +
      theme_classic() +
      theme(axis.text.x = element_text(color="black",size=14)) +
      theme(axis.text.y = element_text(color="black",size=14)) +
      theme(axis.title=element_text(size=18))+
      theme(plot.title = element_text(size=20))
    print(p)
    
    plot_name <- paste("G",speciesCodes[i], "plot",number[i],".pdf", sep="_")
    
    ggsave(paste(plot_name), p,
           width=4,height=4,units = c("in"),
           dpi = 600)
    
    speciesData$speciesCode <-NULL
  }
  
  if (dAICc_m2_spring6mo_lag == 0 | dAICc_m2_spring6mo_lag_AR1 == 0 | dAICc_m2_spring6mo_lag_AR2 == 0) {
    speciesData <- black_original[speciesCode]
    speciesData <-cbind(speciesData,black_climate)
    names(speciesData)[1] <- "speciesCode"
    
    p <- ggplot(speciesData,aes(x=spring6SPEI_prioryear,y=speciesCode^2))+
      geom_point(size=3,color="black")+
      geom_smooth(formula=y~poly(x,2),method="lm",color="black")+
      xlab("Spring SPEI \n(previous year)")+
      ylab("Maximum abundance") +       ggtitle(speciesCode) +
      theme(legend.position = "none") +
      theme_classic() +
      theme(axis.text.x = element_text(color="black",size=14)) +
      theme(axis.text.y = element_text(color="black",size=14)) +
      theme(axis.title=element_text(size=18))+
      theme(plot.title = element_text(size=20))
    print(p)
    
    plot_name <- paste("G",speciesCodes[i], "plot",number[i],".pdf", sep="_")
    
    ggsave(paste(plot_name), p,
           width=4,height=4,units = c("in"),
           dpi = 600)
    
    speciesData$speciesCode <-NULL
  }
  
  if (dAICc_m3_spring6mo_lag == 0 | dAICc_m3_spring6mo_lag_AR1 == 0 | dAICc_m3_spring6mo_lag_AR2 == 0) {
    speciesData <- black_original[speciesCode]
    speciesData <-cbind(speciesData,black_climate)
    names(speciesData)[1] <- "speciesCode"
    
    p <- ggplot(speciesData,aes(x=spring6SPEI_prioryear,y=speciesCode^2))+
      geom_point(size=3,color="black")+
      geom_smooth(formula=y~poly(x,3),method="lm",color="black")+
      xlab("Spring SPEI \n(previous year)")+
      ylab("Maximum abundance") +       ggtitle(speciesCode) +
      theme(legend.position = "none") +
      theme_classic() +
      theme(axis.text.x = element_text(color="black",size=14)) +
      theme(axis.text.y = element_text(color="black",size=14)) +
      theme(axis.title=element_text(size=18))+
      theme(plot.title = element_text(size=20))
    print(p)
    
    plot_name <- paste("G",speciesCodes[i], "plot",number[i],".pdf", sep="_")
    
    ggsave(paste(plot_name), p,
           width=4,height=4,units = c("in"),
           dpi = 600)
    
    speciesData$speciesCode <-NULL
  }
  
  if (dAICc_m1_monsoon6mo_lag == 0 | dAICc_m1_monsoon6mo_lag_AR1 == 0 | dAICc_m1_monsoon6mo_lag_AR2 == 0) {
    speciesData <- black_original[speciesCode]
    speciesData <-cbind(speciesData,black_climate)
    names(speciesData)[1] <- "speciesCode"
    
    p <- ggplot(speciesData,aes(x=monsoon6SPEI_prioryear,y=speciesCode^2))+
      geom_point(size=3,color="black")+
      geom_smooth(formula=y~poly(x,1),method="lm",color="black")+
      xlab("Monsoon SPEI \n(previous year)")+
      ylab("Maximum abundance") +       ggtitle(speciesCode) +
      theme(legend.position = "none") +
      theme_classic() +
      theme(axis.text.x = element_text(color="black",size=14)) +
      theme(axis.text.y = element_text(color="black",size=14)) +
      theme(axis.title=element_text(size=18))+
      theme(plot.title = element_text(size=20))
    print(p)
    
    plot_name <- paste("G",speciesCodes[i], "plot",number[i],".pdf", sep="_")
    
    ggsave(paste(plot_name), p,
           width=4,height=4,units = c("in"),
           dpi = 600)
    
    speciesData$speciesCode <-NULL
  }
  
  if (dAICc_m2_monsoon6mo_lag == 0 | dAICc_m2_monsoon6mo_lag_AR1 == 0 | dAICc_m2_monsoon6mo_lag_AR2 == 0) {
    speciesData <- black_original[speciesCode]
    speciesData <-cbind(speciesData,black_climate)
    names(speciesData)[1] <- "speciesCode"
    
    p <- ggplot(speciesData,aes(x=monsoon6SPEI_prioryear,y=speciesCode^2))+
      geom_point(size=3,color="black")+
      geom_smooth(formula=y~poly(x,2),method="lm",color="black")+
      xlab("Monsoon SPEI \n(previous year)")+
      ylab("Maximum abundance") +       ggtitle(speciesCode) +
      theme(legend.position = "none") +
      theme_classic() +
      theme(axis.text.x = element_text(color="black",size=14)) +
      theme(axis.text.y = element_text(color="black",size=14)) +
      theme(axis.title=element_text(size=18))+
      theme(plot.title = element_text(size=20))
    print(p)
    
    plot_name <- paste("G",speciesCodes[i], "plot",number[i],".pdf", sep="_")
    
    ggsave(paste(plot_name), p,
           width=4,height=4,units = c("in"),
           dpi = 600)
    
    speciesData$speciesCode <-NULL
  }
  
  if (dAICc_m3_monsoon6mo_lag == 0 | dAICc_m3_monsoon6mo_lag_AR1 == 0 | dAICc_m3_monsoon6mo_lag_AR2 == 0) {
    speciesData <- black_original[speciesCode]
    speciesData <-cbind(speciesData,black_climate)
    names(speciesData)[1] <- "speciesCode"
    
    p<-ggplot(speciesData,aes(x=monsoon6SPEI_prioryear,y=speciesCode^2))+
      geom_point(size=3,color="black")+
      geom_smooth(formula=y~poly(x,3),method="lm",color="black")+
      xlab("Monsoon SPEI \n(previous year)")+
      ylab("Maximum abundance") +       ggtitle(speciesCode) +
      theme(legend.position = "none") +
      theme_classic() +
      theme(axis.text.x = element_text(color="black",size=14)) +
      theme(axis.text.y = element_text(color="black",size=14)) +
      theme(axis.title=element_text(size=18))+
      theme(plot.title = element_text(size=20))
    print(p)
    
    plot_name <- paste("G",speciesCodes[i], "plot",number[i],".pdf", sep="_")
    
    ggsave(paste(plot_name), p,
           width=4,height=4,units = c("in"),
           dpi = 600)
    
    speciesData$speciesCode <-NULL
  }
  
  
  # create objects containing statistical values related to each model
  
  numDF_lin_m1_spring6mo<-anova(m1_spring6mo,type="marginal")[2,1]
  denDF_lin_m1_spring6mo<-anova(m1_spring6mo,type="marginal")[2,2]
  F_lin_m1_spring6mo<-anova(m1_spring6mo,type="marginal")[2,3]
  P_lin_m1_spring6mo<-anova(m1_spring6mo,type="marginal")[2,4]
  ParamEst_lin_m1_spring6mo<-coef(summary(m1_spring6mo))[2,1]
  SE_lin_m1_spring6mo<-coef(summary(m1_spring6mo))[2,2]
  
  numDF_lin_m2_spring6mo<-anova(m2_spring6mo,type="marginal")[2,1]
  denDF_lin_m2_spring6mo<-anova(m2_spring6mo,type="marginal")[2,2]
  F_lin_m2_spring6mo<-anova(m2_spring6mo,type="marginal")[2,3]
  P_lin_m2_spring6mo<-anova(m2_spring6mo,type="marginal")[2,4]
  ParamEst_lin_m2_spring6mo<-coef(summary(m2_spring6mo))[2,1]
  SE_lin_m2_spring6mo<-coef(summary(m2_spring6mo))[2,2]
  
  numDF_lin_m3_spring6mo<-anova(m3_spring6mo,type="marginal")[2,1]
  denDF_lin_m3_spring6mo<-anova(m3_spring6mo,type="marginal")[2,2]
  F_lin_m3_spring6mo<-anova(m3_spring6mo,type="marginal")[2,3]
  P_lin_m3_spring6mo<-anova(m3_spring6mo,type="marginal")[2,4]
  ParamEst_lin_m3_spring6mo<-coef(summary(m3_spring6mo))[2,1]
  SE_lin_m3_spring6mo<-coef(summary(m3_spring6mo))[2,2]
  
  numDF_lin_m1_monsoon6mo<-anova(m1_monsoon6mo,type="marginal")[2,1]
  denDF_lin_m1_monsoon6mo<-anova(m1_monsoon6mo,type="marginal")[2,2]
  F_lin_m1_monsoon6mo<-anova(m1_monsoon6mo,type="marginal")[2,3]
  P_lin_m1_monsoon6mo<-anova(m1_monsoon6mo,type="marginal")[2,4]
  ParamEst_lin_m1_monsoon6mo<-coef(summary(m1_monsoon6mo))[2,1]
  SE_lin_m1_monsoon6mo<-coef(summary(m1_monsoon6mo))[2,2]
  
  numDF_lin_m2_monsoon6mo<-anova(m2_monsoon6mo,type="marginal")[2,1]
  denDF_lin_m2_monsoon6mo<-anova(m2_monsoon6mo,type="marginal")[2,2]
  F_lin_m2_monsoon6mo<-anova(m2_monsoon6mo,type="marginal")[2,3]
  P_lin_m2_monsoon6mo<-anova(m2_monsoon6mo,type="marginal")[2,4]
  ParamEst_lin_m2_monsoon6mo<-coef(summary(m2_monsoon6mo))[2,1]
  SE_lin_m2_monsoon6mo<-coef(summary(m2_monsoon6mo))[2,2]
  
  numDF_lin_m3_monsoon6mo<-anova(m3_monsoon6mo,type="marginal")[2,1]
  denDF_lin_m3_monsoon6mo<-anova(m3_monsoon6mo,type="marginal")[2,2]
  F_lin_m3_monsoon6mo<-anova(m3_monsoon6mo,type="marginal")[2,3]
  P_lin_m3_monsoon6mo<-anova(m3_monsoon6mo,type="marginal")[2,4]
  ParamEst_lin_m3_monsoon6mo<-coef(summary(m3_monsoon6mo))[2,1]
  SE_lin_m3_monsoon6mo<-coef(summary(m3_monsoon6mo))[2,2]
  
  numDF_lin_m1_spring6mo_lag<-anova(m1_spring6mo_lag,type="marginal")[2,1]
  denDF_lin_m1_spring6mo_lag<-anova(m1_spring6mo_lag,type="marginal")[2,2]
  F_lin_m1_spring6mo_lag<-anova(m1_spring6mo_lag,type="marginal")[2,3]
  P_lin_m1_spring6mo_lag<-anova(m1_spring6mo_lag,type="marginal")[2,4]
  ParamEst_lin_m1_spring6mo_lag<-coef(summary(m1_spring6mo_lag))[2,1]
  SE_lin_m1_spring6mo_lag<-coef(summary(m1_spring6mo_lag))[2,2]
  
  numDF_lin_m2_spring6mo_lag<-anova(m2_spring6mo_lag,type="marginal")[2,1]
  denDF_lin_m2_spring6mo_lag<-anova(m2_spring6mo_lag,type="marginal")[2,2]
  F_lin_m2_spring6mo_lag<-anova(m2_spring6mo_lag,type="marginal")[2,3]
  P_lin_m2_spring6mo_lag<-anova(m2_spring6mo_lag,type="marginal")[2,4]
  ParamEst_lin_m2_spring6mo_lag<-coef(summary(m2_spring6mo_lag))[2,1]
  SE_lin_m2_spring6mo_lag<-coef(summary(m2_spring6mo_lag))[2,2]
  
  numDF_lin_m3_spring6mo_lag<-anova(m3_spring6mo_lag,type="marginal")[2,1]
  denDF_lin_m3_spring6mo_lag<-anova(m3_spring6mo_lag,type="marginal")[2,2]
  F_lin_m3_spring6mo_lag<-anova(m3_spring6mo_lag,type="marginal")[2,3]
  P_lin_m3_spring6mo_lag<-anova(m3_spring6mo_lag,type="marginal")[2,4]
  ParamEst_lin_m3_spring6mo_lag<-coef(summary(m3_spring6mo_lag))[2,1]
  SE_lin_m3_spring6mo_lag<-coef(summary(m3_spring6mo_lag))[2,2]
  
  numDF_lin_m1_monsoon6mo_lag<-anova(m1_monsoon6mo_lag,type="marginal")[2,1]
  denDF_lin_m1_monsoon6mo_lag<-anova(m1_monsoon6mo_lag,type="marginal")[2,2]
  F_lin_m1_monsoon6mo_lag<-anova(m1_monsoon6mo_lag,type="marginal")[2,3]
  P_lin_m1_monsoon6mo_lag<-anova(m1_monsoon6mo_lag,type="marginal")[2,4]
  ParamEst_lin_m1_monsoon6mo_lag<-coef(summary(m1_monsoon6mo_lag))[2,1]
  SE_lin_m1_monsoon6mo_lag<-coef(summary(m1_monsoon6mo_lag))[2,2]
  
  numDF_lin_m2_monsoon6mo_lag<-anova(m2_monsoon6mo_lag,type="marginal")[2,1]
  denDF_lin_m2_monsoon6mo_lag<-anova(m2_monsoon6mo_lag,type="marginal")[2,2]
  F_lin_m2_monsoon6mo_lag<-anova(m2_monsoon6mo_lag,type="marginal")[2,3]
  P_lin_m2_monsoon6mo_lag<-anova(m2_monsoon6mo_lag,type="marginal")[2,4]
  ParamEst_lin_m2_monsoon6mo_lag<-coef(summary(m2_monsoon6mo_lag))[2,1]
  SE_lin_m2_monsoon6mo_lag<-coef(summary(m2_monsoon6mo_lag))[2,2]
  
  numDF_lin_m3_monsoon6mo_lag<-anova(m3_monsoon6mo_lag,type="marginal")[2,1]
  denDF_lin_m3_monsoon6mo_lag<-anova(m3_monsoon6mo_lag,type="marginal")[2,2]
  F_lin_m3_monsoon6mo_lag<-anova(m3_monsoon6mo_lag,type="marginal")[2,3]
  P_lin_m3_monsoon6mo_lag<-anova(m3_monsoon6mo_lag,type="marginal")[2,4]
  ParamEst_lin_m3_monsoon6mo_lag<-coef(summary(m3_monsoon6mo_lag))[2,1]
  SE_lin_m3_monsoon6mo_lag<-coef(summary(m3_monsoon6mo_lag))[2,2]
  
  numDF_lin_m1_spring6mo_AR1<-anova(m1_spring6mo_AR1,type="marginal")[2,1]
  denDF_lin_m1_spring6mo_AR1<-anova(m1_spring6mo_AR1,type="marginal")[2,2]
  F_lin_m1_spring6mo_AR1<-anova(m1_spring6mo_AR1,type="marginal")[2,3]
  P_lin_m1_spring6mo_AR1<-anova(m1_spring6mo_AR1,type="marginal")[2,4]
  ParamEst_lin_m1_spring6mo_AR1<-coef(summary(m1_spring6mo_AR1))[2,1]
  SE_lin_m1_spring6mo_AR1<-coef(summary(m1_spring6mo_AR1))[2,2]
  
  numDF_lin_m1_spring6mo_AR2<-anova(m1_spring6mo_AR2,type="marginal")[2,1]
  denDF_lin_m1_spring6mo_AR2<-anova(m1_spring6mo_AR2,type="marginal")[2,2]
  F_lin_m1_spring6mo_AR2<-anova(m1_spring6mo_AR2,type="marginal")[2,3]
  P_lin_m1_spring6mo_AR2<-anova(m1_spring6mo_AR2,type="marginal")[2,4]
  ParamEst_lin_m1_spring6mo_AR2<-coef(summary(m1_spring6mo_AR2))[2,1]
  SE_lin_m1_spring6mo_AR2<-coef(summary(m1_spring6mo_AR2))[2,2]
  
  numDF_lin_m2_spring6mo_AR1<-anova(m2_spring6mo_AR1,type="marginal")[2,1]
  denDF_lin_m2_spring6mo_AR1<-anova(m2_spring6mo_AR1,type="marginal")[2,2]
  F_lin_m2_spring6mo_AR1<-anova(m2_spring6mo_AR1,type="marginal")[2,3]
  P_lin_m2_spring6mo_AR1<-anova(m2_spring6mo_AR1,type="marginal")[2,4]
  ParamEst_lin_m2_spring6mo_AR1<-coef(summary(m2_spring6mo_AR1))[2,1]
  SE_lin_m2_spring6mo_AR1<-coef(summary(m2_spring6mo_AR1))[2,2]
  
  numDF_lin_m2_spring6mo_AR2<-anova(m2_spring6mo_AR2,type="marginal")[2,1]
  denDF_lin_m2_spring6mo_AR2<-anova(m2_spring6mo_AR2,type="marginal")[2,2]
  F_lin_m2_spring6mo_AR2<-anova(m2_spring6mo_AR2,type="marginal")[2,3]
  P_lin_m2_spring6mo_AR2<-anova(m2_spring6mo_AR2,type="marginal")[2,4]
  ParamEst_lin_m2_spring6mo_AR2<-coef(summary(m2_spring6mo_AR2))[2,1]
  SE_lin_m2_spring6mo_AR2<-coef(summary(m2_spring6mo_AR2))[2,2]
  
  numDF_lin_m3_spring6mo_AR1<-anova(m3_spring6mo_AR1,type="marginal")[2,1]
  denDF_lin_m3_spring6mo_AR1<-anova(m3_spring6mo_AR1,type="marginal")[2,2]
  F_lin_m3_spring6mo_AR1<-anova(m3_spring6mo_AR1,type="marginal")[2,3]
  P_lin_m3_spring6mo_AR1<-anova(m3_spring6mo_AR1,type="marginal")[2,4]
  ParamEst_lin_m3_spring6mo_AR1<-coef(summary(m3_spring6mo_AR1))[2,1]
  SE_lin_m3_spring6mo_AR1<-coef(summary(m3_spring6mo_AR1))[2,2]
  
  numDF_lin_m3_spring6mo_AR2<-anova(m3_spring6mo_AR2,type="marginal")[2,1]
  denDF_lin_m3_spring6mo_AR2<-anova(m3_spring6mo_AR2,type="marginal")[2,2]
  F_lin_m3_spring6mo_AR2<-anova(m3_spring6mo_AR2,type="marginal")[2,3]
  P_lin_m3_spring6mo_AR2<-anova(m3_spring6mo_AR2,type="marginal")[2,4]
  ParamEst_lin_m3_spring6mo_AR2<-coef(summary(m3_spring6mo_AR2))[2,1]
  SE_lin_m3_spring6mo_AR2<-coef(summary(m3_spring6mo_AR2))[2,2]
  
  numDF_lin_m1_monsoon6mo_AR1<-anova(m1_monsoon6mo_AR1,type="marginal")[2,1]
  denDF_lin_m1_monsoon6mo_AR1<-anova(m1_monsoon6mo_AR1,type="marginal")[2,2]
  F_lin_m1_monsoon6mo_AR1<-anova(m1_monsoon6mo_AR1,type="marginal")[2,3]
  P_lin_m1_monsoon6mo_AR1<-anova(m1_monsoon6mo_AR1,type="marginal")[2,4]
  ParamEst_lin_m1_monsoon6mo_AR1<-coef(summary(m1_monsoon6mo_AR1))[2,1]
  SE_lin_m1_monsoon6mo_AR1<-coef(summary(m1_monsoon6mo_AR1))[2,2]
  
  numDF_lin_m1_monsoon6mo_AR2<-anova(m1_monsoon6mo_AR2,type="marginal")[2,1]
  denDF_lin_m1_monsoon6mo_AR2<-anova(m1_monsoon6mo_AR2,type="marginal")[2,2]
  F_lin_m1_monsoon6mo_AR2<-anova(m1_monsoon6mo_AR2,type="marginal")[2,3]
  P_lin_m1_monsoon6mo_AR2<-anova(m1_monsoon6mo_AR2,type="marginal")[2,4]
  ParamEst_lin_m1_monsoon6mo_AR2<-coef(summary(m1_monsoon6mo_AR2))[2,1]
  SE_lin_m1_monsoon6mo_AR2<-coef(summary(m1_monsoon6mo_AR2))[2,2]
  
  numDF_lin_m2_monsoon6mo_AR1<-anova(m2_monsoon6mo_AR1,type="marginal")[2,1]
  denDF_lin_m2_monsoon6mo_AR1<-anova(m2_monsoon6mo_AR1,type="marginal")[2,2]
  F_lin_m2_monsoon6mo_AR1<-anova(m2_monsoon6mo_AR1,type="marginal")[2,3]
  P_lin_m2_monsoon6mo_AR1<-anova(m2_monsoon6mo_AR1,type="marginal")[2,4]
  ParamEst_lin_m2_monsoon6mo_AR1<-coef(summary(m2_monsoon6mo_AR1))[2,1]
  SE_lin_m2_monsoon6mo_AR1<-coef(summary(m2_monsoon6mo_AR1))[2,2]
  
  numDF_lin_m2_monsoon6mo_AR2<-anova(m2_monsoon6mo_AR2,type="marginal")[2,1]
  denDF_lin_m2_monsoon6mo_AR2<-anova(m2_monsoon6mo_AR2,type="marginal")[2,2]
  F_lin_m2_monsoon6mo_AR2<-anova(m2_monsoon6mo_AR2,type="marginal")[2,3]
  P_lin_m2_monsoon6mo_AR2<-anova(m2_monsoon6mo_AR2,type="marginal")[2,4]
  ParamEst_lin_m2_monsoon6mo_AR2<-coef(summary(m2_monsoon6mo_AR2))[2,1]
  SE_lin_m2_monsoon6mo_AR2<-coef(summary(m2_monsoon6mo_AR2))[2,2]
  
  numDF_lin_m3_monsoon6mo_AR1<-anova(m3_monsoon6mo_AR1,type="marginal")[2,1]
  denDF_lin_m3_monsoon6mo_AR1<-anova(m3_monsoon6mo_AR1,type="marginal")[2,2]
  F_lin_m3_monsoon6mo_AR1<-anova(m3_monsoon6mo_AR1,type="marginal")[2,3]
  P_lin_m3_monsoon6mo_AR1<-anova(m3_monsoon6mo_AR1,type="marginal")[2,4]
  ParamEst_lin_m3_monsoon6mo_AR1<-coef(summary(m3_monsoon6mo_AR1))[2,1]
  SE_lin_m3_monsoon6mo_AR1<-coef(summary(m3_monsoon6mo_AR1))[2,2]
  
  numDF_lin_m3_monsoon6mo_AR2<-anova(m3_monsoon6mo_AR2,type="marginal")[2,1]
  denDF_lin_m3_monsoon6mo_AR2<-anova(m3_monsoon6mo_AR2,type="marginal")[2,2]
  F_lin_m3_monsoon6mo_AR2<-anova(m3_monsoon6mo_AR2,type="marginal")[2,3]
  P_lin_m3_monsoon6mo_AR2<-anova(m3_monsoon6mo_AR2,type="marginal")[2,4]
  ParamEst_lin_m3_monsoon6mo_AR2<-coef(summary(m3_monsoon6mo_AR2))[2,1]
  SE_lin_m3_monsoon6mo_AR2<-coef(summary(m3_monsoon6mo_AR2))[2,2]
  
  numDF_lin_m1_spring6mo_lag_AR1<-anova(m1_spring6mo_lag_AR1,type="marginal")[2,1]
  denDF_lin_m1_spring6mo_lag_AR1<-anova(m1_spring6mo_lag_AR1,type="marginal")[2,2]
  F_lin_m1_spring6mo_lag_AR1<-anova(m1_spring6mo_lag_AR1,type="marginal")[2,3]
  P_lin_m1_spring6mo_lag_AR1<-anova(m1_spring6mo_lag_AR1,type="marginal")[2,4]
  ParamEst_lin_m1_spring6mo_lag_AR1<-coef(summary(m1_spring6mo_lag_AR1))[2,1]
  SE_lin_m1_spring6mo_lag_AR1<-coef(summary(m1_spring6mo_lag_AR1))[2,2]
  
  numDF_lin_m1_spring6mo_lag_AR2<-anova(m1_spring6mo_lag_AR2,type="marginal")[2,1]
  denDF_lin_m1_spring6mo_lag_AR2<-anova(m1_spring6mo_lag_AR2,type="marginal")[2,2]
  F_lin_m1_spring6mo_lag_AR2<-anova(m1_spring6mo_lag_AR2,type="marginal")[2,3]
  P_lin_m1_spring6mo_lag_AR2<-anova(m1_spring6mo_lag_AR2,type="marginal")[2,4]
  ParamEst_lin_m1_spring6mo_lag_AR2<-coef(summary(m1_spring6mo_lag_AR2))[2,1]
  SE_lin_m1_spring6mo_lag_AR2<-coef(summary(m1_spring6mo_lag_AR2))[2,2]
  
  numDF_lin_m2_spring6mo_lag_AR1<-anova(m2_spring6mo_lag_AR1,type="marginal")[2,1]
  denDF_lin_m2_spring6mo_lag_AR1<-anova(m2_spring6mo_lag_AR1,type="marginal")[2,2]
  F_lin_m2_spring6mo_lag_AR1<-anova(m2_spring6mo_lag_AR1,type="marginal")[2,3]
  P_lin_m2_spring6mo_lag_AR1<-anova(m2_spring6mo_lag_AR1,type="marginal")[2,4]
  ParamEst_lin_m2_spring6mo_lag_AR1<-coef(summary(m2_spring6mo_lag_AR1))[2,1]
  SE_lin_m2_spring6mo_lag_AR1<-coef(summary(m2_spring6mo_lag_AR1))[2,2]
  
  numDF_lin_m2_spring6mo_lag_AR2<-anova(m2_spring6mo_lag_AR2,type="marginal")[2,1]
  denDF_lin_m2_spring6mo_lag_AR2<-anova(m2_spring6mo_lag_AR2,type="marginal")[2,2]
  F_lin_m2_spring6mo_lag_AR2<-anova(m2_spring6mo_lag_AR2,type="marginal")[2,3]
  P_lin_m2_spring6mo_lag_AR2<-anova(m2_spring6mo_lag_AR2,type="marginal")[2,4]
  ParamEst_lin_m2_spring6mo_lag_AR2<-coef(summary(m2_spring6mo_lag_AR2))[2,1]
  SE_lin_m2_spring6mo_lag_AR2<-coef(summary(m2_spring6mo_lag_AR2))[2,2]
  
  numDF_lin_m3_spring6mo_lag_AR1<-anova(m3_spring6mo_lag_AR1,type="marginal")[2,1]
  denDF_lin_m3_spring6mo_lag_AR1<-anova(m3_spring6mo_lag_AR1,type="marginal")[2,2]
  F_lin_m3_spring6mo_lag_AR1<-anova(m3_spring6mo_lag_AR1,type="marginal")[2,3]
  P_lin_m3_spring6mo_lag_AR1<-anova(m3_spring6mo_lag_AR1,type="marginal")[2,4]
  ParamEst_lin_m3_spring6mo_lag_AR1<-coef(summary(m3_spring6mo_lag_AR1))[2,1]
  SE_lin_m3_spring6mo_lag_AR1<-coef(summary(m3_spring6mo_lag_AR1))[2,2]
  
  numDF_lin_m3_spring6mo_lag_AR2<-anova(m3_spring6mo_lag_AR2,type="marginal")[2,1]
  denDF_lin_m3_spring6mo_lag_AR2<-anova(m3_spring6mo_lag_AR2,type="marginal")[2,2]
  F_lin_m3_spring6mo_lag_AR2<-anova(m3_spring6mo_lag_AR2,type="marginal")[2,3]
  P_lin_m3_spring6mo_lag_AR2<-anova(m3_spring6mo_lag_AR2,type="marginal")[2,4]
  ParamEst_lin_m3_spring6mo_lag_AR2<-coef(summary(m3_spring6mo_lag_AR2))[2,1]
  SE_lin_m3_spring6mo_lag_AR2<-coef(summary(m3_spring6mo_lag_AR2))[2,2]
  
  numDF_lin_m1_monsoon6mo_lag_AR1<-anova(m1_monsoon6mo_lag_AR1,type="marginal")[2,1]
  denDF_lin_m1_monsoon6mo_lag_AR1<-anova(m1_monsoon6mo_lag_AR1,type="marginal")[2,2]
  F_lin_m1_monsoon6mo_lag_AR1<-anova(m1_monsoon6mo_lag_AR1,type="marginal")[2,3]
  P_lin_m1_monsoon6mo_lag_AR1<-anova(m1_monsoon6mo_lag_AR1,type="marginal")[2,4]
  ParamEst_lin_m1_monsoon6mo_lag_AR1<-coef(summary(m1_monsoon6mo_lag_AR1))[2,1]
  SE_lin_m1_monsoon6mo_lag_AR1<-coef(summary(m1_monsoon6mo_lag_AR1))[2,2]
  
  numDF_lin_m1_monsoon6mo_lag_AR2<-anova(m1_monsoon6mo_lag_AR2,type="marginal")[2,1]
  denDF_lin_m1_monsoon6mo_lag_AR2<-anova(m1_monsoon6mo_lag_AR2,type="marginal")[2,2]
  F_lin_m1_monsoon6mo_lag_AR2<-anova(m1_monsoon6mo_lag_AR2,type="marginal")[2,3]
  P_lin_m1_monsoon6mo_lag_AR2<-anova(m1_monsoon6mo_lag_AR2,type="marginal")[2,4]
  ParamEst_lin_m1_monsoon6mo_lag_AR2<-coef(summary(m1_monsoon6mo_lag_AR2))[2,1]
  SE_lin_m1_monsoon6mo_lag_AR2<-coef(summary(m1_monsoon6mo_lag_AR2))[2,2]
  
  numDF_lin_m2_monsoon6mo_lag_AR1<-anova(m2_monsoon6mo_lag_AR1,type="marginal")[2,1]
  denDF_lin_m2_monsoon6mo_lag_AR1<-anova(m2_monsoon6mo_lag_AR1,type="marginal")[2,2]
  F_lin_m2_monsoon6mo_lag_AR1<-anova(m2_monsoon6mo_lag_AR1,type="marginal")[2,3]
  P_lin_m2_monsoon6mo_lag_AR1<-anova(m2_monsoon6mo_lag_AR1,type="marginal")[2,4]
  ParamEst_lin_m2_monsoon6mo_lag_AR1<-coef(summary(m2_monsoon6mo_lag_AR1))[2,1]
  SE_lin_m2_monsoon6mo_lag_AR1<-coef(summary(m2_monsoon6mo_lag_AR1))[2,2]
  
  numDF_lin_m2_monsoon6mo_lag_AR2<-anova(m2_monsoon6mo_lag_AR2,type="marginal")[2,1]
  denDF_lin_m2_monsoon6mo_lag_AR2<-anova(m2_monsoon6mo_lag_AR2,type="marginal")[2,2]
  F_lin_m2_monsoon6mo_lag_AR2<-anova(m2_monsoon6mo_lag_AR2,type="marginal")[2,3]
  P_lin_m2_monsoon6mo_lag_AR2<-anova(m2_monsoon6mo_lag_AR2,type="marginal")[2,4]
  ParamEst_lin_m2_monsoon6mo_lag_AR2<-coef(summary(m2_monsoon6mo_lag_AR2))[2,1]
  SE_lin_m2_monsoon6mo_lag_AR2<-coef(summary(m2_monsoon6mo_lag_AR2))[2,2]
  
  numDF_lin_m3_monsoon6mo_lag_AR1<-anova(m3_monsoon6mo_lag_AR1,type="marginal")[2,1]
  denDF_lin_m3_monsoon6mo_lag_AR1<-anova(m3_monsoon6mo_lag_AR1,type="marginal")[2,2]
  F_lin_m3_monsoon6mo_lag_AR1<-anova(m3_monsoon6mo_lag_AR1,type="marginal")[2,3]
  P_lin_m3_monsoon6mo_lag_AR1<-anova(m3_monsoon6mo_lag_AR1,type="marginal")[2,4]
  ParamEst_lin_m3_monsoon6mo_lag_AR1<-coef(summary(m3_monsoon6mo_lag_AR1))[2,1]
  SE_lin_m3_monsoon6mo_lag_AR1<-coef(summary(m3_monsoon6mo_lag_AR1))[2,2]
  
  numDF_lin_m3_monsoon6mo_lag_AR2<-anova(m3_monsoon6mo_lag_AR2,type="marginal")[2,1]
  denDF_lin_m3_monsoon6mo_lag_AR2<-anova(m3_monsoon6mo_lag_AR2,type="marginal")[2,2]
  F_lin_m3_monsoon6mo_lag_AR2<-anova(m3_monsoon6mo_lag_AR2,type="marginal")[2,3]
  P_lin_m3_monsoon6mo_lag_AR2<-anova(m3_monsoon6mo_lag_AR2,type="marginal")[2,4]
  ParamEst_lin_m3_monsoon6mo_lag_AR2<-coef(summary(m3_monsoon6mo_lag_AR2))[2,1]
  SE_lin_m3_monsoon6mo_lag_AR2<-coef(summary(m3_monsoon6mo_lag_AR2))[2,2]
  
  Rsquared_marginal_m_null<-rsquared(m_null)[1,5]
  Rsquared_marginal_m1_spring6mo<-rsquared(m1_spring6mo)[1,5]
  Rsquared_marginal_m2_spring6mo<-rsquared(m2_spring6mo)[1,5]
  Rsquared_marginal_m3_spring6mo<-rsquared(m3_spring6mo)[1,5]
  Rsquared_marginal_m1_spring6mo_AR1<-rsquared(m1_spring6mo_AR1)[1,5]
  Rsquared_marginal_m1_spring6mo_AR2<-rsquared(m1_spring6mo_AR2)[1,5]
  Rsquared_marginal_m2_spring6mo_AR1<-rsquared(m2_spring6mo_AR1)[1,5]
  Rsquared_marginal_m2_spring6mo_AR2<-rsquared(m2_spring6mo_AR2)[1,5]
  Rsquared_marginal_m3_spring6mo_AR1<-rsquared(m3_spring6mo_AR1)[1,5]
  Rsquared_marginal_m3_spring6mo_AR2<-rsquared(m3_spring6mo_AR2)[1,5]
  Rsquared_marginal_m1_monsoon6mo<-rsquared(m1_monsoon6mo)[1,5]
  Rsquared_marginal_m2_monsoon6mo<-rsquared(m2_monsoon6mo)[1,5]
  Rsquared_marginal_m3_monsoon6mo<-rsquared(m3_monsoon6mo)[1,5]
  Rsquared_marginal_m1_monsoon6mo_AR1<-rsquared(m1_monsoon6mo_AR1)[1,5]
  Rsquared_marginal_m1_monsoon6mo_AR2<-rsquared(m1_monsoon6mo_AR2)[1,5]
  Rsquared_marginal_m2_monsoon6mo_AR1<-rsquared(m2_monsoon6mo_AR1)[1,5]
  Rsquared_marginal_m2_monsoon6mo_AR2<-rsquared(m2_monsoon6mo_AR2)[1,5]
  Rsquared_marginal_m3_monsoon6mo_AR1<-rsquared(m3_monsoon6mo_AR1)[1,5]
  Rsquared_marginal_m3_monsoon6mo_AR2<-rsquared(m3_monsoon6mo_AR2)[1,5]
  Rsquared_marginal_m1_spring6mo_lag<-rsquared(m1_spring6mo_lag)[1,5]
  Rsquared_marginal_m2_spring6mo_lag<-rsquared(m2_spring6mo_lag)[1,5]
  Rsquared_marginal_m3_spring6mo_lag<-rsquared(m3_spring6mo_lag)[1,5]
  Rsquared_marginal_m1_spring6mo_lag_AR1<-rsquared(m1_spring6mo_lag_AR1)[1,5]
  Rsquared_marginal_m1_spring6mo_lag_AR2<-rsquared(m1_spring6mo_lag_AR2)[1,5]
  Rsquared_marginal_m2_spring6mo_lag_AR1<-rsquared(m2_spring6mo_lag_AR1)[1,5]
  Rsquared_marginal_m2_spring6mo_lag_AR2<-rsquared(m2_spring6mo_lag_AR2)[1,5]
  Rsquared_marginal_m3_spring6mo_lag_AR1<-rsquared(m3_spring6mo_lag_AR1)[1,5]
  Rsquared_marginal_m3_spring6mo_lag_AR2<-rsquared(m3_spring6mo_lag_AR2)[1,5]
  Rsquared_marginal_m1_monsoon6mo_lag<-rsquared(m1_monsoon6mo_lag)[1,5]
  Rsquared_marginal_m2_monsoon6mo_lag<-rsquared(m2_monsoon6mo_lag)[1,5]
  Rsquared_marginal_m3_monsoon6mo_lag<-rsquared(m3_monsoon6mo_lag)[1,5]
  Rsquared_marginal_m1_monsoon6mo_lag_AR1<-rsquared(m1_monsoon6mo_lag_AR1)[1,5]
  Rsquared_marginal_m1_monsoon6mo_lag_AR2<-rsquared(m1_monsoon6mo_lag_AR2)[1,5]
  Rsquared_marginal_m2_monsoon6mo_lag_AR1<-rsquared(m2_monsoon6mo_lag_AR1)[1,5]
  Rsquared_marginal_m2_monsoon6mo_lag_AR2<-rsquared(m2_monsoon6mo_lag_AR2)[1,5]
  Rsquared_marginal_m3_monsoon6mo_lag_AR1<-rsquared(m3_monsoon6mo_lag_AR1)[1,5]
  Rsquared_marginal_m3_monsoon6mo_lag_AR2<-rsquared(m3_monsoon6mo_lag_AR2)[1,5]
  
  Rsquared_conditional_m_null<-rsquared(m_null)[1,6]
  Rsquared_conditional_m1_spring6mo<-rsquared(m1_spring6mo)[1,6]
  Rsquared_conditional_m2_spring6mo<-rsquared(m2_spring6mo)[1,6]
  Rsquared_conditional_m3_spring6mo<-rsquared(m3_spring6mo)[1,6]
  Rsquared_conditional_m1_spring6mo_AR1<-rsquared(m1_spring6mo_AR1)[1,6]
  Rsquared_conditional_m1_spring6mo_AR2<-rsquared(m1_spring6mo_AR2)[1,6]
  Rsquared_conditional_m2_spring6mo_AR1<-rsquared(m2_spring6mo_AR1)[1,6]
  Rsquared_conditional_m2_spring6mo_AR2<-rsquared(m2_spring6mo_AR2)[1,6]
  Rsquared_conditional_m3_spring6mo_AR1<-rsquared(m3_spring6mo_AR1)[1,6]
  Rsquared_conditional_m3_spring6mo_AR2<-rsquared(m3_spring6mo_AR2)[1,6]
  Rsquared_conditional_m1_monsoon6mo<-rsquared(m1_monsoon6mo)[1,6]
  Rsquared_conditional_m2_monsoon6mo<-rsquared(m2_monsoon6mo)[1,6]
  Rsquared_conditional_m3_monsoon6mo<-rsquared(m3_monsoon6mo)[1,6]
  Rsquared_conditional_m1_monsoon6mo_AR1<-rsquared(m1_monsoon6mo_AR1)[1,6]
  Rsquared_conditional_m1_monsoon6mo_AR2<-rsquared(m1_monsoon6mo_AR2)[1,6]
  Rsquared_conditional_m2_monsoon6mo_AR1<-rsquared(m2_monsoon6mo_AR1)[1,6]
  Rsquared_conditional_m2_monsoon6mo_AR2<-rsquared(m2_monsoon6mo_AR2)[1,6]
  Rsquared_conditional_m3_monsoon6mo_AR1<-rsquared(m3_monsoon6mo_AR1)[1,6]
  Rsquared_conditional_m3_monsoon6mo_AR2<-rsquared(m3_monsoon6mo_AR2)[1,6]
  Rsquared_conditional_m1_spring6mo_lag<-rsquared(m1_spring6mo_lag)[1,6]
  Rsquared_conditional_m2_spring6mo_lag<-rsquared(m2_spring6mo_lag)[1,6]
  Rsquared_conditional_m3_spring6mo_lag<-rsquared(m3_spring6mo_lag)[1,6]
  Rsquared_conditional_m1_spring6mo_lag_AR1<-rsquared(m1_spring6mo_lag_AR1)[1,6]
  Rsquared_conditional_m1_spring6mo_lag_AR2<-rsquared(m1_spring6mo_lag_AR2)[1,6]
  Rsquared_conditional_m2_spring6mo_lag_AR1<-rsquared(m2_spring6mo_lag_AR1)[1,6]
  Rsquared_conditional_m2_spring6mo_lag_AR2<-rsquared(m2_spring6mo_lag_AR2)[1,6]
  Rsquared_conditional_m3_spring6mo_lag_AR1<-rsquared(m3_spring6mo_lag_AR1)[1,6]
  Rsquared_conditional_m3_spring6mo_lag_AR2<-rsquared(m3_spring6mo_lag_AR2)[1,6]
  Rsquared_conditional_m1_monsoon6mo_lag<-rsquared(m1_monsoon6mo_lag)[1,6]
  Rsquared_conditional_m2_monsoon6mo_lag<-rsquared(m2_monsoon6mo_lag)[1,6]
  Rsquared_conditional_m3_monsoon6mo_lag<-rsquared(m3_monsoon6mo_lag)[1,6]
  Rsquared_conditional_m1_monsoon6mo_lag_AR1<-rsquared(m1_monsoon6mo_lag_AR1)[1,6]
  Rsquared_conditional_m1_monsoon6mo_lag_AR2<-rsquared(m1_monsoon6mo_lag_AR2)[1,6]
  Rsquared_conditional_m2_monsoon6mo_lag_AR1<-rsquared(m2_monsoon6mo_lag_AR1)[1,6]
  Rsquared_conditional_m2_monsoon6mo_lag_AR2<-rsquared(m2_monsoon6mo_lag_AR2)[1,6]
  Rsquared_conditional_m3_monsoon6mo_lag_AR1<-rsquared(m3_monsoon6mo_lag_AR1)[1,6]
  Rsquared_conditional_m3_monsoon6mo_lag_AR2<-rsquared(m3_monsoon6mo_lag_AR2)[1,6]
  
  numDF_quad_m2_spring6mo<-anova(m2_spring6mo,type="marginal")[3,1]
  numDF_quad_m2_spring6mo_AR1<-anova(m2_spring6mo_AR1,type="marginal")[3,1]
  numDF_quad_m2_spring6mo_AR2<-anova(m2_spring6mo_AR2,type="marginal")[3,1]
  numDF_quad_m2_monsoon6mo<-anova(m2_monsoon6mo,type="marginal")[3,1]
  numDF_quad_m2_monsoon6mo_AR1<-anova(m2_monsoon6mo_AR1,type="marginal")[3,1]
  numDF_quad_m2_monsoon6mo_AR2<-anova(m2_monsoon6mo_AR2,type="marginal")[3,1]
  numDF_quad_m2_spring6mo_lag<-anova(m2_spring6mo_lag,type="marginal")[3,1]
  numDF_quad_m2_spring6mo_lag_AR1<-anova(m2_spring6mo_lag_AR1,type="marginal")[3,1]
  numDF_quad_m2_spring6mo_lag_AR2<-anova(m2_spring6mo_lag_AR2,type="marginal")[3,1]
  numDF_quad_m2_monsoon6mo_lag<-anova(m2_monsoon6mo_lag,type="marginal")[3,1]
  numDF_quad_m2_monsoon6mo_lag_AR1<-anova(m2_monsoon6mo_lag_AR1,type="marginal")[3,1]
  numDF_quad_m2_monsoon6mo_lag_AR2<-anova(m2_monsoon6mo_lag_AR2,type="marginal")[3,1]
  
  numDF_quad_m3_spring6mo<-anova(m3_spring6mo,type="marginal")[3,1]
  numDF_quad_m3_spring6mo_AR1<-anova(m3_spring6mo_AR1,type="marginal")[3,1]
  numDF_quad_m3_spring6mo_AR2<-anova(m3_spring6mo_AR2,type="marginal")[3,1]
  numDF_quad_m3_monsoon6mo<-anova(m3_monsoon6mo,type="marginal")[3,1]
  numDF_quad_m3_monsoon6mo_AR1<-anova(m3_monsoon6mo_AR1,type="marginal")[3,1]
  numDF_quad_m3_monsoon6mo_AR2<-anova(m3_monsoon6mo_AR2,type="marginal")[3,1]
  numDF_quad_m3_spring6mo_lag<-anova(m3_spring6mo_lag,type="marginal")[3,1]
  numDF_quad_m3_spring6mo_lag_AR1<-anova(m3_spring6mo_lag_AR1,type="marginal")[3,1]
  numDF_quad_m3_spring6mo_lag_AR2<-anova(m3_spring6mo_lag_AR2,type="marginal")[3,1]
  numDF_quad_m3_monsoon6mo_lag<-anova(m3_monsoon6mo_lag,type="marginal")[3,1]
  numDF_quad_m3_monsoon6mo_lag_AR1<-anova(m3_monsoon6mo_lag_AR1,type="marginal")[3,1]
  numDF_quad_m3_monsoon6mo_lag_AR2<-anova(m3_monsoon6mo_lag_AR2,type="marginal")[3,1]
  
  denDF_quad_m2_spring6mo<-anova(m2_spring6mo,type="marginal")[3,2]
  denDF_quad_m2_spring6mo_AR1<-anova(m2_spring6mo_AR1,type="marginal")[3,2]
  denDF_quad_m2_spring6mo_AR2<-anova(m2_spring6mo_AR2,type="marginal")[3,2]
  denDF_quad_m2_monsoon6mo<-anova(m2_monsoon6mo,type="marginal")[3,2]
  denDF_quad_m2_monsoon6mo_AR1<-anova(m2_monsoon6mo_AR1,type="marginal")[3,2]
  denDF_quad_m2_monsoon6mo_AR2<-anova(m2_monsoon6mo_AR2,type="marginal")[3,2]
  denDF_quad_m2_spring6mo_lag<-anova(m2_spring6mo_lag,type="marginal")[3,2]
  denDF_quad_m2_spring6mo_lag_AR1<-anova(m2_spring6mo_lag_AR1,type="marginal")[3,2]
  denDF_quad_m2_spring6mo_lag_AR2<-anova(m2_spring6mo_lag_AR2,type="marginal")[3,2]
  denDF_quad_m2_monsoon6mo_lag<-anova(m2_monsoon6mo_lag,type="marginal")[3,2]
  denDF_quad_m2_monsoon6mo_lag_AR1<-anova(m2_monsoon6mo_lag_AR1,type="marginal")[3,2]
  denDF_quad_m2_monsoon6mo_lag_AR2<-anova(m2_monsoon6mo_lag_AR2,type="marginal")[3,2]
  
  denDF_quad_m3_spring6mo<-anova(m3_spring6mo,type="marginal")[3,2]
  denDF_quad_m3_spring6mo_AR1<-anova(m3_spring6mo_AR1,type="marginal")[3,2]
  denDF_quad_m3_spring6mo_AR2<-anova(m3_spring6mo_AR2,type="marginal")[3,2]
  denDF_quad_m3_monsoon6mo<-anova(m3_monsoon6mo,type="marginal")[3,2]
  denDF_quad_m3_monsoon6mo_AR1<-anova(m3_monsoon6mo_AR1,type="marginal")[3,2]
  denDF_quad_m3_monsoon6mo_AR2<-anova(m3_monsoon6mo_AR2,type="marginal")[3,2]
  denDF_quad_m3_spring6mo_lag<-anova(m3_spring6mo_lag,type="marginal")[3,2]
  denDF_quad_m3_spring6mo_lag_AR1<-anova(m3_spring6mo_lag_AR1,type="marginal")[3,2]
  denDF_quad_m3_spring6mo_lag_AR2<-anova(m3_spring6mo_lag_AR2,type="marginal")[3,2]
  denDF_quad_m3_monsoon6mo_lag<-anova(m3_monsoon6mo_lag,type="marginal")[3,2]
  denDF_quad_m3_monsoon6mo_lag_AR1<-anova(m3_monsoon6mo_lag_AR1,type="marginal")[3,2]
  denDF_quad_m3_monsoon6mo_lag_AR2<-anova(m3_monsoon6mo_lag_AR2,type="marginal")[3,2]
  
  F_quad_m2_spring6mo<-anova(m2_spring6mo,type="marginal")[3,3]
  F_quad_m2_spring6mo_AR1<-anova(m2_spring6mo_AR1,type="marginal")[3,3]
  F_quad_m2_spring6mo_AR2<-anova(m2_spring6mo_AR2,type="marginal")[3,3]
  F_quad_m2_monsoon6mo<-anova(m2_monsoon6mo,type="marginal")[3,3]
  F_quad_m2_monsoon6mo_AR1<-anova(m2_monsoon6mo_AR1,type="marginal")[3,3]
  F_quad_m2_monsoon6mo_AR2<-anova(m2_monsoon6mo_AR2,type="marginal")[3,3]
  F_quad_m2_spring6mo_lag<-anova(m2_spring6mo_lag,type="marginal")[3,3]
  F_quad_m2_spring6mo_lag_AR1<-anova(m2_spring6mo_lag_AR1,type="marginal")[3,3]
  F_quad_m2_spring6mo_lag_AR2<-anova(m2_spring6mo_lag_AR2,type="marginal")[3,3]
  F_quad_m2_monsoon6mo_lag<-anova(m2_monsoon6mo_lag,type="marginal")[3,3]
  F_quad_m2_monsoon6mo_lag_AR1<-anova(m2_monsoon6mo_lag_AR1,type="marginal")[3,3]
  F_quad_m2_monsoon6mo_lag_AR2<-anova(m2_monsoon6mo_lag_AR2,type="marginal")[3,3]
  
  F_quad_m3_spring6mo<-anova(m3_spring6mo,type="marginal")[3,3]
  F_quad_m3_spring6mo_AR1<-anova(m3_spring6mo_AR1,type="marginal")[3,3]
  F_quad_m3_spring6mo_AR2<-anova(m3_spring6mo_AR2,type="marginal")[3,3]
  F_quad_m3_monsoon6mo<-anova(m3_monsoon6mo,type="marginal")[3,3]
  F_quad_m3_monsoon6mo_AR1<-anova(m3_monsoon6mo_AR1,type="marginal")[3,3]
  F_quad_m3_monsoon6mo_AR2<-anova(m3_monsoon6mo_AR2,type="marginal")[3,3]
  F_quad_m3_spring6mo_lag<-anova(m3_spring6mo_lag,type="marginal")[3,3]
  F_quad_m3_spring6mo_lag_AR1<-anova(m3_spring6mo_lag_AR1,type="marginal")[3,3]
  F_quad_m3_spring6mo_lag_AR2<-anova(m3_spring6mo_lag_AR2,type="marginal")[3,3]
  F_quad_m3_monsoon6mo_lag<-anova(m3_monsoon6mo_lag,type="marginal")[3,3]
  F_quad_m3_monsoon6mo_lag_AR1<-anova(m3_monsoon6mo_lag_AR1,type="marginal")[3,3]
  F_quad_m3_monsoon6mo_lag_AR2<-anova(m3_monsoon6mo_lag_AR2,type="marginal")[3,3]
  
  P_quad_m2_spring6mo<-anova(m2_spring6mo,type="marginal")[3,4]
  P_quad_m2_spring6mo_AR1<-anova(m2_spring6mo_AR1,type="marginal")[3,4]
  P_quad_m2_spring6mo_AR2<-anova(m2_spring6mo_AR2,type="marginal")[3,4]
  P_quad_m2_monsoon6mo<-anova(m2_monsoon6mo,type="marginal")[3,4]
  P_quad_m2_monsoon6mo_AR1<-anova(m2_monsoon6mo_AR1,type="marginal")[3,4]
  P_quad_m2_monsoon6mo_AR2<-anova(m2_monsoon6mo_AR2,type="marginal")[3,4]
  P_quad_m2_spring6mo_lag<-anova(m2_spring6mo_lag,type="marginal")[3,4]
  P_quad_m2_spring6mo_lag_AR1<-anova(m2_spring6mo_lag_AR1,type="marginal")[3,4]
  P_quad_m2_spring6mo_lag_AR2<-anova(m2_spring6mo_lag_AR2,type="marginal")[3,4]
  P_quad_m2_monsoon6mo_lag<-anova(m2_monsoon6mo_lag,type="marginal")[3,4]
  P_quad_m2_monsoon6mo_lag_AR1<-anova(m2_monsoon6mo_lag_AR1,type="marginal")[3,4]
  P_quad_m2_monsoon6mo_lag_AR2<-anova(m2_monsoon6mo_lag_AR2,type="marginal")[3,4]
  
  P_quad_m3_spring6mo<-anova(m3_spring6mo,type="marginal")[3,4]
  P_quad_m3_spring6mo_AR1<-anova(m3_spring6mo_AR1,type="marginal")[3,4]
  P_quad_m3_spring6mo_AR2<-anova(m3_spring6mo_AR2,type="marginal")[3,4]
  P_quad_m3_monsoon6mo<-anova(m3_monsoon6mo,type="marginal")[3,4]
  P_quad_m3_monsoon6mo_AR1<-anova(m3_monsoon6mo_AR1,type="marginal")[3,4]
  P_quad_m3_monsoon6mo_AR2<-anova(m3_monsoon6mo_AR2,type="marginal")[3,4]
  P_quad_m3_spring6mo_lag<-anova(m3_spring6mo_lag,type="marginal")[3,4]
  P_quad_m3_spring6mo_lag_AR1<-anova(m3_spring6mo_lag_AR1,type="marginal")[3,4]
  P_quad_m3_spring6mo_lag_AR2<-anova(m3_spring6mo_lag_AR2,type="marginal")[3,4]
  P_quad_m3_monsoon6mo_lag<-anova(m3_monsoon6mo_lag,type="marginal")[3,4]
  P_quad_m3_monsoon6mo_lag_AR1<-anova(m3_monsoon6mo_lag_AR1,type="marginal")[3,4]
  P_quad_m3_monsoon6mo_lag_AR2<-anova(m3_monsoon6mo_lag_AR2,type="marginal")[3,4]
  
  ParamEst_quad_m2_spring6mo<-coef(summary(m2_spring6mo))[3,1]
  ParamEst_quad_m2_spring6mo_AR1<-coef(summary(m2_spring6mo_AR1))[3,1]
  ParamEst_quad_m2_spring6mo_AR2<-coef(summary(m2_spring6mo_AR2))[3,1]
  ParamEst_quad_m2_monsoon6mo<-coef(summary(m2_monsoon6mo))[3,1]
  ParamEst_quad_m2_monsoon6mo_AR1<-coef(summary(m2_monsoon6mo_AR1))[3,1]
  ParamEst_quad_m2_monsoon6mo_AR2<-coef(summary(m2_monsoon6mo_AR2))[3,1]
  ParamEst_quad_m2_spring6mo_lag<-coef(summary(m2_spring6mo_lag))[3,1]
  ParamEst_quad_m2_spring6mo_lag_AR1<-coef(summary(m2_spring6mo_lag_AR1))[3,1]
  ParamEst_quad_m2_spring6mo_lag_AR2<-coef(summary(m2_spring6mo_lag_AR2))[3,1]
  ParamEst_quad_m2_monsoon6mo_lag<-coef(summary(m2_monsoon6mo_lag))[3,1]
  ParamEst_quad_m2_monsoon6mo_lag_AR1<-coef(summary(m2_monsoon6mo_lag_AR1))[3,1]
  ParamEst_quad_m2_monsoon6mo_lag_AR2<-coef(summary(m2_monsoon6mo_lag_AR2))[3,1]
  
  ParamEst_quad_m3_spring6mo<-coef(summary(m3_spring6mo))[3,1]
  ParamEst_quad_m3_spring6mo_AR1<-coef(summary(m3_spring6mo_AR1))[3,1]
  ParamEst_quad_m3_spring6mo_AR2<-coef(summary(m3_spring6mo_AR2))[3,1]
  ParamEst_quad_m3_monsoon6mo<-coef(summary(m3_monsoon6mo))[3,1]
  ParamEst_quad_m3_monsoon6mo_AR1<-coef(summary(m3_monsoon6mo_AR1))[3,1]
  ParamEst_quad_m3_monsoon6mo_AR2<-coef(summary(m3_monsoon6mo_AR2))[3,1]
  ParamEst_quad_m3_spring6mo_lag<-coef(summary(m3_spring6mo_lag))[3,1]
  ParamEst_quad_m3_spring6mo_lag_AR1<-coef(summary(m3_spring6mo_lag_AR1))[3,1]
  ParamEst_quad_m3_spring6mo_lag_AR2<-coef(summary(m3_spring6mo_lag_AR2))[3,1]
  ParamEst_quad_m3_monsoon6mo_lag<-coef(summary(m3_monsoon6mo_lag))[3,1]
  ParamEst_quad_m3_monsoon6mo_lag_AR1<-coef(summary(m3_monsoon6mo_lag_AR1))[3,1]
  ParamEst_quad_m3_monsoon6mo_lag_AR2<-coef(summary(m3_monsoon6mo_lag_AR2))[3,1]
  
  SE_quad_m2_spring6mo<-coef(summary(m2_spring6mo))[3,2]
  SE_quad_m2_spring6mo_AR1<-coef(summary(m2_spring6mo_AR1))[3,2]
  SE_quad_m2_spring6mo_AR2<-coef(summary(m2_spring6mo_AR2))[3,2]
  SE_quad_m2_monsoon6mo<-coef(summary(m2_monsoon6mo))[3,2]
  SE_quad_m2_monsoon6mo_AR1<-coef(summary(m2_monsoon6mo_AR1))[3,2]
  SE_quad_m2_monsoon6mo_AR2<-coef(summary(m2_monsoon6mo_AR2))[3,2]
  SE_quad_m2_spring6mo_lag<-coef(summary(m2_spring6mo_lag))[3,2]
  SE_quad_m2_spring6mo_lag_AR1<-coef(summary(m2_spring6mo_lag_AR1))[3,2]
  SE_quad_m2_spring6mo_lag_AR2<-coef(summary(m2_spring6mo_lag_AR2))[3,2]
  SE_quad_m2_monsoon6mo_lag<-coef(summary(m2_monsoon6mo_lag))[3,2]
  SE_quad_m2_monsoon6mo_lag_AR1<-coef(summary(m2_monsoon6mo_lag_AR1))[3,2]
  SE_quad_m2_monsoon6mo_lag_AR2<-coef(summary(m2_monsoon6mo_lag_AR2))[3,2]
  
  SE_quad_m3_spring6mo<-coef(summary(m3_spring6mo))[3,2]
  SE_quad_m3_spring6mo_AR1<-coef(summary(m3_spring6mo_AR1))[3,2]
  SE_quad_m3_spring6mo_AR2<-coef(summary(m3_spring6mo_AR2))[3,2]
  SE_quad_m3_monsoon6mo<-coef(summary(m3_monsoon6mo))[3,2]
  SE_quad_m3_monsoon6mo_AR1<-coef(summary(m3_monsoon6mo_AR1))[3,2]
  SE_quad_m3_monsoon6mo_AR2<-coef(summary(m3_monsoon6mo_AR2))[3,2]
  SE_quad_m3_spring6mo_lag<-coef(summary(m3_spring6mo_lag))[3,2]
  SE_quad_m3_spring6mo_lag_AR1<-coef(summary(m3_spring6mo_lag_AR1))[3,2]
  SE_quad_m3_spring6mo_lag_AR2<-coef(summary(m3_spring6mo_lag_AR2))[3,2]
  SE_quad_m3_monsoon6mo_lag<-coef(summary(m3_monsoon6mo_lag))[3,2]
  SE_quad_m3_monsoon6mo_lag_AR1<-coef(summary(m3_monsoon6mo_lag_AR1))[3,2]
  SE_quad_m3_monsoon6mo_lag_AR2<-coef(summary(m3_monsoon6mo_lag_AR2))[3,2]
  
  numDF_cub_m3_spring6mo<-anova(m3_spring6mo,type="marginal")[4,1]
  numDF_cub_m3_spring6mo_AR1<-anova(m3_spring6mo_AR1,type="marginal")[4,1]
  numDF_cub_m3_spring6mo_AR2<-anova(m3_spring6mo_AR2,type="marginal")[4,1]
  numDF_cub_m3_monsoon6mo<-anova(m3_monsoon6mo,type="marginal")[4,1]
  numDF_cub_m3_monsoon6mo_AR1<-anova(m3_monsoon6mo_AR1,type="marginal")[4,1]
  numDF_cub_m3_monsoon6mo_AR2<-anova(m3_monsoon6mo_AR2,type="marginal")[4,1]
  numDF_cub_m3_spring6mo_lag<-anova(m3_spring6mo_lag,type="marginal")[4,1]
  numDF_cub_m3_spring6mo_lag_AR1<-anova(m3_spring6mo_lag_AR1,type="marginal")[4,1]
  numDF_cub_m3_spring6mo_lag_AR2<-anova(m3_spring6mo_lag_AR2,type="marginal")[4,1]
  numDF_cub_m3_monsoon6mo_lag<-anova(m3_monsoon6mo_lag,type="marginal")[4,1]
  numDF_cub_m3_monsoon6mo_lag_AR1<-anova(m3_monsoon6mo_lag_AR1,type="marginal")[4,1]
  numDF_cub_m3_monsoon6mo_lag_AR2<-anova(m3_monsoon6mo_lag_AR2,type="marginal")[4,1]
  
  denDF_cub_m3_spring6mo<-anova(m3_spring6mo,type="marginal")[4,2]
  denDF_cub_m3_spring6mo_AR1<-anova(m3_spring6mo_AR1,type="marginal")[4,2]
  denDF_cub_m3_spring6mo_AR2<-anova(m3_spring6mo_AR2,type="marginal")[4,2]
  denDF_cub_m3_monsoon6mo<-anova(m3_monsoon6mo,type="marginal")[4,2]
  denDF_cub_m3_monsoon6mo_AR1<-anova(m3_monsoon6mo_AR1,type="marginal")[4,2]
  denDF_cub_m3_monsoon6mo_AR2<-anova(m3_monsoon6mo_AR2,type="marginal")[4,2]
  denDF_cub_m3_spring6mo_lag<-anova(m3_spring6mo_lag,type="marginal")[4,2]
  denDF_cub_m3_spring6mo_lag_AR1<-anova(m3_spring6mo_lag_AR1,type="marginal")[4,2]
  denDF_cub_m3_spring6mo_lag_AR2<-anova(m3_spring6mo_lag_AR2,type="marginal")[4,2]
  denDF_cub_m3_monsoon6mo_lag<-anova(m3_monsoon6mo_lag,type="marginal")[4,2]
  denDF_cub_m3_monsoon6mo_lag_AR1<-anova(m3_monsoon6mo_lag_AR1,type="marginal")[4,2]
  denDF_cub_m3_monsoon6mo_lag_AR2<-anova(m3_monsoon6mo_lag_AR2,type="marginal")[4,2]
  
  F_cub_m3_spring6mo<-anova(m3_spring6mo,type="marginal")[4,3]
  F_cub_m3_spring6mo_AR1<-anova(m3_spring6mo_AR1,type="marginal")[4,3]
  F_cub_m3_spring6mo_AR2<-anova(m3_spring6mo_AR2,type="marginal")[4,3]
  F_cub_m3_monsoon6mo<-anova(m3_monsoon6mo,type="marginal")[4,3]
  F_cub_m3_monsoon6mo_AR1<-anova(m3_monsoon6mo_AR1,type="marginal")[4,3]
  F_cub_m3_monsoon6mo_AR2<-anova(m3_monsoon6mo_AR2,type="marginal")[4,3]
  F_cub_m3_spring6mo_lag<-anova(m3_spring6mo_lag,type="marginal")[4,3]
  F_cub_m3_spring6mo_lag_AR1<-anova(m3_spring6mo_lag_AR1,type="marginal")[4,3]
  F_cub_m3_spring6mo_lag_AR2<-anova(m3_spring6mo_lag_AR2,type="marginal")[4,3]
  F_cub_m3_monsoon6mo_lag<-anova(m3_monsoon6mo_lag,type="marginal")[4,3]
  F_cub_m3_monsoon6mo_lag_AR1<-anova(m3_monsoon6mo_lag_AR1,type="marginal")[4,3]
  F_cub_m3_monsoon6mo_lag_AR2<-anova(m3_monsoon6mo_lag_AR2,type="marginal")[4,3]
  
  P_cub_m3_spring6mo<-anova(m3_spring6mo,type="marginal")[4,4]
  P_cub_m3_spring6mo_AR1<-anova(m3_spring6mo_AR1,type="marginal")[4,4]
  P_cub_m3_spring6mo_AR2<-anova(m3_spring6mo_AR2,type="marginal")[4,4]
  P_cub_m3_monsoon6mo<-anova(m3_monsoon6mo,type="marginal")[4,4]
  P_cub_m3_monsoon6mo_AR1<-anova(m3_monsoon6mo_AR1,type="marginal")[4,4]
  P_cub_m3_monsoon6mo_AR2<-anova(m3_monsoon6mo_AR2,type="marginal")[4,4]
  P_cub_m3_spring6mo_lag<-anova(m3_spring6mo_lag,type="marginal")[4,4]
  P_cub_m3_spring6mo_lag_AR1<-anova(m3_spring6mo_lag_AR1,type="marginal")[4,4]
  P_cub_m3_spring6mo_lag_AR2<-anova(m3_spring6mo_lag_AR2,type="marginal")[4,4]
  P_cub_m3_monsoon6mo_lag<-anova(m3_monsoon6mo_lag,type="marginal")[4,4]
  P_cub_m3_monsoon6mo_lag_AR1<-anova(m3_monsoon6mo_lag_AR1,type="marginal")[4,4]
  P_cub_m3_monsoon6mo_lag_AR2<-anova(m3_monsoon6mo_lag_AR2,type="marginal")[4,4]
  
  ParamEst_cub_m3_spring6mo<-coef(summary(m3_spring6mo))[4,1]
  ParamEst_cub_m3_spring6mo_AR1<-coef(summary(m3_spring6mo_AR1))[4,1]
  ParamEst_cub_m3_spring6mo_AR2<-coef(summary(m3_spring6mo_AR2))[4,1]
  ParamEst_cub_m3_monsoon6mo<-coef(summary(m3_monsoon6mo))[4,1]
  ParamEst_cub_m3_monsoon6mo_AR1<-coef(summary(m3_monsoon6mo_AR1))[4,1]
  ParamEst_cub_m3_monsoon6mo_AR2<-coef(summary(m3_monsoon6mo_AR2))[4,1]
  ParamEst_cub_m3_spring6mo_lag<-coef(summary(m3_spring6mo_lag))[4,1]
  ParamEst_cub_m3_spring6mo_lag_AR1<-coef(summary(m3_spring6mo_lag_AR1))[4,1]
  ParamEst_cub_m3_spring6mo_lag_AR2<-coef(summary(m3_spring6mo_lag_AR2))[4,1]
  ParamEst_cub_m3_monsoon6mo_lag<-coef(summary(m3_monsoon6mo_lag))[4,1]
  ParamEst_cub_m3_monsoon6mo_lag_AR1<-coef(summary(m3_monsoon6mo_lag_AR1))[4,1]
  ParamEst_cub_m3_monsoon6mo_lag_AR2<-coef(summary(m3_monsoon6mo_lag_AR2))[4,1]
  
  SE_cub_m3_spring6mo<-coef(summary(m3_spring6mo))[4,2]
  SE_cub_m3_spring6mo_AR1<-coef(summary(m3_spring6mo_AR1))[4,2]
  SE_cub_m3_spring6mo_AR2<-coef(summary(m3_spring6mo_AR2))[4,2]
  SE_cub_m3_monsoon6mo<-coef(summary(m3_monsoon6mo))[4,2]
  SE_cub_m3_monsoon6mo_AR1<-coef(summary(m3_monsoon6mo_AR1))[4,2]
  SE_cub_m3_monsoon6mo_AR2<-coef(summary(m3_monsoon6mo_AR2))[4,2]
  SE_cub_m3_spring6mo_lag<-coef(summary(m3_spring6mo_lag))[4,2]
  SE_cub_m3_spring6mo_lag_AR1<-coef(summary(m3_spring6mo_lag_AR1))[4,2]
  SE_cub_m3_spring6mo_lag_AR2<-coef(summary(m3_spring6mo_lag_AR2))[4,2]
  SE_cub_m3_monsoon6mo_lag<-coef(summary(m3_monsoon6mo_lag))[4,2]
  SE_cub_m3_monsoon6mo_lag_AR1<-coef(summary(m3_monsoon6mo_lag_AR1))[4,2]
  SE_cub_m3_monsoon6mo_lag_AR2<-coef(summary(m3_monsoon6mo_lag_AR2))[4,2]
  
  # bind all target output values together
  output_id<-cbind(speciesCode, 
                   ecosystem, 
                   P_shapiro_null,
                   dAICc_m1_monsoon6mo,
                   dAICc_m1_monsoon6mo_lag,
                   dAICc_m1_spring6mo,
                   dAICc_m1_spring6mo_lag,
                   dAICc_m2_monsoon6mo,
                   dAICc_m2_monsoon6mo_lag,
                   dAICc_m2_spring6mo,
                   dAICc_m2_spring6mo_lag,
                   dAICc_m3_monsoon6mo,
                   dAICc_m3_monsoon6mo_lag,
                   dAICc_m3_spring6mo,
                   dAICc_m3_spring6mo_lag,
                   dAICc_null,
                   dAICc_m1_spring6mo_AR1,
                   dAICc_m1_spring6mo_AR2,
                   dAICc_m2_spring6mo_AR1,
                   dAICc_m2_spring6mo_AR2,
                   dAICc_m3_spring6mo_AR1,
                   dAICc_m3_spring6mo_AR2,
                   dAICc_m1_monsoon6mo_AR1,
                   dAICc_m1_monsoon6mo_AR2,
                   dAICc_m2_monsoon6mo_AR1,
                   dAICc_m2_monsoon6mo_AR2,
                   dAICc_m3_monsoon6mo_AR1,
                   dAICc_m3_monsoon6mo_AR2,
                   dAICc_m1_spring6mo_lag_AR1,
                   dAICc_m1_spring6mo_lag_AR2,
                   dAICc_m2_spring6mo_lag_AR1,
                   dAICc_m2_spring6mo_lag_AR2,
                   dAICc_m3_spring6mo_lag_AR1,
                   dAICc_m3_spring6mo_lag_AR2,
                   dAICc_m1_monsoon6mo_lag_AR1,
                   dAICc_m1_monsoon6mo_lag_AR2,
                   dAICc_m2_monsoon6mo_lag_AR1,
                   dAICc_m2_monsoon6mo_lag_AR2,
                   dAICc_m3_monsoon6mo_lag_AR1,
                   dAICc_m3_monsoon6mo_lag_AR2,
                   numDF_lin_m1_spring6mo,
                   denDF_lin_m1_spring6mo,
                   F_lin_m1_spring6mo,
                   P_lin_m1_spring6mo,
                   numDF_lin_m2_spring6mo,
                   denDF_lin_m2_spring6mo,
                   F_lin_m2_spring6mo,
                   P_lin_m2_spring6mo,
                   numDF_lin_m3_spring6mo,
                   denDF_lin_m3_spring6mo,
                   F_lin_m3_spring6mo,
                   P_lin_m3_spring6mo,
                   numDF_lin_m1_monsoon6mo,
                   denDF_lin_m1_monsoon6mo,
                   F_lin_m1_monsoon6mo,
                   P_lin_m1_monsoon6mo,
                   numDF_lin_m2_monsoon6mo,
                   denDF_lin_m2_monsoon6mo,
                   F_lin_m2_monsoon6mo,
                   P_lin_m2_monsoon6mo,
                   numDF_lin_m3_monsoon6mo,
                   denDF_lin_m3_monsoon6mo,
                   F_lin_m3_monsoon6mo,
                   P_lin_m3_monsoon6mo,
                   numDF_lin_m1_spring6mo_lag,
                   denDF_lin_m1_spring6mo_lag,
                   F_lin_m1_spring6mo_lag,
                   P_lin_m1_spring6mo_lag,
                   numDF_lin_m2_spring6mo_lag,
                   denDF_lin_m2_spring6mo_lag,
                   F_lin_m2_spring6mo_lag,
                   P_lin_m2_spring6mo_lag,
                   numDF_lin_m3_spring6mo_lag,
                   denDF_lin_m3_spring6mo_lag,
                   F_lin_m3_spring6mo_lag,
                   P_lin_m3_spring6mo_lag,
                   numDF_lin_m1_monsoon6mo_lag,
                   denDF_lin_m1_monsoon6mo_lag,
                   F_lin_m1_monsoon6mo_lag,
                   P_lin_m1_monsoon6mo_lag,
                   numDF_lin_m2_monsoon6mo_lag,
                   denDF_lin_m2_monsoon6mo_lag,
                   F_lin_m2_monsoon6mo_lag,
                   P_lin_m2_monsoon6mo_lag,
                   numDF_lin_m3_monsoon6mo_lag,
                   denDF_lin_m3_monsoon6mo_lag,
                   F_lin_m3_monsoon6mo_lag,
                   P_lin_m3_monsoon6mo_lag,
                   ParamEst_lin_m1_monsoon6mo,
                   ParamEst_lin_m1_monsoon6mo_lag,
                   ParamEst_lin_m1_spring6mo,
                   ParamEst_lin_m1_spring6mo_lag,
                   ParamEst_lin_m2_monsoon6mo,
                   ParamEst_lin_m2_monsoon6mo_lag,
                   ParamEst_lin_m2_spring6mo,
                   ParamEst_lin_m2_spring6mo_lag,
                   ParamEst_lin_m3_monsoon6mo,
                   ParamEst_lin_m3_monsoon6mo_lag,
                   ParamEst_lin_m3_spring6mo,
                   ParamEst_lin_m3_spring6mo_lag,
                   SE_lin_m1_monsoon6mo,
                   SE_lin_m1_monsoon6mo_lag,
                   SE_lin_m1_spring6mo,
                   SE_lin_m1_spring6mo_lag,
                   SE_lin_m2_monsoon6mo,
                   SE_lin_m2_monsoon6mo_lag,
                   SE_lin_m2_spring6mo,
                   SE_lin_m2_spring6mo_lag,
                   SE_lin_m3_monsoon6mo,
                   SE_lin_m3_monsoon6mo_lag,
                   SE_lin_m3_spring6mo,
                   SE_lin_m3_spring6mo_lag,
                   numDF_lin_m1_spring6mo_AR1,
                   denDF_lin_m1_spring6mo_AR1,
                   F_lin_m1_spring6mo_AR1,
                   P_lin_m1_spring6mo_AR1,
                   ParamEst_lin_m1_spring6mo_AR1,
                   SE_lin_m1_spring6mo_AR1,
                   numDF_lin_m1_spring6mo_AR2,
                   denDF_lin_m1_spring6mo_AR2,
                   F_lin_m1_spring6mo_AR2,
                   P_lin_m1_spring6mo_AR2,
                   ParamEst_lin_m1_spring6mo_AR2,
                   SE_lin_m1_spring6mo_AR2,
                   numDF_lin_m2_spring6mo_AR1,
                   denDF_lin_m2_spring6mo_AR1,
                   F_lin_m2_spring6mo_AR1,
                   P_lin_m2_spring6mo_AR1,
                   ParamEst_lin_m2_spring6mo_AR1,
                   SE_lin_m2_spring6mo_AR1,
                   numDF_lin_m2_spring6mo_AR2,
                   denDF_lin_m2_spring6mo_AR2,
                   F_lin_m2_spring6mo_AR2,
                   P_lin_m2_spring6mo_AR2,
                   ParamEst_lin_m2_spring6mo_AR2,
                   SE_lin_m2_spring6mo_AR2,
                   numDF_lin_m3_spring6mo_AR1,
                   denDF_lin_m3_spring6mo_AR1,
                   F_lin_m3_spring6mo_AR1,
                   P_lin_m3_spring6mo_AR1,
                   ParamEst_lin_m3_spring6mo_AR1,
                   SE_lin_m3_spring6mo_AR1,
                   numDF_lin_m3_spring6mo_AR2,
                   denDF_lin_m3_spring6mo_AR2,
                   F_lin_m3_spring6mo_AR2,
                   P_lin_m3_spring6mo_AR2,
                   ParamEst_lin_m3_spring6mo_AR2,
                   SE_lin_m3_spring6mo_AR2,
                   numDF_lin_m1_monsoon6mo_AR1,
                   denDF_lin_m1_monsoon6mo_AR1,
                   F_lin_m1_monsoon6mo_AR1,
                   P_lin_m1_monsoon6mo_AR1,
                   ParamEst_lin_m1_monsoon6mo_AR1,
                   SE_lin_m1_monsoon6mo_AR1,
                   numDF_lin_m1_monsoon6mo_AR2,
                   denDF_lin_m1_monsoon6mo_AR2,
                   F_lin_m1_monsoon6mo_AR2,
                   P_lin_m1_monsoon6mo_AR2,
                   ParamEst_lin_m1_monsoon6mo_AR2,
                   SE_lin_m1_monsoon6mo_AR2,
                   numDF_lin_m2_monsoon6mo_AR1,
                   denDF_lin_m2_monsoon6mo_AR1,
                   F_lin_m2_monsoon6mo_AR1,
                   P_lin_m2_monsoon6mo_AR1,
                   ParamEst_lin_m2_monsoon6mo_AR1,
                   SE_lin_m2_monsoon6mo_AR1,
                   numDF_lin_m2_monsoon6mo_AR2,
                   denDF_lin_m2_monsoon6mo_AR2,
                   F_lin_m2_monsoon6mo_AR2,
                   P_lin_m2_monsoon6mo_AR2,
                   ParamEst_lin_m2_monsoon6mo_AR2,
                   SE_lin_m2_monsoon6mo_AR2,
                   numDF_lin_m3_monsoon6mo_AR1,
                   denDF_lin_m3_monsoon6mo_AR1,
                   F_lin_m3_monsoon6mo_AR1,
                   P_lin_m3_monsoon6mo_AR1,
                   ParamEst_lin_m3_monsoon6mo_AR1,
                   SE_lin_m3_monsoon6mo_AR1,
                   numDF_lin_m3_monsoon6mo_AR2,
                   denDF_lin_m3_monsoon6mo_AR2,
                   F_lin_m3_monsoon6mo_AR2,
                   P_lin_m3_monsoon6mo_AR2,
                   ParamEst_lin_m3_monsoon6mo_AR2,
                   SE_lin_m3_monsoon6mo_AR2,
                   numDF_lin_m1_spring6mo_lag_AR1,
                   denDF_lin_m1_spring6mo_lag_AR1,
                   F_lin_m1_spring6mo_lag_AR1,
                   P_lin_m1_spring6mo_lag_AR1,
                   ParamEst_lin_m1_spring6mo_lag_AR1,
                   SE_lin_m1_spring6mo_lag_AR1,
                   numDF_lin_m1_spring6mo_lag_AR2,
                   denDF_lin_m1_spring6mo_lag_AR2,
                   F_lin_m1_spring6mo_lag_AR2,
                   P_lin_m1_spring6mo_lag_AR2,
                   ParamEst_lin_m1_spring6mo_lag_AR2,
                   SE_lin_m1_spring6mo_lag_AR2,
                   numDF_lin_m2_spring6mo_lag_AR1,
                   denDF_lin_m2_spring6mo_lag_AR1,
                   F_lin_m2_spring6mo_lag_AR1,
                   P_lin_m2_spring6mo_lag_AR1,
                   ParamEst_lin_m2_spring6mo_lag_AR1,
                   SE_lin_m2_spring6mo_lag_AR1,
                   numDF_lin_m2_spring6mo_lag_AR2,
                   denDF_lin_m2_spring6mo_lag_AR2,
                   F_lin_m2_spring6mo_lag_AR2,
                   P_lin_m2_spring6mo_lag_AR2,
                   ParamEst_lin_m2_spring6mo_lag_AR2,
                   SE_lin_m2_spring6mo_lag_AR2,
                   numDF_lin_m3_spring6mo_lag_AR1,
                   denDF_lin_m3_spring6mo_lag_AR1,
                   F_lin_m3_spring6mo_lag_AR1,
                   P_lin_m3_spring6mo_lag_AR1,
                   ParamEst_lin_m3_spring6mo_lag_AR1,
                   SE_lin_m3_spring6mo_lag_AR1,
                   numDF_lin_m3_spring6mo_lag_AR2,
                   denDF_lin_m3_spring6mo_lag_AR2,
                   F_lin_m3_spring6mo_lag_AR2,
                   P_lin_m3_spring6mo_lag_AR2,
                   ParamEst_lin_m3_spring6mo_lag_AR2,
                   SE_lin_m3_spring6mo_lag_AR2,
                   numDF_lin_m1_monsoon6mo_lag_AR1,
                   denDF_lin_m1_monsoon6mo_lag_AR1,
                   F_lin_m1_monsoon6mo_lag_AR1,
                   P_lin_m1_monsoon6mo_lag_AR1,
                   ParamEst_lin_m1_monsoon6mo_lag_AR1,
                   SE_lin_m1_monsoon6mo_lag_AR1,
                   numDF_lin_m1_monsoon6mo_lag_AR2,
                   denDF_lin_m1_monsoon6mo_lag_AR2,
                   F_lin_m1_monsoon6mo_lag_AR2,
                   P_lin_m1_monsoon6mo_lag_AR2,
                   ParamEst_lin_m1_monsoon6mo_lag_AR2,
                   SE_lin_m1_monsoon6mo_lag_AR2,
                   numDF_lin_m2_monsoon6mo_lag_AR1,
                   denDF_lin_m2_monsoon6mo_lag_AR1,
                   F_lin_m2_monsoon6mo_lag_AR1,
                   P_lin_m2_monsoon6mo_lag_AR1,
                   ParamEst_lin_m2_monsoon6mo_lag_AR1,
                   SE_lin_m2_monsoon6mo_lag_AR1,
                   numDF_lin_m2_monsoon6mo_lag_AR2,
                   denDF_lin_m2_monsoon6mo_lag_AR2,
                   F_lin_m2_monsoon6mo_lag_AR2,
                   P_lin_m2_monsoon6mo_lag_AR2,
                   ParamEst_lin_m2_monsoon6mo_lag_AR2,
                   SE_lin_m2_monsoon6mo_lag_AR2,
                   numDF_lin_m3_monsoon6mo_lag_AR1,
                   denDF_lin_m3_monsoon6mo_lag_AR1,
                   F_lin_m3_monsoon6mo_lag_AR1,
                   P_lin_m3_monsoon6mo_lag_AR1,
                   ParamEst_lin_m3_monsoon6mo_lag_AR1,
                   SE_lin_m3_monsoon6mo_lag_AR1,
                   numDF_lin_m3_monsoon6mo_lag_AR2,
                   denDF_lin_m3_monsoon6mo_lag_AR2,
                   F_lin_m3_monsoon6mo_lag_AR2,
                   P_lin_m3_monsoon6mo_lag_AR2,
                   ParamEst_lin_m3_monsoon6mo_lag_AR2,
                   SE_lin_m3_monsoon6mo_lag_AR2,
                   Rsquared_marginal_m_null,
                   Rsquared_marginal_m1_spring6mo,
                   Rsquared_marginal_m2_spring6mo,
                   Rsquared_marginal_m3_spring6mo,
                   Rsquared_marginal_m1_spring6mo_AR1,
                   Rsquared_marginal_m1_spring6mo_AR2,
                   Rsquared_marginal_m2_spring6mo_AR1,
                   Rsquared_marginal_m2_spring6mo_AR2,
                   Rsquared_marginal_m3_spring6mo_AR1,
                   Rsquared_marginal_m3_spring6mo_AR2,
                   Rsquared_marginal_m1_monsoon6mo,
                   Rsquared_marginal_m2_monsoon6mo,
                   Rsquared_marginal_m3_monsoon6mo,
                   Rsquared_marginal_m1_monsoon6mo_AR1,
                   Rsquared_marginal_m1_monsoon6mo_AR2,
                   Rsquared_marginal_m2_monsoon6mo_AR1,
                   Rsquared_marginal_m2_monsoon6mo_AR2,
                   Rsquared_marginal_m3_monsoon6mo_AR1,
                   Rsquared_marginal_m3_monsoon6mo_AR2,
                   Rsquared_marginal_m1_spring6mo_lag,
                   Rsquared_marginal_m2_spring6mo_lag,
                   Rsquared_marginal_m3_spring6mo_lag,
                   Rsquared_marginal_m1_spring6mo_lag_AR1,
                   Rsquared_marginal_m1_spring6mo_lag_AR2,
                   Rsquared_marginal_m2_spring6mo_lag_AR1,
                   Rsquared_marginal_m2_spring6mo_lag_AR2,
                   Rsquared_marginal_m3_spring6mo_lag_AR1,
                   Rsquared_marginal_m3_spring6mo_lag_AR2,
                   Rsquared_marginal_m1_monsoon6mo_lag,
                   Rsquared_marginal_m2_monsoon6mo_lag,
                   Rsquared_marginal_m3_monsoon6mo_lag,
                   Rsquared_marginal_m1_monsoon6mo_lag_AR1,
                   Rsquared_marginal_m1_monsoon6mo_lag_AR2,
                   Rsquared_marginal_m2_monsoon6mo_lag_AR1,
                   Rsquared_marginal_m2_monsoon6mo_lag_AR2,
                   Rsquared_marginal_m3_monsoon6mo_lag_AR1,
                   Rsquared_marginal_m3_monsoon6mo_lag_AR2,
                   Rsquared_conditional_m_null,
                   Rsquared_conditional_m1_spring6mo,
                   Rsquared_conditional_m2_spring6mo,
                   Rsquared_conditional_m3_spring6mo,
                   Rsquared_conditional_m1_spring6mo_AR1,
                   Rsquared_conditional_m1_spring6mo_AR2,
                   Rsquared_conditional_m2_spring6mo_AR1,
                   Rsquared_conditional_m2_spring6mo_AR2,
                   Rsquared_conditional_m3_spring6mo_AR1,
                   Rsquared_conditional_m3_spring6mo_AR2,
                   Rsquared_conditional_m1_monsoon6mo,
                   Rsquared_conditional_m2_monsoon6mo,
                   Rsquared_conditional_m3_monsoon6mo,
                   Rsquared_conditional_m1_monsoon6mo_AR1,
                   Rsquared_conditional_m1_monsoon6mo_AR2,
                   Rsquared_conditional_m2_monsoon6mo_AR1,
                   Rsquared_conditional_m2_monsoon6mo_AR2,
                   Rsquared_conditional_m3_monsoon6mo_AR1,
                   Rsquared_conditional_m3_monsoon6mo_AR2,
                   Rsquared_conditional_m1_spring6mo_lag,
                   Rsquared_conditional_m2_spring6mo_lag,
                   Rsquared_conditional_m3_spring6mo_lag,
                   Rsquared_conditional_m1_spring6mo_lag_AR1,
                   Rsquared_conditional_m1_spring6mo_lag_AR2,
                   Rsquared_conditional_m2_spring6mo_lag_AR1,
                   Rsquared_conditional_m2_spring6mo_lag_AR2,
                   Rsquared_conditional_m3_spring6mo_lag_AR1,
                   Rsquared_conditional_m3_spring6mo_lag_AR2,
                   Rsquared_conditional_m1_monsoon6mo_lag,
                   Rsquared_conditional_m2_monsoon6mo_lag,
                   Rsquared_conditional_m3_monsoon6mo_lag,
                   Rsquared_conditional_m1_monsoon6mo_lag_AR1,
                   Rsquared_conditional_m1_monsoon6mo_lag_AR2,
                   Rsquared_conditional_m2_monsoon6mo_lag_AR1,
                   Rsquared_conditional_m2_monsoon6mo_lag_AR2,
                   Rsquared_conditional_m3_monsoon6mo_lag_AR1,
                   Rsquared_conditional_m3_monsoon6mo_lag_AR2,
                   numDF_quad_m2_spring6mo,
                   numDF_quad_m2_spring6mo_AR1,
                   numDF_quad_m2_spring6mo_AR2,
                   numDF_quad_m2_monsoon6mo,
                   numDF_quad_m2_monsoon6mo_AR1,
                   numDF_quad_m2_monsoon6mo_AR2,
                   numDF_quad_m2_spring6mo_lag,
                   numDF_quad_m2_spring6mo_lag_AR1,
                   numDF_quad_m2_spring6mo_lag_AR2,
                   numDF_quad_m2_monsoon6mo_lag,
                   numDF_quad_m2_monsoon6mo_lag_AR1,
                   numDF_quad_m2_monsoon6mo_lag_AR2,
                   numDF_quad_m3_spring6mo,
                   numDF_quad_m3_spring6mo_AR1,
                   numDF_quad_m3_spring6mo_AR2,
                   numDF_quad_m3_monsoon6mo,
                   numDF_quad_m3_monsoon6mo_AR1,
                   numDF_quad_m3_monsoon6mo_AR2,
                   numDF_quad_m3_spring6mo_lag,
                   numDF_quad_m3_spring6mo_lag_AR1,
                   numDF_quad_m3_spring6mo_lag_AR2,
                   numDF_quad_m3_monsoon6mo_lag,
                   numDF_quad_m3_monsoon6mo_lag_AR1,
                   numDF_quad_m3_monsoon6mo_lag_AR2,
                   denDF_quad_m2_spring6mo,
                   denDF_quad_m2_spring6mo_AR1,
                   denDF_quad_m2_spring6mo_AR2,
                   denDF_quad_m2_monsoon6mo,
                   denDF_quad_m2_monsoon6mo_AR1,
                   denDF_quad_m2_monsoon6mo_AR2,
                   denDF_quad_m2_spring6mo_lag,
                   denDF_quad_m2_spring6mo_lag_AR1,
                   denDF_quad_m2_spring6mo_lag_AR2,
                   denDF_quad_m2_monsoon6mo_lag,
                   denDF_quad_m2_monsoon6mo_lag_AR1,
                   denDF_quad_m2_monsoon6mo_lag_AR2,
                   denDF_quad_m3_spring6mo,
                   denDF_quad_m3_spring6mo_AR1,
                   denDF_quad_m3_spring6mo_AR2,
                   denDF_quad_m3_monsoon6mo,
                   denDF_quad_m3_monsoon6mo_AR1,
                   denDF_quad_m3_monsoon6mo_AR2,
                   denDF_quad_m3_spring6mo_lag,
                   denDF_quad_m3_spring6mo_lag_AR1,
                   denDF_quad_m3_spring6mo_lag_AR2,
                   denDF_quad_m3_monsoon6mo_lag,
                   denDF_quad_m3_monsoon6mo_lag_AR1,
                   denDF_quad_m3_monsoon6mo_lag_AR2,
                   F_quad_m2_spring6mo,
                   F_quad_m2_spring6mo_AR1,
                   F_quad_m2_spring6mo_AR2,
                   F_quad_m2_monsoon6mo,
                   F_quad_m2_monsoon6mo_AR1,
                   F_quad_m2_monsoon6mo_AR2,
                   F_quad_m2_spring6mo_lag,
                   F_quad_m2_spring6mo_lag_AR1,
                   F_quad_m2_spring6mo_lag_AR2,
                   F_quad_m2_monsoon6mo_lag,
                   F_quad_m2_monsoon6mo_lag_AR1,
                   F_quad_m2_monsoon6mo_lag_AR2,
                   F_quad_m3_spring6mo,
                   F_quad_m3_spring6mo_AR1,
                   F_quad_m3_spring6mo_AR2,
                   F_quad_m3_monsoon6mo,
                   F_quad_m3_monsoon6mo_AR1,
                   F_quad_m3_monsoon6mo_AR2,
                   F_quad_m3_spring6mo_lag,
                   F_quad_m3_spring6mo_lag_AR1,
                   F_quad_m3_spring6mo_lag_AR2,
                   F_quad_m3_monsoon6mo_lag,
                   F_quad_m3_monsoon6mo_lag_AR1,
                   F_quad_m3_monsoon6mo_lag_AR2,
                   P_quad_m2_spring6mo,
                   P_quad_m2_spring6mo_AR1,
                   P_quad_m2_spring6mo_AR2,
                   P_quad_m2_monsoon6mo,
                   P_quad_m2_monsoon6mo_AR1,
                   P_quad_m2_monsoon6mo_AR2,
                   P_quad_m2_spring6mo_lag,
                   P_quad_m2_spring6mo_lag_AR1,
                   P_quad_m2_spring6mo_lag_AR2,
                   P_quad_m2_monsoon6mo_lag,
                   P_quad_m2_monsoon6mo_lag_AR1,
                   P_quad_m2_monsoon6mo_lag_AR2,
                   P_quad_m3_spring6mo,
                   P_quad_m3_spring6mo_AR1,
                   P_quad_m3_spring6mo_AR2,
                   P_quad_m3_monsoon6mo,
                   P_quad_m3_monsoon6mo_AR1,
                   P_quad_m3_monsoon6mo_AR2,
                   P_quad_m3_spring6mo_lag,
                   P_quad_m3_spring6mo_lag_AR1,
                   P_quad_m3_spring6mo_lag_AR2,
                   P_quad_m3_monsoon6mo_lag,
                   P_quad_m3_monsoon6mo_lag_AR1,
                   P_quad_m3_monsoon6mo_lag_AR2,
                   ParamEst_quad_m2_spring6mo,
                   ParamEst_quad_m2_spring6mo_AR1,
                   ParamEst_quad_m2_spring6mo_AR2,
                   ParamEst_quad_m2_monsoon6mo,
                   ParamEst_quad_m2_monsoon6mo_AR1,
                   ParamEst_quad_m2_monsoon6mo_AR2,
                   ParamEst_quad_m2_spring6mo_lag,
                   ParamEst_quad_m2_spring6mo_lag_AR1,
                   ParamEst_quad_m2_spring6mo_lag_AR2,
                   ParamEst_quad_m2_monsoon6mo_lag,
                   ParamEst_quad_m2_monsoon6mo_lag_AR1,
                   ParamEst_quad_m2_monsoon6mo_lag_AR2,
                   ParamEst_quad_m3_spring6mo,
                   ParamEst_quad_m3_spring6mo_AR1,
                   ParamEst_quad_m3_spring6mo_AR2,
                   ParamEst_quad_m3_monsoon6mo,
                   ParamEst_quad_m3_monsoon6mo_AR1,
                   ParamEst_quad_m3_monsoon6mo_AR2,
                   ParamEst_quad_m3_spring6mo_lag,
                   ParamEst_quad_m3_spring6mo_lag_AR1,
                   ParamEst_quad_m3_spring6mo_lag_AR2,
                   ParamEst_quad_m3_monsoon6mo_lag,
                   ParamEst_quad_m3_monsoon6mo_lag_AR1,
                   ParamEst_quad_m3_monsoon6mo_lag_AR2,
                   SE_quad_m2_spring6mo,
                   SE_quad_m2_spring6mo_AR1,
                   SE_quad_m2_spring6mo_AR2,
                   SE_quad_m2_monsoon6mo,
                   SE_quad_m2_monsoon6mo_AR1,
                   SE_quad_m2_monsoon6mo_AR2,
                   SE_quad_m2_spring6mo_lag,
                   SE_quad_m2_spring6mo_lag_AR1,
                   SE_quad_m2_spring6mo_lag_AR2,
                   SE_quad_m2_monsoon6mo_lag,
                   SE_quad_m2_monsoon6mo_lag_AR1,
                   SE_quad_m2_monsoon6mo_lag_AR2,
                   SE_quad_m3_spring6mo,
                   SE_quad_m3_spring6mo_AR1,
                   SE_quad_m3_spring6mo_AR2,
                   SE_quad_m3_monsoon6mo,
                   SE_quad_m3_monsoon6mo_AR1,
                   SE_quad_m3_monsoon6mo_AR2,
                   SE_quad_m3_spring6mo_lag,
                   SE_quad_m3_spring6mo_lag_AR1,
                   SE_quad_m3_spring6mo_lag_AR2,
                   SE_quad_m3_monsoon6mo_lag,
                   SE_quad_m3_monsoon6mo_lag_AR1,
                   SE_quad_m3_monsoon6mo_lag_AR2,
                   numDF_cub_m3_spring6mo,
                   numDF_cub_m3_spring6mo_AR1,
                   numDF_cub_m3_spring6mo_AR2,
                   numDF_cub_m3_monsoon6mo,
                   numDF_cub_m3_monsoon6mo_AR1,
                   numDF_cub_m3_monsoon6mo_AR2,
                   numDF_cub_m3_spring6mo_lag,
                   numDF_cub_m3_spring6mo_lag_AR1,
                   numDF_cub_m3_spring6mo_lag_AR2,
                   numDF_cub_m3_monsoon6mo_lag,
                   numDF_cub_m3_monsoon6mo_lag_AR1,
                   numDF_cub_m3_monsoon6mo_lag_AR2,
                   denDF_cub_m3_spring6mo,
                   denDF_cub_m3_spring6mo_AR1,
                   denDF_cub_m3_spring6mo_AR2,
                   denDF_cub_m3_monsoon6mo,
                   denDF_cub_m3_monsoon6mo_AR1,
                   denDF_cub_m3_monsoon6mo_AR2,
                   denDF_cub_m3_spring6mo_lag,
                   denDF_cub_m3_spring6mo_lag_AR1,
                   denDF_cub_m3_spring6mo_lag_AR2,
                   denDF_cub_m3_monsoon6mo_lag,
                   denDF_cub_m3_monsoon6mo_lag_AR1,
                   denDF_cub_m3_monsoon6mo_lag_AR2,
                   F_cub_m3_spring6mo,
                   F_cub_m3_spring6mo_AR1,
                   F_cub_m3_spring6mo_AR2,
                   F_cub_m3_monsoon6mo,
                   F_cub_m3_monsoon6mo_AR1,
                   F_cub_m3_monsoon6mo_AR2,
                   F_cub_m3_spring6mo_lag,
                   F_cub_m3_spring6mo_lag_AR1,
                   F_cub_m3_spring6mo_lag_AR2,
                   F_cub_m3_monsoon6mo_lag,
                   F_cub_m3_monsoon6mo_lag_AR1,
                   F_cub_m3_monsoon6mo_lag_AR2,
                   P_cub_m3_spring6mo,
                   P_cub_m3_spring6mo_AR1,
                   P_cub_m3_spring6mo_AR2,
                   P_cub_m3_monsoon6mo,
                   P_cub_m3_monsoon6mo_AR1,
                   P_cub_m3_monsoon6mo_AR2,
                   P_cub_m3_spring6mo_lag,
                   P_cub_m3_spring6mo_lag_AR1,
                   P_cub_m3_spring6mo_lag_AR2,
                   P_cub_m3_monsoon6mo_lag,
                   P_cub_m3_monsoon6mo_lag_AR1,
                   P_cub_m3_monsoon6mo_lag_AR2,
                   ParamEst_cub_m3_spring6mo,
                   ParamEst_cub_m3_spring6mo_AR1,
                   ParamEst_cub_m3_spring6mo_AR2,
                   ParamEst_cub_m3_monsoon6mo,
                   ParamEst_cub_m3_monsoon6mo_AR1,
                   ParamEst_cub_m3_monsoon6mo_AR2,
                   ParamEst_cub_m3_spring6mo_lag,
                   ParamEst_cub_m3_spring6mo_lag_AR1,
                   ParamEst_cub_m3_spring6mo_lag_AR2,
                   ParamEst_cub_m3_monsoon6mo_lag,
                   ParamEst_cub_m3_monsoon6mo_lag_AR1,
                   ParamEst_cub_m3_monsoon6mo_lag_AR2,
                   SE_cub_m3_spring6mo,
                   SE_cub_m3_spring6mo_AR1,
                   SE_cub_m3_spring6mo_AR2,
                   SE_cub_m3_monsoon6mo,
                   SE_cub_m3_monsoon6mo_AR1,
                   SE_cub_m3_monsoon6mo_AR2,
                   SE_cub_m3_spring6mo_lag,
                   SE_cub_m3_spring6mo_lag_AR1,
                   SE_cub_m3_spring6mo_lag_AR2,
                   SE_cub_m3_monsoon6mo_lag,
                   SE_cub_m3_monsoon6mo_lag_AR1,
                   SE_cub_m3_monsoon6mo_lag_AR2)
  
  # append results for the bee species to the output dataset
  beeCSF_output<-rbind(beeCSF_output,output_id)
  
}



##### CSFs: Plains Grassland #####

# Create a new data frame of the original data
blue_original <- blue

# Create a data frame of just climate data
blue_climate <- blue[,5:8]

# Create a data frame of just the bee abundance matrix (descriptor variables removed)
speciesMatrix <- blue[,9:232]

# Create vector of species codes
speciesCodes <- colnames(speciesMatrix)

# Create vector of number with a number for each species
number <-1:length(speciesCodes)


### Loop through each column of speciesMatrix, running CSFs and putting the model output in the beeCSF_output matrix ###

for (i in 1:length(speciesMatrix[1,])) {
  
  # save the species code for column i
  speciesCode <- speciesCodes[i]
  
  # create an object with the name of the ecosystem type
  ecosystem <-"B"
  
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
  
  
  # normality test for null model
  shapiro_null<-shapiro.test(resid(m_null))
  P_shapiro_null<-shapiro_null$p.value
  print(P_shapiro_null)
  
  
  # calculate delta AICc for each model
  
  # delta AICc for AR models
  dAICc_m1_spring6mo_AR1<-(AICc(m1_spring6mo_AR1))-min(AICc(m_null),AICc(m1_spring6mo),AICc(m2_spring6mo),AICc(m3_spring6mo), AICc(m1_monsoon6mo), AICc(m2_monsoon6mo),AICc(m3_monsoon6mo),AICc(m1_spring6mo_lag),AICc(m2_spring6mo_lag),AICc(m3_spring6mo_lag),AICc(m1_monsoon6mo_lag),AICc(m2_monsoon6mo_lag),AICc(m3_monsoon6mo_lag),AICc(m1_spring6mo_AR1),AICc(m1_spring6mo_AR2),AICc(m2_spring6mo_AR1),AICc(m2_spring6mo_AR2),AICc(m3_spring6mo_AR1),AICc(m3_spring6mo_AR2),AICc(m1_monsoon6mo_AR1),AICc(m1_monsoon6mo_AR2),AICc(m2_monsoon6mo_AR1),AICc(m2_monsoon6mo_AR2),AICc(m3_monsoon6mo_AR1),AICc(m3_monsoon6mo_AR2),AICc(m1_spring6mo_lag_AR1),AICc(m1_spring6mo_lag_AR2),AICc(m2_spring6mo_lag_AR1),AICc(m2_spring6mo_lag_AR2),AICc(m3_spring6mo_lag_AR1),AICc(m3_spring6mo_lag_AR2),AICc(m1_monsoon6mo_lag_AR1),AICc(m1_monsoon6mo_lag_AR2),AICc(m2_monsoon6mo_lag_AR1),AICc(m2_monsoon6mo_lag_AR2),AICc(m3_monsoon6mo_lag_AR1),AICc(m3_monsoon6mo_lag_AR2))
  
  dAICc_m1_spring6mo_AR2<-(AICc(m1_spring6mo_AR2))-min(AICc(m_null),AICc(m1_spring6mo),AICc(m2_spring6mo),AICc(m3_spring6mo), AICc(m1_monsoon6mo), AICc(m2_monsoon6mo),AICc(m3_monsoon6mo),AICc(m1_spring6mo_lag),AICc(m2_spring6mo_lag),AICc(m3_spring6mo_lag),AICc(m1_monsoon6mo_lag),AICc(m2_monsoon6mo_lag),AICc(m3_monsoon6mo_lag),AICc(m1_spring6mo_AR1),AICc(m1_spring6mo_AR2),AICc(m2_spring6mo_AR1),AICc(m2_spring6mo_AR2),AICc(m3_spring6mo_AR1),AICc(m3_spring6mo_AR2),AICc(m1_monsoon6mo_AR1),AICc(m1_monsoon6mo_AR2),AICc(m2_monsoon6mo_AR1),AICc(m2_monsoon6mo_AR2),AICc(m3_monsoon6mo_AR1),AICc(m3_monsoon6mo_AR2),AICc(m1_spring6mo_lag_AR1),AICc(m1_spring6mo_lag_AR2),AICc(m2_spring6mo_lag_AR1),AICc(m2_spring6mo_lag_AR2),AICc(m3_spring6mo_lag_AR1),AICc(m3_spring6mo_lag_AR2),AICc(m1_monsoon6mo_lag_AR1),AICc(m1_monsoon6mo_lag_AR2),AICc(m2_monsoon6mo_lag_AR1),AICc(m2_monsoon6mo_lag_AR2),AICc(m3_monsoon6mo_lag_AR1),AICc(m3_monsoon6mo_lag_AR2))
  
  dAICc_m2_spring6mo_AR1<-(AICc(m2_spring6mo_AR1))-min(AICc(m_null),AICc(m1_spring6mo),AICc(m2_spring6mo),AICc(m3_spring6mo), AICc(m1_monsoon6mo), AICc(m2_monsoon6mo),AICc(m3_monsoon6mo),AICc(m1_spring6mo_lag),AICc(m2_spring6mo_lag),AICc(m3_spring6mo_lag),AICc(m1_monsoon6mo_lag),AICc(m2_monsoon6mo_lag),AICc(m3_monsoon6mo_lag),AICc(m1_spring6mo_AR1),AICc(m1_spring6mo_AR2),AICc(m2_spring6mo_AR1),AICc(m2_spring6mo_AR2),AICc(m3_spring6mo_AR1),AICc(m3_spring6mo_AR2),AICc(m1_monsoon6mo_AR1),AICc(m1_monsoon6mo_AR2),AICc(m2_monsoon6mo_AR1),AICc(m2_monsoon6mo_AR2),AICc(m3_monsoon6mo_AR1),AICc(m3_monsoon6mo_AR2),AICc(m1_spring6mo_lag_AR1),AICc(m1_spring6mo_lag_AR2),AICc(m2_spring6mo_lag_AR1),AICc(m2_spring6mo_lag_AR2),AICc(m3_spring6mo_lag_AR1),AICc(m3_spring6mo_lag_AR2),AICc(m1_monsoon6mo_lag_AR1),AICc(m1_monsoon6mo_lag_AR2),AICc(m2_monsoon6mo_lag_AR1),AICc(m2_monsoon6mo_lag_AR2),AICc(m3_monsoon6mo_lag_AR1),AICc(m3_monsoon6mo_lag_AR2))
  
  dAICc_m2_spring6mo_AR2<-(AICc(m2_spring6mo_AR2))-min(AICc(m_null),AICc(m1_spring6mo),AICc(m2_spring6mo),AICc(m3_spring6mo), AICc(m1_monsoon6mo), AICc(m2_monsoon6mo),AICc(m3_monsoon6mo),AICc(m1_spring6mo_lag),AICc(m2_spring6mo_lag),AICc(m3_spring6mo_lag),AICc(m1_monsoon6mo_lag),AICc(m2_monsoon6mo_lag),AICc(m3_monsoon6mo_lag),AICc(m1_spring6mo_AR1),AICc(m1_spring6mo_AR2),AICc(m2_spring6mo_AR1),AICc(m2_spring6mo_AR2),AICc(m3_spring6mo_AR1),AICc(m3_spring6mo_AR2),AICc(m1_monsoon6mo_AR1),AICc(m1_monsoon6mo_AR2),AICc(m2_monsoon6mo_AR1),AICc(m2_monsoon6mo_AR2),AICc(m3_monsoon6mo_AR1),AICc(m3_monsoon6mo_AR2),AICc(m1_spring6mo_lag_AR1),AICc(m1_spring6mo_lag_AR2),AICc(m2_spring6mo_lag_AR1),AICc(m2_spring6mo_lag_AR2),AICc(m3_spring6mo_lag_AR1),AICc(m3_spring6mo_lag_AR2),AICc(m1_monsoon6mo_lag_AR1),AICc(m1_monsoon6mo_lag_AR2),AICc(m2_monsoon6mo_lag_AR1),AICc(m2_monsoon6mo_lag_AR2),AICc(m3_monsoon6mo_lag_AR1),AICc(m3_monsoon6mo_lag_AR2))
  
  dAICc_m3_spring6mo_AR1<-(AICc(m3_spring6mo_AR1))-min(AICc(m_null),AICc(m1_spring6mo),AICc(m2_spring6mo),AICc(m3_spring6mo), AICc(m1_monsoon6mo), AICc(m2_monsoon6mo),AICc(m3_monsoon6mo),AICc(m1_spring6mo_lag),AICc(m2_spring6mo_lag),AICc(m3_spring6mo_lag),AICc(m1_monsoon6mo_lag),AICc(m2_monsoon6mo_lag),AICc(m3_monsoon6mo_lag),AICc(m1_spring6mo_AR1),AICc(m1_spring6mo_AR2),AICc(m2_spring6mo_AR1),AICc(m2_spring6mo_AR2),AICc(m3_spring6mo_AR1),AICc(m3_spring6mo_AR2),AICc(m1_monsoon6mo_AR1),AICc(m1_monsoon6mo_AR2),AICc(m2_monsoon6mo_AR1),AICc(m2_monsoon6mo_AR2),AICc(m3_monsoon6mo_AR1),AICc(m3_monsoon6mo_AR2),AICc(m1_spring6mo_lag_AR1),AICc(m1_spring6mo_lag_AR2),AICc(m2_spring6mo_lag_AR1),AICc(m2_spring6mo_lag_AR2),AICc(m3_spring6mo_lag_AR1),AICc(m3_spring6mo_lag_AR2),AICc(m1_monsoon6mo_lag_AR1),AICc(m1_monsoon6mo_lag_AR2),AICc(m2_monsoon6mo_lag_AR1),AICc(m2_monsoon6mo_lag_AR2),AICc(m3_monsoon6mo_lag_AR1),AICc(m3_monsoon6mo_lag_AR2))
  
  dAICc_m3_spring6mo_AR2<-(AICc(m3_spring6mo_AR2))-min(AICc(m_null),AICc(m1_spring6mo),AICc(m2_spring6mo),AICc(m3_spring6mo), AICc(m1_monsoon6mo), AICc(m2_monsoon6mo),AICc(m3_monsoon6mo),AICc(m1_spring6mo_lag),AICc(m2_spring6mo_lag),AICc(m3_spring6mo_lag),AICc(m1_monsoon6mo_lag),AICc(m2_monsoon6mo_lag),AICc(m3_monsoon6mo_lag),AICc(m1_spring6mo_AR1),AICc(m1_spring6mo_AR2),AICc(m2_spring6mo_AR1),AICc(m2_spring6mo_AR2),AICc(m3_spring6mo_AR1),AICc(m3_spring6mo_AR2),AICc(m1_monsoon6mo_AR1),AICc(m1_monsoon6mo_AR2),AICc(m2_monsoon6mo_AR1),AICc(m2_monsoon6mo_AR2),AICc(m3_monsoon6mo_AR1),AICc(m3_monsoon6mo_AR2),AICc(m1_spring6mo_lag_AR1),AICc(m1_spring6mo_lag_AR2),AICc(m2_spring6mo_lag_AR1),AICc(m2_spring6mo_lag_AR2),AICc(m3_spring6mo_lag_AR1),AICc(m3_spring6mo_lag_AR2),AICc(m1_monsoon6mo_lag_AR1),AICc(m1_monsoon6mo_lag_AR2),AICc(m2_monsoon6mo_lag_AR1),AICc(m2_monsoon6mo_lag_AR2),AICc(m3_monsoon6mo_lag_AR1),AICc(m3_monsoon6mo_lag_AR2))
  
  dAICc_m1_monsoon6mo_AR1<-(AICc(m1_monsoon6mo_AR1))-min(AICc(m_null),AICc(m1_spring6mo),AICc(m2_spring6mo),AICc(m3_spring6mo), AICc(m1_monsoon6mo), AICc(m2_monsoon6mo),AICc(m3_monsoon6mo),AICc(m1_spring6mo_lag),AICc(m2_spring6mo_lag),AICc(m3_spring6mo_lag),AICc(m1_monsoon6mo_lag),AICc(m2_monsoon6mo_lag),AICc(m3_monsoon6mo_lag),AICc(m1_spring6mo_AR1),AICc(m1_spring6mo_AR2),AICc(m2_spring6mo_AR1),AICc(m2_spring6mo_AR2),AICc(m3_spring6mo_AR1),AICc(m3_spring6mo_AR2),AICc(m1_monsoon6mo_AR1),AICc(m1_monsoon6mo_AR2),AICc(m2_monsoon6mo_AR1),AICc(m2_monsoon6mo_AR2),AICc(m3_monsoon6mo_AR1),AICc(m3_monsoon6mo_AR2),AICc(m1_spring6mo_lag_AR1),AICc(m1_spring6mo_lag_AR2),AICc(m2_spring6mo_lag_AR1),AICc(m2_spring6mo_lag_AR2),AICc(m3_spring6mo_lag_AR1),AICc(m3_spring6mo_lag_AR2),AICc(m1_monsoon6mo_lag_AR1),AICc(m1_monsoon6mo_lag_AR2),AICc(m2_monsoon6mo_lag_AR1),AICc(m2_monsoon6mo_lag_AR2),AICc(m3_monsoon6mo_lag_AR1),AICc(m3_monsoon6mo_lag_AR2))
  
  dAICc_m1_monsoon6mo_AR2<-(AICc(m1_monsoon6mo_AR2))-min(AICc(m_null),AICc(m1_spring6mo),AICc(m2_spring6mo),AICc(m3_spring6mo), AICc(m1_monsoon6mo), AICc(m2_monsoon6mo),AICc(m3_monsoon6mo),AICc(m1_spring6mo_lag),AICc(m2_spring6mo_lag),AICc(m3_spring6mo_lag),AICc(m1_monsoon6mo_lag),AICc(m2_monsoon6mo_lag),AICc(m3_monsoon6mo_lag),AICc(m1_spring6mo_AR1),AICc(m1_spring6mo_AR2),AICc(m2_spring6mo_AR1),AICc(m2_spring6mo_AR2),AICc(m3_spring6mo_AR1),AICc(m3_spring6mo_AR2),AICc(m1_monsoon6mo_AR1),AICc(m1_monsoon6mo_AR2),AICc(m2_monsoon6mo_AR1),AICc(m2_monsoon6mo_AR2),AICc(m3_monsoon6mo_AR1),AICc(m3_monsoon6mo_AR2),AICc(m1_spring6mo_lag_AR1),AICc(m1_spring6mo_lag_AR2),AICc(m2_spring6mo_lag_AR1),AICc(m2_spring6mo_lag_AR2),AICc(m3_spring6mo_lag_AR1),AICc(m3_spring6mo_lag_AR2),AICc(m1_monsoon6mo_lag_AR1),AICc(m1_monsoon6mo_lag_AR2),AICc(m2_monsoon6mo_lag_AR1),AICc(m2_monsoon6mo_lag_AR2),AICc(m3_monsoon6mo_lag_AR1),AICc(m3_monsoon6mo_lag_AR2))
  
  dAICc_m2_monsoon6mo_AR1<-(AICc(m2_monsoon6mo_AR1))-min(AICc(m_null),AICc(m1_spring6mo),AICc(m2_spring6mo),AICc(m3_spring6mo), AICc(m1_monsoon6mo), AICc(m2_monsoon6mo),AICc(m3_monsoon6mo),AICc(m1_spring6mo_lag),AICc(m2_spring6mo_lag),AICc(m3_spring6mo_lag),AICc(m1_monsoon6mo_lag),AICc(m2_monsoon6mo_lag),AICc(m3_monsoon6mo_lag),AICc(m1_spring6mo_AR1),AICc(m1_spring6mo_AR2),AICc(m2_spring6mo_AR1),AICc(m2_spring6mo_AR2),AICc(m3_spring6mo_AR1),AICc(m3_spring6mo_AR2),AICc(m1_monsoon6mo_AR1),AICc(m1_monsoon6mo_AR2),AICc(m2_monsoon6mo_AR1),AICc(m2_monsoon6mo_AR2),AICc(m3_monsoon6mo_AR1),AICc(m3_monsoon6mo_AR2),AICc(m1_spring6mo_lag_AR1),AICc(m1_spring6mo_lag_AR2),AICc(m2_spring6mo_lag_AR1),AICc(m2_spring6mo_lag_AR2),AICc(m3_spring6mo_lag_AR1),AICc(m3_spring6mo_lag_AR2),AICc(m1_monsoon6mo_lag_AR1),AICc(m1_monsoon6mo_lag_AR2),AICc(m2_monsoon6mo_lag_AR1),AICc(m2_monsoon6mo_lag_AR2),AICc(m3_monsoon6mo_lag_AR1),AICc(m3_monsoon6mo_lag_AR2))
  
  dAICc_m2_monsoon6mo_AR2<-(AICc(m2_monsoon6mo_AR2))-min(AICc(m_null),AICc(m1_spring6mo),AICc(m2_spring6mo),AICc(m3_spring6mo), AICc(m1_monsoon6mo), AICc(m2_monsoon6mo),AICc(m3_monsoon6mo),AICc(m1_spring6mo_lag),AICc(m2_spring6mo_lag),AICc(m3_spring6mo_lag),AICc(m1_monsoon6mo_lag),AICc(m2_monsoon6mo_lag),AICc(m3_monsoon6mo_lag),AICc(m1_spring6mo_AR1),AICc(m1_spring6mo_AR2),AICc(m2_spring6mo_AR1),AICc(m2_spring6mo_AR2),AICc(m3_spring6mo_AR1),AICc(m3_spring6mo_AR2),AICc(m1_monsoon6mo_AR1),AICc(m1_monsoon6mo_AR2),AICc(m2_monsoon6mo_AR1),AICc(m2_monsoon6mo_AR2),AICc(m3_monsoon6mo_AR1),AICc(m3_monsoon6mo_AR2),AICc(m1_spring6mo_lag_AR1),AICc(m1_spring6mo_lag_AR2),AICc(m2_spring6mo_lag_AR1),AICc(m2_spring6mo_lag_AR2),AICc(m3_spring6mo_lag_AR1),AICc(m3_spring6mo_lag_AR2),AICc(m1_monsoon6mo_lag_AR1),AICc(m1_monsoon6mo_lag_AR2),AICc(m2_monsoon6mo_lag_AR1),AICc(m2_monsoon6mo_lag_AR2),AICc(m3_monsoon6mo_lag_AR1),AICc(m3_monsoon6mo_lag_AR2))
  
  dAICc_m3_monsoon6mo_AR1<-(AICc(m3_monsoon6mo_AR1))-min(AICc(m_null),AICc(m1_spring6mo),AICc(m2_spring6mo),AICc(m3_spring6mo), AICc(m1_monsoon6mo), AICc(m2_monsoon6mo),AICc(m3_monsoon6mo),AICc(m1_spring6mo_lag),AICc(m2_spring6mo_lag),AICc(m3_spring6mo_lag),AICc(m1_monsoon6mo_lag),AICc(m2_monsoon6mo_lag),AICc(m3_monsoon6mo_lag),AICc(m1_spring6mo_AR1),AICc(m1_spring6mo_AR2),AICc(m2_spring6mo_AR1),AICc(m2_spring6mo_AR2),AICc(m3_spring6mo_AR1),AICc(m3_spring6mo_AR2),AICc(m1_monsoon6mo_AR1),AICc(m1_monsoon6mo_AR2),AICc(m2_monsoon6mo_AR1),AICc(m2_monsoon6mo_AR2),AICc(m3_monsoon6mo_AR1),AICc(m3_monsoon6mo_AR2),AICc(m1_spring6mo_lag_AR1),AICc(m1_spring6mo_lag_AR2),AICc(m2_spring6mo_lag_AR1),AICc(m2_spring6mo_lag_AR2),AICc(m3_spring6mo_lag_AR1),AICc(m3_spring6mo_lag_AR2),AICc(m1_monsoon6mo_lag_AR1),AICc(m1_monsoon6mo_lag_AR2),AICc(m2_monsoon6mo_lag_AR1),AICc(m2_monsoon6mo_lag_AR2),AICc(m3_monsoon6mo_lag_AR1),AICc(m3_monsoon6mo_lag_AR2))
  
  dAICc_m3_monsoon6mo_AR2<-(AICc(m3_monsoon6mo_AR2))-min(AICc(m_null),AICc(m1_spring6mo),AICc(m2_spring6mo),AICc(m3_spring6mo), AICc(m1_monsoon6mo), AICc(m2_monsoon6mo),AICc(m3_monsoon6mo),AICc(m1_spring6mo_lag),AICc(m2_spring6mo_lag),AICc(m3_spring6mo_lag),AICc(m1_monsoon6mo_lag),AICc(m2_monsoon6mo_lag),AICc(m3_monsoon6mo_lag),AICc(m1_spring6mo_AR1),AICc(m1_spring6mo_AR2),AICc(m2_spring6mo_AR1),AICc(m2_spring6mo_AR2),AICc(m3_spring6mo_AR1),AICc(m3_spring6mo_AR2),AICc(m1_monsoon6mo_AR1),AICc(m1_monsoon6mo_AR2),AICc(m2_monsoon6mo_AR1),AICc(m2_monsoon6mo_AR2),AICc(m3_monsoon6mo_AR1),AICc(m3_monsoon6mo_AR2),AICc(m1_spring6mo_lag_AR1),AICc(m1_spring6mo_lag_AR2),AICc(m2_spring6mo_lag_AR1),AICc(m2_spring6mo_lag_AR2),AICc(m3_spring6mo_lag_AR1),AICc(m3_spring6mo_lag_AR2),AICc(m1_monsoon6mo_lag_AR1),AICc(m1_monsoon6mo_lag_AR2),AICc(m2_monsoon6mo_lag_AR1),AICc(m2_monsoon6mo_lag_AR2),AICc(m3_monsoon6mo_lag_AR1),AICc(m3_monsoon6mo_lag_AR2))
  
  dAICc_m1_spring6mo_lag_AR1<-(AICc(m1_spring6mo_lag_AR1))-min(AICc(m_null),AICc(m1_spring6mo),AICc(m2_spring6mo),AICc(m3_spring6mo), AICc(m1_monsoon6mo), AICc(m2_monsoon6mo),AICc(m3_monsoon6mo),AICc(m1_spring6mo_lag),AICc(m2_spring6mo_lag),AICc(m3_spring6mo_lag),AICc(m1_monsoon6mo_lag),AICc(m2_monsoon6mo_lag),AICc(m3_monsoon6mo_lag),AICc(m1_spring6mo_AR1),AICc(m1_spring6mo_AR2),AICc(m2_spring6mo_AR1),AICc(m2_spring6mo_AR2),AICc(m3_spring6mo_AR1),AICc(m3_spring6mo_AR2),AICc(m1_monsoon6mo_AR1),AICc(m1_monsoon6mo_AR2),AICc(m2_monsoon6mo_AR1),AICc(m2_monsoon6mo_AR2),AICc(m3_monsoon6mo_AR1),AICc(m3_monsoon6mo_AR2),AICc(m1_spring6mo_lag_AR1),AICc(m1_spring6mo_lag_AR2),AICc(m2_spring6mo_lag_AR1),AICc(m2_spring6mo_lag_AR2),AICc(m3_spring6mo_lag_AR1),AICc(m3_spring6mo_lag_AR2),AICc(m1_monsoon6mo_lag_AR1),AICc(m1_monsoon6mo_lag_AR2),AICc(m2_monsoon6mo_lag_AR1),AICc(m2_monsoon6mo_lag_AR2),AICc(m3_monsoon6mo_lag_AR1),AICc(m3_monsoon6mo_lag_AR2))
  
  dAICc_m1_spring6mo_lag_AR2<-(AICc(m1_spring6mo_lag_AR2))-min(AICc(m_null),AICc(m1_spring6mo),AICc(m2_spring6mo),AICc(m3_spring6mo), AICc(m1_monsoon6mo), AICc(m2_monsoon6mo),AICc(m3_monsoon6mo),AICc(m1_spring6mo_lag),AICc(m2_spring6mo_lag),AICc(m3_spring6mo_lag),AICc(m1_monsoon6mo_lag),AICc(m2_monsoon6mo_lag),AICc(m3_monsoon6mo_lag),AICc(m1_spring6mo_AR1),AICc(m1_spring6mo_AR2),AICc(m2_spring6mo_AR1),AICc(m2_spring6mo_AR2),AICc(m3_spring6mo_AR1),AICc(m3_spring6mo_AR2),AICc(m1_monsoon6mo_AR1),AICc(m1_monsoon6mo_AR2),AICc(m2_monsoon6mo_AR1),AICc(m2_monsoon6mo_AR2),AICc(m3_monsoon6mo_AR1),AICc(m3_monsoon6mo_AR2),AICc(m1_spring6mo_lag_AR1),AICc(m1_spring6mo_lag_AR2),AICc(m2_spring6mo_lag_AR1),AICc(m2_spring6mo_lag_AR2),AICc(m3_spring6mo_lag_AR1),AICc(m3_spring6mo_lag_AR2),AICc(m1_monsoon6mo_lag_AR1),AICc(m1_monsoon6mo_lag_AR2),AICc(m2_monsoon6mo_lag_AR1),AICc(m2_monsoon6mo_lag_AR2),AICc(m3_monsoon6mo_lag_AR1),AICc(m3_monsoon6mo_lag_AR2))
  
  dAICc_m2_spring6mo_lag_AR1<-(AICc(m2_spring6mo_lag_AR1))-min(AICc(m_null),AICc(m1_spring6mo),AICc(m2_spring6mo),AICc(m3_spring6mo), AICc(m1_monsoon6mo), AICc(m2_monsoon6mo),AICc(m3_monsoon6mo),AICc(m1_spring6mo_lag),AICc(m2_spring6mo_lag),AICc(m3_spring6mo_lag),AICc(m1_monsoon6mo_lag),AICc(m2_monsoon6mo_lag),AICc(m3_monsoon6mo_lag),AICc(m1_spring6mo_AR1),AICc(m1_spring6mo_AR2),AICc(m2_spring6mo_AR1),AICc(m2_spring6mo_AR2),AICc(m3_spring6mo_AR1),AICc(m3_spring6mo_AR2),AICc(m1_monsoon6mo_AR1),AICc(m1_monsoon6mo_AR2),AICc(m2_monsoon6mo_AR1),AICc(m2_monsoon6mo_AR2),AICc(m3_monsoon6mo_AR1),AICc(m3_monsoon6mo_AR2),AICc(m1_spring6mo_lag_AR1),AICc(m1_spring6mo_lag_AR2),AICc(m2_spring6mo_lag_AR1),AICc(m2_spring6mo_lag_AR2),AICc(m3_spring6mo_lag_AR1),AICc(m3_spring6mo_lag_AR2),AICc(m1_monsoon6mo_lag_AR1),AICc(m1_monsoon6mo_lag_AR2),AICc(m2_monsoon6mo_lag_AR1),AICc(m2_monsoon6mo_lag_AR2),AICc(m3_monsoon6mo_lag_AR1),AICc(m3_monsoon6mo_lag_AR2))
  
  dAICc_m2_spring6mo_lag_AR2<-(AICc(m2_spring6mo_lag_AR2))-min(AICc(m_null),AICc(m1_spring6mo),AICc(m2_spring6mo),AICc(m3_spring6mo), AICc(m1_monsoon6mo), AICc(m2_monsoon6mo),AICc(m3_monsoon6mo),AICc(m1_spring6mo_lag),AICc(m2_spring6mo_lag),AICc(m3_spring6mo_lag),AICc(m1_monsoon6mo_lag),AICc(m2_monsoon6mo_lag),AICc(m3_monsoon6mo_lag),AICc(m1_spring6mo_AR1),AICc(m1_spring6mo_AR2),AICc(m2_spring6mo_AR1),AICc(m2_spring6mo_AR2),AICc(m3_spring6mo_AR1),AICc(m3_spring6mo_AR2),AICc(m1_monsoon6mo_AR1),AICc(m1_monsoon6mo_AR2),AICc(m2_monsoon6mo_AR1),AICc(m2_monsoon6mo_AR2),AICc(m3_monsoon6mo_AR1),AICc(m3_monsoon6mo_AR2),AICc(m1_spring6mo_lag_AR1),AICc(m1_spring6mo_lag_AR2),AICc(m2_spring6mo_lag_AR1),AICc(m2_spring6mo_lag_AR2),AICc(m3_spring6mo_lag_AR1),AICc(m3_spring6mo_lag_AR2),AICc(m1_monsoon6mo_lag_AR1),AICc(m1_monsoon6mo_lag_AR2),AICc(m2_monsoon6mo_lag_AR1),AICc(m2_monsoon6mo_lag_AR2),AICc(m3_monsoon6mo_lag_AR1),AICc(m3_monsoon6mo_lag_AR2))
  
  dAICc_m3_spring6mo_lag_AR1<-(AICc(m3_spring6mo_lag_AR1))-min(AICc(m_null),AICc(m1_spring6mo),AICc(m2_spring6mo),AICc(m3_spring6mo), AICc(m1_monsoon6mo), AICc(m2_monsoon6mo),AICc(m3_monsoon6mo),AICc(m1_spring6mo_lag),AICc(m2_spring6mo_lag),AICc(m3_spring6mo_lag),AICc(m1_monsoon6mo_lag),AICc(m2_monsoon6mo_lag),AICc(m3_monsoon6mo_lag),AICc(m1_spring6mo_AR1),AICc(m1_spring6mo_AR2),AICc(m2_spring6mo_AR1),AICc(m2_spring6mo_AR2),AICc(m3_spring6mo_AR1),AICc(m3_spring6mo_AR2),AICc(m1_monsoon6mo_AR1),AICc(m1_monsoon6mo_AR2),AICc(m2_monsoon6mo_AR1),AICc(m2_monsoon6mo_AR2),AICc(m3_monsoon6mo_AR1),AICc(m3_monsoon6mo_AR2),AICc(m1_spring6mo_lag_AR1),AICc(m1_spring6mo_lag_AR2),AICc(m2_spring6mo_lag_AR1),AICc(m2_spring6mo_lag_AR2),AICc(m3_spring6mo_lag_AR1),AICc(m3_spring6mo_lag_AR2),AICc(m1_monsoon6mo_lag_AR1),AICc(m1_monsoon6mo_lag_AR2),AICc(m2_monsoon6mo_lag_AR1),AICc(m2_monsoon6mo_lag_AR2),AICc(m3_monsoon6mo_lag_AR1),AICc(m3_monsoon6mo_lag_AR2))
  
  dAICc_m3_spring6mo_lag_AR2<-(AICc(m3_spring6mo_lag_AR2))-min(AICc(m_null),AICc(m1_spring6mo),AICc(m2_spring6mo),AICc(m3_spring6mo), AICc(m1_monsoon6mo), AICc(m2_monsoon6mo),AICc(m3_monsoon6mo),AICc(m1_spring6mo_lag),AICc(m2_spring6mo_lag),AICc(m3_spring6mo_lag),AICc(m1_monsoon6mo_lag),AICc(m2_monsoon6mo_lag),AICc(m3_monsoon6mo_lag),AICc(m1_spring6mo_AR1),AICc(m1_spring6mo_AR2),AICc(m2_spring6mo_AR1),AICc(m2_spring6mo_AR2),AICc(m3_spring6mo_AR1),AICc(m3_spring6mo_AR2),AICc(m1_monsoon6mo_AR1),AICc(m1_monsoon6mo_AR2),AICc(m2_monsoon6mo_AR1),AICc(m2_monsoon6mo_AR2),AICc(m3_monsoon6mo_AR1),AICc(m3_monsoon6mo_AR2),AICc(m1_spring6mo_lag_AR1),AICc(m1_spring6mo_lag_AR2),AICc(m2_spring6mo_lag_AR1),AICc(m2_spring6mo_lag_AR2),AICc(m3_spring6mo_lag_AR1),AICc(m3_spring6mo_lag_AR2),AICc(m1_monsoon6mo_lag_AR1),AICc(m1_monsoon6mo_lag_AR2),AICc(m2_monsoon6mo_lag_AR1),AICc(m2_monsoon6mo_lag_AR2),AICc(m3_monsoon6mo_lag_AR1),AICc(m3_monsoon6mo_lag_AR2))
  
  dAICc_m1_monsoon6mo_lag_AR1<-(AICc(m1_monsoon6mo_lag_AR1))-min(AICc(m_null),AICc(m1_spring6mo),AICc(m2_spring6mo),AICc(m3_spring6mo), AICc(m1_monsoon6mo), AICc(m2_monsoon6mo),AICc(m3_monsoon6mo),AICc(m1_spring6mo_lag),AICc(m2_spring6mo_lag),AICc(m3_spring6mo_lag),AICc(m1_monsoon6mo_lag),AICc(m2_monsoon6mo_lag),AICc(m3_monsoon6mo_lag),AICc(m1_spring6mo_AR1),AICc(m1_spring6mo_AR2),AICc(m2_spring6mo_AR1),AICc(m2_spring6mo_AR2),AICc(m3_spring6mo_AR1),AICc(m3_spring6mo_AR2),AICc(m1_monsoon6mo_AR1),AICc(m1_monsoon6mo_AR2),AICc(m2_monsoon6mo_AR1),AICc(m2_monsoon6mo_AR2),AICc(m3_monsoon6mo_AR1),AICc(m3_monsoon6mo_AR2),AICc(m1_spring6mo_lag_AR1),AICc(m1_spring6mo_lag_AR2),AICc(m2_spring6mo_lag_AR1),AICc(m2_spring6mo_lag_AR2),AICc(m3_spring6mo_lag_AR1),AICc(m3_spring6mo_lag_AR2),AICc(m1_monsoon6mo_lag_AR1),AICc(m1_monsoon6mo_lag_AR2),AICc(m2_monsoon6mo_lag_AR1),AICc(m2_monsoon6mo_lag_AR2),AICc(m3_monsoon6mo_lag_AR1),AICc(m3_monsoon6mo_lag_AR2))
  
  dAICc_m1_monsoon6mo_lag_AR2<-(AICc(m1_monsoon6mo_lag_AR2))-min(AICc(m_null),AICc(m1_spring6mo),AICc(m2_spring6mo),AICc(m3_spring6mo), AICc(m1_monsoon6mo), AICc(m2_monsoon6mo),AICc(m3_monsoon6mo),AICc(m1_spring6mo_lag),AICc(m2_spring6mo_lag),AICc(m3_spring6mo_lag),AICc(m1_monsoon6mo_lag),AICc(m2_monsoon6mo_lag),AICc(m3_monsoon6mo_lag),AICc(m1_spring6mo_AR1),AICc(m1_spring6mo_AR2),AICc(m2_spring6mo_AR1),AICc(m2_spring6mo_AR2),AICc(m3_spring6mo_AR1),AICc(m3_spring6mo_AR2),AICc(m1_monsoon6mo_AR1),AICc(m1_monsoon6mo_AR2),AICc(m2_monsoon6mo_AR1),AICc(m2_monsoon6mo_AR2),AICc(m3_monsoon6mo_AR1),AICc(m3_monsoon6mo_AR2),AICc(m1_spring6mo_lag_AR1),AICc(m1_spring6mo_lag_AR2),AICc(m2_spring6mo_lag_AR1),AICc(m2_spring6mo_lag_AR2),AICc(m3_spring6mo_lag_AR1),AICc(m3_spring6mo_lag_AR2),AICc(m1_monsoon6mo_lag_AR1),AICc(m1_monsoon6mo_lag_AR2),AICc(m2_monsoon6mo_lag_AR1),AICc(m2_monsoon6mo_lag_AR2),AICc(m3_monsoon6mo_lag_AR1),AICc(m3_monsoon6mo_lag_AR2))
  
  dAICc_m2_monsoon6mo_lag_AR1<-(AICc(m2_monsoon6mo_lag_AR1))-min(AICc(m_null),AICc(m1_spring6mo),AICc(m2_spring6mo),AICc(m3_spring6mo), AICc(m1_monsoon6mo), AICc(m2_monsoon6mo),AICc(m3_monsoon6mo),AICc(m1_spring6mo_lag),AICc(m2_spring6mo_lag),AICc(m3_spring6mo_lag),AICc(m1_monsoon6mo_lag),AICc(m2_monsoon6mo_lag),AICc(m3_monsoon6mo_lag),AICc(m1_spring6mo_AR1),AICc(m1_spring6mo_AR2),AICc(m2_spring6mo_AR1),AICc(m2_spring6mo_AR2),AICc(m3_spring6mo_AR1),AICc(m3_spring6mo_AR2),AICc(m1_monsoon6mo_AR1),AICc(m1_monsoon6mo_AR2),AICc(m2_monsoon6mo_AR1),AICc(m2_monsoon6mo_AR2),AICc(m3_monsoon6mo_AR1),AICc(m3_monsoon6mo_AR2),AICc(m1_spring6mo_lag_AR1),AICc(m1_spring6mo_lag_AR2),AICc(m2_spring6mo_lag_AR1),AICc(m2_spring6mo_lag_AR2),AICc(m3_spring6mo_lag_AR1),AICc(m3_spring6mo_lag_AR2),AICc(m1_monsoon6mo_lag_AR1),AICc(m1_monsoon6mo_lag_AR2),AICc(m2_monsoon6mo_lag_AR1),AICc(m2_monsoon6mo_lag_AR2),AICc(m3_monsoon6mo_lag_AR1),AICc(m3_monsoon6mo_lag_AR2))
  
  dAICc_m2_monsoon6mo_lag_AR2<-(AICc(m2_monsoon6mo_lag_AR2))-min(AICc(m_null),AICc(m1_spring6mo),AICc(m2_spring6mo),AICc(m3_spring6mo), AICc(m1_monsoon6mo), AICc(m2_monsoon6mo),AICc(m3_monsoon6mo),AICc(m1_spring6mo_lag),AICc(m2_spring6mo_lag),AICc(m3_spring6mo_lag),AICc(m1_monsoon6mo_lag),AICc(m2_monsoon6mo_lag),AICc(m3_monsoon6mo_lag),AICc(m1_spring6mo_AR1),AICc(m1_spring6mo_AR2),AICc(m2_spring6mo_AR1),AICc(m2_spring6mo_AR2),AICc(m3_spring6mo_AR1),AICc(m3_spring6mo_AR2),AICc(m1_monsoon6mo_AR1),AICc(m1_monsoon6mo_AR2),AICc(m2_monsoon6mo_AR1),AICc(m2_monsoon6mo_AR2),AICc(m3_monsoon6mo_AR1),AICc(m3_monsoon6mo_AR2),AICc(m1_spring6mo_lag_AR1),AICc(m1_spring6mo_lag_AR2),AICc(m2_spring6mo_lag_AR1),AICc(m2_spring6mo_lag_AR2),AICc(m3_spring6mo_lag_AR1),AICc(m3_spring6mo_lag_AR2),AICc(m1_monsoon6mo_lag_AR1),AICc(m1_monsoon6mo_lag_AR2),AICc(m2_monsoon6mo_lag_AR1),AICc(m2_monsoon6mo_lag_AR2),AICc(m3_monsoon6mo_lag_AR1),AICc(m3_monsoon6mo_lag_AR2))
  
  dAICc_m3_monsoon6mo_lag_AR1<-(AICc(m3_monsoon6mo_lag_AR1))-min(AICc(m_null),AICc(m1_spring6mo),AICc(m2_spring6mo),AICc(m3_spring6mo), AICc(m1_monsoon6mo), AICc(m2_monsoon6mo),AICc(m3_monsoon6mo),AICc(m1_spring6mo_lag),AICc(m2_spring6mo_lag),AICc(m3_spring6mo_lag),AICc(m1_monsoon6mo_lag),AICc(m2_monsoon6mo_lag),AICc(m3_monsoon6mo_lag),AICc(m1_spring6mo_AR1),AICc(m1_spring6mo_AR2),AICc(m2_spring6mo_AR1),AICc(m2_spring6mo_AR2),AICc(m3_spring6mo_AR1),AICc(m3_spring6mo_AR2),AICc(m1_monsoon6mo_AR1),AICc(m1_monsoon6mo_AR2),AICc(m2_monsoon6mo_AR1),AICc(m2_monsoon6mo_AR2),AICc(m3_monsoon6mo_AR1),AICc(m3_monsoon6mo_AR2),AICc(m1_spring6mo_lag_AR1),AICc(m1_spring6mo_lag_AR2),AICc(m2_spring6mo_lag_AR1),AICc(m2_spring6mo_lag_AR2),AICc(m3_spring6mo_lag_AR1),AICc(m3_spring6mo_lag_AR2),AICc(m1_monsoon6mo_lag_AR1),AICc(m1_monsoon6mo_lag_AR2),AICc(m2_monsoon6mo_lag_AR1),AICc(m2_monsoon6mo_lag_AR2),AICc(m3_monsoon6mo_lag_AR1),AICc(m3_monsoon6mo_lag_AR2))
  
  dAICc_m3_monsoon6mo_lag_AR2<-(AICc(m3_monsoon6mo_lag_AR2))-min(AICc(m_null),AICc(m1_spring6mo),AICc(m2_spring6mo),AICc(m3_spring6mo), AICc(m1_monsoon6mo), AICc(m2_monsoon6mo),AICc(m3_monsoon6mo),AICc(m1_spring6mo_lag),AICc(m2_spring6mo_lag),AICc(m3_spring6mo_lag),AICc(m1_monsoon6mo_lag),AICc(m2_monsoon6mo_lag),AICc(m3_monsoon6mo_lag),AICc(m1_spring6mo_AR1),AICc(m1_spring6mo_AR2),AICc(m2_spring6mo_AR1),AICc(m2_spring6mo_AR2),AICc(m3_spring6mo_AR1),AICc(m3_spring6mo_AR2),AICc(m1_monsoon6mo_AR1),AICc(m1_monsoon6mo_AR2),AICc(m2_monsoon6mo_AR1),AICc(m2_monsoon6mo_AR2),AICc(m3_monsoon6mo_AR1),AICc(m3_monsoon6mo_AR2),AICc(m1_spring6mo_lag_AR1),AICc(m1_spring6mo_lag_AR2),AICc(m2_spring6mo_lag_AR1),AICc(m2_spring6mo_lag_AR2),AICc(m3_spring6mo_lag_AR1),AICc(m3_spring6mo_lag_AR2),AICc(m1_monsoon6mo_lag_AR1),AICc(m1_monsoon6mo_lag_AR2),AICc(m2_monsoon6mo_lag_AR1),AICc(m2_monsoon6mo_lag_AR2),AICc(m3_monsoon6mo_lag_AR1),AICc(m3_monsoon6mo_lag_AR2))
  
  # delta AICc for other models
  dAICc_null<-(AICc(m_null))-min(AICc(m_null),AICc(m1_spring6mo),AICc(m2_spring6mo),AICc(m3_spring6mo), AICc(m1_monsoon6mo), AICc(m2_monsoon6mo),AICc(m3_monsoon6mo),AICc(m1_spring6mo_lag),AICc(m2_spring6mo_lag),AICc(m3_spring6mo_lag),AICc(m1_monsoon6mo_lag),AICc(m2_monsoon6mo_lag),AICc(m3_monsoon6mo_lag),AICc(m1_spring6mo_AR1),AICc(m1_spring6mo_AR2),AICc(m2_spring6mo_AR1),AICc(m2_spring6mo_AR2),AICc(m3_spring6mo_AR1),AICc(m3_spring6mo_AR2),AICc(m1_monsoon6mo_AR1),AICc(m1_monsoon6mo_AR2),AICc(m2_monsoon6mo_AR1),AICc(m2_monsoon6mo_AR2),AICc(m3_monsoon6mo_AR1),AICc(m3_monsoon6mo_AR2),AICc(m1_spring6mo_lag_AR1),AICc(m1_spring6mo_lag_AR2),AICc(m2_spring6mo_lag_AR1),AICc(m2_spring6mo_lag_AR2),AICc(m3_spring6mo_lag_AR1),AICc(m3_spring6mo_lag_AR2),AICc(m1_monsoon6mo_lag_AR1),AICc(m1_monsoon6mo_lag_AR2),AICc(m2_monsoon6mo_lag_AR1),AICc(m2_monsoon6mo_lag_AR2),AICc(m3_monsoon6mo_lag_AR1),AICc(m3_monsoon6mo_lag_AR2))
  
  dAICc_m1_spring6mo<-(AICc(m1_spring6mo))-min(AICc(m_null),AICc(m1_spring6mo),AICc(m2_spring6mo),AICc(m3_spring6mo), AICc(m1_monsoon6mo), AICc(m2_monsoon6mo),AICc(m3_monsoon6mo),AICc(m1_spring6mo_lag),AICc(m2_spring6mo_lag),AICc(m3_spring6mo_lag),AICc(m1_monsoon6mo_lag),AICc(m2_monsoon6mo_lag),AICc(m3_monsoon6mo_lag),AICc(m1_spring6mo_AR1),AICc(m1_spring6mo_AR2),AICc(m2_spring6mo_AR1),AICc(m2_spring6mo_AR2),AICc(m3_spring6mo_AR1),AICc(m3_spring6mo_AR2),AICc(m1_monsoon6mo_AR1),AICc(m1_monsoon6mo_AR2),AICc(m2_monsoon6mo_AR1),AICc(m2_monsoon6mo_AR2),AICc(m3_monsoon6mo_AR1),AICc(m3_monsoon6mo_AR2),AICc(m1_spring6mo_lag_AR1),AICc(m1_spring6mo_lag_AR2),AICc(m2_spring6mo_lag_AR1),AICc(m2_spring6mo_lag_AR2),AICc(m3_spring6mo_lag_AR1),AICc(m3_spring6mo_lag_AR2),AICc(m1_monsoon6mo_lag_AR1),AICc(m1_monsoon6mo_lag_AR2),AICc(m2_monsoon6mo_lag_AR1),AICc(m2_monsoon6mo_lag_AR2),AICc(m3_monsoon6mo_lag_AR1),AICc(m3_monsoon6mo_lag_AR2))
  
  dAICc_m2_spring6mo<-(AICc(m2_spring6mo))-min(AICc(m_null),AICc(m1_spring6mo),AICc(m2_spring6mo),AICc(m3_spring6mo), AICc(m1_monsoon6mo), AICc(m2_monsoon6mo),AICc(m3_monsoon6mo),AICc(m1_spring6mo_lag),AICc(m2_spring6mo_lag),AICc(m3_spring6mo_lag),AICc(m1_monsoon6mo_lag),AICc(m2_monsoon6mo_lag),AICc(m3_monsoon6mo_lag),AICc(m1_spring6mo_AR1),AICc(m1_spring6mo_AR2),AICc(m2_spring6mo_AR1),AICc(m2_spring6mo_AR2),AICc(m3_spring6mo_AR1),AICc(m3_spring6mo_AR2),AICc(m1_monsoon6mo_AR1),AICc(m1_monsoon6mo_AR2),AICc(m2_monsoon6mo_AR1),AICc(m2_monsoon6mo_AR2),AICc(m3_monsoon6mo_AR1),AICc(m3_monsoon6mo_AR2),AICc(m1_spring6mo_lag_AR1),AICc(m1_spring6mo_lag_AR2),AICc(m2_spring6mo_lag_AR1),AICc(m2_spring6mo_lag_AR2),AICc(m3_spring6mo_lag_AR1),AICc(m3_spring6mo_lag_AR2),AICc(m1_monsoon6mo_lag_AR1),AICc(m1_monsoon6mo_lag_AR2),AICc(m2_monsoon6mo_lag_AR1),AICc(m2_monsoon6mo_lag_AR2),AICc(m3_monsoon6mo_lag_AR1),AICc(m3_monsoon6mo_lag_AR2))
  
  dAICc_m3_spring6mo<-(AICc(m3_spring6mo))-min(AICc(m_null),AICc(m1_spring6mo),AICc(m2_spring6mo),AICc(m3_spring6mo), AICc(m1_monsoon6mo), AICc(m2_monsoon6mo),AICc(m3_monsoon6mo),AICc(m1_spring6mo_lag),AICc(m2_spring6mo_lag),AICc(m3_spring6mo_lag),AICc(m1_monsoon6mo_lag),AICc(m2_monsoon6mo_lag),AICc(m3_monsoon6mo_lag),AICc(m1_spring6mo_AR1),AICc(m1_spring6mo_AR2),AICc(m2_spring6mo_AR1),AICc(m2_spring6mo_AR2),AICc(m3_spring6mo_AR1),AICc(m3_spring6mo_AR2),AICc(m1_monsoon6mo_AR1),AICc(m1_monsoon6mo_AR2),AICc(m2_monsoon6mo_AR1),AICc(m2_monsoon6mo_AR2),AICc(m3_monsoon6mo_AR1),AICc(m3_monsoon6mo_AR2),AICc(m1_spring6mo_lag_AR1),AICc(m1_spring6mo_lag_AR2),AICc(m2_spring6mo_lag_AR1),AICc(m2_spring6mo_lag_AR2),AICc(m3_spring6mo_lag_AR1),AICc(m3_spring6mo_lag_AR2),AICc(m1_monsoon6mo_lag_AR1),AICc(m1_monsoon6mo_lag_AR2),AICc(m2_monsoon6mo_lag_AR1),AICc(m2_monsoon6mo_lag_AR2),AICc(m3_monsoon6mo_lag_AR1),AICc(m3_monsoon6mo_lag_AR2))
  
  dAICc_m1_monsoon6mo<-(AICc(m1_monsoon6mo))-min(AICc(m_null),AICc(m1_spring6mo),AICc(m2_spring6mo),AICc(m3_spring6mo), AICc(m1_monsoon6mo), AICc(m2_monsoon6mo),AICc(m3_monsoon6mo),AICc(m1_spring6mo_lag),AICc(m2_spring6mo_lag),AICc(m3_spring6mo_lag),AICc(m1_monsoon6mo_lag),AICc(m2_monsoon6mo_lag),AICc(m3_monsoon6mo_lag),AICc(m1_spring6mo_AR1),AICc(m1_spring6mo_AR2),AICc(m2_spring6mo_AR1),AICc(m2_spring6mo_AR2),AICc(m3_spring6mo_AR1),AICc(m3_spring6mo_AR2),AICc(m1_monsoon6mo_AR1),AICc(m1_monsoon6mo_AR2),AICc(m2_monsoon6mo_AR1),AICc(m2_monsoon6mo_AR2),AICc(m3_monsoon6mo_AR1),AICc(m3_monsoon6mo_AR2),AICc(m1_spring6mo_lag_AR1),AICc(m1_spring6mo_lag_AR2),AICc(m2_spring6mo_lag_AR1),AICc(m2_spring6mo_lag_AR2),AICc(m3_spring6mo_lag_AR1),AICc(m3_spring6mo_lag_AR2),AICc(m1_monsoon6mo_lag_AR1),AICc(m1_monsoon6mo_lag_AR2),AICc(m2_monsoon6mo_lag_AR1),AICc(m2_monsoon6mo_lag_AR2),AICc(m3_monsoon6mo_lag_AR1),AICc(m3_monsoon6mo_lag_AR2))
  
  dAICc_m2_monsoon6mo<-(AICc(m2_monsoon6mo))-min(AICc(m_null),AICc(m1_spring6mo),AICc(m2_spring6mo),AICc(m3_spring6mo), AICc(m1_monsoon6mo), AICc(m2_monsoon6mo),AICc(m3_monsoon6mo),AICc(m1_spring6mo_lag),AICc(m2_spring6mo_lag),AICc(m3_spring6mo_lag),AICc(m1_monsoon6mo_lag),AICc(m2_monsoon6mo_lag),AICc(m3_monsoon6mo_lag),AICc(m1_spring6mo_AR1),AICc(m1_spring6mo_AR2),AICc(m2_spring6mo_AR1),AICc(m2_spring6mo_AR2),AICc(m3_spring6mo_AR1),AICc(m3_spring6mo_AR2),AICc(m1_monsoon6mo_AR1),AICc(m1_monsoon6mo_AR2),AICc(m2_monsoon6mo_AR1),AICc(m2_monsoon6mo_AR2),AICc(m3_monsoon6mo_AR1),AICc(m3_monsoon6mo_AR2),AICc(m1_spring6mo_lag_AR1),AICc(m1_spring6mo_lag_AR2),AICc(m2_spring6mo_lag_AR1),AICc(m2_spring6mo_lag_AR2),AICc(m3_spring6mo_lag_AR1),AICc(m3_spring6mo_lag_AR2),AICc(m1_monsoon6mo_lag_AR1),AICc(m1_monsoon6mo_lag_AR2),AICc(m2_monsoon6mo_lag_AR1),AICc(m2_monsoon6mo_lag_AR2),AICc(m3_monsoon6mo_lag_AR1),AICc(m3_monsoon6mo_lag_AR2))
  
  dAICc_m3_monsoon6mo<-(AICc(m3_monsoon6mo))-min(AICc(m_null),AICc(m1_spring6mo),AICc(m2_spring6mo),AICc(m3_spring6mo), AICc(m1_monsoon6mo), AICc(m2_monsoon6mo),AICc(m3_monsoon6mo),AICc(m1_spring6mo_lag),AICc(m2_spring6mo_lag),AICc(m3_spring6mo_lag),AICc(m1_monsoon6mo_lag),AICc(m2_monsoon6mo_lag),AICc(m3_monsoon6mo_lag),AICc(m1_spring6mo_AR1),AICc(m1_spring6mo_AR2),AICc(m2_spring6mo_AR1),AICc(m2_spring6mo_AR2),AICc(m3_spring6mo_AR1),AICc(m3_spring6mo_AR2),AICc(m1_monsoon6mo_AR1),AICc(m1_monsoon6mo_AR2),AICc(m2_monsoon6mo_AR1),AICc(m2_monsoon6mo_AR2),AICc(m3_monsoon6mo_AR1),AICc(m3_monsoon6mo_AR2),AICc(m1_spring6mo_lag_AR1),AICc(m1_spring6mo_lag_AR2),AICc(m2_spring6mo_lag_AR1),AICc(m2_spring6mo_lag_AR2),AICc(m3_spring6mo_lag_AR1),AICc(m3_spring6mo_lag_AR2),AICc(m1_monsoon6mo_lag_AR1),AICc(m1_monsoon6mo_lag_AR2),AICc(m2_monsoon6mo_lag_AR1),AICc(m2_monsoon6mo_lag_AR2),AICc(m3_monsoon6mo_lag_AR1),AICc(m3_monsoon6mo_lag_AR2))
  
  dAICc_m1_spring6mo_lag<-(AICc(m1_spring6mo_lag))-min(AICc(m_null),AICc(m1_spring6mo),AICc(m2_spring6mo),AICc(m3_spring6mo), AICc(m1_monsoon6mo), AICc(m2_monsoon6mo),AICc(m3_monsoon6mo),AICc(m1_spring6mo_lag),AICc(m2_spring6mo_lag),AICc(m3_spring6mo_lag),AICc(m1_monsoon6mo_lag),AICc(m2_monsoon6mo_lag),AICc(m3_monsoon6mo_lag),AICc(m1_spring6mo_AR1),AICc(m1_spring6mo_AR2),AICc(m2_spring6mo_AR1),AICc(m2_spring6mo_AR2),AICc(m3_spring6mo_AR1),AICc(m3_spring6mo_AR2),AICc(m1_monsoon6mo_AR1),AICc(m1_monsoon6mo_AR2),AICc(m2_monsoon6mo_AR1),AICc(m2_monsoon6mo_AR2),AICc(m3_monsoon6mo_AR1),AICc(m3_monsoon6mo_AR2),AICc(m1_spring6mo_lag_AR1),AICc(m1_spring6mo_lag_AR2),AICc(m2_spring6mo_lag_AR1),AICc(m2_spring6mo_lag_AR2),AICc(m3_spring6mo_lag_AR1),AICc(m3_spring6mo_lag_AR2),AICc(m1_monsoon6mo_lag_AR1),AICc(m1_monsoon6mo_lag_AR2),AICc(m2_monsoon6mo_lag_AR1),AICc(m2_monsoon6mo_lag_AR2),AICc(m3_monsoon6mo_lag_AR1),AICc(m3_monsoon6mo_lag_AR2))
  
  dAICc_m2_spring6mo_lag<-(AICc(m2_spring6mo_lag))-min(AICc(m_null),AICc(m1_spring6mo),AICc(m2_spring6mo),AICc(m3_spring6mo), AICc(m1_monsoon6mo), AICc(m2_monsoon6mo),AICc(m3_monsoon6mo),AICc(m1_spring6mo_lag),AICc(m2_spring6mo_lag),AICc(m3_spring6mo_lag),AICc(m1_monsoon6mo_lag),AICc(m2_monsoon6mo_lag),AICc(m3_monsoon6mo_lag),AICc(m1_spring6mo_AR1),AICc(m1_spring6mo_AR2),AICc(m2_spring6mo_AR1),AICc(m2_spring6mo_AR2),AICc(m3_spring6mo_AR1),AICc(m3_spring6mo_AR2),AICc(m1_monsoon6mo_AR1),AICc(m1_monsoon6mo_AR2),AICc(m2_monsoon6mo_AR1),AICc(m2_monsoon6mo_AR2),AICc(m3_monsoon6mo_AR1),AICc(m3_monsoon6mo_AR2),AICc(m1_spring6mo_lag_AR1),AICc(m1_spring6mo_lag_AR2),AICc(m2_spring6mo_lag_AR1),AICc(m2_spring6mo_lag_AR2),AICc(m3_spring6mo_lag_AR1),AICc(m3_spring6mo_lag_AR2),AICc(m1_monsoon6mo_lag_AR1),AICc(m1_monsoon6mo_lag_AR2),AICc(m2_monsoon6mo_lag_AR1),AICc(m2_monsoon6mo_lag_AR2),AICc(m3_monsoon6mo_lag_AR1),AICc(m3_monsoon6mo_lag_AR2))
  
  dAICc_m3_spring6mo_lag<-(AICc(m3_spring6mo_lag))-min(AICc(m_null),AICc(m1_spring6mo),AICc(m2_spring6mo),AICc(m3_spring6mo), AICc(m1_monsoon6mo), AICc(m2_monsoon6mo),AICc(m3_monsoon6mo),AICc(m1_spring6mo_lag),AICc(m2_spring6mo_lag),AICc(m3_spring6mo_lag),AICc(m1_monsoon6mo_lag),AICc(m2_monsoon6mo_lag),AICc(m3_monsoon6mo_lag),AICc(m1_spring6mo_AR1),AICc(m1_spring6mo_AR2),AICc(m2_spring6mo_AR1),AICc(m2_spring6mo_AR2),AICc(m3_spring6mo_AR1),AICc(m3_spring6mo_AR2),AICc(m1_monsoon6mo_AR1),AICc(m1_monsoon6mo_AR2),AICc(m2_monsoon6mo_AR1),AICc(m2_monsoon6mo_AR2),AICc(m3_monsoon6mo_AR1),AICc(m3_monsoon6mo_AR2),AICc(m1_spring6mo_lag_AR1),AICc(m1_spring6mo_lag_AR2),AICc(m2_spring6mo_lag_AR1),AICc(m2_spring6mo_lag_AR2),AICc(m3_spring6mo_lag_AR1),AICc(m3_spring6mo_lag_AR2),AICc(m1_monsoon6mo_lag_AR1),AICc(m1_monsoon6mo_lag_AR2),AICc(m2_monsoon6mo_lag_AR1),AICc(m2_monsoon6mo_lag_AR2),AICc(m3_monsoon6mo_lag_AR1),AICc(m3_monsoon6mo_lag_AR2))
  
  dAICc_m1_monsoon6mo_lag<-(AICc(m1_monsoon6mo_lag))-min(AICc(m_null),AICc(m1_spring6mo),AICc(m2_spring6mo),AICc(m3_spring6mo), AICc(m1_monsoon6mo), AICc(m2_monsoon6mo),AICc(m3_monsoon6mo),AICc(m1_spring6mo_lag),AICc(m2_spring6mo_lag),AICc(m3_spring6mo_lag),AICc(m1_monsoon6mo_lag),AICc(m2_monsoon6mo_lag),AICc(m3_monsoon6mo_lag),AICc(m1_spring6mo_AR1),AICc(m1_spring6mo_AR2),AICc(m2_spring6mo_AR1),AICc(m2_spring6mo_AR2),AICc(m3_spring6mo_AR1),AICc(m3_spring6mo_AR2),AICc(m1_monsoon6mo_AR1),AICc(m1_monsoon6mo_AR2),AICc(m2_monsoon6mo_AR1),AICc(m2_monsoon6mo_AR2),AICc(m3_monsoon6mo_AR1),AICc(m3_monsoon6mo_AR2),AICc(m1_spring6mo_lag_AR1),AICc(m1_spring6mo_lag_AR2),AICc(m2_spring6mo_lag_AR1),AICc(m2_spring6mo_lag_AR2),AICc(m3_spring6mo_lag_AR1),AICc(m3_spring6mo_lag_AR2),AICc(m1_monsoon6mo_lag_AR1),AICc(m1_monsoon6mo_lag_AR2),AICc(m2_monsoon6mo_lag_AR1),AICc(m2_monsoon6mo_lag_AR2),AICc(m3_monsoon6mo_lag_AR1),AICc(m3_monsoon6mo_lag_AR2))
  
  dAICc_m2_monsoon6mo_lag<-(AICc(m2_monsoon6mo_lag))-min(AICc(m_null),AICc(m1_spring6mo),AICc(m2_spring6mo),AICc(m3_spring6mo), AICc(m1_monsoon6mo), AICc(m2_monsoon6mo),AICc(m3_monsoon6mo),AICc(m1_spring6mo_lag),AICc(m2_spring6mo_lag),AICc(m3_spring6mo_lag),AICc(m1_monsoon6mo_lag),AICc(m2_monsoon6mo_lag),AICc(m3_monsoon6mo_lag),AICc(m1_spring6mo_AR1),AICc(m1_spring6mo_AR2),AICc(m2_spring6mo_AR1),AICc(m2_spring6mo_AR2),AICc(m3_spring6mo_AR1),AICc(m3_spring6mo_AR2),AICc(m1_monsoon6mo_AR1),AICc(m1_monsoon6mo_AR2),AICc(m2_monsoon6mo_AR1),AICc(m2_monsoon6mo_AR2),AICc(m3_monsoon6mo_AR1),AICc(m3_monsoon6mo_AR2),AICc(m1_spring6mo_lag_AR1),AICc(m1_spring6mo_lag_AR2),AICc(m2_spring6mo_lag_AR1),AICc(m2_spring6mo_lag_AR2),AICc(m3_spring6mo_lag_AR1),AICc(m3_spring6mo_lag_AR2),AICc(m1_monsoon6mo_lag_AR1),AICc(m1_monsoon6mo_lag_AR2),AICc(m2_monsoon6mo_lag_AR1),AICc(m2_monsoon6mo_lag_AR2),AICc(m3_monsoon6mo_lag_AR1),AICc(m3_monsoon6mo_lag_AR2))
  
  dAICc_m3_monsoon6mo_lag<-(AICc(m3_monsoon6mo_lag))-min(AICc(m_null),AICc(m1_spring6mo),AICc(m2_spring6mo),AICc(m3_spring6mo), AICc(m1_monsoon6mo), AICc(m2_monsoon6mo),AICc(m3_monsoon6mo),AICc(m1_spring6mo_lag),AICc(m2_spring6mo_lag),AICc(m3_spring6mo_lag),AICc(m1_monsoon6mo_lag),AICc(m2_monsoon6mo_lag),AICc(m3_monsoon6mo_lag),AICc(m1_spring6mo_AR1),AICc(m1_spring6mo_AR2),AICc(m2_spring6mo_AR1),AICc(m2_spring6mo_AR2),AICc(m3_spring6mo_AR1),AICc(m3_spring6mo_AR2),AICc(m1_monsoon6mo_AR1),AICc(m1_monsoon6mo_AR2),AICc(m2_monsoon6mo_AR1),AICc(m2_monsoon6mo_AR2),AICc(m3_monsoon6mo_AR1),AICc(m3_monsoon6mo_AR2),AICc(m1_spring6mo_lag_AR1),AICc(m1_spring6mo_lag_AR2),AICc(m2_spring6mo_lag_AR1),AICc(m2_spring6mo_lag_AR2),AICc(m3_spring6mo_lag_AR1),AICc(m3_spring6mo_lag_AR2),AICc(m1_monsoon6mo_lag_AR1),AICc(m1_monsoon6mo_lag_AR2),AICc(m2_monsoon6mo_lag_AR1),AICc(m2_monsoon6mo_lag_AR2),AICc(m3_monsoon6mo_lag_AR1),AICc(m3_monsoon6mo_lag_AR2))
  
  # graph the best CSF for the species (maximum abundance as a function of SPEI), if the best model is not the null, and save the graph
  
  if (dAICc_m1_spring6mo == 0 | dAICc_m1_spring6mo_AR1 == 0 | dAICc_m1_spring6mo_AR2 == 0) {
    
    speciesData <- blue_original[speciesCode]
    speciesData <-cbind(speciesData,blue_climate)
    names(speciesData)[1] <- "speciesCode"
    
    p <- ggplot(speciesData,aes(x=spring6SPEI,y=speciesCode^2))+
      geom_point(size=3,color="#3182bd")+
      geom_smooth(formula=y~poly(x,1),method="lm",color="black")+
      xlab("Spring SPEI\n")+
      ylab("Maximum abundance") +       ggtitle(speciesCode) +
      ggtitle(speciesCode) +
      theme(legend.position = "none") +
      theme_classic() +
      theme(axis.text.x = element_text(color="black",size=14)) +
      theme(axis.text.y = element_text(color="black",size=14)) +
      theme(axis.title=element_text(size=18))+
      theme(plot.title = element_text(size=20))
    print(p)
    
    plot_name <- paste("B",speciesCodes[i], "plot",number[i],".pdf", sep="_")
    
    ggsave(paste(plot_name), p,
           width=4,height=4,units = c("in"),
           dpi = 600)
    
    speciesData$speciesCode <-NULL
  }
  
  if (dAICc_m2_spring6mo == 0 | dAICc_m2_spring6mo_AR1 == 0 | dAICc_m2_spring6mo_AR2 == 0) {
    speciesData <- blue_original[speciesCode]
    speciesData <-cbind(speciesData,blue_climate)
    names(speciesData)[1] <- "speciesCode"
    
    p <- ggplot(speciesData,aes(x=spring6SPEI,y=speciesCode^2))+
      geom_point(size=3,color="#3182bd")+
      geom_smooth(formula=y~poly(x,2),method="lm",color="black")+
      xlab("Spring SPEI\n")+
      ylab("Maximum abundance") +       ggtitle(speciesCode) +
      theme(legend.position = "none") +
      theme_classic() +
      theme(axis.text.x = element_text(color="black",size=14)) +
      theme(axis.text.y = element_text(color="black",size=14)) +
      theme(axis.title=element_text(size=18))+
      theme(plot.title = element_text(size=20))
    print(p)
    
    plot_name <- paste("B",speciesCodes[i], "plot",number[i],".pdf", sep="_")
    
    ggsave(paste(plot_name), p,
           width=4,height=4,units = c("in"),
           dpi = 600)
    
    speciesData$speciesCode <-NULL
  }
  
  if (dAICc_m3_spring6mo == 0 | dAICc_m3_spring6mo_AR1 == 0 | dAICc_m3_spring6mo_AR2 == 0) {
    speciesData <- blue_original[speciesCode]
    speciesData <-cbind(speciesData,blue_climate)
    names(speciesData)[1] <- "speciesCode"
    
    p <- ggplot(speciesData,aes(x=spring6SPEI,y=speciesCode^2))+
      geom_point(size=3,color="#3182bd")+
      geom_smooth(formula=y~poly(x,3),method="lm",color="black")+
      xlab("Spring SPEI\n")+
      ylab("Maximum abundance") +       ggtitle(speciesCode) +
      theme(legend.position = "none") +
      theme_classic() +
      theme(axis.text.x = element_text(color="black",size=14)) +
      theme(axis.text.y = element_text(color="black",size=14)) +
      theme(axis.title=element_text(size=18))+
      theme(plot.title = element_text(size=20))
    print(p)
    
    plot_name <- paste("B",speciesCodes[i], "plot",number[i],".pdf", sep="_")
    
    ggsave(paste(plot_name), p,
           width=4,height=4,units = c("in"),
           dpi = 600)    
    speciesData$speciesCode <-NULL
  }
  
  if (dAICc_m1_monsoon6mo == 0 | dAICc_m1_monsoon6mo_AR1 == 0 | dAICc_m1_monsoon6mo_AR2 == 0) {
    speciesData <- blue_original[speciesCode]
    speciesData <-cbind(speciesData,blue_climate)
    names(speciesData)[1] <- "speciesCode"
    
    p <- ggplot(speciesData,aes(x=monsoon6SPEI,y=speciesCode^2))+
      geom_point(size=3,color="#3182bd")+
      geom_smooth(formula=y~poly(x,1),method="lm",color="black")+
      xlab("Monsoon SPEI\n")+
      ylab("Maximum abundance") +       ggtitle(speciesCode) +
      theme(legend.position = "none") +
      theme_classic() +
      theme(axis.text.x = element_text(color="black",size=14)) +
      theme(axis.text.y = element_text(color="black",size=14)) +
      theme(axis.title=element_text(size=18))+
      theme(plot.title = element_text(size=20))
    print(p)
    
    plot_name <- paste("B",speciesCodes[i], "plot",number[i],".pdf", sep="_")
    
    ggsave(paste(plot_name), p,
           width=4,height=4,units = c("in"),
           dpi = 600)
    
    speciesData$speciesCode <-NULL
  }
  
  if (dAICc_m2_monsoon6mo == 0 | dAICc_m2_monsoon6mo_AR1 == 0 | dAICc_m2_monsoon6mo_AR2 == 0) {
    speciesData <- blue_original[speciesCode]
    speciesData <-cbind(speciesData,blue_climate)
    names(speciesData)[1] <- "speciesCode"
    
    p <- ggplot(speciesData,aes(x=monsoon6SPEI,y=speciesCode^2))+
      geom_point(size=3,color="#3182bd")+
      geom_smooth(formula=y~poly(x,2),method="lm",color="black")+
      xlab("Monsoon SPEI\n")+
      ylab("Maximum abundance") +       ggtitle(speciesCode) +
      theme(legend.position = "none") +
      theme_classic() +
      theme(axis.text.x = element_text(color="black",size=14)) +
      theme(axis.text.y = element_text(color="black",size=14)) +
      theme(axis.title=element_text(size=18))+
      theme(plot.title = element_text(size=20))
    print(p)
    
    plot_name <- paste("B",speciesCodes[i], "plot",number[i],".pdf", sep="_")
    
    ggsave(paste(plot_name), p,
           width=4,height=4,units = c("in"),
           dpi = 600)    
    speciesData$speciesCode <-NULL
  }
  
  if (dAICc_m3_monsoon6mo == 0 | dAICc_m3_monsoon6mo_AR1 == 0 | dAICc_m3_monsoon6mo_AR2 == 0) {
    speciesData <- blue_original[speciesCode]
    speciesData <-cbind(speciesData,blue_climate)
    names(speciesData)[1] <- "speciesCode"
    
    p <- ggplot(speciesData,aes(x=monsoon6SPEI,y=speciesCode^2))+
      geom_point(size=3,color="#3182bd")+
      geom_smooth(formula=y~poly(x,3),method="lm",color="black")+
      xlab("Monsoon SPEI\n")+
      ylab("Maximum abundance") +       ggtitle(speciesCode) +
      theme(legend.position = "none") +
      theme_classic() +
      theme(axis.text.x = element_text(color="black",size=14)) +
      theme(axis.text.y = element_text(color="black",size=14)) +
      theme(axis.title=element_text(size=18))+
      theme(plot.title = element_text(size=20))
    print(p)
    
    plot_name <- paste("B",speciesCodes[i], "plot",number[i],".pdf", sep="_")
    
    ggsave(paste(plot_name), p,
           width=4,height=4,units = c("in"),
           dpi = 600)
    
    speciesData$speciesCode <-NULL
  }
  
  if (dAICc_m1_spring6mo_lag == 0 | dAICc_m1_spring6mo_lag_AR1 == 0 | dAICc_m1_spring6mo_lag_AR2 == 0) {
    speciesData <- blue_original[speciesCode]
    speciesData <-cbind(speciesData,blue_climate)
    names(speciesData)[1] <- "speciesCode"
    
    p <- ggplot(speciesData,aes(x=spring6SPEI_prioryear,y=speciesCode^2))+
      geom_point(size=3,color="#3182bd")+
      geom_smooth(formula=y~poly(x,1),method="lm",color="black")+
      xlab("Spring SPEI \n(previous year)")+
      ylab("Maximum abundance") +       ggtitle(speciesCode) +
      theme(legend.position = "none") +
      theme_classic() +
      theme(axis.text.x = element_text(color="black",size=14)) +
      theme(axis.text.y = element_text(color="black",size=14)) +
      theme(axis.title=element_text(size=18))+
      theme(plot.title = element_text(size=20))
    print(p)
    
    plot_name <- paste("B",speciesCodes[i], "plot",number[i],".pdf", sep="_")
    
    ggsave(paste(plot_name), p,
           width=4,height=4,units = c("in"),
           dpi = 600)
    
    speciesData$speciesCode <-NULL
  }
  
  if (dAICc_m2_spring6mo_lag == 0 | dAICc_m2_spring6mo_lag_AR1 == 0 | dAICc_m2_spring6mo_lag_AR2 == 0) {
    speciesData <- blue_original[speciesCode]
    speciesData <-cbind(speciesData,blue_climate)
    names(speciesData)[1] <- "speciesCode"
    
    p <- ggplot(speciesData,aes(x=spring6SPEI_prioryear,y=speciesCode^2))+
      geom_point(size=3,color="#3182bd")+
      geom_smooth(formula=y~poly(x,2),method="lm",color="black")+
      xlab("Spring SPEI \n(previous year)")+
      ylab("Maximum abundance") +       ggtitle(speciesCode) +
      theme(legend.position = "none") +
      theme_classic() +
      theme(axis.text.x = element_text(color="black",size=14)) +
      theme(axis.text.y = element_text(color="black",size=14)) +
      theme(axis.title=element_text(size=18))+
      theme(plot.title = element_text(size=20))
    print(p)
    
    plot_name <- paste("B",speciesCodes[i], "plot",number[i],".pdf", sep="_")
    
    ggsave(paste(plot_name), p,
           width=4,height=4,units = c("in"),
           dpi = 600)
    
    speciesData$speciesCode <-NULL
  }
  
  if (dAICc_m3_spring6mo_lag == 0 | dAICc_m3_spring6mo_lag_AR1 == 0 | dAICc_m3_spring6mo_lag_AR2 == 0) {
    speciesData <- blue_original[speciesCode]
    speciesData <-cbind(speciesData,blue_climate)
    names(speciesData)[1] <- "speciesCode"
    
    p <- ggplot(speciesData,aes(x=spring6SPEI_prioryear,y=speciesCode^2))+
      geom_point(size=3,color="#3182bd")+
      geom_smooth(formula=y~poly(x,3),method="lm",color="black")+
      xlab("Spring SPEI \n(previous year)")+
      ylab("Maximum abundance") +       ggtitle(speciesCode) +
      theme(legend.position = "none") +
      theme_classic() +
      theme(axis.text.x = element_text(color="black",size=14)) +
      theme(axis.text.y = element_text(color="black",size=14)) +
      theme(axis.title=element_text(size=18))+
      theme(plot.title = element_text(size=20))
    print(p)
    
    plot_name <- paste("B",speciesCodes[i], "plot",number[i],".pdf", sep="_")
    
    ggsave(paste(plot_name), p,
           width=4,height=4,units = c("in"),
           dpi = 600)
    
    speciesData$speciesCode <-NULL
  }
  
  if (dAICc_m1_monsoon6mo_lag == 0 | dAICc_m1_monsoon6mo_lag_AR1 == 0 | dAICc_m1_monsoon6mo_lag_AR2 == 0) {
    speciesData <- blue_original[speciesCode]
    speciesData <-cbind(speciesData,blue_climate)
    names(speciesData)[1] <- "speciesCode"
    
    p <- ggplot(speciesData,aes(x=monsoon6SPEI_prioryear,y=speciesCode^2))+
      geom_point(size=3,color="#3182bd")+
      geom_smooth(formula=y~poly(x,1),method="lm",color="black")+
      xlab("Monsoon SPEI \n(previous year)")+
      ylab("Maximum abundance") +       ggtitle(speciesCode) +
      theme(legend.position = "none") +
      theme_classic() +
      theme(axis.text.x = element_text(color="black",size=14)) +
      theme(axis.text.y = element_text(color="black",size=14)) +
      theme(axis.title=element_text(size=18))+
      theme(plot.title = element_text(size=20))
    print(p)
    
    plot_name <- paste("B",speciesCodes[i], "plot",number[i],".pdf", sep="_")
    
    ggsave(paste(plot_name), p,
           width=4,height=4,units = c("in"),
           dpi = 600)
    
    speciesData$speciesCode <-NULL
  }
  
  if (dAICc_m2_monsoon6mo_lag == 0 | dAICc_m2_monsoon6mo_lag_AR1 == 0 | dAICc_m2_monsoon6mo_lag_AR2 == 0) {
    speciesData <- blue_original[speciesCode]
    speciesData <-cbind(speciesData,blue_climate)
    names(speciesData)[1] <- "speciesCode"
    
    p <- ggplot(speciesData,aes(x=monsoon6SPEI_prioryear,y=speciesCode^2))+
      geom_point(size=3,color="#3182bd")+
      geom_smooth(formula=y~poly(x,2),method="lm",color="black")+
      xlab("Monsoon SPEI \n(previous year)")+
      ylab("Maximum abundance") +       ggtitle(speciesCode) +
      theme(legend.position = "none") +
      theme_classic() +
      theme(axis.text.x = element_text(color="black",size=14)) +
      theme(axis.text.y = element_text(color="black",size=14)) +
      theme(axis.title=element_text(size=18))+
      theme(plot.title = element_text(size=20))
    print(p)
    
    plot_name <- paste("B",speciesCodes[i], "plot",number[i],".pdf", sep="_")
    
    ggsave(paste(plot_name), p,
           width=4,height=4,units = c("in"),
           dpi = 600)
    
    speciesData$speciesCode <-NULL
  }
  
  if (dAICc_m3_monsoon6mo_lag == 0 | dAICc_m3_monsoon6mo_lag_AR1 == 0 | dAICc_m3_monsoon6mo_lag_AR2 == 0) {
    speciesData <- blue_original[speciesCode]
    speciesData <-cbind(speciesData,blue_climate)
    names(speciesData)[1] <- "speciesCode"
    
    p<-ggplot(speciesData,aes(x=monsoon6SPEI_prioryear,y=speciesCode^2))+
      geom_point(size=3,color="#3182bd")+
      geom_smooth(formula=y~poly(x,3),method="lm",color="black")+
      xlab("Monsoon SPEI \n(previous year)")+
      ylab("Maximum abundance") +       ggtitle(speciesCode) +
      theme(legend.position = "none") +
      theme_classic() +
      theme(axis.text.x = element_text(color="black",size=14)) +
      theme(axis.text.y = element_text(color="black",size=14)) +
      theme(axis.title=element_text(size=18))+
      theme(plot.title = element_text(size=20))
    print(p)
    
    plot_name <- paste("B",speciesCodes[i], "plot",number[i],".pdf", sep="_")
    
    ggsave(paste(plot_name), p,
           width=4,height=4,units = c("in"),
           dpi = 600)
    
    speciesData$speciesCode <-NULL
  }
  
  
  # create objects containing statistical values related to each model
  
  numDF_lin_m1_spring6mo<-anova(m1_spring6mo,type="marginal")[2,1]
  denDF_lin_m1_spring6mo<-anova(m1_spring6mo,type="marginal")[2,2]
  F_lin_m1_spring6mo<-anova(m1_spring6mo,type="marginal")[2,3]
  P_lin_m1_spring6mo<-anova(m1_spring6mo,type="marginal")[2,4]
  ParamEst_lin_m1_spring6mo<-coef(summary(m1_spring6mo))[2,1]
  SE_lin_m1_spring6mo<-coef(summary(m1_spring6mo))[2,2]
  
  numDF_lin_m2_spring6mo<-anova(m2_spring6mo,type="marginal")[2,1]
  denDF_lin_m2_spring6mo<-anova(m2_spring6mo,type="marginal")[2,2]
  F_lin_m2_spring6mo<-anova(m2_spring6mo,type="marginal")[2,3]
  P_lin_m2_spring6mo<-anova(m2_spring6mo,type="marginal")[2,4]
  ParamEst_lin_m2_spring6mo<-coef(summary(m2_spring6mo))[2,1]
  SE_lin_m2_spring6mo<-coef(summary(m2_spring6mo))[2,2]
  
  numDF_lin_m3_spring6mo<-anova(m3_spring6mo,type="marginal")[2,1]
  denDF_lin_m3_spring6mo<-anova(m3_spring6mo,type="marginal")[2,2]
  F_lin_m3_spring6mo<-anova(m3_spring6mo,type="marginal")[2,3]
  P_lin_m3_spring6mo<-anova(m3_spring6mo,type="marginal")[2,4]
  ParamEst_lin_m3_spring6mo<-coef(summary(m3_spring6mo))[2,1]
  SE_lin_m3_spring6mo<-coef(summary(m3_spring6mo))[2,2]
  
  numDF_lin_m1_monsoon6mo<-anova(m1_monsoon6mo,type="marginal")[2,1]
  denDF_lin_m1_monsoon6mo<-anova(m1_monsoon6mo,type="marginal")[2,2]
  F_lin_m1_monsoon6mo<-anova(m1_monsoon6mo,type="marginal")[2,3]
  P_lin_m1_monsoon6mo<-anova(m1_monsoon6mo,type="marginal")[2,4]
  ParamEst_lin_m1_monsoon6mo<-coef(summary(m1_monsoon6mo))[2,1]
  SE_lin_m1_monsoon6mo<-coef(summary(m1_monsoon6mo))[2,2]
  
  numDF_lin_m2_monsoon6mo<-anova(m2_monsoon6mo,type="marginal")[2,1]
  denDF_lin_m2_monsoon6mo<-anova(m2_monsoon6mo,type="marginal")[2,2]
  F_lin_m2_monsoon6mo<-anova(m2_monsoon6mo,type="marginal")[2,3]
  P_lin_m2_monsoon6mo<-anova(m2_monsoon6mo,type="marginal")[2,4]
  ParamEst_lin_m2_monsoon6mo<-coef(summary(m2_monsoon6mo))[2,1]
  SE_lin_m2_monsoon6mo<-coef(summary(m2_monsoon6mo))[2,2]
  
  numDF_lin_m3_monsoon6mo<-anova(m3_monsoon6mo,type="marginal")[2,1]
  denDF_lin_m3_monsoon6mo<-anova(m3_monsoon6mo,type="marginal")[2,2]
  F_lin_m3_monsoon6mo<-anova(m3_monsoon6mo,type="marginal")[2,3]
  P_lin_m3_monsoon6mo<-anova(m3_monsoon6mo,type="marginal")[2,4]
  ParamEst_lin_m3_monsoon6mo<-coef(summary(m3_monsoon6mo))[2,1]
  SE_lin_m3_monsoon6mo<-coef(summary(m3_monsoon6mo))[2,2]
  
  numDF_lin_m1_spring6mo_lag<-anova(m1_spring6mo_lag,type="marginal")[2,1]
  denDF_lin_m1_spring6mo_lag<-anova(m1_spring6mo_lag,type="marginal")[2,2]
  F_lin_m1_spring6mo_lag<-anova(m1_spring6mo_lag,type="marginal")[2,3]
  P_lin_m1_spring6mo_lag<-anova(m1_spring6mo_lag,type="marginal")[2,4]
  ParamEst_lin_m1_spring6mo_lag<-coef(summary(m1_spring6mo_lag))[2,1]
  SE_lin_m1_spring6mo_lag<-coef(summary(m1_spring6mo_lag))[2,2]
  
  numDF_lin_m2_spring6mo_lag<-anova(m2_spring6mo_lag,type="marginal")[2,1]
  denDF_lin_m2_spring6mo_lag<-anova(m2_spring6mo_lag,type="marginal")[2,2]
  F_lin_m2_spring6mo_lag<-anova(m2_spring6mo_lag,type="marginal")[2,3]
  P_lin_m2_spring6mo_lag<-anova(m2_spring6mo_lag,type="marginal")[2,4]
  ParamEst_lin_m2_spring6mo_lag<-coef(summary(m2_spring6mo_lag))[2,1]
  SE_lin_m2_spring6mo_lag<-coef(summary(m2_spring6mo_lag))[2,2]
  
  numDF_lin_m3_spring6mo_lag<-anova(m3_spring6mo_lag,type="marginal")[2,1]
  denDF_lin_m3_spring6mo_lag<-anova(m3_spring6mo_lag,type="marginal")[2,2]
  F_lin_m3_spring6mo_lag<-anova(m3_spring6mo_lag,type="marginal")[2,3]
  P_lin_m3_spring6mo_lag<-anova(m3_spring6mo_lag,type="marginal")[2,4]
  ParamEst_lin_m3_spring6mo_lag<-coef(summary(m3_spring6mo_lag))[2,1]
  SE_lin_m3_spring6mo_lag<-coef(summary(m3_spring6mo_lag))[2,2]
  
  numDF_lin_m1_monsoon6mo_lag<-anova(m1_monsoon6mo_lag,type="marginal")[2,1]
  denDF_lin_m1_monsoon6mo_lag<-anova(m1_monsoon6mo_lag,type="marginal")[2,2]
  F_lin_m1_monsoon6mo_lag<-anova(m1_monsoon6mo_lag,type="marginal")[2,3]
  P_lin_m1_monsoon6mo_lag<-anova(m1_monsoon6mo_lag,type="marginal")[2,4]
  ParamEst_lin_m1_monsoon6mo_lag<-coef(summary(m1_monsoon6mo_lag))[2,1]
  SE_lin_m1_monsoon6mo_lag<-coef(summary(m1_monsoon6mo_lag))[2,2]
  
  numDF_lin_m2_monsoon6mo_lag<-anova(m2_monsoon6mo_lag,type="marginal")[2,1]
  denDF_lin_m2_monsoon6mo_lag<-anova(m2_monsoon6mo_lag,type="marginal")[2,2]
  F_lin_m2_monsoon6mo_lag<-anova(m2_monsoon6mo_lag,type="marginal")[2,3]
  P_lin_m2_monsoon6mo_lag<-anova(m2_monsoon6mo_lag,type="marginal")[2,4]
  ParamEst_lin_m2_monsoon6mo_lag<-coef(summary(m2_monsoon6mo_lag))[2,1]
  SE_lin_m2_monsoon6mo_lag<-coef(summary(m2_monsoon6mo_lag))[2,2]
  
  numDF_lin_m3_monsoon6mo_lag<-anova(m3_monsoon6mo_lag,type="marginal")[2,1]
  denDF_lin_m3_monsoon6mo_lag<-anova(m3_monsoon6mo_lag,type="marginal")[2,2]
  F_lin_m3_monsoon6mo_lag<-anova(m3_monsoon6mo_lag,type="marginal")[2,3]
  P_lin_m3_monsoon6mo_lag<-anova(m3_monsoon6mo_lag,type="marginal")[2,4]
  ParamEst_lin_m3_monsoon6mo_lag<-coef(summary(m3_monsoon6mo_lag))[2,1]
  SE_lin_m3_monsoon6mo_lag<-coef(summary(m3_monsoon6mo_lag))[2,2]
  
  numDF_lin_m1_spring6mo_AR1<-anova(m1_spring6mo_AR1,type="marginal")[2,1]
  denDF_lin_m1_spring6mo_AR1<-anova(m1_spring6mo_AR1,type="marginal")[2,2]
  F_lin_m1_spring6mo_AR1<-anova(m1_spring6mo_AR1,type="marginal")[2,3]
  P_lin_m1_spring6mo_AR1<-anova(m1_spring6mo_AR1,type="marginal")[2,4]
  ParamEst_lin_m1_spring6mo_AR1<-coef(summary(m1_spring6mo_AR1))[2,1]
  SE_lin_m1_spring6mo_AR1<-coef(summary(m1_spring6mo_AR1))[2,2]
  
  numDF_lin_m1_spring6mo_AR2<-anova(m1_spring6mo_AR2,type="marginal")[2,1]
  denDF_lin_m1_spring6mo_AR2<-anova(m1_spring6mo_AR2,type="marginal")[2,2]
  F_lin_m1_spring6mo_AR2<-anova(m1_spring6mo_AR2,type="marginal")[2,3]
  P_lin_m1_spring6mo_AR2<-anova(m1_spring6mo_AR2,type="marginal")[2,4]
  ParamEst_lin_m1_spring6mo_AR2<-coef(summary(m1_spring6mo_AR2))[2,1]
  SE_lin_m1_spring6mo_AR2<-coef(summary(m1_spring6mo_AR2))[2,2]
  
  numDF_lin_m2_spring6mo_AR1<-anova(m2_spring6mo_AR1,type="marginal")[2,1]
  denDF_lin_m2_spring6mo_AR1<-anova(m2_spring6mo_AR1,type="marginal")[2,2]
  F_lin_m2_spring6mo_AR1<-anova(m2_spring6mo_AR1,type="marginal")[2,3]
  P_lin_m2_spring6mo_AR1<-anova(m2_spring6mo_AR1,type="marginal")[2,4]
  ParamEst_lin_m2_spring6mo_AR1<-coef(summary(m2_spring6mo_AR1))[2,1]
  SE_lin_m2_spring6mo_AR1<-coef(summary(m2_spring6mo_AR1))[2,2]
  
  numDF_lin_m2_spring6mo_AR2<-anova(m2_spring6mo_AR2,type="marginal")[2,1]
  denDF_lin_m2_spring6mo_AR2<-anova(m2_spring6mo_AR2,type="marginal")[2,2]
  F_lin_m2_spring6mo_AR2<-anova(m2_spring6mo_AR2,type="marginal")[2,3]
  P_lin_m2_spring6mo_AR2<-anova(m2_spring6mo_AR2,type="marginal")[2,4]
  ParamEst_lin_m2_spring6mo_AR2<-coef(summary(m2_spring6mo_AR2))[2,1]
  SE_lin_m2_spring6mo_AR2<-coef(summary(m2_spring6mo_AR2))[2,2]
  
  numDF_lin_m3_spring6mo_AR1<-anova(m3_spring6mo_AR1,type="marginal")[2,1]
  denDF_lin_m3_spring6mo_AR1<-anova(m3_spring6mo_AR1,type="marginal")[2,2]
  F_lin_m3_spring6mo_AR1<-anova(m3_spring6mo_AR1,type="marginal")[2,3]
  P_lin_m3_spring6mo_AR1<-anova(m3_spring6mo_AR1,type="marginal")[2,4]
  ParamEst_lin_m3_spring6mo_AR1<-coef(summary(m3_spring6mo_AR1))[2,1]
  SE_lin_m3_spring6mo_AR1<-coef(summary(m3_spring6mo_AR1))[2,2]
  
  numDF_lin_m3_spring6mo_AR2<-anova(m3_spring6mo_AR2,type="marginal")[2,1]
  denDF_lin_m3_spring6mo_AR2<-anova(m3_spring6mo_AR2,type="marginal")[2,2]
  F_lin_m3_spring6mo_AR2<-anova(m3_spring6mo_AR2,type="marginal")[2,3]
  P_lin_m3_spring6mo_AR2<-anova(m3_spring6mo_AR2,type="marginal")[2,4]
  ParamEst_lin_m3_spring6mo_AR2<-coef(summary(m3_spring6mo_AR2))[2,1]
  SE_lin_m3_spring6mo_AR2<-coef(summary(m3_spring6mo_AR2))[2,2]
  
  numDF_lin_m1_monsoon6mo_AR1<-anova(m1_monsoon6mo_AR1,type="marginal")[2,1]
  denDF_lin_m1_monsoon6mo_AR1<-anova(m1_monsoon6mo_AR1,type="marginal")[2,2]
  F_lin_m1_monsoon6mo_AR1<-anova(m1_monsoon6mo_AR1,type="marginal")[2,3]
  P_lin_m1_monsoon6mo_AR1<-anova(m1_monsoon6mo_AR1,type="marginal")[2,4]
  ParamEst_lin_m1_monsoon6mo_AR1<-coef(summary(m1_monsoon6mo_AR1))[2,1]
  SE_lin_m1_monsoon6mo_AR1<-coef(summary(m1_monsoon6mo_AR1))[2,2]
  
  numDF_lin_m1_monsoon6mo_AR2<-anova(m1_monsoon6mo_AR2,type="marginal")[2,1]
  denDF_lin_m1_monsoon6mo_AR2<-anova(m1_monsoon6mo_AR2,type="marginal")[2,2]
  F_lin_m1_monsoon6mo_AR2<-anova(m1_monsoon6mo_AR2,type="marginal")[2,3]
  P_lin_m1_monsoon6mo_AR2<-anova(m1_monsoon6mo_AR2,type="marginal")[2,4]
  ParamEst_lin_m1_monsoon6mo_AR2<-coef(summary(m1_monsoon6mo_AR2))[2,1]
  SE_lin_m1_monsoon6mo_AR2<-coef(summary(m1_monsoon6mo_AR2))[2,2]
  
  numDF_lin_m2_monsoon6mo_AR1<-anova(m2_monsoon6mo_AR1,type="marginal")[2,1]
  denDF_lin_m2_monsoon6mo_AR1<-anova(m2_monsoon6mo_AR1,type="marginal")[2,2]
  F_lin_m2_monsoon6mo_AR1<-anova(m2_monsoon6mo_AR1,type="marginal")[2,3]
  P_lin_m2_monsoon6mo_AR1<-anova(m2_monsoon6mo_AR1,type="marginal")[2,4]
  ParamEst_lin_m2_monsoon6mo_AR1<-coef(summary(m2_monsoon6mo_AR1))[2,1]
  SE_lin_m2_monsoon6mo_AR1<-coef(summary(m2_monsoon6mo_AR1))[2,2]
  
  numDF_lin_m2_monsoon6mo_AR2<-anova(m2_monsoon6mo_AR2,type="marginal")[2,1]
  denDF_lin_m2_monsoon6mo_AR2<-anova(m2_monsoon6mo_AR2,type="marginal")[2,2]
  F_lin_m2_monsoon6mo_AR2<-anova(m2_monsoon6mo_AR2,type="marginal")[2,3]
  P_lin_m2_monsoon6mo_AR2<-anova(m2_monsoon6mo_AR2,type="marginal")[2,4]
  ParamEst_lin_m2_monsoon6mo_AR2<-coef(summary(m2_monsoon6mo_AR2))[2,1]
  SE_lin_m2_monsoon6mo_AR2<-coef(summary(m2_monsoon6mo_AR2))[2,2]
  
  numDF_lin_m3_monsoon6mo_AR1<-anova(m3_monsoon6mo_AR1,type="marginal")[2,1]
  denDF_lin_m3_monsoon6mo_AR1<-anova(m3_monsoon6mo_AR1,type="marginal")[2,2]
  F_lin_m3_monsoon6mo_AR1<-anova(m3_monsoon6mo_AR1,type="marginal")[2,3]
  P_lin_m3_monsoon6mo_AR1<-anova(m3_monsoon6mo_AR1,type="marginal")[2,4]
  ParamEst_lin_m3_monsoon6mo_AR1<-coef(summary(m3_monsoon6mo_AR1))[2,1]
  SE_lin_m3_monsoon6mo_AR1<-coef(summary(m3_monsoon6mo_AR1))[2,2]
  
  numDF_lin_m3_monsoon6mo_AR2<-anova(m3_monsoon6mo_AR2,type="marginal")[2,1]
  denDF_lin_m3_monsoon6mo_AR2<-anova(m3_monsoon6mo_AR2,type="marginal")[2,2]
  F_lin_m3_monsoon6mo_AR2<-anova(m3_monsoon6mo_AR2,type="marginal")[2,3]
  P_lin_m3_monsoon6mo_AR2<-anova(m3_monsoon6mo_AR2,type="marginal")[2,4]
  ParamEst_lin_m3_monsoon6mo_AR2<-coef(summary(m3_monsoon6mo_AR2))[2,1]
  SE_lin_m3_monsoon6mo_AR2<-coef(summary(m3_monsoon6mo_AR2))[2,2]
  
  numDF_lin_m1_spring6mo_lag_AR1<-anova(m1_spring6mo_lag_AR1,type="marginal")[2,1]
  denDF_lin_m1_spring6mo_lag_AR1<-anova(m1_spring6mo_lag_AR1,type="marginal")[2,2]
  F_lin_m1_spring6mo_lag_AR1<-anova(m1_spring6mo_lag_AR1,type="marginal")[2,3]
  P_lin_m1_spring6mo_lag_AR1<-anova(m1_spring6mo_lag_AR1,type="marginal")[2,4]
  ParamEst_lin_m1_spring6mo_lag_AR1<-coef(summary(m1_spring6mo_lag_AR1))[2,1]
  SE_lin_m1_spring6mo_lag_AR1<-coef(summary(m1_spring6mo_lag_AR1))[2,2]
  
  numDF_lin_m1_spring6mo_lag_AR2<-anova(m1_spring6mo_lag_AR2,type="marginal")[2,1]
  denDF_lin_m1_spring6mo_lag_AR2<-anova(m1_spring6mo_lag_AR2,type="marginal")[2,2]
  F_lin_m1_spring6mo_lag_AR2<-anova(m1_spring6mo_lag_AR2,type="marginal")[2,3]
  P_lin_m1_spring6mo_lag_AR2<-anova(m1_spring6mo_lag_AR2,type="marginal")[2,4]
  ParamEst_lin_m1_spring6mo_lag_AR2<-coef(summary(m1_spring6mo_lag_AR2))[2,1]
  SE_lin_m1_spring6mo_lag_AR2<-coef(summary(m1_spring6mo_lag_AR2))[2,2]
  
  numDF_lin_m2_spring6mo_lag_AR1<-anova(m2_spring6mo_lag_AR1,type="marginal")[2,1]
  denDF_lin_m2_spring6mo_lag_AR1<-anova(m2_spring6mo_lag_AR1,type="marginal")[2,2]
  F_lin_m2_spring6mo_lag_AR1<-anova(m2_spring6mo_lag_AR1,type="marginal")[2,3]
  P_lin_m2_spring6mo_lag_AR1<-anova(m2_spring6mo_lag_AR1,type="marginal")[2,4]
  ParamEst_lin_m2_spring6mo_lag_AR1<-coef(summary(m2_spring6mo_lag_AR1))[2,1]
  SE_lin_m2_spring6mo_lag_AR1<-coef(summary(m2_spring6mo_lag_AR1))[2,2]
  
  numDF_lin_m2_spring6mo_lag_AR2<-anova(m2_spring6mo_lag_AR2,type="marginal")[2,1]
  denDF_lin_m2_spring6mo_lag_AR2<-anova(m2_spring6mo_lag_AR2,type="marginal")[2,2]
  F_lin_m2_spring6mo_lag_AR2<-anova(m2_spring6mo_lag_AR2,type="marginal")[2,3]
  P_lin_m2_spring6mo_lag_AR2<-anova(m2_spring6mo_lag_AR2,type="marginal")[2,4]
  ParamEst_lin_m2_spring6mo_lag_AR2<-coef(summary(m2_spring6mo_lag_AR2))[2,1]
  SE_lin_m2_spring6mo_lag_AR2<-coef(summary(m2_spring6mo_lag_AR2))[2,2]
  
  numDF_lin_m3_spring6mo_lag_AR1<-anova(m3_spring6mo_lag_AR1,type="marginal")[2,1]
  denDF_lin_m3_spring6mo_lag_AR1<-anova(m3_spring6mo_lag_AR1,type="marginal")[2,2]
  F_lin_m3_spring6mo_lag_AR1<-anova(m3_spring6mo_lag_AR1,type="marginal")[2,3]
  P_lin_m3_spring6mo_lag_AR1<-anova(m3_spring6mo_lag_AR1,type="marginal")[2,4]
  ParamEst_lin_m3_spring6mo_lag_AR1<-coef(summary(m3_spring6mo_lag_AR1))[2,1]
  SE_lin_m3_spring6mo_lag_AR1<-coef(summary(m3_spring6mo_lag_AR1))[2,2]
  
  numDF_lin_m3_spring6mo_lag_AR2<-anova(m3_spring6mo_lag_AR2,type="marginal")[2,1]
  denDF_lin_m3_spring6mo_lag_AR2<-anova(m3_spring6mo_lag_AR2,type="marginal")[2,2]
  F_lin_m3_spring6mo_lag_AR2<-anova(m3_spring6mo_lag_AR2,type="marginal")[2,3]
  P_lin_m3_spring6mo_lag_AR2<-anova(m3_spring6mo_lag_AR2,type="marginal")[2,4]
  ParamEst_lin_m3_spring6mo_lag_AR2<-coef(summary(m3_spring6mo_lag_AR2))[2,1]
  SE_lin_m3_spring6mo_lag_AR2<-coef(summary(m3_spring6mo_lag_AR2))[2,2]
  
  numDF_lin_m1_monsoon6mo_lag_AR1<-anova(m1_monsoon6mo_lag_AR1,type="marginal")[2,1]
  denDF_lin_m1_monsoon6mo_lag_AR1<-anova(m1_monsoon6mo_lag_AR1,type="marginal")[2,2]
  F_lin_m1_monsoon6mo_lag_AR1<-anova(m1_monsoon6mo_lag_AR1,type="marginal")[2,3]
  P_lin_m1_monsoon6mo_lag_AR1<-anova(m1_monsoon6mo_lag_AR1,type="marginal")[2,4]
  ParamEst_lin_m1_monsoon6mo_lag_AR1<-coef(summary(m1_monsoon6mo_lag_AR1))[2,1]
  SE_lin_m1_monsoon6mo_lag_AR1<-coef(summary(m1_monsoon6mo_lag_AR1))[2,2]
  
  numDF_lin_m1_monsoon6mo_lag_AR2<-anova(m1_monsoon6mo_lag_AR2,type="marginal")[2,1]
  denDF_lin_m1_monsoon6mo_lag_AR2<-anova(m1_monsoon6mo_lag_AR2,type="marginal")[2,2]
  F_lin_m1_monsoon6mo_lag_AR2<-anova(m1_monsoon6mo_lag_AR2,type="marginal")[2,3]
  P_lin_m1_monsoon6mo_lag_AR2<-anova(m1_monsoon6mo_lag_AR2,type="marginal")[2,4]
  ParamEst_lin_m1_monsoon6mo_lag_AR2<-coef(summary(m1_monsoon6mo_lag_AR2))[2,1]
  SE_lin_m1_monsoon6mo_lag_AR2<-coef(summary(m1_monsoon6mo_lag_AR2))[2,2]
  
  numDF_lin_m2_monsoon6mo_lag_AR1<-anova(m2_monsoon6mo_lag_AR1,type="marginal")[2,1]
  denDF_lin_m2_monsoon6mo_lag_AR1<-anova(m2_monsoon6mo_lag_AR1,type="marginal")[2,2]
  F_lin_m2_monsoon6mo_lag_AR1<-anova(m2_monsoon6mo_lag_AR1,type="marginal")[2,3]
  P_lin_m2_monsoon6mo_lag_AR1<-anova(m2_monsoon6mo_lag_AR1,type="marginal")[2,4]
  ParamEst_lin_m2_monsoon6mo_lag_AR1<-coef(summary(m2_monsoon6mo_lag_AR1))[2,1]
  SE_lin_m2_monsoon6mo_lag_AR1<-coef(summary(m2_monsoon6mo_lag_AR1))[2,2]
  
  numDF_lin_m2_monsoon6mo_lag_AR2<-anova(m2_monsoon6mo_lag_AR2,type="marginal")[2,1]
  denDF_lin_m2_monsoon6mo_lag_AR2<-anova(m2_monsoon6mo_lag_AR2,type="marginal")[2,2]
  F_lin_m2_monsoon6mo_lag_AR2<-anova(m2_monsoon6mo_lag_AR2,type="marginal")[2,3]
  P_lin_m2_monsoon6mo_lag_AR2<-anova(m2_monsoon6mo_lag_AR2,type="marginal")[2,4]
  ParamEst_lin_m2_monsoon6mo_lag_AR2<-coef(summary(m2_monsoon6mo_lag_AR2))[2,1]
  SE_lin_m2_monsoon6mo_lag_AR2<-coef(summary(m2_monsoon6mo_lag_AR2))[2,2]
  
  numDF_lin_m3_monsoon6mo_lag_AR1<-anova(m3_monsoon6mo_lag_AR1,type="marginal")[2,1]
  denDF_lin_m3_monsoon6mo_lag_AR1<-anova(m3_monsoon6mo_lag_AR1,type="marginal")[2,2]
  F_lin_m3_monsoon6mo_lag_AR1<-anova(m3_monsoon6mo_lag_AR1,type="marginal")[2,3]
  P_lin_m3_monsoon6mo_lag_AR1<-anova(m3_monsoon6mo_lag_AR1,type="marginal")[2,4]
  ParamEst_lin_m3_monsoon6mo_lag_AR1<-coef(summary(m3_monsoon6mo_lag_AR1))[2,1]
  SE_lin_m3_monsoon6mo_lag_AR1<-coef(summary(m3_monsoon6mo_lag_AR1))[2,2]
  
  numDF_lin_m3_monsoon6mo_lag_AR2<-anova(m3_monsoon6mo_lag_AR2,type="marginal")[2,1]
  denDF_lin_m3_monsoon6mo_lag_AR2<-anova(m3_monsoon6mo_lag_AR2,type="marginal")[2,2]
  F_lin_m3_monsoon6mo_lag_AR2<-anova(m3_monsoon6mo_lag_AR2,type="marginal")[2,3]
  P_lin_m3_monsoon6mo_lag_AR2<-anova(m3_monsoon6mo_lag_AR2,type="marginal")[2,4]
  ParamEst_lin_m3_monsoon6mo_lag_AR2<-coef(summary(m3_monsoon6mo_lag_AR2))[2,1]
  SE_lin_m3_monsoon6mo_lag_AR2<-coef(summary(m3_monsoon6mo_lag_AR2))[2,2]
  
  Rsquared_marginal_m_null<-rsquared(m_null)[1,5]
  Rsquared_marginal_m1_spring6mo<-rsquared(m1_spring6mo)[1,5]
  Rsquared_marginal_m2_spring6mo<-rsquared(m2_spring6mo)[1,5]
  Rsquared_marginal_m3_spring6mo<-rsquared(m3_spring6mo)[1,5]
  Rsquared_marginal_m1_spring6mo_AR1<-rsquared(m1_spring6mo_AR1)[1,5]
  Rsquared_marginal_m1_spring6mo_AR2<-rsquared(m1_spring6mo_AR2)[1,5]
  Rsquared_marginal_m2_spring6mo_AR1<-rsquared(m2_spring6mo_AR1)[1,5]
  Rsquared_marginal_m2_spring6mo_AR2<-rsquared(m2_spring6mo_AR2)[1,5]
  Rsquared_marginal_m3_spring6mo_AR1<-rsquared(m3_spring6mo_AR1)[1,5]
  Rsquared_marginal_m3_spring6mo_AR2<-rsquared(m3_spring6mo_AR2)[1,5]
  Rsquared_marginal_m1_monsoon6mo<-rsquared(m1_monsoon6mo)[1,5]
  Rsquared_marginal_m2_monsoon6mo<-rsquared(m2_monsoon6mo)[1,5]
  Rsquared_marginal_m3_monsoon6mo<-rsquared(m3_monsoon6mo)[1,5]
  Rsquared_marginal_m1_monsoon6mo_AR1<-rsquared(m1_monsoon6mo_AR1)[1,5]
  Rsquared_marginal_m1_monsoon6mo_AR2<-rsquared(m1_monsoon6mo_AR2)[1,5]
  Rsquared_marginal_m2_monsoon6mo_AR1<-rsquared(m2_monsoon6mo_AR1)[1,5]
  Rsquared_marginal_m2_monsoon6mo_AR2<-rsquared(m2_monsoon6mo_AR2)[1,5]
  Rsquared_marginal_m3_monsoon6mo_AR1<-rsquared(m3_monsoon6mo_AR1)[1,5]
  Rsquared_marginal_m3_monsoon6mo_AR2<-rsquared(m3_monsoon6mo_AR2)[1,5]
  Rsquared_marginal_m1_spring6mo_lag<-rsquared(m1_spring6mo_lag)[1,5]
  Rsquared_marginal_m2_spring6mo_lag<-rsquared(m2_spring6mo_lag)[1,5]
  Rsquared_marginal_m3_spring6mo_lag<-rsquared(m3_spring6mo_lag)[1,5]
  Rsquared_marginal_m1_spring6mo_lag_AR1<-rsquared(m1_spring6mo_lag_AR1)[1,5]
  Rsquared_marginal_m1_spring6mo_lag_AR2<-rsquared(m1_spring6mo_lag_AR2)[1,5]
  Rsquared_marginal_m2_spring6mo_lag_AR1<-rsquared(m2_spring6mo_lag_AR1)[1,5]
  Rsquared_marginal_m2_spring6mo_lag_AR2<-rsquared(m2_spring6mo_lag_AR2)[1,5]
  Rsquared_marginal_m3_spring6mo_lag_AR1<-rsquared(m3_spring6mo_lag_AR1)[1,5]
  Rsquared_marginal_m3_spring6mo_lag_AR2<-rsquared(m3_spring6mo_lag_AR2)[1,5]
  Rsquared_marginal_m1_monsoon6mo_lag<-rsquared(m1_monsoon6mo_lag)[1,5]
  Rsquared_marginal_m2_monsoon6mo_lag<-rsquared(m2_monsoon6mo_lag)[1,5]
  Rsquared_marginal_m3_monsoon6mo_lag<-rsquared(m3_monsoon6mo_lag)[1,5]
  Rsquared_marginal_m1_monsoon6mo_lag_AR1<-rsquared(m1_monsoon6mo_lag_AR1)[1,5]
  Rsquared_marginal_m1_monsoon6mo_lag_AR2<-rsquared(m1_monsoon6mo_lag_AR2)[1,5]
  Rsquared_marginal_m2_monsoon6mo_lag_AR1<-rsquared(m2_monsoon6mo_lag_AR1)[1,5]
  Rsquared_marginal_m2_monsoon6mo_lag_AR2<-rsquared(m2_monsoon6mo_lag_AR2)[1,5]
  Rsquared_marginal_m3_monsoon6mo_lag_AR1<-rsquared(m3_monsoon6mo_lag_AR1)[1,5]
  Rsquared_marginal_m3_monsoon6mo_lag_AR2<-rsquared(m3_monsoon6mo_lag_AR2)[1,5]
  
  Rsquared_conditional_m_null<-rsquared(m_null)[1,6]
  Rsquared_conditional_m1_spring6mo<-rsquared(m1_spring6mo)[1,6]
  Rsquared_conditional_m2_spring6mo<-rsquared(m2_spring6mo)[1,6]
  Rsquared_conditional_m3_spring6mo<-rsquared(m3_spring6mo)[1,6]
  Rsquared_conditional_m1_spring6mo_AR1<-rsquared(m1_spring6mo_AR1)[1,6]
  Rsquared_conditional_m1_spring6mo_AR2<-rsquared(m1_spring6mo_AR2)[1,6]
  Rsquared_conditional_m2_spring6mo_AR1<-rsquared(m2_spring6mo_AR1)[1,6]
  Rsquared_conditional_m2_spring6mo_AR2<-rsquared(m2_spring6mo_AR2)[1,6]
  Rsquared_conditional_m3_spring6mo_AR1<-rsquared(m3_spring6mo_AR1)[1,6]
  Rsquared_conditional_m3_spring6mo_AR2<-rsquared(m3_spring6mo_AR2)[1,6]
  Rsquared_conditional_m1_monsoon6mo<-rsquared(m1_monsoon6mo)[1,6]
  Rsquared_conditional_m2_monsoon6mo<-rsquared(m2_monsoon6mo)[1,6]
  Rsquared_conditional_m3_monsoon6mo<-rsquared(m3_monsoon6mo)[1,6]
  Rsquared_conditional_m1_monsoon6mo_AR1<-rsquared(m1_monsoon6mo_AR1)[1,6]
  Rsquared_conditional_m1_monsoon6mo_AR2<-rsquared(m1_monsoon6mo_AR2)[1,6]
  Rsquared_conditional_m2_monsoon6mo_AR1<-rsquared(m2_monsoon6mo_AR1)[1,6]
  Rsquared_conditional_m2_monsoon6mo_AR2<-rsquared(m2_monsoon6mo_AR2)[1,6]
  Rsquared_conditional_m3_monsoon6mo_AR1<-rsquared(m3_monsoon6mo_AR1)[1,6]
  Rsquared_conditional_m3_monsoon6mo_AR2<-rsquared(m3_monsoon6mo_AR2)[1,6]
  Rsquared_conditional_m1_spring6mo_lag<-rsquared(m1_spring6mo_lag)[1,6]
  Rsquared_conditional_m2_spring6mo_lag<-rsquared(m2_spring6mo_lag)[1,6]
  Rsquared_conditional_m3_spring6mo_lag<-rsquared(m3_spring6mo_lag)[1,6]
  Rsquared_conditional_m1_spring6mo_lag_AR1<-rsquared(m1_spring6mo_lag_AR1)[1,6]
  Rsquared_conditional_m1_spring6mo_lag_AR2<-rsquared(m1_spring6mo_lag_AR2)[1,6]
  Rsquared_conditional_m2_spring6mo_lag_AR1<-rsquared(m2_spring6mo_lag_AR1)[1,6]
  Rsquared_conditional_m2_spring6mo_lag_AR2<-rsquared(m2_spring6mo_lag_AR2)[1,6]
  Rsquared_conditional_m3_spring6mo_lag_AR1<-rsquared(m3_spring6mo_lag_AR1)[1,6]
  Rsquared_conditional_m3_spring6mo_lag_AR2<-rsquared(m3_spring6mo_lag_AR2)[1,6]
  Rsquared_conditional_m1_monsoon6mo_lag<-rsquared(m1_monsoon6mo_lag)[1,6]
  Rsquared_conditional_m2_monsoon6mo_lag<-rsquared(m2_monsoon6mo_lag)[1,6]
  Rsquared_conditional_m3_monsoon6mo_lag<-rsquared(m3_monsoon6mo_lag)[1,6]
  Rsquared_conditional_m1_monsoon6mo_lag_AR1<-rsquared(m1_monsoon6mo_lag_AR1)[1,6]
  Rsquared_conditional_m1_monsoon6mo_lag_AR2<-rsquared(m1_monsoon6mo_lag_AR2)[1,6]
  Rsquared_conditional_m2_monsoon6mo_lag_AR1<-rsquared(m2_monsoon6mo_lag_AR1)[1,6]
  Rsquared_conditional_m2_monsoon6mo_lag_AR2<-rsquared(m2_monsoon6mo_lag_AR2)[1,6]
  Rsquared_conditional_m3_monsoon6mo_lag_AR1<-rsquared(m3_monsoon6mo_lag_AR1)[1,6]
  Rsquared_conditional_m3_monsoon6mo_lag_AR2<-rsquared(m3_monsoon6mo_lag_AR2)[1,6]
  
  numDF_quad_m2_spring6mo<-anova(m2_spring6mo,type="marginal")[3,1]
  numDF_quad_m2_spring6mo_AR1<-anova(m2_spring6mo_AR1,type="marginal")[3,1]
  numDF_quad_m2_spring6mo_AR2<-anova(m2_spring6mo_AR2,type="marginal")[3,1]
  numDF_quad_m2_monsoon6mo<-anova(m2_monsoon6mo,type="marginal")[3,1]
  numDF_quad_m2_monsoon6mo_AR1<-anova(m2_monsoon6mo_AR1,type="marginal")[3,1]
  numDF_quad_m2_monsoon6mo_AR2<-anova(m2_monsoon6mo_AR2,type="marginal")[3,1]
  numDF_quad_m2_spring6mo_lag<-anova(m2_spring6mo_lag,type="marginal")[3,1]
  numDF_quad_m2_spring6mo_lag_AR1<-anova(m2_spring6mo_lag_AR1,type="marginal")[3,1]
  numDF_quad_m2_spring6mo_lag_AR2<-anova(m2_spring6mo_lag_AR2,type="marginal")[3,1]
  numDF_quad_m2_monsoon6mo_lag<-anova(m2_monsoon6mo_lag,type="marginal")[3,1]
  numDF_quad_m2_monsoon6mo_lag_AR1<-anova(m2_monsoon6mo_lag_AR1,type="marginal")[3,1]
  numDF_quad_m2_monsoon6mo_lag_AR2<-anova(m2_monsoon6mo_lag_AR2,type="marginal")[3,1]
  
  numDF_quad_m3_spring6mo<-anova(m3_spring6mo,type="marginal")[3,1]
  numDF_quad_m3_spring6mo_AR1<-anova(m3_spring6mo_AR1,type="marginal")[3,1]
  numDF_quad_m3_spring6mo_AR2<-anova(m3_spring6mo_AR2,type="marginal")[3,1]
  numDF_quad_m3_monsoon6mo<-anova(m3_monsoon6mo,type="marginal")[3,1]
  numDF_quad_m3_monsoon6mo_AR1<-anova(m3_monsoon6mo_AR1,type="marginal")[3,1]
  numDF_quad_m3_monsoon6mo_AR2<-anova(m3_monsoon6mo_AR2,type="marginal")[3,1]
  numDF_quad_m3_spring6mo_lag<-anova(m3_spring6mo_lag,type="marginal")[3,1]
  numDF_quad_m3_spring6mo_lag_AR1<-anova(m3_spring6mo_lag_AR1,type="marginal")[3,1]
  numDF_quad_m3_spring6mo_lag_AR2<-anova(m3_spring6mo_lag_AR2,type="marginal")[3,1]
  numDF_quad_m3_monsoon6mo_lag<-anova(m3_monsoon6mo_lag,type="marginal")[3,1]
  numDF_quad_m3_monsoon6mo_lag_AR1<-anova(m3_monsoon6mo_lag_AR1,type="marginal")[3,1]
  numDF_quad_m3_monsoon6mo_lag_AR2<-anova(m3_monsoon6mo_lag_AR2,type="marginal")[3,1]
  
  denDF_quad_m2_spring6mo<-anova(m2_spring6mo,type="marginal")[3,2]
  denDF_quad_m2_spring6mo_AR1<-anova(m2_spring6mo_AR1,type="marginal")[3,2]
  denDF_quad_m2_spring6mo_AR2<-anova(m2_spring6mo_AR2,type="marginal")[3,2]
  denDF_quad_m2_monsoon6mo<-anova(m2_monsoon6mo,type="marginal")[3,2]
  denDF_quad_m2_monsoon6mo_AR1<-anova(m2_monsoon6mo_AR1,type="marginal")[3,2]
  denDF_quad_m2_monsoon6mo_AR2<-anova(m2_monsoon6mo_AR2,type="marginal")[3,2]
  denDF_quad_m2_spring6mo_lag<-anova(m2_spring6mo_lag,type="marginal")[3,2]
  denDF_quad_m2_spring6mo_lag_AR1<-anova(m2_spring6mo_lag_AR1,type="marginal")[3,2]
  denDF_quad_m2_spring6mo_lag_AR2<-anova(m2_spring6mo_lag_AR2,type="marginal")[3,2]
  denDF_quad_m2_monsoon6mo_lag<-anova(m2_monsoon6mo_lag,type="marginal")[3,2]
  denDF_quad_m2_monsoon6mo_lag_AR1<-anova(m2_monsoon6mo_lag_AR1,type="marginal")[3,2]
  denDF_quad_m2_monsoon6mo_lag_AR2<-anova(m2_monsoon6mo_lag_AR2,type="marginal")[3,2]
  
  denDF_quad_m3_spring6mo<-anova(m3_spring6mo,type="marginal")[3,2]
  denDF_quad_m3_spring6mo_AR1<-anova(m3_spring6mo_AR1,type="marginal")[3,2]
  denDF_quad_m3_spring6mo_AR2<-anova(m3_spring6mo_AR2,type="marginal")[3,2]
  denDF_quad_m3_monsoon6mo<-anova(m3_monsoon6mo,type="marginal")[3,2]
  denDF_quad_m3_monsoon6mo_AR1<-anova(m3_monsoon6mo_AR1,type="marginal")[3,2]
  denDF_quad_m3_monsoon6mo_AR2<-anova(m3_monsoon6mo_AR2,type="marginal")[3,2]
  denDF_quad_m3_spring6mo_lag<-anova(m3_spring6mo_lag,type="marginal")[3,2]
  denDF_quad_m3_spring6mo_lag_AR1<-anova(m3_spring6mo_lag_AR1,type="marginal")[3,2]
  denDF_quad_m3_spring6mo_lag_AR2<-anova(m3_spring6mo_lag_AR2,type="marginal")[3,2]
  denDF_quad_m3_monsoon6mo_lag<-anova(m3_monsoon6mo_lag,type="marginal")[3,2]
  denDF_quad_m3_monsoon6mo_lag_AR1<-anova(m3_monsoon6mo_lag_AR1,type="marginal")[3,2]
  denDF_quad_m3_monsoon6mo_lag_AR2<-anova(m3_monsoon6mo_lag_AR2,type="marginal")[3,2]
  
  F_quad_m2_spring6mo<-anova(m2_spring6mo,type="marginal")[3,3]
  F_quad_m2_spring6mo_AR1<-anova(m2_spring6mo_AR1,type="marginal")[3,3]
  F_quad_m2_spring6mo_AR2<-anova(m2_spring6mo_AR2,type="marginal")[3,3]
  F_quad_m2_monsoon6mo<-anova(m2_monsoon6mo,type="marginal")[3,3]
  F_quad_m2_monsoon6mo_AR1<-anova(m2_monsoon6mo_AR1,type="marginal")[3,3]
  F_quad_m2_monsoon6mo_AR2<-anova(m2_monsoon6mo_AR2,type="marginal")[3,3]
  F_quad_m2_spring6mo_lag<-anova(m2_spring6mo_lag,type="marginal")[3,3]
  F_quad_m2_spring6mo_lag_AR1<-anova(m2_spring6mo_lag_AR1,type="marginal")[3,3]
  F_quad_m2_spring6mo_lag_AR2<-anova(m2_spring6mo_lag_AR2,type="marginal")[3,3]
  F_quad_m2_monsoon6mo_lag<-anova(m2_monsoon6mo_lag,type="marginal")[3,3]
  F_quad_m2_monsoon6mo_lag_AR1<-anova(m2_monsoon6mo_lag_AR1,type="marginal")[3,3]
  F_quad_m2_monsoon6mo_lag_AR2<-anova(m2_monsoon6mo_lag_AR2,type="marginal")[3,3]
  
  F_quad_m3_spring6mo<-anova(m3_spring6mo,type="marginal")[3,3]
  F_quad_m3_spring6mo_AR1<-anova(m3_spring6mo_AR1,type="marginal")[3,3]
  F_quad_m3_spring6mo_AR2<-anova(m3_spring6mo_AR2,type="marginal")[3,3]
  F_quad_m3_monsoon6mo<-anova(m3_monsoon6mo,type="marginal")[3,3]
  F_quad_m3_monsoon6mo_AR1<-anova(m3_monsoon6mo_AR1,type="marginal")[3,3]
  F_quad_m3_monsoon6mo_AR2<-anova(m3_monsoon6mo_AR2,type="marginal")[3,3]
  F_quad_m3_spring6mo_lag<-anova(m3_spring6mo_lag,type="marginal")[3,3]
  F_quad_m3_spring6mo_lag_AR1<-anova(m3_spring6mo_lag_AR1,type="marginal")[3,3]
  F_quad_m3_spring6mo_lag_AR2<-anova(m3_spring6mo_lag_AR2,type="marginal")[3,3]
  F_quad_m3_monsoon6mo_lag<-anova(m3_monsoon6mo_lag,type="marginal")[3,3]
  F_quad_m3_monsoon6mo_lag_AR1<-anova(m3_monsoon6mo_lag_AR1,type="marginal")[3,3]
  F_quad_m3_monsoon6mo_lag_AR2<-anova(m3_monsoon6mo_lag_AR2,type="marginal")[3,3]
  
  P_quad_m2_spring6mo<-anova(m2_spring6mo,type="marginal")[3,4]
  P_quad_m2_spring6mo_AR1<-anova(m2_spring6mo_AR1,type="marginal")[3,4]
  P_quad_m2_spring6mo_AR2<-anova(m2_spring6mo_AR2,type="marginal")[3,4]
  P_quad_m2_monsoon6mo<-anova(m2_monsoon6mo,type="marginal")[3,4]
  P_quad_m2_monsoon6mo_AR1<-anova(m2_monsoon6mo_AR1,type="marginal")[3,4]
  P_quad_m2_monsoon6mo_AR2<-anova(m2_monsoon6mo_AR2,type="marginal")[3,4]
  P_quad_m2_spring6mo_lag<-anova(m2_spring6mo_lag,type="marginal")[3,4]
  P_quad_m2_spring6mo_lag_AR1<-anova(m2_spring6mo_lag_AR1,type="marginal")[3,4]
  P_quad_m2_spring6mo_lag_AR2<-anova(m2_spring6mo_lag_AR2,type="marginal")[3,4]
  P_quad_m2_monsoon6mo_lag<-anova(m2_monsoon6mo_lag,type="marginal")[3,4]
  P_quad_m2_monsoon6mo_lag_AR1<-anova(m2_monsoon6mo_lag_AR1,type="marginal")[3,4]
  P_quad_m2_monsoon6mo_lag_AR2<-anova(m2_monsoon6mo_lag_AR2,type="marginal")[3,4]
  
  P_quad_m3_spring6mo<-anova(m3_spring6mo,type="marginal")[3,4]
  P_quad_m3_spring6mo_AR1<-anova(m3_spring6mo_AR1,type="marginal")[3,4]
  P_quad_m3_spring6mo_AR2<-anova(m3_spring6mo_AR2,type="marginal")[3,4]
  P_quad_m3_monsoon6mo<-anova(m3_monsoon6mo,type="marginal")[3,4]
  P_quad_m3_monsoon6mo_AR1<-anova(m3_monsoon6mo_AR1,type="marginal")[3,4]
  P_quad_m3_monsoon6mo_AR2<-anova(m3_monsoon6mo_AR2,type="marginal")[3,4]
  P_quad_m3_spring6mo_lag<-anova(m3_spring6mo_lag,type="marginal")[3,4]
  P_quad_m3_spring6mo_lag_AR1<-anova(m3_spring6mo_lag_AR1,type="marginal")[3,4]
  P_quad_m3_spring6mo_lag_AR2<-anova(m3_spring6mo_lag_AR2,type="marginal")[3,4]
  P_quad_m3_monsoon6mo_lag<-anova(m3_monsoon6mo_lag,type="marginal")[3,4]
  P_quad_m3_monsoon6mo_lag_AR1<-anova(m3_monsoon6mo_lag_AR1,type="marginal")[3,4]
  P_quad_m3_monsoon6mo_lag_AR2<-anova(m3_monsoon6mo_lag_AR2,type="marginal")[3,4]
  
  ParamEst_quad_m2_spring6mo<-coef(summary(m2_spring6mo))[3,1]
  ParamEst_quad_m2_spring6mo_AR1<-coef(summary(m2_spring6mo_AR1))[3,1]
  ParamEst_quad_m2_spring6mo_AR2<-coef(summary(m2_spring6mo_AR2))[3,1]
  ParamEst_quad_m2_monsoon6mo<-coef(summary(m2_monsoon6mo))[3,1]
  ParamEst_quad_m2_monsoon6mo_AR1<-coef(summary(m2_monsoon6mo_AR1))[3,1]
  ParamEst_quad_m2_monsoon6mo_AR2<-coef(summary(m2_monsoon6mo_AR2))[3,1]
  ParamEst_quad_m2_spring6mo_lag<-coef(summary(m2_spring6mo_lag))[3,1]
  ParamEst_quad_m2_spring6mo_lag_AR1<-coef(summary(m2_spring6mo_lag_AR1))[3,1]
  ParamEst_quad_m2_spring6mo_lag_AR2<-coef(summary(m2_spring6mo_lag_AR2))[3,1]
  ParamEst_quad_m2_monsoon6mo_lag<-coef(summary(m2_monsoon6mo_lag))[3,1]
  ParamEst_quad_m2_monsoon6mo_lag_AR1<-coef(summary(m2_monsoon6mo_lag_AR1))[3,1]
  ParamEst_quad_m2_monsoon6mo_lag_AR2<-coef(summary(m2_monsoon6mo_lag_AR2))[3,1]
  
  ParamEst_quad_m3_spring6mo<-coef(summary(m3_spring6mo))[3,1]
  ParamEst_quad_m3_spring6mo_AR1<-coef(summary(m3_spring6mo_AR1))[3,1]
  ParamEst_quad_m3_spring6mo_AR2<-coef(summary(m3_spring6mo_AR2))[3,1]
  ParamEst_quad_m3_monsoon6mo<-coef(summary(m3_monsoon6mo))[3,1]
  ParamEst_quad_m3_monsoon6mo_AR1<-coef(summary(m3_monsoon6mo_AR1))[3,1]
  ParamEst_quad_m3_monsoon6mo_AR2<-coef(summary(m3_monsoon6mo_AR2))[3,1]
  ParamEst_quad_m3_spring6mo_lag<-coef(summary(m3_spring6mo_lag))[3,1]
  ParamEst_quad_m3_spring6mo_lag_AR1<-coef(summary(m3_spring6mo_lag_AR1))[3,1]
  ParamEst_quad_m3_spring6mo_lag_AR2<-coef(summary(m3_spring6mo_lag_AR2))[3,1]
  ParamEst_quad_m3_monsoon6mo_lag<-coef(summary(m3_monsoon6mo_lag))[3,1]
  ParamEst_quad_m3_monsoon6mo_lag_AR1<-coef(summary(m3_monsoon6mo_lag_AR1))[3,1]
  ParamEst_quad_m3_monsoon6mo_lag_AR2<-coef(summary(m3_monsoon6mo_lag_AR2))[3,1]
  
  SE_quad_m2_spring6mo<-coef(summary(m2_spring6mo))[3,2]
  SE_quad_m2_spring6mo_AR1<-coef(summary(m2_spring6mo_AR1))[3,2]
  SE_quad_m2_spring6mo_AR2<-coef(summary(m2_spring6mo_AR2))[3,2]
  SE_quad_m2_monsoon6mo<-coef(summary(m2_monsoon6mo))[3,2]
  SE_quad_m2_monsoon6mo_AR1<-coef(summary(m2_monsoon6mo_AR1))[3,2]
  SE_quad_m2_monsoon6mo_AR2<-coef(summary(m2_monsoon6mo_AR2))[3,2]
  SE_quad_m2_spring6mo_lag<-coef(summary(m2_spring6mo_lag))[3,2]
  SE_quad_m2_spring6mo_lag_AR1<-coef(summary(m2_spring6mo_lag_AR1))[3,2]
  SE_quad_m2_spring6mo_lag_AR2<-coef(summary(m2_spring6mo_lag_AR2))[3,2]
  SE_quad_m2_monsoon6mo_lag<-coef(summary(m2_monsoon6mo_lag))[3,2]
  SE_quad_m2_monsoon6mo_lag_AR1<-coef(summary(m2_monsoon6mo_lag_AR1))[3,2]
  SE_quad_m2_monsoon6mo_lag_AR2<-coef(summary(m2_monsoon6mo_lag_AR2))[3,2]
  
  SE_quad_m3_spring6mo<-coef(summary(m3_spring6mo))[3,2]
  SE_quad_m3_spring6mo_AR1<-coef(summary(m3_spring6mo_AR1))[3,2]
  SE_quad_m3_spring6mo_AR2<-coef(summary(m3_spring6mo_AR2))[3,2]
  SE_quad_m3_monsoon6mo<-coef(summary(m3_monsoon6mo))[3,2]
  SE_quad_m3_monsoon6mo_AR1<-coef(summary(m3_monsoon6mo_AR1))[3,2]
  SE_quad_m3_monsoon6mo_AR2<-coef(summary(m3_monsoon6mo_AR2))[3,2]
  SE_quad_m3_spring6mo_lag<-coef(summary(m3_spring6mo_lag))[3,2]
  SE_quad_m3_spring6mo_lag_AR1<-coef(summary(m3_spring6mo_lag_AR1))[3,2]
  SE_quad_m3_spring6mo_lag_AR2<-coef(summary(m3_spring6mo_lag_AR2))[3,2]
  SE_quad_m3_monsoon6mo_lag<-coef(summary(m3_monsoon6mo_lag))[3,2]
  SE_quad_m3_monsoon6mo_lag_AR1<-coef(summary(m3_monsoon6mo_lag_AR1))[3,2]
  SE_quad_m3_monsoon6mo_lag_AR2<-coef(summary(m3_monsoon6mo_lag_AR2))[3,2]
  
  numDF_cub_m3_spring6mo<-anova(m3_spring6mo,type="marginal")[4,1]
  numDF_cub_m3_spring6mo_AR1<-anova(m3_spring6mo_AR1,type="marginal")[4,1]
  numDF_cub_m3_spring6mo_AR2<-anova(m3_spring6mo_AR2,type="marginal")[4,1]
  numDF_cub_m3_monsoon6mo<-anova(m3_monsoon6mo,type="marginal")[4,1]
  numDF_cub_m3_monsoon6mo_AR1<-anova(m3_monsoon6mo_AR1,type="marginal")[4,1]
  numDF_cub_m3_monsoon6mo_AR2<-anova(m3_monsoon6mo_AR2,type="marginal")[4,1]
  numDF_cub_m3_spring6mo_lag<-anova(m3_spring6mo_lag,type="marginal")[4,1]
  numDF_cub_m3_spring6mo_lag_AR1<-anova(m3_spring6mo_lag_AR1,type="marginal")[4,1]
  numDF_cub_m3_spring6mo_lag_AR2<-anova(m3_spring6mo_lag_AR2,type="marginal")[4,1]
  numDF_cub_m3_monsoon6mo_lag<-anova(m3_monsoon6mo_lag,type="marginal")[4,1]
  numDF_cub_m3_monsoon6mo_lag_AR1<-anova(m3_monsoon6mo_lag_AR1,type="marginal")[4,1]
  numDF_cub_m3_monsoon6mo_lag_AR2<-anova(m3_monsoon6mo_lag_AR2,type="marginal")[4,1]
  
  denDF_cub_m3_spring6mo<-anova(m3_spring6mo,type="marginal")[4,2]
  denDF_cub_m3_spring6mo_AR1<-anova(m3_spring6mo_AR1,type="marginal")[4,2]
  denDF_cub_m3_spring6mo_AR2<-anova(m3_spring6mo_AR2,type="marginal")[4,2]
  denDF_cub_m3_monsoon6mo<-anova(m3_monsoon6mo,type="marginal")[4,2]
  denDF_cub_m3_monsoon6mo_AR1<-anova(m3_monsoon6mo_AR1,type="marginal")[4,2]
  denDF_cub_m3_monsoon6mo_AR2<-anova(m3_monsoon6mo_AR2,type="marginal")[4,2]
  denDF_cub_m3_spring6mo_lag<-anova(m3_spring6mo_lag,type="marginal")[4,2]
  denDF_cub_m3_spring6mo_lag_AR1<-anova(m3_spring6mo_lag_AR1,type="marginal")[4,2]
  denDF_cub_m3_spring6mo_lag_AR2<-anova(m3_spring6mo_lag_AR2,type="marginal")[4,2]
  denDF_cub_m3_monsoon6mo_lag<-anova(m3_monsoon6mo_lag,type="marginal")[4,2]
  denDF_cub_m3_monsoon6mo_lag_AR1<-anova(m3_monsoon6mo_lag_AR1,type="marginal")[4,2]
  denDF_cub_m3_monsoon6mo_lag_AR2<-anova(m3_monsoon6mo_lag_AR2,type="marginal")[4,2]
  
  F_cub_m3_spring6mo<-anova(m3_spring6mo,type="marginal")[4,3]
  F_cub_m3_spring6mo_AR1<-anova(m3_spring6mo_AR1,type="marginal")[4,3]
  F_cub_m3_spring6mo_AR2<-anova(m3_spring6mo_AR2,type="marginal")[4,3]
  F_cub_m3_monsoon6mo<-anova(m3_monsoon6mo,type="marginal")[4,3]
  F_cub_m3_monsoon6mo_AR1<-anova(m3_monsoon6mo_AR1,type="marginal")[4,3]
  F_cub_m3_monsoon6mo_AR2<-anova(m3_monsoon6mo_AR2,type="marginal")[4,3]
  F_cub_m3_spring6mo_lag<-anova(m3_spring6mo_lag,type="marginal")[4,3]
  F_cub_m3_spring6mo_lag_AR1<-anova(m3_spring6mo_lag_AR1,type="marginal")[4,3]
  F_cub_m3_spring6mo_lag_AR2<-anova(m3_spring6mo_lag_AR2,type="marginal")[4,3]
  F_cub_m3_monsoon6mo_lag<-anova(m3_monsoon6mo_lag,type="marginal")[4,3]
  F_cub_m3_monsoon6mo_lag_AR1<-anova(m3_monsoon6mo_lag_AR1,type="marginal")[4,3]
  F_cub_m3_monsoon6mo_lag_AR2<-anova(m3_monsoon6mo_lag_AR2,type="marginal")[4,3]
  
  P_cub_m3_spring6mo<-anova(m3_spring6mo,type="marginal")[4,4]
  P_cub_m3_spring6mo_AR1<-anova(m3_spring6mo_AR1,type="marginal")[4,4]
  P_cub_m3_spring6mo_AR2<-anova(m3_spring6mo_AR2,type="marginal")[4,4]
  P_cub_m3_monsoon6mo<-anova(m3_monsoon6mo,type="marginal")[4,4]
  P_cub_m3_monsoon6mo_AR1<-anova(m3_monsoon6mo_AR1,type="marginal")[4,4]
  P_cub_m3_monsoon6mo_AR2<-anova(m3_monsoon6mo_AR2,type="marginal")[4,4]
  P_cub_m3_spring6mo_lag<-anova(m3_spring6mo_lag,type="marginal")[4,4]
  P_cub_m3_spring6mo_lag_AR1<-anova(m3_spring6mo_lag_AR1,type="marginal")[4,4]
  P_cub_m3_spring6mo_lag_AR2<-anova(m3_spring6mo_lag_AR2,type="marginal")[4,4]
  P_cub_m3_monsoon6mo_lag<-anova(m3_monsoon6mo_lag,type="marginal")[4,4]
  P_cub_m3_monsoon6mo_lag_AR1<-anova(m3_monsoon6mo_lag_AR1,type="marginal")[4,4]
  P_cub_m3_monsoon6mo_lag_AR2<-anova(m3_monsoon6mo_lag_AR2,type="marginal")[4,4]
  
  ParamEst_cub_m3_spring6mo<-coef(summary(m3_spring6mo))[4,1]
  ParamEst_cub_m3_spring6mo_AR1<-coef(summary(m3_spring6mo_AR1))[4,1]
  ParamEst_cub_m3_spring6mo_AR2<-coef(summary(m3_spring6mo_AR2))[4,1]
  ParamEst_cub_m3_monsoon6mo<-coef(summary(m3_monsoon6mo))[4,1]
  ParamEst_cub_m3_monsoon6mo_AR1<-coef(summary(m3_monsoon6mo_AR1))[4,1]
  ParamEst_cub_m3_monsoon6mo_AR2<-coef(summary(m3_monsoon6mo_AR2))[4,1]
  ParamEst_cub_m3_spring6mo_lag<-coef(summary(m3_spring6mo_lag))[4,1]
  ParamEst_cub_m3_spring6mo_lag_AR1<-coef(summary(m3_spring6mo_lag_AR1))[4,1]
  ParamEst_cub_m3_spring6mo_lag_AR2<-coef(summary(m3_spring6mo_lag_AR2))[4,1]
  ParamEst_cub_m3_monsoon6mo_lag<-coef(summary(m3_monsoon6mo_lag))[4,1]
  ParamEst_cub_m3_monsoon6mo_lag_AR1<-coef(summary(m3_monsoon6mo_lag_AR1))[4,1]
  ParamEst_cub_m3_monsoon6mo_lag_AR2<-coef(summary(m3_monsoon6mo_lag_AR2))[4,1]
  
  SE_cub_m3_spring6mo<-coef(summary(m3_spring6mo))[4,2]
  SE_cub_m3_spring6mo_AR1<-coef(summary(m3_spring6mo_AR1))[4,2]
  SE_cub_m3_spring6mo_AR2<-coef(summary(m3_spring6mo_AR2))[4,2]
  SE_cub_m3_monsoon6mo<-coef(summary(m3_monsoon6mo))[4,2]
  SE_cub_m3_monsoon6mo_AR1<-coef(summary(m3_monsoon6mo_AR1))[4,2]
  SE_cub_m3_monsoon6mo_AR2<-coef(summary(m3_monsoon6mo_AR2))[4,2]
  SE_cub_m3_spring6mo_lag<-coef(summary(m3_spring6mo_lag))[4,2]
  SE_cub_m3_spring6mo_lag_AR1<-coef(summary(m3_spring6mo_lag_AR1))[4,2]
  SE_cub_m3_spring6mo_lag_AR2<-coef(summary(m3_spring6mo_lag_AR2))[4,2]
  SE_cub_m3_monsoon6mo_lag<-coef(summary(m3_monsoon6mo_lag))[4,2]
  SE_cub_m3_monsoon6mo_lag_AR1<-coef(summary(m3_monsoon6mo_lag_AR1))[4,2]
  SE_cub_m3_monsoon6mo_lag_AR2<-coef(summary(m3_monsoon6mo_lag_AR2))[4,2]
  
  # bind all target output values together
  output_id<-cbind(speciesCode, 
                   ecosystem, 
                   P_shapiro_null,
                   dAICc_m1_monsoon6mo,
                   dAICc_m1_monsoon6mo_lag,
                   dAICc_m1_spring6mo,
                   dAICc_m1_spring6mo_lag,
                   dAICc_m2_monsoon6mo,
                   dAICc_m2_monsoon6mo_lag,
                   dAICc_m2_spring6mo,
                   dAICc_m2_spring6mo_lag,
                   dAICc_m3_monsoon6mo,
                   dAICc_m3_monsoon6mo_lag,
                   dAICc_m3_spring6mo,
                   dAICc_m3_spring6mo_lag,
                   dAICc_null,
                   dAICc_m1_spring6mo_AR1,
                   dAICc_m1_spring6mo_AR2,
                   dAICc_m2_spring6mo_AR1,
                   dAICc_m2_spring6mo_AR2,
                   dAICc_m3_spring6mo_AR1,
                   dAICc_m3_spring6mo_AR2,
                   dAICc_m1_monsoon6mo_AR1,
                   dAICc_m1_monsoon6mo_AR2,
                   dAICc_m2_monsoon6mo_AR1,
                   dAICc_m2_monsoon6mo_AR2,
                   dAICc_m3_monsoon6mo_AR1,
                   dAICc_m3_monsoon6mo_AR2,
                   dAICc_m1_spring6mo_lag_AR1,
                   dAICc_m1_spring6mo_lag_AR2,
                   dAICc_m2_spring6mo_lag_AR1,
                   dAICc_m2_spring6mo_lag_AR2,
                   dAICc_m3_spring6mo_lag_AR1,
                   dAICc_m3_spring6mo_lag_AR2,
                   dAICc_m1_monsoon6mo_lag_AR1,
                   dAICc_m1_monsoon6mo_lag_AR2,
                   dAICc_m2_monsoon6mo_lag_AR1,
                   dAICc_m2_monsoon6mo_lag_AR2,
                   dAICc_m3_monsoon6mo_lag_AR1,
                   dAICc_m3_monsoon6mo_lag_AR2,
                   numDF_lin_m1_spring6mo,
                   denDF_lin_m1_spring6mo,
                   F_lin_m1_spring6mo,
                   P_lin_m1_spring6mo,
                   numDF_lin_m2_spring6mo,
                   denDF_lin_m2_spring6mo,
                   F_lin_m2_spring6mo,
                   P_lin_m2_spring6mo,
                   numDF_lin_m3_spring6mo,
                   denDF_lin_m3_spring6mo,
                   F_lin_m3_spring6mo,
                   P_lin_m3_spring6mo,
                   numDF_lin_m1_monsoon6mo,
                   denDF_lin_m1_monsoon6mo,
                   F_lin_m1_monsoon6mo,
                   P_lin_m1_monsoon6mo,
                   numDF_lin_m2_monsoon6mo,
                   denDF_lin_m2_monsoon6mo,
                   F_lin_m2_monsoon6mo,
                   P_lin_m2_monsoon6mo,
                   numDF_lin_m3_monsoon6mo,
                   denDF_lin_m3_monsoon6mo,
                   F_lin_m3_monsoon6mo,
                   P_lin_m3_monsoon6mo,
                   numDF_lin_m1_spring6mo_lag,
                   denDF_lin_m1_spring6mo_lag,
                   F_lin_m1_spring6mo_lag,
                   P_lin_m1_spring6mo_lag,
                   numDF_lin_m2_spring6mo_lag,
                   denDF_lin_m2_spring6mo_lag,
                   F_lin_m2_spring6mo_lag,
                   P_lin_m2_spring6mo_lag,
                   numDF_lin_m3_spring6mo_lag,
                   denDF_lin_m3_spring6mo_lag,
                   F_lin_m3_spring6mo_lag,
                   P_lin_m3_spring6mo_lag,
                   numDF_lin_m1_monsoon6mo_lag,
                   denDF_lin_m1_monsoon6mo_lag,
                   F_lin_m1_monsoon6mo_lag,
                   P_lin_m1_monsoon6mo_lag,
                   numDF_lin_m2_monsoon6mo_lag,
                   denDF_lin_m2_monsoon6mo_lag,
                   F_lin_m2_monsoon6mo_lag,
                   P_lin_m2_monsoon6mo_lag,
                   numDF_lin_m3_monsoon6mo_lag,
                   denDF_lin_m3_monsoon6mo_lag,
                   F_lin_m3_monsoon6mo_lag,
                   P_lin_m3_monsoon6mo_lag,
                   ParamEst_lin_m1_monsoon6mo,
                   ParamEst_lin_m1_monsoon6mo_lag,
                   ParamEst_lin_m1_spring6mo,
                   ParamEst_lin_m1_spring6mo_lag,
                   ParamEst_lin_m2_monsoon6mo,
                   ParamEst_lin_m2_monsoon6mo_lag,
                   ParamEst_lin_m2_spring6mo,
                   ParamEst_lin_m2_spring6mo_lag,
                   ParamEst_lin_m3_monsoon6mo,
                   ParamEst_lin_m3_monsoon6mo_lag,
                   ParamEst_lin_m3_spring6mo,
                   ParamEst_lin_m3_spring6mo_lag,
                   SE_lin_m1_monsoon6mo,
                   SE_lin_m1_monsoon6mo_lag,
                   SE_lin_m1_spring6mo,
                   SE_lin_m1_spring6mo_lag,
                   SE_lin_m2_monsoon6mo,
                   SE_lin_m2_monsoon6mo_lag,
                   SE_lin_m2_spring6mo,
                   SE_lin_m2_spring6mo_lag,
                   SE_lin_m3_monsoon6mo,
                   SE_lin_m3_monsoon6mo_lag,
                   SE_lin_m3_spring6mo,
                   SE_lin_m3_spring6mo_lag,
                   numDF_lin_m1_spring6mo_AR1,
                   denDF_lin_m1_spring6mo_AR1,
                   F_lin_m1_spring6mo_AR1,
                   P_lin_m1_spring6mo_AR1,
                   ParamEst_lin_m1_spring6mo_AR1,
                   SE_lin_m1_spring6mo_AR1,
                   numDF_lin_m1_spring6mo_AR2,
                   denDF_lin_m1_spring6mo_AR2,
                   F_lin_m1_spring6mo_AR2,
                   P_lin_m1_spring6mo_AR2,
                   ParamEst_lin_m1_spring6mo_AR2,
                   SE_lin_m1_spring6mo_AR2,
                   numDF_lin_m2_spring6mo_AR1,
                   denDF_lin_m2_spring6mo_AR1,
                   F_lin_m2_spring6mo_AR1,
                   P_lin_m2_spring6mo_AR1,
                   ParamEst_lin_m2_spring6mo_AR1,
                   SE_lin_m2_spring6mo_AR1,
                   numDF_lin_m2_spring6mo_AR2,
                   denDF_lin_m2_spring6mo_AR2,
                   F_lin_m2_spring6mo_AR2,
                   P_lin_m2_spring6mo_AR2,
                   ParamEst_lin_m2_spring6mo_AR2,
                   SE_lin_m2_spring6mo_AR2,
                   numDF_lin_m3_spring6mo_AR1,
                   denDF_lin_m3_spring6mo_AR1,
                   F_lin_m3_spring6mo_AR1,
                   P_lin_m3_spring6mo_AR1,
                   ParamEst_lin_m3_spring6mo_AR1,
                   SE_lin_m3_spring6mo_AR1,
                   numDF_lin_m3_spring6mo_AR2,
                   denDF_lin_m3_spring6mo_AR2,
                   F_lin_m3_spring6mo_AR2,
                   P_lin_m3_spring6mo_AR2,
                   ParamEst_lin_m3_spring6mo_AR2,
                   SE_lin_m3_spring6mo_AR2,
                   numDF_lin_m1_monsoon6mo_AR1,
                   denDF_lin_m1_monsoon6mo_AR1,
                   F_lin_m1_monsoon6mo_AR1,
                   P_lin_m1_monsoon6mo_AR1,
                   ParamEst_lin_m1_monsoon6mo_AR1,
                   SE_lin_m1_monsoon6mo_AR1,
                   numDF_lin_m1_monsoon6mo_AR2,
                   denDF_lin_m1_monsoon6mo_AR2,
                   F_lin_m1_monsoon6mo_AR2,
                   P_lin_m1_monsoon6mo_AR2,
                   ParamEst_lin_m1_monsoon6mo_AR2,
                   SE_lin_m1_monsoon6mo_AR2,
                   numDF_lin_m2_monsoon6mo_AR1,
                   denDF_lin_m2_monsoon6mo_AR1,
                   F_lin_m2_monsoon6mo_AR1,
                   P_lin_m2_monsoon6mo_AR1,
                   ParamEst_lin_m2_monsoon6mo_AR1,
                   SE_lin_m2_monsoon6mo_AR1,
                   numDF_lin_m2_monsoon6mo_AR2,
                   denDF_lin_m2_monsoon6mo_AR2,
                   F_lin_m2_monsoon6mo_AR2,
                   P_lin_m2_monsoon6mo_AR2,
                   ParamEst_lin_m2_monsoon6mo_AR2,
                   SE_lin_m2_monsoon6mo_AR2,
                   numDF_lin_m3_monsoon6mo_AR1,
                   denDF_lin_m3_monsoon6mo_AR1,
                   F_lin_m3_monsoon6mo_AR1,
                   P_lin_m3_monsoon6mo_AR1,
                   ParamEst_lin_m3_monsoon6mo_AR1,
                   SE_lin_m3_monsoon6mo_AR1,
                   numDF_lin_m3_monsoon6mo_AR2,
                   denDF_lin_m3_monsoon6mo_AR2,
                   F_lin_m3_monsoon6mo_AR2,
                   P_lin_m3_monsoon6mo_AR2,
                   ParamEst_lin_m3_monsoon6mo_AR2,
                   SE_lin_m3_monsoon6mo_AR2,
                   numDF_lin_m1_spring6mo_lag_AR1,
                   denDF_lin_m1_spring6mo_lag_AR1,
                   F_lin_m1_spring6mo_lag_AR1,
                   P_lin_m1_spring6mo_lag_AR1,
                   ParamEst_lin_m1_spring6mo_lag_AR1,
                   SE_lin_m1_spring6mo_lag_AR1,
                   numDF_lin_m1_spring6mo_lag_AR2,
                   denDF_lin_m1_spring6mo_lag_AR2,
                   F_lin_m1_spring6mo_lag_AR2,
                   P_lin_m1_spring6mo_lag_AR2,
                   ParamEst_lin_m1_spring6mo_lag_AR2,
                   SE_lin_m1_spring6mo_lag_AR2,
                   numDF_lin_m2_spring6mo_lag_AR1,
                   denDF_lin_m2_spring6mo_lag_AR1,
                   F_lin_m2_spring6mo_lag_AR1,
                   P_lin_m2_spring6mo_lag_AR1,
                   ParamEst_lin_m2_spring6mo_lag_AR1,
                   SE_lin_m2_spring6mo_lag_AR1,
                   numDF_lin_m2_spring6mo_lag_AR2,
                   denDF_lin_m2_spring6mo_lag_AR2,
                   F_lin_m2_spring6mo_lag_AR2,
                   P_lin_m2_spring6mo_lag_AR2,
                   ParamEst_lin_m2_spring6mo_lag_AR2,
                   SE_lin_m2_spring6mo_lag_AR2,
                   numDF_lin_m3_spring6mo_lag_AR1,
                   denDF_lin_m3_spring6mo_lag_AR1,
                   F_lin_m3_spring6mo_lag_AR1,
                   P_lin_m3_spring6mo_lag_AR1,
                   ParamEst_lin_m3_spring6mo_lag_AR1,
                   SE_lin_m3_spring6mo_lag_AR1,
                   numDF_lin_m3_spring6mo_lag_AR2,
                   denDF_lin_m3_spring6mo_lag_AR2,
                   F_lin_m3_spring6mo_lag_AR2,
                   P_lin_m3_spring6mo_lag_AR2,
                   ParamEst_lin_m3_spring6mo_lag_AR2,
                   SE_lin_m3_spring6mo_lag_AR2,
                   numDF_lin_m1_monsoon6mo_lag_AR1,
                   denDF_lin_m1_monsoon6mo_lag_AR1,
                   F_lin_m1_monsoon6mo_lag_AR1,
                   P_lin_m1_monsoon6mo_lag_AR1,
                   ParamEst_lin_m1_monsoon6mo_lag_AR1,
                   SE_lin_m1_monsoon6mo_lag_AR1,
                   numDF_lin_m1_monsoon6mo_lag_AR2,
                   denDF_lin_m1_monsoon6mo_lag_AR2,
                   F_lin_m1_monsoon6mo_lag_AR2,
                   P_lin_m1_monsoon6mo_lag_AR2,
                   ParamEst_lin_m1_monsoon6mo_lag_AR2,
                   SE_lin_m1_monsoon6mo_lag_AR2,
                   numDF_lin_m2_monsoon6mo_lag_AR1,
                   denDF_lin_m2_monsoon6mo_lag_AR1,
                   F_lin_m2_monsoon6mo_lag_AR1,
                   P_lin_m2_monsoon6mo_lag_AR1,
                   ParamEst_lin_m2_monsoon6mo_lag_AR1,
                   SE_lin_m2_monsoon6mo_lag_AR1,
                   numDF_lin_m2_monsoon6mo_lag_AR2,
                   denDF_lin_m2_monsoon6mo_lag_AR2,
                   F_lin_m2_monsoon6mo_lag_AR2,
                   P_lin_m2_monsoon6mo_lag_AR2,
                   ParamEst_lin_m2_monsoon6mo_lag_AR2,
                   SE_lin_m2_monsoon6mo_lag_AR2,
                   numDF_lin_m3_monsoon6mo_lag_AR1,
                   denDF_lin_m3_monsoon6mo_lag_AR1,
                   F_lin_m3_monsoon6mo_lag_AR1,
                   P_lin_m3_monsoon6mo_lag_AR1,
                   ParamEst_lin_m3_monsoon6mo_lag_AR1,
                   SE_lin_m3_monsoon6mo_lag_AR1,
                   numDF_lin_m3_monsoon6mo_lag_AR2,
                   denDF_lin_m3_monsoon6mo_lag_AR2,
                   F_lin_m3_monsoon6mo_lag_AR2,
                   P_lin_m3_monsoon6mo_lag_AR2,
                   ParamEst_lin_m3_monsoon6mo_lag_AR2,
                   SE_lin_m3_monsoon6mo_lag_AR2,
                   Rsquared_marginal_m_null,
                   Rsquared_marginal_m1_spring6mo,
                   Rsquared_marginal_m2_spring6mo,
                   Rsquared_marginal_m3_spring6mo,
                   Rsquared_marginal_m1_spring6mo_AR1,
                   Rsquared_marginal_m1_spring6mo_AR2,
                   Rsquared_marginal_m2_spring6mo_AR1,
                   Rsquared_marginal_m2_spring6mo_AR2,
                   Rsquared_marginal_m3_spring6mo_AR1,
                   Rsquared_marginal_m3_spring6mo_AR2,
                   Rsquared_marginal_m1_monsoon6mo,
                   Rsquared_marginal_m2_monsoon6mo,
                   Rsquared_marginal_m3_monsoon6mo,
                   Rsquared_marginal_m1_monsoon6mo_AR1,
                   Rsquared_marginal_m1_monsoon6mo_AR2,
                   Rsquared_marginal_m2_monsoon6mo_AR1,
                   Rsquared_marginal_m2_monsoon6mo_AR2,
                   Rsquared_marginal_m3_monsoon6mo_AR1,
                   Rsquared_marginal_m3_monsoon6mo_AR2,
                   Rsquared_marginal_m1_spring6mo_lag,
                   Rsquared_marginal_m2_spring6mo_lag,
                   Rsquared_marginal_m3_spring6mo_lag,
                   Rsquared_marginal_m1_spring6mo_lag_AR1,
                   Rsquared_marginal_m1_spring6mo_lag_AR2,
                   Rsquared_marginal_m2_spring6mo_lag_AR1,
                   Rsquared_marginal_m2_spring6mo_lag_AR2,
                   Rsquared_marginal_m3_spring6mo_lag_AR1,
                   Rsquared_marginal_m3_spring6mo_lag_AR2,
                   Rsquared_marginal_m1_monsoon6mo_lag,
                   Rsquared_marginal_m2_monsoon6mo_lag,
                   Rsquared_marginal_m3_monsoon6mo_lag,
                   Rsquared_marginal_m1_monsoon6mo_lag_AR1,
                   Rsquared_marginal_m1_monsoon6mo_lag_AR2,
                   Rsquared_marginal_m2_monsoon6mo_lag_AR1,
                   Rsquared_marginal_m2_monsoon6mo_lag_AR2,
                   Rsquared_marginal_m3_monsoon6mo_lag_AR1,
                   Rsquared_marginal_m3_monsoon6mo_lag_AR2,
                   Rsquared_conditional_m_null,
                   Rsquared_conditional_m1_spring6mo,
                   Rsquared_conditional_m2_spring6mo,
                   Rsquared_conditional_m3_spring6mo,
                   Rsquared_conditional_m1_spring6mo_AR1,
                   Rsquared_conditional_m1_spring6mo_AR2,
                   Rsquared_conditional_m2_spring6mo_AR1,
                   Rsquared_conditional_m2_spring6mo_AR2,
                   Rsquared_conditional_m3_spring6mo_AR1,
                   Rsquared_conditional_m3_spring6mo_AR2,
                   Rsquared_conditional_m1_monsoon6mo,
                   Rsquared_conditional_m2_monsoon6mo,
                   Rsquared_conditional_m3_monsoon6mo,
                   Rsquared_conditional_m1_monsoon6mo_AR1,
                   Rsquared_conditional_m1_monsoon6mo_AR2,
                   Rsquared_conditional_m2_monsoon6mo_AR1,
                   Rsquared_conditional_m2_monsoon6mo_AR2,
                   Rsquared_conditional_m3_monsoon6mo_AR1,
                   Rsquared_conditional_m3_monsoon6mo_AR2,
                   Rsquared_conditional_m1_spring6mo_lag,
                   Rsquared_conditional_m2_spring6mo_lag,
                   Rsquared_conditional_m3_spring6mo_lag,
                   Rsquared_conditional_m1_spring6mo_lag_AR1,
                   Rsquared_conditional_m1_spring6mo_lag_AR2,
                   Rsquared_conditional_m2_spring6mo_lag_AR1,
                   Rsquared_conditional_m2_spring6mo_lag_AR2,
                   Rsquared_conditional_m3_spring6mo_lag_AR1,
                   Rsquared_conditional_m3_spring6mo_lag_AR2,
                   Rsquared_conditional_m1_monsoon6mo_lag,
                   Rsquared_conditional_m2_monsoon6mo_lag,
                   Rsquared_conditional_m3_monsoon6mo_lag,
                   Rsquared_conditional_m1_monsoon6mo_lag_AR1,
                   Rsquared_conditional_m1_monsoon6mo_lag_AR2,
                   Rsquared_conditional_m2_monsoon6mo_lag_AR1,
                   Rsquared_conditional_m2_monsoon6mo_lag_AR2,
                   Rsquared_conditional_m3_monsoon6mo_lag_AR1,
                   Rsquared_conditional_m3_monsoon6mo_lag_AR2,
                   numDF_quad_m2_spring6mo,
                   numDF_quad_m2_spring6mo_AR1,
                   numDF_quad_m2_spring6mo_AR2,
                   numDF_quad_m2_monsoon6mo,
                   numDF_quad_m2_monsoon6mo_AR1,
                   numDF_quad_m2_monsoon6mo_AR2,
                   numDF_quad_m2_spring6mo_lag,
                   numDF_quad_m2_spring6mo_lag_AR1,
                   numDF_quad_m2_spring6mo_lag_AR2,
                   numDF_quad_m2_monsoon6mo_lag,
                   numDF_quad_m2_monsoon6mo_lag_AR1,
                   numDF_quad_m2_monsoon6mo_lag_AR2,
                   numDF_quad_m3_spring6mo,
                   numDF_quad_m3_spring6mo_AR1,
                   numDF_quad_m3_spring6mo_AR2,
                   numDF_quad_m3_monsoon6mo,
                   numDF_quad_m3_monsoon6mo_AR1,
                   numDF_quad_m3_monsoon6mo_AR2,
                   numDF_quad_m3_spring6mo_lag,
                   numDF_quad_m3_spring6mo_lag_AR1,
                   numDF_quad_m3_spring6mo_lag_AR2,
                   numDF_quad_m3_monsoon6mo_lag,
                   numDF_quad_m3_monsoon6mo_lag_AR1,
                   numDF_quad_m3_monsoon6mo_lag_AR2,
                   denDF_quad_m2_spring6mo,
                   denDF_quad_m2_spring6mo_AR1,
                   denDF_quad_m2_spring6mo_AR2,
                   denDF_quad_m2_monsoon6mo,
                   denDF_quad_m2_monsoon6mo_AR1,
                   denDF_quad_m2_monsoon6mo_AR2,
                   denDF_quad_m2_spring6mo_lag,
                   denDF_quad_m2_spring6mo_lag_AR1,
                   denDF_quad_m2_spring6mo_lag_AR2,
                   denDF_quad_m2_monsoon6mo_lag,
                   denDF_quad_m2_monsoon6mo_lag_AR1,
                   denDF_quad_m2_monsoon6mo_lag_AR2,
                   denDF_quad_m3_spring6mo,
                   denDF_quad_m3_spring6mo_AR1,
                   denDF_quad_m3_spring6mo_AR2,
                   denDF_quad_m3_monsoon6mo,
                   denDF_quad_m3_monsoon6mo_AR1,
                   denDF_quad_m3_monsoon6mo_AR2,
                   denDF_quad_m3_spring6mo_lag,
                   denDF_quad_m3_spring6mo_lag_AR1,
                   denDF_quad_m3_spring6mo_lag_AR2,
                   denDF_quad_m3_monsoon6mo_lag,
                   denDF_quad_m3_monsoon6mo_lag_AR1,
                   denDF_quad_m3_monsoon6mo_lag_AR2,
                   F_quad_m2_spring6mo,
                   F_quad_m2_spring6mo_AR1,
                   F_quad_m2_spring6mo_AR2,
                   F_quad_m2_monsoon6mo,
                   F_quad_m2_monsoon6mo_AR1,
                   F_quad_m2_monsoon6mo_AR2,
                   F_quad_m2_spring6mo_lag,
                   F_quad_m2_spring6mo_lag_AR1,
                   F_quad_m2_spring6mo_lag_AR2,
                   F_quad_m2_monsoon6mo_lag,
                   F_quad_m2_monsoon6mo_lag_AR1,
                   F_quad_m2_monsoon6mo_lag_AR2,
                   F_quad_m3_spring6mo,
                   F_quad_m3_spring6mo_AR1,
                   F_quad_m3_spring6mo_AR2,
                   F_quad_m3_monsoon6mo,
                   F_quad_m3_monsoon6mo_AR1,
                   F_quad_m3_monsoon6mo_AR2,
                   F_quad_m3_spring6mo_lag,
                   F_quad_m3_spring6mo_lag_AR1,
                   F_quad_m3_spring6mo_lag_AR2,
                   F_quad_m3_monsoon6mo_lag,
                   F_quad_m3_monsoon6mo_lag_AR1,
                   F_quad_m3_monsoon6mo_lag_AR2,
                   P_quad_m2_spring6mo,
                   P_quad_m2_spring6mo_AR1,
                   P_quad_m2_spring6mo_AR2,
                   P_quad_m2_monsoon6mo,
                   P_quad_m2_monsoon6mo_AR1,
                   P_quad_m2_monsoon6mo_AR2,
                   P_quad_m2_spring6mo_lag,
                   P_quad_m2_spring6mo_lag_AR1,
                   P_quad_m2_spring6mo_lag_AR2,
                   P_quad_m2_monsoon6mo_lag,
                   P_quad_m2_monsoon6mo_lag_AR1,
                   P_quad_m2_monsoon6mo_lag_AR2,
                   P_quad_m3_spring6mo,
                   P_quad_m3_spring6mo_AR1,
                   P_quad_m3_spring6mo_AR2,
                   P_quad_m3_monsoon6mo,
                   P_quad_m3_monsoon6mo_AR1,
                   P_quad_m3_monsoon6mo_AR2,
                   P_quad_m3_spring6mo_lag,
                   P_quad_m3_spring6mo_lag_AR1,
                   P_quad_m3_spring6mo_lag_AR2,
                   P_quad_m3_monsoon6mo_lag,
                   P_quad_m3_monsoon6mo_lag_AR1,
                   P_quad_m3_monsoon6mo_lag_AR2,
                   ParamEst_quad_m2_spring6mo,
                   ParamEst_quad_m2_spring6mo_AR1,
                   ParamEst_quad_m2_spring6mo_AR2,
                   ParamEst_quad_m2_monsoon6mo,
                   ParamEst_quad_m2_monsoon6mo_AR1,
                   ParamEst_quad_m2_monsoon6mo_AR2,
                   ParamEst_quad_m2_spring6mo_lag,
                   ParamEst_quad_m2_spring6mo_lag_AR1,
                   ParamEst_quad_m2_spring6mo_lag_AR2,
                   ParamEst_quad_m2_monsoon6mo_lag,
                   ParamEst_quad_m2_monsoon6mo_lag_AR1,
                   ParamEst_quad_m2_monsoon6mo_lag_AR2,
                   ParamEst_quad_m3_spring6mo,
                   ParamEst_quad_m3_spring6mo_AR1,
                   ParamEst_quad_m3_spring6mo_AR2,
                   ParamEst_quad_m3_monsoon6mo,
                   ParamEst_quad_m3_monsoon6mo_AR1,
                   ParamEst_quad_m3_monsoon6mo_AR2,
                   ParamEst_quad_m3_spring6mo_lag,
                   ParamEst_quad_m3_spring6mo_lag_AR1,
                   ParamEst_quad_m3_spring6mo_lag_AR2,
                   ParamEst_quad_m3_monsoon6mo_lag,
                   ParamEst_quad_m3_monsoon6mo_lag_AR1,
                   ParamEst_quad_m3_monsoon6mo_lag_AR2,
                   SE_quad_m2_spring6mo,
                   SE_quad_m2_spring6mo_AR1,
                   SE_quad_m2_spring6mo_AR2,
                   SE_quad_m2_monsoon6mo,
                   SE_quad_m2_monsoon6mo_AR1,
                   SE_quad_m2_monsoon6mo_AR2,
                   SE_quad_m2_spring6mo_lag,
                   SE_quad_m2_spring6mo_lag_AR1,
                   SE_quad_m2_spring6mo_lag_AR2,
                   SE_quad_m2_monsoon6mo_lag,
                   SE_quad_m2_monsoon6mo_lag_AR1,
                   SE_quad_m2_monsoon6mo_lag_AR2,
                   SE_quad_m3_spring6mo,
                   SE_quad_m3_spring6mo_AR1,
                   SE_quad_m3_spring6mo_AR2,
                   SE_quad_m3_monsoon6mo,
                   SE_quad_m3_monsoon6mo_AR1,
                   SE_quad_m3_monsoon6mo_AR2,
                   SE_quad_m3_spring6mo_lag,
                   SE_quad_m3_spring6mo_lag_AR1,
                   SE_quad_m3_spring6mo_lag_AR2,
                   SE_quad_m3_monsoon6mo_lag,
                   SE_quad_m3_monsoon6mo_lag_AR1,
                   SE_quad_m3_monsoon6mo_lag_AR2,
                   numDF_cub_m3_spring6mo,
                   numDF_cub_m3_spring6mo_AR1,
                   numDF_cub_m3_spring6mo_AR2,
                   numDF_cub_m3_monsoon6mo,
                   numDF_cub_m3_monsoon6mo_AR1,
                   numDF_cub_m3_monsoon6mo_AR2,
                   numDF_cub_m3_spring6mo_lag,
                   numDF_cub_m3_spring6mo_lag_AR1,
                   numDF_cub_m3_spring6mo_lag_AR2,
                   numDF_cub_m3_monsoon6mo_lag,
                   numDF_cub_m3_monsoon6mo_lag_AR1,
                   numDF_cub_m3_monsoon6mo_lag_AR2,
                   denDF_cub_m3_spring6mo,
                   denDF_cub_m3_spring6mo_AR1,
                   denDF_cub_m3_spring6mo_AR2,
                   denDF_cub_m3_monsoon6mo,
                   denDF_cub_m3_monsoon6mo_AR1,
                   denDF_cub_m3_monsoon6mo_AR2,
                   denDF_cub_m3_spring6mo_lag,
                   denDF_cub_m3_spring6mo_lag_AR1,
                   denDF_cub_m3_spring6mo_lag_AR2,
                   denDF_cub_m3_monsoon6mo_lag,
                   denDF_cub_m3_monsoon6mo_lag_AR1,
                   denDF_cub_m3_monsoon6mo_lag_AR2,
                   F_cub_m3_spring6mo,
                   F_cub_m3_spring6mo_AR1,
                   F_cub_m3_spring6mo_AR2,
                   F_cub_m3_monsoon6mo,
                   F_cub_m3_monsoon6mo_AR1,
                   F_cub_m3_monsoon6mo_AR2,
                   F_cub_m3_spring6mo_lag,
                   F_cub_m3_spring6mo_lag_AR1,
                   F_cub_m3_spring6mo_lag_AR2,
                   F_cub_m3_monsoon6mo_lag,
                   F_cub_m3_monsoon6mo_lag_AR1,
                   F_cub_m3_monsoon6mo_lag_AR2,
                   P_cub_m3_spring6mo,
                   P_cub_m3_spring6mo_AR1,
                   P_cub_m3_spring6mo_AR2,
                   P_cub_m3_monsoon6mo,
                   P_cub_m3_monsoon6mo_AR1,
                   P_cub_m3_monsoon6mo_AR2,
                   P_cub_m3_spring6mo_lag,
                   P_cub_m3_spring6mo_lag_AR1,
                   P_cub_m3_spring6mo_lag_AR2,
                   P_cub_m3_monsoon6mo_lag,
                   P_cub_m3_monsoon6mo_lag_AR1,
                   P_cub_m3_monsoon6mo_lag_AR2,
                   ParamEst_cub_m3_spring6mo,
                   ParamEst_cub_m3_spring6mo_AR1,
                   ParamEst_cub_m3_spring6mo_AR2,
                   ParamEst_cub_m3_monsoon6mo,
                   ParamEst_cub_m3_monsoon6mo_AR1,
                   ParamEst_cub_m3_monsoon6mo_AR2,
                   ParamEst_cub_m3_spring6mo_lag,
                   ParamEst_cub_m3_spring6mo_lag_AR1,
                   ParamEst_cub_m3_spring6mo_lag_AR2,
                   ParamEst_cub_m3_monsoon6mo_lag,
                   ParamEst_cub_m3_monsoon6mo_lag_AR1,
                   ParamEst_cub_m3_monsoon6mo_lag_AR2,
                   SE_cub_m3_spring6mo,
                   SE_cub_m3_spring6mo_AR1,
                   SE_cub_m3_spring6mo_AR2,
                   SE_cub_m3_monsoon6mo,
                   SE_cub_m3_monsoon6mo_AR1,
                   SE_cub_m3_monsoon6mo_AR2,
                   SE_cub_m3_spring6mo_lag,
                   SE_cub_m3_spring6mo_lag_AR1,
                   SE_cub_m3_spring6mo_lag_AR2,
                   SE_cub_m3_monsoon6mo_lag,
                   SE_cub_m3_monsoon6mo_lag_AR1,
                   SE_cub_m3_monsoon6mo_lag_AR2)
  
  # append results for the bee species to the output dataset
  beeCSF_output<-rbind(beeCSF_output,output_id)
  
}


# Write .csv file of output
beeCSF_output <- data.frame(beeCSF_output) # make matrix into data frame
beeCSF_output2 <- beeCSF_output[-1,] # remove empty first row

#write.csv(beeCSF_output2, "bee_CSFs_2023-08-29.csv",row.names=FALSE)
