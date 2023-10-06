################################################################################### 
# Predicting future bee abundances based on historic climate sensitivities and projected future aridity values
# Using projected future climate data from the CanESM2 GCM

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


##### Chihuahuan Desert Shrubland #####

# Read in data frame of historic and projected future SPEI data
climate <- read.csv("spei_historic_and_future_byscenario_CanESM2_2023-08-29.csv")

# Create data frame just containing projected future data (2021-2100)
climate_forecast<-subset(climate,source!="SEV-met" & station=="Five Points")

# Create a data frame of just 2020 data
clim_2020<-subset(climate,year==2020 & station=="Five Points")

# Create a copy of the 2020 data to be associated with the RCP 2.6 data series
clim_2020_2.6 <- clim_2020
clim_2020_2.6$source <- "rcp2.6"

# Create a copy of the 2020 data to be associated with the RCP 4.5 data series
clim_2020_4.5 <- clim_2020
clim_2020_4.5$source <- "rcp4.5"

# Create a copy of the 2020 data to be associated with the RCP 8.5 data series
clim_2020_8.5 <- clim_2020
clim_2020_8.5$source <- "rcp8.5"

# Create a merged climate data frame for 2020-2100
climate_forecast<-bind_rows(climate_forecast,clim_2020_2.6,clim_2020_4.5,clim_2020_8.5)

# Rename column
names(climate_forecast)[2]<-"scenario"

# Create different data frames for each transect, and add transect labels
climate_forecast_c1<-climate_forecast
climate_forecast_c1$transect <- "C1"

climate_forecast_c2<-climate_forecast
climate_forecast_c2$transect <- "C2"

climate_forecast_c3<-climate_forecast
climate_forecast_c3$transect <- "C3"

climate_forecast_c4<-climate_forecast
climate_forecast_c4$transect <- "C4"

climate_forecast_c5<-climate_forecast
climate_forecast_c5$transect <- "C5"

# Combine the transect data frames into a final data frame
climate_forecast<-bind_rows(climate_forecast_c1,climate_forecast_c2,climate_forecast_c3,climate_forecast_c4,climate_forecast_c5)

# Create matrix to hold results (predicted future abundance values)
beeCSF_output<-matrix(nrow=1,ncol=9,byrow=TRUE,dimnames=list(c("row1"),c("code","ecosystem","transect","predicted_max_abundance_per_transect","scenario","year","model_id","spei_integ","spei_value")))

# Create an ecosystem code object
ecosystem <-"C"

# Create a data frame of just the species abundance matrix
speciesMatrixC <- creo[,9:229] 

# Create vector of species codes
speciesCodesC <- colnames(speciesMatrixC)

# Create vector of transect names
transectsC <- c("C1","C2","C3","C4","C5")

# Treat SPEI variables as numeric
climate_forecast$monsoon6SPEI<-as.numeric(climate_forecast$monsoon6SPEI)
climate_forecast$spring6SPEI<-as.numeric(climate_forecast$spring6SPEI)
climate_forecast$spring6SPEI_prioryear<-as.numeric(climate_forecast$spring6SPEI_prioryear)

creo$spring6SPEI<-as.numeric(creo$spring6SPEI)
creo$spring6SPEI_prioryear<-as.numeric(creo$spring6SPEI_prioryear)


### Loop: if a species is found in a given ecosystem, run CSFs, predict future bee abundances under different climate scenarios, and put the results in the beeCSF_output matrix
for (i in 1:length(speciesCodesC)) {
  
  # save the species code for column i
  speciesCode <- speciesCodesC[i]
  
  ### build CSFs
  
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
  
  # create list of models
  model_list <- list(m_null,m1_spring6mo,m2_spring6mo,m3_spring6mo,m1_spring6mo_AR1,m1_spring6mo_AR2,m2_spring6mo_AR1,m2_spring6mo_AR2,m3_spring6mo_AR1,m3_spring6mo_AR2,m1_monsoon6mo,m2_monsoon6mo,m3_monsoon6mo,m1_monsoon6mo_AR1,m1_monsoon6mo_AR2,m2_monsoon6mo_AR1,m2_monsoon6mo_AR2,m3_monsoon6mo_AR1,m3_monsoon6mo_AR2,m1_spring6mo_lag,m2_spring6mo_lag,m3_spring6mo_lag,m1_spring6mo_lag_AR1,m1_spring6mo_lag_AR2,m2_spring6mo_lag_AR1,m2_spring6mo_lag_AR2,m3_spring6mo_lag_AR1,m3_spring6mo_lag_AR2,m1_monsoon6mo_lag,m2_monsoon6mo_lag,m3_monsoon6mo_lag,m1_monsoon6mo_lag_AR1,m1_monsoon6mo_lag_AR2,m2_monsoon6mo_lag_AR1,m2_monsoon6mo_lag_AR2,m3_monsoon6mo_lag_AR1,m3_monsoon6mo_lag_AR2)
  
  # create data frame of model names and AICc values
  model_names <- c("m_null","m1_spring6mo","m2_spring6mo","m3_spring6mo","m1_spring6mo_AR1","m1_spring6mo_AR2","m2_spring6mo_AR1","m2_spring6mo_AR2","m3_spring6mo_AR1","m3_spring6mo_AR2","m1_monsoon6mo","m2_monsoon6mo","m3_monsoon6mo","m1_monsoon6mo_AR1","m1_monsoon6mo_AR2","m2_monsoon6mo_AR1","m2_monsoon6mo_AR2","m3_monsoon6mo_AR1","m3_monsoon6mo_AR2","m1_spring6mo_lag","m2_spring6mo_lag","m3_spring6mo_lag","m1_spring6mo_lag_AR1","m1_spring6mo_lag_AR2","m2_spring6mo_lag_AR1","m2_spring6mo_lag_AR2","m3_spring6mo_lag_AR1","m3_spring6mo_lag_AR2","m1_monsoon6mo_lag","m2_monsoon6mo_lag","m3_monsoon6mo_lag","m1_monsoon6mo_lag_AR1","m1_monsoon6mo_lag_AR2","m2_monsoon6mo_lag_AR1","m2_monsoon6mo_lag_AR2","m3_monsoon6mo_lag_AR1","m3_monsoon6mo_lag_AR2")
  
  model_aicc<-c(AICc(m_null),AICc(m1_spring6mo),AICc(m2_spring6mo),AICc(m3_spring6mo),AICc(m1_spring6mo_AR1),AICc(m1_spring6mo_AR2),AICc(m2_spring6mo_AR1),AICc(m2_spring6mo_AR2),AICc(m3_spring6mo_AR1),AICc(m3_spring6mo_AR2),AICc(m1_monsoon6mo),AICc(m2_monsoon6mo),AICc(m3_monsoon6mo),AICc(m1_monsoon6mo_AR1),AICc(m1_monsoon6mo_AR2),AICc(m2_monsoon6mo_AR1),AICc(m2_monsoon6mo_AR2),AICc(m3_monsoon6mo_AR1),AICc(m3_monsoon6mo_AR2),AICc(m1_spring6mo_lag),AICc(m2_spring6mo_lag),AICc(m3_spring6mo_lag),AICc(m1_spring6mo_lag_AR1),AICc(m1_spring6mo_lag_AR2),AICc(m2_spring6mo_lag_AR1),AICc(m2_spring6mo_lag_AR2),AICc(m3_spring6mo_lag_AR1),AICc(m3_spring6mo_lag_AR2),AICc(m1_monsoon6mo_lag),AICc(m2_monsoon6mo_lag),AICc(m3_monsoon6mo_lag),AICc(m1_monsoon6mo_lag_AR1),AICc(m1_monsoon6mo_lag_AR2),AICc(m2_monsoon6mo_lag_AR1),AICc(m2_monsoon6mo_lag_AR2),AICc(m3_monsoon6mo_lag_AR1),AICc(m3_monsoon6mo_lag_AR2))
  
  model_df<-data.frame(model_names,model_aicc)
  
  
  # calculate minimum AICc value
  min_aicc<-min(model_aicc)
  
  # create vectors of scenarios and years
  scenarios <- levels(as.factor(climate_forecast$scenario))
  years <- levels(as.factor(climate_forecast$year))
  
  # find the best CSF for the bee species, and use it to predict future abundance in different years under different climate scenarios
  for(j in 1:length(model_aicc)) {
    for(k in 1:length(scenarios)) {
      for (l in 1:length(years)) {
        for(m in 1:length(transectsC)) {
          rcp<-subset(climate_forecast,scenario==scenarios[k] & year == years[l] & transect == transectsC[m])
          if (model_aicc[j]==min_aicc) {
            best_model<-subset(model_df,model_aicc==min_aicc)
            if (best_model$model_names != "m_null") {
              model<-model_list[[j]]
              predicted_abundance_sqrt<-predict(model,rcp)[1]
              predicted_abundance<-predicted_abundance_sqrt^2
              scenario <- scenarios[k]
              year <- years[l]
              transect <- transectsC[m]
              
              model_id<-best_model[1,1]
              spei_integ<-names(model_list[[j]][4]$coefficients$fixed[2])
              spei_value<-rcp[,spei_integ]
              
              #bind relevant data columns together
              output_id<-cbind(speciesCode, ecosystem, transect, predicted_abundance, scenario, year,model_id,spei_integ,spei_value)
              
              #append new results to output dataset
              beeCSF_output<-rbind(beeCSF_output,output_id)
            }
          }
        }
      }
    }
  }
  
}



##### Chihuahuan Desert Grassland #####

# Read in data frame of historic and projected future SPEI data
climate <- read.csv("spei_historic_and_future_byscenario_CanESM2_2023-08-29.csv")

# Create data frame just containing projected future data (2021-2100)
climate_forecast<-subset(climate,source!="SEV-met" & station=="Five Points")

# Create a data frame of just 2020 data
clim_2020<-subset(climate,year==2020 & station=="Five Points")

# Create a copy of the 2020 data to be associated with the RCP 2.6 data series
clim_2020_2.6 <- clim_2020
clim_2020_2.6$source <- "rcp2.6"

# Create a copy of the 2020 data to be associated with the RCP 4.5 data series
clim_2020_4.5 <- clim_2020
clim_2020_4.5$source <- "rcp4.5"

# Create a copy of the 2020 data to be associated with the RCP 8.5 data series
clim_2020_8.5 <- clim_2020
clim_2020_8.5$source <- "rcp8.5"

# Create a merged climate data frame for 2020-2100
climate_forecast<-bind_rows(climate_forecast,clim_2020_2.6,clim_2020_4.5,clim_2020_8.5)

# Rename column
names(climate_forecast)[2]<-"scenario"

# Create different data frames for each transect, and add transect labels
climate_forecast_g1<-climate_forecast
climate_forecast_g1$transect <- "G1"

climate_forecast_g2<-climate_forecast
climate_forecast_g2$transect <- "G2"

climate_forecast_g3<-climate_forecast
climate_forecast_g3$transect <- "G3"

climate_forecast_g4<-climate_forecast
climate_forecast_g4$transect <- "G4"

climate_forecast_g5<-climate_forecast
climate_forecast_g5$transect <- "G5"

# Combine the transect data frames into a final data frame
climate_forecast<-bind_rows(climate_forecast_g1,climate_forecast_g2,climate_forecast_g3,climate_forecast_g4,climate_forecast_g5)

# Create an ecosystem code object
ecosystem <-"G"

# Create a data frame of just the species abundance matrix
speciesMatrixG <- black[,9:224] 

# Create vector of species codes
speciesCodesG <- colnames(speciesMatrixG)

# Create vector of transect names
transectsG <- c("G1","G2","G3","G4","G5")

# Treat SPEI variables as numeric
climate_forecast$monsoon6SPEI<-as.numeric(climate_forecast$monsoon6SPEI)
climate_forecast$spring6SPEI<-as.numeric(climate_forecast$spring6SPEI)
climate_forecast$spring6SPEI_prioryear<-as.numeric(climate_forecast$spring6SPEI_prioryear)

### Loop: if a species is found in a given ecosystem, run CSFs, predict future bee abundances under different climate scenarios, and put the results in the beeCSF_output matrix
for (i in 1:length(speciesCodesG)) {
  
  # save the species code for column i
  speciesCode <- speciesCodesG[i]
  
  ### build CSFs
  
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
  
  # create list of models
  model_list <- list(m_null,m1_spring6mo,m2_spring6mo,m3_spring6mo,m1_spring6mo_AR1,m1_spring6mo_AR2,m2_spring6mo_AR1,m2_spring6mo_AR2,m3_spring6mo_AR1,m3_spring6mo_AR2,m1_monsoon6mo,m2_monsoon6mo,m3_monsoon6mo,m1_monsoon6mo_AR1,m1_monsoon6mo_AR2,m2_monsoon6mo_AR1,m2_monsoon6mo_AR2,m3_monsoon6mo_AR1,m3_monsoon6mo_AR2,m1_spring6mo_lag,m2_spring6mo_lag,m3_spring6mo_lag,m1_spring6mo_lag_AR1,m1_spring6mo_lag_AR2,m2_spring6mo_lag_AR1,m2_spring6mo_lag_AR2,m3_spring6mo_lag_AR1,m3_spring6mo_lag_AR2,m1_monsoon6mo_lag,m2_monsoon6mo_lag,m3_monsoon6mo_lag,m1_monsoon6mo_lag_AR1,m1_monsoon6mo_lag_AR2,m2_monsoon6mo_lag_AR1,m2_monsoon6mo_lag_AR2,m3_monsoon6mo_lag_AR1,m3_monsoon6mo_lag_AR2)
  
  # create data frame of model names and AICc values
  model_names <- c("m_null","m1_spring6mo","m2_spring6mo","m3_spring6mo","m1_spring6mo_AR1","m1_spring6mo_AR2","m2_spring6mo_AR1","m2_spring6mo_AR2","m3_spring6mo_AR1","m3_spring6mo_AR2","m1_monsoon6mo","m2_monsoon6mo","m3_monsoon6mo","m1_monsoon6mo_AR1","m1_monsoon6mo_AR2","m2_monsoon6mo_AR1","m2_monsoon6mo_AR2","m3_monsoon6mo_AR1","m3_monsoon6mo_AR2","m1_spring6mo_lag","m2_spring6mo_lag","m3_spring6mo_lag","m1_spring6mo_lag_AR1","m1_spring6mo_lag_AR2","m2_spring6mo_lag_AR1","m2_spring6mo_lag_AR2","m3_spring6mo_lag_AR1","m3_spring6mo_lag_AR2","m1_monsoon6mo_lag","m2_monsoon6mo_lag","m3_monsoon6mo_lag","m1_monsoon6mo_lag_AR1","m1_monsoon6mo_lag_AR2","m2_monsoon6mo_lag_AR1","m2_monsoon6mo_lag_AR2","m3_monsoon6mo_lag_AR1","m3_monsoon6mo_lag_AR2")
  
  model_aicc<-c(AICc(m_null),AICc(m1_spring6mo),AICc(m2_spring6mo),AICc(m3_spring6mo),AICc(m1_spring6mo_AR1),AICc(m1_spring6mo_AR2),AICc(m2_spring6mo_AR1),AICc(m2_spring6mo_AR2),AICc(m3_spring6mo_AR1),AICc(m3_spring6mo_AR2),AICc(m1_monsoon6mo),AICc(m2_monsoon6mo),AICc(m3_monsoon6mo),AICc(m1_monsoon6mo_AR1),AICc(m1_monsoon6mo_AR2),AICc(m2_monsoon6mo_AR1),AICc(m2_monsoon6mo_AR2),AICc(m3_monsoon6mo_AR1),AICc(m3_monsoon6mo_AR2),AICc(m1_spring6mo_lag),AICc(m2_spring6mo_lag),AICc(m3_spring6mo_lag),AICc(m1_spring6mo_lag_AR1),AICc(m1_spring6mo_lag_AR2),AICc(m2_spring6mo_lag_AR1),AICc(m2_spring6mo_lag_AR2),AICc(m3_spring6mo_lag_AR1),AICc(m3_spring6mo_lag_AR2),AICc(m1_monsoon6mo_lag),AICc(m2_monsoon6mo_lag),AICc(m3_monsoon6mo_lag),AICc(m1_monsoon6mo_lag_AR1),AICc(m1_monsoon6mo_lag_AR2),AICc(m2_monsoon6mo_lag_AR1),AICc(m2_monsoon6mo_lag_AR2),AICc(m3_monsoon6mo_lag_AR1),AICc(m3_monsoon6mo_lag_AR2))
  
  model_df<-data.frame(model_names,model_aicc)
  
  # calculate minimum AICc value
  min_aicc<-min(model_aicc)
  
  # create vectors of scenarios and years
  scenarios <- levels(as.factor(climate_forecast$scenario))
  years <- levels(as.factor(climate_forecast$year))
  
  # find the best CSF for the bee species, and use it to predict future abundance in different years under different climate scenarios
  for(j in 1:length(model_aicc)) {
    for(k in 1:length(scenarios)) {
      for (l in 1:length(years)) {
        for(m in 1:length(transectsG)) {
          rcp<-subset(climate_forecast,scenario==scenarios[k] & year == years[l] & transect == transectsG[m])
          if (model_aicc[j]==min_aicc) {
            best_model<-subset(model_df,model_aicc==min_aicc)
            if (best_model$model_names != "m_null") {
              model<-model_list[[j]]
              predicted_abundance_sqrt<-predict(model,rcp)[1]
              predicted_abundance<-predicted_abundance_sqrt^2
              scenario <- scenarios[k]
              year <- years[l]
              transect <- transectsG[m]
              
              model_id<-best_model[1,1]
              spei_integ<-names(model_list[[j]][4]$coefficients$fixed[2])
              spei_value<-rcp[,spei_integ]
              
              #bind relevant data columns together
              output_id<-cbind(speciesCode, ecosystem, transect, predicted_abundance, scenario, year,model_id,spei_integ,spei_value)
              
              #append new results to output dataset
              beeCSF_output<-rbind(beeCSF_output,output_id)
            }
          }
        }
      }
    }
  }
  
}



##### Plains Grassland #####

# Read in data frame of historic and projected future SPEI data
climate <- read.csv("spei_historic_and_future_byscenario_CanESM2_2023-08-29.csv")

# Create data frame just containing projected future data (2021-2100)
climate_forecast<-subset(climate,source!="SEV-met" & station=="Blue Grama")

# Create a data frame of just 2020 data
clim_2020<-subset(climate,year==2020 & station=="Blue Grama")

# Create a copy of the 2020 data to be associated with the RCP 2.6 data series
clim_2020_2.6 <- clim_2020
clim_2020_2.6$source <- "rcp2.6"

# Create a copy of the 2020 data to be associated with the RCP 4.5 data series
clim_2020_4.5 <- clim_2020
clim_2020_4.5$source <- "rcp4.5"

# Create a copy of the 2020 data to be associated with the RCP 8.5 data series
clim_2020_8.5 <- clim_2020
clim_2020_8.5$source <- "rcp8.5"

# Create a merged climate data frame for 2020-2100
climate_forecast<-bind_rows(climate_forecast,clim_2020_2.6,clim_2020_4.5,clim_2020_8.5)

# Rename column
names(climate_forecast)[2]<-"scenario"

# Create different data frames for each transect, and add transect labels
climate_forecast_b1<-climate_forecast
climate_forecast_b1$transect <- "B1"

climate_forecast_b2<-climate_forecast
climate_forecast_b2$transect <- "B2"

climate_forecast_b3<-climate_forecast
climate_forecast_b3$transect <- "B3"

climate_forecast_b4<-climate_forecast
climate_forecast_b4$transect <- "B4"

climate_forecast_b5<-climate_forecast
climate_forecast_b5$transect <- "B5"

# Combine the transect data frames into a final data frame
climate_forecast<-bind_rows(climate_forecast_b1,climate_forecast_b2,climate_forecast_b3,climate_forecast_b4,climate_forecast_b5)

# Create an ecosystem code object
ecosystem <-"B"

# Create a data frame of just the species abundance matrix
speciesMatrixB <- blue[,9:232] 

# Create vector of species codes
speciesCodesB <- colnames(speciesMatrixB)

# Create vector of transect names
transectsB <- c("B1","B2","B3","B4","B5")

# Treat SPEI variables as numeric
climate_forecast$monsoon6SPEI<-as.numeric(climate_forecast$monsoon6SPEI)
climate_forecast$spring6SPEI<-as.numeric(climate_forecast$spring6SPEI)
climate_forecast$spring6SPEI_prioryear<-as.numeric(climate_forecast$spring6SPEI_prioryear)


### Loop: if a species is found in a given ecosystem, run CSFs, predict future bee abundances under different climate scenarios, and put the results in the beeCSF_output matrix
for (i in 1:length(speciesCodesB)) {
  
  # save the species code for column i
  speciesCode <- speciesCodesB[i]
  
  ### build CSFs
  
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
  
  # create list of models
  model_list <- list(m_null,m1_spring6mo,m2_spring6mo,m3_spring6mo,m1_spring6mo_AR1,m1_spring6mo_AR2,m2_spring6mo_AR1,m2_spring6mo_AR2,m3_spring6mo_AR1,m3_spring6mo_AR2,m1_monsoon6mo,m2_monsoon6mo,m3_monsoon6mo,m1_monsoon6mo_AR1,m1_monsoon6mo_AR2,m2_monsoon6mo_AR1,m2_monsoon6mo_AR2,m3_monsoon6mo_AR1,m3_monsoon6mo_AR2,m1_spring6mo_lag,m2_spring6mo_lag,m3_spring6mo_lag,m1_spring6mo_lag_AR1,m1_spring6mo_lag_AR2,m2_spring6mo_lag_AR1,m2_spring6mo_lag_AR2,m3_spring6mo_lag_AR1,m3_spring6mo_lag_AR2,m1_monsoon6mo_lag,m2_monsoon6mo_lag,m3_monsoon6mo_lag,m1_monsoon6mo_lag_AR1,m1_monsoon6mo_lag_AR2,m2_monsoon6mo_lag_AR1,m2_monsoon6mo_lag_AR2,m3_monsoon6mo_lag_AR1,m3_monsoon6mo_lag_AR2)
  
  # create data frame of model names and AICc values
  model_names <- c("m_null","m1_spring6mo","m2_spring6mo","m3_spring6mo","m1_spring6mo_AR1","m1_spring6mo_AR2","m2_spring6mo_AR1","m2_spring6mo_AR2","m3_spring6mo_AR1","m3_spring6mo_AR2","m1_monsoon6mo","m2_monsoon6mo","m3_monsoon6mo","m1_monsoon6mo_AR1","m1_monsoon6mo_AR2","m2_monsoon6mo_AR1","m2_monsoon6mo_AR2","m3_monsoon6mo_AR1","m3_monsoon6mo_AR2","m1_spring6mo_lag","m2_spring6mo_lag","m3_spring6mo_lag","m1_spring6mo_lag_AR1","m1_spring6mo_lag_AR2","m2_spring6mo_lag_AR1","m2_spring6mo_lag_AR2","m3_spring6mo_lag_AR1","m3_spring6mo_lag_AR2","m1_monsoon6mo_lag","m2_monsoon6mo_lag","m3_monsoon6mo_lag","m1_monsoon6mo_lag_AR1","m1_monsoon6mo_lag_AR2","m2_monsoon6mo_lag_AR1","m2_monsoon6mo_lag_AR2","m3_monsoon6mo_lag_AR1","m3_monsoon6mo_lag_AR2")
  
  model_aicc<-c(AICc(m_null),AICc(m1_spring6mo),AICc(m2_spring6mo),AICc(m3_spring6mo),AICc(m1_spring6mo_AR1),AICc(m1_spring6mo_AR2),AICc(m2_spring6mo_AR1),AICc(m2_spring6mo_AR2),AICc(m3_spring6mo_AR1),AICc(m3_spring6mo_AR2),AICc(m1_monsoon6mo),AICc(m2_monsoon6mo),AICc(m3_monsoon6mo),AICc(m1_monsoon6mo_AR1),AICc(m1_monsoon6mo_AR2),AICc(m2_monsoon6mo_AR1),AICc(m2_monsoon6mo_AR2),AICc(m3_monsoon6mo_AR1),AICc(m3_monsoon6mo_AR2),AICc(m1_spring6mo_lag),AICc(m2_spring6mo_lag),AICc(m3_spring6mo_lag),AICc(m1_spring6mo_lag_AR1),AICc(m1_spring6mo_lag_AR2),AICc(m2_spring6mo_lag_AR1),AICc(m2_spring6mo_lag_AR2),AICc(m3_spring6mo_lag_AR1),AICc(m3_spring6mo_lag_AR2),AICc(m1_monsoon6mo_lag),AICc(m2_monsoon6mo_lag),AICc(m3_monsoon6mo_lag),AICc(m1_monsoon6mo_lag_AR1),AICc(m1_monsoon6mo_lag_AR2),AICc(m2_monsoon6mo_lag_AR1),AICc(m2_monsoon6mo_lag_AR2),AICc(m3_monsoon6mo_lag_AR1),AICc(m3_monsoon6mo_lag_AR2))
  
  model_df<-data.frame(model_names,model_aicc)
  
  
  # calculate minimum AICc value
  min_aicc<-min(model_aicc)
  
  # create vectors of scenarios and years
  scenarios <- levels(as.factor(climate_forecast$scenario))
  years <- levels(as.factor(climate_forecast$year))
  
  # find the best CSF for the bee species, and use it to predict future abundance in different years under different climate scenarios
  for(j in 1:length(model_aicc)) {
    for(k in 1:length(scenarios)) {
      for (l in 1:length(years)) {
        for(m in 1:length(transectsB)) {
          rcp<-subset(climate_forecast,scenario==scenarios[k] & year == years[l] & transect == transectsB[m])
          if (model_aicc[j]==min_aicc) {
            best_model<-subset(model_df,model_aicc==min_aicc)
            if (best_model$model_names != "m_null") {
              model<-model_list[[j]]
              predicted_abundance_sqrt<-predict(model,rcp)[1]
              predicted_abundance<-predicted_abundance_sqrt^2
              scenario <- scenarios[k]
              year <- years[l]
              transect <- transectsB[m]
              
              model_id<-best_model[1,1]
              spei_integ<-names(model_list[[j]][4]$coefficients$fixed[2])
              spei_value<-rcp[,spei_integ]
              
              #bind relevant data columns together
              output_id<-cbind(speciesCode, ecosystem, transect, predicted_abundance, scenario, year,model_id,spei_integ,spei_value)
              
              #append new results to output dataset
              beeCSF_output<-rbind(beeCSF_output,output_id)
            }
          }
        }
      }
    }
  }
  
}

# create a data frame of the output
beeCSF_output <- data.frame(beeCSF_output)
beeCSF_output2<-beeCSF_output[-1,] # remove empty row

# Write .csv file of output
# write.csv(beeCSF_output2,"predicted_abundances_CanESM2-2023-08-29.csv", row.names=FALSE)
