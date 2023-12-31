################################################################################### 
# Climate sensitivity functions for bees in three ecosystems at the Sevilleta National Wildlife Refuge
# Code for determining whether CSF models run for each species x ecosystem combination

# Heat and desiccation tolerances predict bee abundance under climate change
# Melanie R. Kazenel, Karen W. Wright, Terry Griswold, Kenneth D. Whitney, and Jennifer A. Rudgers

# Date: 2023-12-01
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
library(tidyr)


##### Format data in preparation for running code #####

# Read in bee count data for 2002-2019
bee_wide<-read.csv("SEVBeeData2002-2019_revised2023-07-19.csv")

# Subset to exclude data from 2016 and 2017 due to incomplete sampling in these years
bee_wide<-subset(bee_wide,year!=2016 & year!=2017)

# Convert the data to long format
bee_melt<-pivot_longer(data=bee_wide,cols=11:351,names_to="code",values_to="count")

# Sum the bee counts for each month x year x transect x bee species combination (summing across the two traps on a given transect)
bee_melt_year <- bee_melt %>% group_by(month,year,ecosystem,transect,code) %>% summarise(count=sum(count))

# Calculate maximum bee abundance for each year x transect x bee species combination
dataset <- bee_melt_year %>% group_by(year,ecosystem,transect,code) %>% summarise(response=max(count))

# Remove bee species only observed in 2016 and/or 2017
dataset <- subset(dataset, code!="APTRILUN" & code!="HALASDI57")

# Create wide version of data frame
dataset_wide <- pivot_wider(data=dataset,id_cols=1:3,names_from="code",values_from="response")

# Add weather station column
dataset_wide$station<-"FivePoints"
dataset_wide$station[dataset_wide$ecosystem=="B"]<-"BlueGrama"

# Add station_year column
dataset_wide$station_year<-paste(dataset_wide$station,dataset_wide$year,sep = "_")

# Reorder columns
dataset_wide<-dataset_wide[,c(1:3,343:344,4:342)]

# Read in climate data
climate <- read.csv("spei_historic_and_future_byscenario_CanESM2_2023-12-01.csv")

# Add new weather station column
climate$station_no_space <- gsub(" ", "", climate$station)

# Add station_year column to climate data
climate$station_year<-paste(climate$station_no_space,climate$year,sep = "_")

# Join bee and climate data frames
dataset_wide<-left_join(dataset_wide,climate[,c(3:6,9)],by="station_year")

# Create new data frame of select columns for analyses
dataset_wide2<-dataset_wide[,c(1:5,345:348,6:344)]

# Write data file
# write.csv(dataset_wide2,"bee_wide_year_2002-2019_no2016or2017_maxabund_2023-12-01.csv",row.names = FALSE)

# Square root transform bee counts prior to analysis
dataset_wide2[,10:348] <- sqrt(dataset_wide2[,10:348])

# Write data file
# write.csv(dataset_wide2,"bee_wide_year_2002-2019_no2016or2017_sqrtmaxabund_2023-12-01.csv",row.names = FALSE)

# Create a separate data frame for each site
creo<-subset(dataset_wide2,ecosystem=="C")
black <- subset(dataset_wide2, ecosystem == "G")
blue <- subset(dataset_wide2, ecosystem == "B" & year != 2002 & year != 2003)


##### Loop to Determine which CSF Models Run: Chihuahuan Desert Shrubland Site #####

# Create matrix to hold results from the loop
beeCSF_output<-matrix(nrow=1,ncol=3,byrow=TRUE,dimnames=list(c("row1"),c("code","ecosystem","skip_to_next")))

# Create an ecosystem code object
ecosystem <-"C"

# Create a data frame of just the bee abundance matrix
speciesMatrixC <- creo[,10:348] 

# Create a vector of species codes
speciesCodesC <- colnames(speciesMatrixC)


# Loop:
# For each bee species, run each CSF model variant
# If there is an error in running the model (i.e., if the model does not run), code the skip_to_next variable as TRUE. Otherwise, this variable is coded as FALSE.
# The output is a data frame indicating, for each species and CSF model variant, whether the model ran (skip_to_next = FALSE) or did not run (skip_to_next = TRUE)

for (i in 1:length(speciesCodesC)) {
  
  # save the species code for column i
  speciesCode <- speciesCodesC[i]
  
  ### build CSFs
  
  skip_to_next <- FALSE
  
  tryCatch(lme(formula(paste(speciesCode, "~1")), random=~1|transect, data=creo, method="ML"),error = function(e) { skip_to_next <<- TRUE})
  #bind relevant data columns together
  output_id<-cbind(speciesCode, ecosystem, skip_to_next)
  #append new results to output dataset
  beeCSF_output<-rbind(beeCSF_output,output_id)
  
  skip_to_next <- FALSE
  
  tryCatch(lme(formula(paste(speciesCode, "~spring6SPEI")), random=~1|transect, data=creo, method="ML"),error = function(e) { skip_to_next <<- TRUE})
  #bind relevant data columns together
  output_id<-cbind(speciesCode, ecosystem, skip_to_next)
  #append new results to output dataset
  beeCSF_output<-rbind(beeCSF_output,output_id)
  
  skip_to_next <- FALSE
  
  tryCatch(lme(formula(paste(speciesCode, "~spring6SPEI +I(spring6SPEI^2)")), random=~1|transect, data=creo, method="ML"),error = function(e) { skip_to_next <<- TRUE})
  #bind relevant data columns together
  output_id<-cbind(speciesCode, ecosystem, skip_to_next)
  #append new results to output dataset
  beeCSF_output<-rbind(beeCSF_output,output_id)
  
  skip_to_next <- FALSE
  
  tryCatch(lme(formula(paste(speciesCode, "~spring6SPEI +I(spring6SPEI^2) + I(spring6SPEI^3)")), random=~1|transect, data=creo, method="ML") ,error = function(e) { skip_to_next <<- TRUE})
  #bind relevant data columns together
  output_id<-cbind(speciesCode, ecosystem, skip_to_next)
  #append new results to output dataset
  beeCSF_output<-rbind(beeCSF_output,output_id)
  
  skip_to_next <- FALSE
  
  tryCatch(lme(formula(paste(speciesCode, "~spring6SPEI")), correlation=corAR1(form=~year|transect), random=~1|transect, data=creo, method="ML"),error = function(e) { skip_to_next <<- TRUE})
  #bind relevant data columns together
  output_id<-cbind(speciesCode, ecosystem, skip_to_next)
  #append new results to output dataset
  beeCSF_output<-rbind(beeCSF_output,output_id)
  
  skip_to_next <- FALSE
  
  tryCatch(lme(formula(paste(speciesCode, "~spring6SPEI")), correlation=corARMA(form=~year|transect,p=2), random=~1|transect, data=creo, method="ML"),error = function(e) { skip_to_next <<- TRUE})
  #bind relevant data columns together
  output_id<-cbind(speciesCode, ecosystem, skip_to_next)
  #append new results to output dataset
  beeCSF_output<-rbind(beeCSF_output,output_id)
  
  skip_to_next <- FALSE
  
  tryCatch(lme(formula(paste(speciesCode, "~spring6SPEI +I(spring6SPEI^2)")), random=~1|transect, correlation=corAR1(form=~year|transect), data=creo, method="ML"),error = function(e) { skip_to_next <<- TRUE})
  #bind relevant data columns together
  output_id<-cbind(speciesCode, ecosystem, skip_to_next)
  #append new results to output dataset
  beeCSF_output<-rbind(beeCSF_output,output_id)
  
  skip_to_next <- FALSE
  
  tryCatch(lme(formula(paste(speciesCode, "~spring6SPEI +I(spring6SPEI^2)")), random=~1|transect, correlation=corARMA(form=~year|transect,p=2), data=creo, method="ML"),error = function(e) { skip_to_next <<- TRUE})
  #bind relevant data columns together
  output_id<-cbind(speciesCode, ecosystem, skip_to_next)
  #append new results to output dataset
  beeCSF_output<-rbind(beeCSF_output,output_id)
  
  skip_to_next <- FALSE
  
  tryCatch(lme(formula(paste(speciesCode, "~spring6SPEI +I(spring6SPEI^2) +I(spring6SPEI^3)")), random=~1|transect,correlation=corAR1(form=~year|transect), data=creo, method="ML"),error = function(e) { skip_to_next <<- TRUE})
  #bind relevant data columns together
  output_id<-cbind(speciesCode, ecosystem, skip_to_next)
  #append new results to output dataset
  beeCSF_output<-rbind(beeCSF_output,output_id)
  
  skip_to_next <- FALSE
  
  tryCatch(lme(formula(paste(speciesCode, "~spring6SPEI +I(spring6SPEI^2) +I(spring6SPEI^3)")), random=~1|transect, correlation=corARMA(form=~year|transect,p=2), data=creo, method="ML"),error = function(e) { skip_to_next <<- TRUE})
  #bind relevant data columns together
  output_id<-cbind(speciesCode, ecosystem, skip_to_next)
  #append new results to output dataset
  beeCSF_output<-rbind(beeCSF_output,output_id)
  
  skip_to_next <- FALSE
  
  tryCatch(lme(formula(paste(speciesCode, "~monsoon6SPEI")), random=~1|transect, data=creo, method="ML"),error = function(e) { skip_to_next <<- TRUE})
  #bind relevant data columns together
  output_id<-cbind(speciesCode, ecosystem, skip_to_next)
  #append new results to output dataset
  beeCSF_output<-rbind(beeCSF_output,output_id)
  
  skip_to_next <- FALSE
  
  tryCatch(lme(formula(paste(speciesCode, "~monsoon6SPEI +I(monsoon6SPEI^2)")), random=~1|transect, data=creo, method="ML"),error = function(e) { skip_to_next <<- TRUE})
  #bind relevant data columns together
  output_id<-cbind(speciesCode, ecosystem, skip_to_next)
  #append new results to output dataset
  beeCSF_output<-rbind(beeCSF_output,output_id)
  
  skip_to_next <- FALSE
  
  
  tryCatch(lme(formula(paste(speciesCode, "~monsoon6SPEI +I(monsoon6SPEI^2) +I(monsoon6SPEI^3)")), random=~1|transect, data=creo, method="ML"),error = function(e) { skip_to_next <<- TRUE})
  #bind relevant data columns together
  output_id<-cbind(speciesCode, ecosystem, skip_to_next)
  #append new results to output dataset
  beeCSF_output<-rbind(beeCSF_output,output_id)
  
  skip_to_next <- FALSE
  
  tryCatch(lme(formula(paste(speciesCode, "~monsoon6SPEI")), correlation=corAR1(form=~year|transect), random=~1|transect, data=creo, method="ML"),error = function(e) { skip_to_next <<- TRUE})
  #bind relevant data columns together
  output_id<-cbind(speciesCode, ecosystem, skip_to_next)
  #append new results to output dataset
  beeCSF_output<-rbind(beeCSF_output,output_id)
  
  skip_to_next <- FALSE
  
  tryCatch(lme(formula(paste(speciesCode, "~monsoon6SPEI")), correlation=corARMA(form=~year|transect,p=2), random=~1|transect, data=creo, method="ML"),error = function(e) { skip_to_next <<- TRUE})
  #bind relevant data columns together
  output_id<-cbind(speciesCode, ecosystem, skip_to_next)
  #append new results to output dataset
  beeCSF_output<-rbind(beeCSF_output,output_id)
  
  skip_to_next <- FALSE
  
  tryCatch(lme(formula(paste(speciesCode, "~monsoon6SPEI + I(monsoon6SPEI^2)")), correlation=corAR1(form=~year|transect), random=~1|transect, data=creo, method="ML"),error = function(e) { skip_to_next <<- TRUE})
  #bind relevant data columns together
  output_id<-cbind(speciesCode, ecosystem, skip_to_next)
  #append new results to output dataset
  beeCSF_output<-rbind(beeCSF_output,output_id)
  
  skip_to_next <- FALSE
  
  tryCatch(lme(formula(paste(speciesCode, "~monsoon6SPEI + I(monsoon6SPEI^2)")), correlation=corARMA(form=~year|transect,p=2), random=~1|transect, data=creo, method="ML"),error = function(e) { skip_to_next <<- TRUE})
  #bind relevant data columns together
  output_id<-cbind(speciesCode, ecosystem, skip_to_next)
  #append new results to output dataset
  beeCSF_output<-rbind(beeCSF_output,output_id)
  
  skip_to_next <- FALSE
  
  tryCatch(lme(formula(paste(speciesCode, "~monsoon6SPEI + I(monsoon6SPEI^2) + I(monsoon6SPEI^3)")), correlation=corAR1(form=~year|transect), random=~1|transect, data=creo, method="ML"),error = function(e) { skip_to_next <<- TRUE})
  #bind relevant data columns together
  output_id<-cbind(speciesCode, ecosystem, skip_to_next)
  #append new results to output dataset
  beeCSF_output<-rbind(beeCSF_output,output_id)
  
  skip_to_next <- FALSE
  
  tryCatch(lme(formula(paste(speciesCode, "~monsoon6SPEI + I(monsoon6SPEI^2) + I(monsoon6SPEI^3)")), correlation=corARMA(form=~year|transect,p=2), random=~1|transect, data=creo, method="ML"),error = function(e) { skip_to_next <<- TRUE})
  #bind relevant data columns together
  output_id<-cbind(speciesCode, ecosystem, skip_to_next)
  #append new results to output dataset
  beeCSF_output<-rbind(beeCSF_output,output_id)
  
  skip_to_next <- FALSE
  
  
  tryCatch(lme(formula(paste(speciesCode, "~spring6SPEI_prioryear")), random=~1|transect, data=creo, method="ML"),error = function(e) { skip_to_next <<- TRUE})
  #bind relevant data columns together
  output_id<-cbind(speciesCode, ecosystem, skip_to_next)
  #append new results to output dataset
  beeCSF_output<-rbind(beeCSF_output,output_id)
  
  skip_to_next <- FALSE
  
  
  tryCatch(lme(formula(paste(speciesCode, "~spring6SPEI_prioryear +I(spring6SPEI_prioryear^2)")), random=~1|transect, data=creo, method="ML"),error = function(e) { skip_to_next <<- TRUE})
  #bind relevant data columns together
  output_id<-cbind(speciesCode, ecosystem, skip_to_next)
  #append new results to output dataset
  beeCSF_output<-rbind(beeCSF_output,output_id)
  
  skip_to_next <- FALSE
  
  tryCatch(lme(formula(paste(speciesCode, "~spring6SPEI_prioryear +I(spring6SPEI_prioryear^2) + I(spring6SPEI_prioryear^3)")), random=~1|transect, data=creo, method="ML") ,error = function(e) { skip_to_next <<- TRUE})
  #bind relevant data columns together
  output_id<-cbind(speciesCode, ecosystem, skip_to_next)
  #append new results to output dataset
  beeCSF_output<-rbind(beeCSF_output,output_id)
  
  skip_to_next <- FALSE
  
  tryCatch(lme(formula(paste(speciesCode, "~spring6SPEI_prioryear")), correlation=corAR1(form=~year|transect), random=~1|transect, data=creo, method="ML"),error = function(e) { skip_to_next <<- TRUE})
  #bind relevant data columns together
  output_id<-cbind(speciesCode, ecosystem, skip_to_next)
  #append new results to output dataset
  beeCSF_output<-rbind(beeCSF_output,output_id)
  
  skip_to_next <- FALSE
  
  tryCatch(lme(formula(paste(speciesCode, "~spring6SPEI_prioryear")), correlation=corARMA(form=~year|transect,p=2), random=~1|transect, data=creo, method="ML"),error = function(e) { skip_to_next <<- TRUE})
  #bind relevant data columns together
  output_id<-cbind(speciesCode, ecosystem, skip_to_next)
  #append new results to output dataset
  beeCSF_output<-rbind(beeCSF_output,output_id)
  
  skip_to_next <- FALSE
  
  tryCatch(lme(formula(paste(speciesCode, "~spring6SPEI_prioryear +I(spring6SPEI_prioryear^2)")), random=~1|transect, correlation=corAR1(form=~year|transect), data=creo, method="ML"),error = function(e) { skip_to_next <<- TRUE})
  #bind relevant data columns together
  output_id<-cbind(speciesCode, ecosystem, skip_to_next)
  #append new results to output dataset
  beeCSF_output<-rbind(beeCSF_output,output_id)
  
  skip_to_next <- FALSE
  
  tryCatch(lme(formula(paste(speciesCode, "~spring6SPEI_prioryear +I(spring6SPEI_prioryear^2)")), random=~1|transect, correlation=corARMA(form=~year|transect,p=2), data=creo, method="ML"),error = function(e) { skip_to_next <<- TRUE})
  #bind relevant data columns together
  output_id<-cbind(speciesCode, ecosystem, skip_to_next)
  #append new results to output dataset
  beeCSF_output<-rbind(beeCSF_output,output_id)
  
  skip_to_next <- FALSE
  
  tryCatch(lme(formula(paste(speciesCode, "~spring6SPEI_prioryear +I(spring6SPEI_prioryear^2) +I(spring6SPEI_prioryear^3)")), random=~1|transect,correlation=corAR1(form=~year|transect), data=creo, method="ML"),error = function(e) { skip_to_next <<- TRUE})
  #bind relevant data columns together
  output_id<-cbind(speciesCode, ecosystem, skip_to_next)
  #append new results to output dataset
  beeCSF_output<-rbind(beeCSF_output,output_id)
  
  skip_to_next <- FALSE
  
  tryCatch(lme(formula(paste(speciesCode, "~spring6SPEI_prioryear +I(spring6SPEI_prioryear^2) +I(spring6SPEI_prioryear^3)")), random=~1|transect, correlation=corARMA(form=~year|transect,p=2), data=creo, method="ML"),error = function(e) { skip_to_next <<- TRUE})
  #bind relevant data columns together
  output_id<-cbind(speciesCode, ecosystem, skip_to_next)
  #append new results to output dataset
  beeCSF_output<-rbind(beeCSF_output,output_id)
  
  skip_to_next <- FALSE
  
  tryCatch(lme(formula(paste(speciesCode, "~monsoon6SPEI_prioryear")), random=~1|transect, data=creo, method="ML"),error = function(e) { skip_to_next <<- TRUE})
  #bind relevant data columns together
  output_id<-cbind(speciesCode, ecosystem, skip_to_next)
  #append new results to output dataset
  beeCSF_output<-rbind(beeCSF_output,output_id)
  
  skip_to_next <- FALSE
  
  tryCatch(lme(formula(paste(speciesCode, "~monsoon6SPEI_prioryear +I(monsoon6SPEI_prioryear^2)")), random=~1|transect, data=creo, method="ML"),error = function(e) { skip_to_next <<- TRUE})
  #bind relevant data columns together
  output_id<-cbind(speciesCode, ecosystem, skip_to_next)
  #append new results to output dataset
  beeCSF_output<-rbind(beeCSF_output,output_id)
  
  skip_to_next <- FALSE
  
  tryCatch(lme(formula(paste(speciesCode, "~monsoon6SPEI_prioryear +I(monsoon6SPEI_prioryear^2) +I(monsoon6SPEI_prioryear^3)")), random=~1|transect, data=creo, method="ML"),error = function(e) { skip_to_next <<- TRUE})
  #bind relevant data columns together
  output_id<-cbind(speciesCode, ecosystem, skip_to_next)
  #append new results to output dataset
  beeCSF_output<-rbind(beeCSF_output,output_id)
  
  skip_to_next <- FALSE
  
  tryCatch(lme(formula(paste(speciesCode, "~monsoon6SPEI_prioryear")), correlation=corAR1(form=~year|transect), random=~1|transect, data=creo, method="ML"),error = function(e) { skip_to_next <<- TRUE})
  #bind relevant data columns together
  output_id<-cbind(speciesCode, ecosystem, skip_to_next)
  #append new results to output dataset
  beeCSF_output<-rbind(beeCSF_output,output_id)
  
  skip_to_next <- FALSE
  
  tryCatch(lme(formula(paste(speciesCode, "~monsoon6SPEI_prioryear")), correlation=corARMA(form=~year|transect,p=2), random=~1|transect, data=creo, method="ML"),error = function(e) { skip_to_next <<- TRUE})
  #bind relevant data columns together
  output_id<-cbind(speciesCode, ecosystem, skip_to_next)
  #append new results to output dataset
  beeCSF_output<-rbind(beeCSF_output,output_id)
  
  skip_to_next <- FALSE
  
  tryCatch(lme(formula(paste(speciesCode, "~monsoon6SPEI_prioryear + I(monsoon6SPEI_prioryear^2)")), correlation=corAR1(form=~year|transect), random=~1|transect, data=creo, method="ML"),error = function(e) { skip_to_next <<- TRUE})
  #bind relevant data columns together
  output_id<-cbind(speciesCode, ecosystem, skip_to_next)
  #append new results to output dataset
  beeCSF_output<-rbind(beeCSF_output,output_id)
  
  skip_to_next <- FALSE
  
  tryCatch(lme(formula(paste(speciesCode, "~monsoon6SPEI_prioryear + I(monsoon6SPEI_prioryear^2)")), correlation=corARMA(form=~year|transect,p=2), random=~1|transect, data=creo, method="ML"),error = function(e) { skip_to_next <<- TRUE})
  #bind relevant data columns together
  output_id<-cbind(speciesCode, ecosystem, skip_to_next)
  #append new results to output dataset
  beeCSF_output<-rbind(beeCSF_output,output_id)
  
  skip_to_next <- FALSE
  
  tryCatch(lme(formula(paste(speciesCode, "~monsoon6SPEI_prioryear + I(monsoon6SPEI_prioryear^2) + I(monsoon6SPEI_prioryear^3)")), correlation=corAR1(form=~year|transect), random=~1|transect, data=creo, method="ML"),error = function(e) { skip_to_next <<- TRUE})
  #bind relevant data columns together
  output_id<-cbind(speciesCode, ecosystem, skip_to_next)
  #append new results to output dataset
  beeCSF_output<-rbind(beeCSF_output,output_id)
  
  skip_to_next <- FALSE
  
  tryCatch(lme(formula(paste(speciesCode, "~monsoon6SPEI_prioryear + I(monsoon6SPEI_prioryear^2) + I(monsoon6SPEI_prioryear^3)")), correlation=corARMA(form=~year|transect,p=2), random=~1|transect, data=creo, method="ML"),error = function(e) { skip_to_next <<- TRUE})
  #bind relevant data columns together
  output_id<-cbind(speciesCode, ecosystem, skip_to_next)
  #append new results to output dataset
  beeCSF_output<-rbind(beeCSF_output,output_id)
  
}


##### Loop to Determine which CSF Models Run: Chihuahuan Desert Grassland Site #####

# Create an ecosystem code object
ecosystem <-"G"

# Create a data frame of just the bee abundance matrix
speciesMatrixG <- black[,10:348] 

### Create a vector of species codes
speciesCodesG <- colnames(speciesMatrixG)


# Loop:
# For each bee species, run each CSF model variant
# If there is an error in running the model (i.e., if the model does not run), code the skip_to_next variable as TRUE. Otherwise, this variable is coded as FALSE.
# The output is a data frame indicating, for each species and CSF model variant, whether the model ran (skip_to_next = FALSE) or did not run (skip_to_next = TRUE)

for (i in 1:length(speciesCodesG)) {
  
  # save the species code for column i
  speciesCode <- speciesCodesG[i]
  
  ### build CSFs
  
  skip_to_next <- FALSE
  
  tryCatch(lme(formula(paste(speciesCode, "~1")), random=~1|transect, data=black, method="ML"),error = function(e) { skip_to_next <<- TRUE})
  #bind relevant data columns together
  output_id<-cbind(speciesCode, ecosystem, skip_to_next)
  #append new results to output dataset
  beeCSF_output<-rbind(beeCSF_output,output_id)
  
  skip_to_next <- FALSE
  
  tryCatch(lme(formula(paste(speciesCode, "~spring6SPEI")), random=~1|transect, data=black, method="ML"),error = function(e) { skip_to_next <<- TRUE})
  #bind relevant data columns together
  output_id<-cbind(speciesCode, ecosystem, skip_to_next)
  #append new results to output dataset
  beeCSF_output<-rbind(beeCSF_output,output_id)
  
  skip_to_next <- FALSE
  
  tryCatch(lme(formula(paste(speciesCode, "~spring6SPEI +I(spring6SPEI^2)")), random=~1|transect, data=black, method="ML"),error = function(e) { skip_to_next <<- TRUE})
  #bind relevant data columns together
  output_id<-cbind(speciesCode, ecosystem, skip_to_next)
  #append new results to output dataset
  beeCSF_output<-rbind(beeCSF_output,output_id)
  
  skip_to_next <- FALSE
  
  tryCatch(lme(formula(paste(speciesCode, "~spring6SPEI +I(spring6SPEI^2) + I(spring6SPEI^3)")), random=~1|transect, data=black, method="ML") ,error = function(e) { skip_to_next <<- TRUE})
  #bind relevant data columns together
  output_id<-cbind(speciesCode, ecosystem, skip_to_next)
  #append new results to output dataset
  beeCSF_output<-rbind(beeCSF_output,output_id)
  
  skip_to_next <- FALSE
  
  tryCatch(lme(formula(paste(speciesCode, "~spring6SPEI")), correlation=corAR1(form=~year|transect), random=~1|transect, data=black, method="ML"),error = function(e) { skip_to_next <<- TRUE})
  #bind relevant data columns together
  output_id<-cbind(speciesCode, ecosystem, skip_to_next)
  #append new results to output dataset
  beeCSF_output<-rbind(beeCSF_output,output_id)
  
  skip_to_next <- FALSE
  
  tryCatch(lme(formula(paste(speciesCode, "~spring6SPEI")), correlation=corARMA(form=~year|transect,p=2), random=~1|transect, data=black, method="ML"),error = function(e) { skip_to_next <<- TRUE})
  #bind relevant data columns together
  output_id<-cbind(speciesCode, ecosystem, skip_to_next)
  #append new results to output dataset
  beeCSF_output<-rbind(beeCSF_output,output_id)
  
  skip_to_next <- FALSE
  
  tryCatch(lme(formula(paste(speciesCode, "~spring6SPEI +I(spring6SPEI^2)")), random=~1|transect, correlation=corAR1(form=~year|transect), data=black, method="ML"),error = function(e) { skip_to_next <<- TRUE})
  #bind relevant data columns together
  output_id<-cbind(speciesCode, ecosystem, skip_to_next)
  #append new results to output dataset
  beeCSF_output<-rbind(beeCSF_output,output_id)
  
  skip_to_next <- FALSE
  
  tryCatch(lme(formula(paste(speciesCode, "~spring6SPEI +I(spring6SPEI^2)")), random=~1|transect, correlation=corARMA(form=~year|transect,p=2), data=black, method="ML"),error = function(e) { skip_to_next <<- TRUE})
  #bind relevant data columns together
  output_id<-cbind(speciesCode, ecosystem, skip_to_next)
  #append new results to output dataset
  beeCSF_output<-rbind(beeCSF_output,output_id)
  
  skip_to_next <- FALSE
  
  tryCatch(lme(formula(paste(speciesCode, "~spring6SPEI +I(spring6SPEI^2) +I(spring6SPEI^3)")), random=~1|transect,correlation=corAR1(form=~year|transect), data=black, method="ML"),error = function(e) { skip_to_next <<- TRUE})
  #bind relevant data columns together
  output_id<-cbind(speciesCode, ecosystem, skip_to_next)
  #append new results to output dataset
  beeCSF_output<-rbind(beeCSF_output,output_id)
  
  skip_to_next <- FALSE
  
  tryCatch(lme(formula(paste(speciesCode, "~spring6SPEI +I(spring6SPEI^2) +I(spring6SPEI^3)")), random=~1|transect, correlation=corARMA(form=~year|transect,p=2), data=black, method="ML"),error = function(e) { skip_to_next <<- TRUE})
  #bind relevant data columns together
  output_id<-cbind(speciesCode, ecosystem, skip_to_next)
  #append new results to output dataset
  beeCSF_output<-rbind(beeCSF_output,output_id)
  
  skip_to_next <- FALSE
  
  tryCatch(lme(formula(paste(speciesCode, "~monsoon6SPEI")), random=~1|transect, data=black, method="ML"),error = function(e) { skip_to_next <<- TRUE})
  #bind relevant data columns together
  output_id<-cbind(speciesCode, ecosystem, skip_to_next)
  #append new results to output dataset
  beeCSF_output<-rbind(beeCSF_output,output_id)
  
  skip_to_next <- FALSE
  
  tryCatch(lme(formula(paste(speciesCode, "~monsoon6SPEI +I(monsoon6SPEI^2)")), random=~1|transect, data=black, method="ML"),error = function(e) { skip_to_next <<- TRUE})
  #bind relevant data columns together
  output_id<-cbind(speciesCode, ecosystem, skip_to_next)
  #append new results to output dataset
  beeCSF_output<-rbind(beeCSF_output,output_id)
  
  skip_to_next <- FALSE
  
  
  tryCatch(lme(formula(paste(speciesCode, "~monsoon6SPEI +I(monsoon6SPEI^2) +I(monsoon6SPEI^3)")), random=~1|transect, data=black, method="ML"),error = function(e) { skip_to_next <<- TRUE})
  #bind relevant data columns together
  output_id<-cbind(speciesCode, ecosystem, skip_to_next)
  #append new results to output dataset
  beeCSF_output<-rbind(beeCSF_output,output_id)
  
  skip_to_next <- FALSE
  
  tryCatch(lme(formula(paste(speciesCode, "~monsoon6SPEI")), correlation=corAR1(form=~year|transect), random=~1|transect, data=black, method="ML"),error = function(e) { skip_to_next <<- TRUE})
  #bind relevant data columns together
  output_id<-cbind(speciesCode, ecosystem, skip_to_next)
  #append new results to output dataset
  beeCSF_output<-rbind(beeCSF_output,output_id)
  
  skip_to_next <- FALSE
  
  tryCatch(lme(formula(paste(speciesCode, "~monsoon6SPEI")), correlation=corARMA(form=~year|transect,p=2), random=~1|transect, data=black, method="ML"),error = function(e) { skip_to_next <<- TRUE})
  #bind relevant data columns together
  output_id<-cbind(speciesCode, ecosystem, skip_to_next)
  #append new results to output dataset
  beeCSF_output<-rbind(beeCSF_output,output_id)
  
  skip_to_next <- FALSE
  
  tryCatch(lme(formula(paste(speciesCode, "~monsoon6SPEI + I(monsoon6SPEI^2)")), correlation=corAR1(form=~year|transect), random=~1|transect, data=black, method="ML"),error = function(e) { skip_to_next <<- TRUE})
  #bind relevant data columns together
  output_id<-cbind(speciesCode, ecosystem, skip_to_next)
  #append new results to output dataset
  beeCSF_output<-rbind(beeCSF_output,output_id)
  
  skip_to_next <- FALSE
  
  tryCatch(lme(formula(paste(speciesCode, "~monsoon6SPEI + I(monsoon6SPEI^2)")), correlation=corARMA(form=~year|transect,p=2), random=~1|transect, data=black, method="ML"),error = function(e) { skip_to_next <<- TRUE})
  #bind relevant data columns together
  output_id<-cbind(speciesCode, ecosystem, skip_to_next)
  #append new results to output dataset
  beeCSF_output<-rbind(beeCSF_output,output_id)
  
  skip_to_next <- FALSE
  
  tryCatch(lme(formula(paste(speciesCode, "~monsoon6SPEI + I(monsoon6SPEI^2) + I(monsoon6SPEI^3)")), correlation=corAR1(form=~year|transect), random=~1|transect, data=black, method="ML"),error = function(e) { skip_to_next <<- TRUE})
  #bind relevant data columns together
  output_id<-cbind(speciesCode, ecosystem, skip_to_next)
  #append new results to output dataset
  beeCSF_output<-rbind(beeCSF_output,output_id)
  
  skip_to_next <- FALSE
  
  tryCatch(lme(formula(paste(speciesCode, "~monsoon6SPEI + I(monsoon6SPEI^2) + I(monsoon6SPEI^3)")), correlation=corARMA(form=~year|transect,p=2), random=~1|transect, data=black, method="ML"),error = function(e) { skip_to_next <<- TRUE})
  #bind relevant data columns together
  output_id<-cbind(speciesCode, ecosystem, skip_to_next)
  #append new results to output dataset
  beeCSF_output<-rbind(beeCSF_output,output_id)
  
  skip_to_next <- FALSE
  
  
  tryCatch(lme(formula(paste(speciesCode, "~spring6SPEI_prioryear")), random=~1|transect, data=black, method="ML"),error = function(e) { skip_to_next <<- TRUE})
  #bind relevant data columns together
  output_id<-cbind(speciesCode, ecosystem, skip_to_next)
  #append new results to output dataset
  beeCSF_output<-rbind(beeCSF_output,output_id)
  
  skip_to_next <- FALSE
  
  
  tryCatch(lme(formula(paste(speciesCode, "~spring6SPEI_prioryear +I(spring6SPEI_prioryear^2)")), random=~1|transect, data=black, method="ML"),error = function(e) { skip_to_next <<- TRUE})
  #bind relevant data columns together
  output_id<-cbind(speciesCode, ecosystem, skip_to_next)
  #append new results to output dataset
  beeCSF_output<-rbind(beeCSF_output,output_id)
  
  skip_to_next <- FALSE
  
  tryCatch(lme(formula(paste(speciesCode, "~spring6SPEI_prioryear +I(spring6SPEI_prioryear^2) + I(spring6SPEI_prioryear^3)")), random=~1|transect, data=black, method="ML") ,error = function(e) { skip_to_next <<- TRUE})
  #bind relevant data columns together
  output_id<-cbind(speciesCode, ecosystem, skip_to_next)
  #append new results to output dataset
  beeCSF_output<-rbind(beeCSF_output,output_id)
  
  skip_to_next <- FALSE
  
  tryCatch(lme(formula(paste(speciesCode, "~spring6SPEI_prioryear")), correlation=corAR1(form=~year|transect), random=~1|transect, data=black, method="ML"),error = function(e) { skip_to_next <<- TRUE})
  #bind relevant data columns together
  output_id<-cbind(speciesCode, ecosystem, skip_to_next)
  #append new results to output dataset
  beeCSF_output<-rbind(beeCSF_output,output_id)
  
  skip_to_next <- FALSE
  
  tryCatch(lme(formula(paste(speciesCode, "~spring6SPEI_prioryear")), correlation=corARMA(form=~year|transect,p=2), random=~1|transect, data=black, method="ML"),error = function(e) { skip_to_next <<- TRUE})
  #bind relevant data columns together
  output_id<-cbind(speciesCode, ecosystem, skip_to_next)
  #append new results to output dataset
  beeCSF_output<-rbind(beeCSF_output,output_id)
  
  skip_to_next <- FALSE
  
  tryCatch(lme(formula(paste(speciesCode, "~spring6SPEI_prioryear +I(spring6SPEI_prioryear^2)")), random=~1|transect, correlation=corAR1(form=~year|transect), data=black, method="ML"),error = function(e) { skip_to_next <<- TRUE})
  #bind relevant data columns together
  output_id<-cbind(speciesCode, ecosystem, skip_to_next)
  #append new results to output dataset
  beeCSF_output<-rbind(beeCSF_output,output_id)
  
  skip_to_next <- FALSE
  
  tryCatch(lme(formula(paste(speciesCode, "~spring6SPEI_prioryear +I(spring6SPEI_prioryear^2)")), random=~1|transect, correlation=corARMA(form=~year|transect,p=2), data=black, method="ML"),error = function(e) { skip_to_next <<- TRUE})
  #bind relevant data columns together
  output_id<-cbind(speciesCode, ecosystem, skip_to_next)
  #append new results to output dataset
  beeCSF_output<-rbind(beeCSF_output,output_id)
  
  skip_to_next <- FALSE
  
  tryCatch(lme(formula(paste(speciesCode, "~spring6SPEI_prioryear +I(spring6SPEI_prioryear^2) +I(spring6SPEI_prioryear^3)")), random=~1|transect,correlation=corAR1(form=~year|transect), data=black, method="ML"),error = function(e) { skip_to_next <<- TRUE})
  #bind relevant data columns together
  output_id<-cbind(speciesCode, ecosystem, skip_to_next)
  #append new results to output dataset
  beeCSF_output<-rbind(beeCSF_output,output_id)
  
  skip_to_next <- FALSE
  
  tryCatch(lme(formula(paste(speciesCode, "~spring6SPEI_prioryear +I(spring6SPEI_prioryear^2) +I(spring6SPEI_prioryear^3)")), random=~1|transect, correlation=corARMA(form=~year|transect,p=2), data=black, method="ML"),error = function(e) { skip_to_next <<- TRUE})
  #bind relevant data columns together
  output_id<-cbind(speciesCode, ecosystem, skip_to_next)
  #append new results to output dataset
  beeCSF_output<-rbind(beeCSF_output,output_id)
  
  skip_to_next <- FALSE
  
  tryCatch(lme(formula(paste(speciesCode, "~monsoon6SPEI_prioryear")), random=~1|transect, data=black, method="ML"),error = function(e) { skip_to_next <<- TRUE})
  #bind relevant data columns together
  output_id<-cbind(speciesCode, ecosystem, skip_to_next)
  #append new results to output dataset
  beeCSF_output<-rbind(beeCSF_output,output_id)
  
  skip_to_next <- FALSE
  
  tryCatch(lme(formula(paste(speciesCode, "~monsoon6SPEI_prioryear +I(monsoon6SPEI_prioryear^2)")), random=~1|transect, data=black, method="ML"),error = function(e) { skip_to_next <<- TRUE})
  #bind relevant data columns together
  output_id<-cbind(speciesCode, ecosystem, skip_to_next)
  #append new results to output dataset
  beeCSF_output<-rbind(beeCSF_output,output_id)
  
  skip_to_next <- FALSE
  
  tryCatch(lme(formula(paste(speciesCode, "~monsoon6SPEI_prioryear +I(monsoon6SPEI_prioryear^2) +I(monsoon6SPEI_prioryear^3)")), random=~1|transect, data=black, method="ML"),error = function(e) { skip_to_next <<- TRUE})
  #bind relevant data columns together
  output_id<-cbind(speciesCode, ecosystem, skip_to_next)
  #append new results to output dataset
  beeCSF_output<-rbind(beeCSF_output,output_id)
  
  skip_to_next <- FALSE
  
  tryCatch(lme(formula(paste(speciesCode, "~monsoon6SPEI_prioryear")), correlation=corAR1(form=~year|transect), random=~1|transect, data=black, method="ML"),error = function(e) { skip_to_next <<- TRUE})
  #bind relevant data columns together
  output_id<-cbind(speciesCode, ecosystem, skip_to_next)
  #append new results to output dataset
  beeCSF_output<-rbind(beeCSF_output,output_id)
  
  skip_to_next <- FALSE
  
  tryCatch(lme(formula(paste(speciesCode, "~monsoon6SPEI_prioryear")), correlation=corARMA(form=~year|transect,p=2), random=~1|transect, data=black, method="ML"),error = function(e) { skip_to_next <<- TRUE})
  #bind relevant data columns together
  output_id<-cbind(speciesCode, ecosystem, skip_to_next)
  #append new results to output dataset
  beeCSF_output<-rbind(beeCSF_output,output_id)
  
  skip_to_next <- FALSE
  
  tryCatch(lme(formula(paste(speciesCode, "~monsoon6SPEI_prioryear + I(monsoon6SPEI_prioryear^2)")), correlation=corAR1(form=~year|transect), random=~1|transect, data=black, method="ML"),error = function(e) { skip_to_next <<- TRUE})
  #bind relevant data columns together
  output_id<-cbind(speciesCode, ecosystem, skip_to_next)
  #append new results to output dataset
  beeCSF_output<-rbind(beeCSF_output,output_id)
  
  skip_to_next <- FALSE
  
  tryCatch(lme(formula(paste(speciesCode, "~monsoon6SPEI_prioryear + I(monsoon6SPEI_prioryear^2)")), correlation=corARMA(form=~year|transect,p=2), random=~1|transect, data=black, method="ML"),error = function(e) { skip_to_next <<- TRUE})
  #bind relevant data columns together
  output_id<-cbind(speciesCode, ecosystem, skip_to_next)
  #append new results to output dataset
  beeCSF_output<-rbind(beeCSF_output,output_id)
  
  skip_to_next <- FALSE
  
  tryCatch(lme(formula(paste(speciesCode, "~monsoon6SPEI_prioryear + I(monsoon6SPEI_prioryear^2) + I(monsoon6SPEI_prioryear^3)")), correlation=corAR1(form=~year|transect), random=~1|transect, data=black, method="ML"),error = function(e) { skip_to_next <<- TRUE})
  #bind relevant data columns together
  output_id<-cbind(speciesCode, ecosystem, skip_to_next)
  #append new results to output dataset
  beeCSF_output<-rbind(beeCSF_output,output_id)
  
  skip_to_next <- FALSE
  
  tryCatch(lme(formula(paste(speciesCode, "~monsoon6SPEI_prioryear + I(monsoon6SPEI_prioryear^2) + I(monsoon6SPEI_prioryear^3)")), correlation=corARMA(form=~year|transect,p=2), random=~1|transect, data=black, method="ML"),error = function(e) { skip_to_next <<- TRUE})
  #bind relevant data columns together
  output_id<-cbind(speciesCode, ecosystem, skip_to_next)
  #append new results to output dataset
  beeCSF_output<-rbind(beeCSF_output,output_id)
  
}


##### Loop to Determine which CSF Models Run: Plains Grassland Site #####

# Create an ecosystem code object
ecosystem <-"B"

# Create a data frame of just the bee abundance matrix
speciesMatrixB <- blue[,10:348] 

# Create vector of species codes
speciesCodesB <- colnames(speciesMatrixB)


# Loop:
# For each bee species, run each CSF model variant
# If there is an error in running the model (i.e., if the model does not run), code the skip_to_next variable as TRUE. Otherwise, this variable is coded as FALSE.
# The output is a data frame indicating, for each species and CSF model variant, whether the model ran (skip_to_next = FALSE) or did not run (skip_to_next = TRUE)

for (i in 1:length(speciesCodesB)) {
  
  # save the species code for column i
  speciesCode <- speciesCodesB[i]
  
  ### build CSFs
  
  skip_to_next <- FALSE
  
  tryCatch(lme(formula(paste(speciesCode, "~1")), random=~1|transect, data=blue, method="ML"),error = function(e) { skip_to_next <<- TRUE})
  #bind relevant data columns together
  output_id<-cbind(speciesCode, ecosystem, skip_to_next)
  #append new results to output dataset
  beeCSF_output<-rbind(beeCSF_output,output_id)
  
  skip_to_next <- FALSE
  
  tryCatch(lme(formula(paste(speciesCode, "~spring6SPEI")), random=~1|transect, data=blue, method="ML"),error = function(e) { skip_to_next <<- TRUE})
  #bind relevant data columns together
  output_id<-cbind(speciesCode, ecosystem, skip_to_next)
  #append new results to output dataset
  beeCSF_output<-rbind(beeCSF_output,output_id)
  
  skip_to_next <- FALSE
  
  tryCatch(lme(formula(paste(speciesCode, "~spring6SPEI +I(spring6SPEI^2)")), random=~1|transect, data=blue, method="ML"),error = function(e) { skip_to_next <<- TRUE})
  #bind relevant data columns together
  output_id<-cbind(speciesCode, ecosystem, skip_to_next)
  #append new results to output dataset
  beeCSF_output<-rbind(beeCSF_output,output_id)
  
  skip_to_next <- FALSE
  
  tryCatch(lme(formula(paste(speciesCode, "~spring6SPEI +I(spring6SPEI^2) + I(spring6SPEI^3)")), random=~1|transect, data=blue, method="ML") ,error = function(e) { skip_to_next <<- TRUE})
  #bind relevant data columns together
  output_id<-cbind(speciesCode, ecosystem, skip_to_next)
  #append new results to output dataset
  beeCSF_output<-rbind(beeCSF_output,output_id)
  
  skip_to_next <- FALSE
  
  tryCatch(lme(formula(paste(speciesCode, "~spring6SPEI")), correlation=corAR1(form=~year|transect), random=~1|transect, data=blue, method="ML"),error = function(e) { skip_to_next <<- TRUE})
  #bind relevant data columns together
  output_id<-cbind(speciesCode, ecosystem, skip_to_next)
  #append new results to output dataset
  beeCSF_output<-rbind(beeCSF_output,output_id)
  
  skip_to_next <- FALSE
  
  tryCatch(lme(formula(paste(speciesCode, "~spring6SPEI")), correlation=corARMA(form=~year|transect,p=2), random=~1|transect, data=blue, method="ML"),error = function(e) { skip_to_next <<- TRUE})
  #bind relevant data columns together
  output_id<-cbind(speciesCode, ecosystem, skip_to_next)
  #append new results to output dataset
  beeCSF_output<-rbind(beeCSF_output,output_id)
  
  skip_to_next <- FALSE
  
  tryCatch(lme(formula(paste(speciesCode, "~spring6SPEI +I(spring6SPEI^2)")), random=~1|transect, correlation=corAR1(form=~year|transect), data=blue, method="ML"),error = function(e) { skip_to_next <<- TRUE})
  #bind relevant data columns together
  output_id<-cbind(speciesCode, ecosystem, skip_to_next)
  #append new results to output dataset
  beeCSF_output<-rbind(beeCSF_output,output_id)
  
  skip_to_next <- FALSE
  
  tryCatch(lme(formula(paste(speciesCode, "~spring6SPEI +I(spring6SPEI^2)")), random=~1|transect, correlation=corARMA(form=~year|transect,p=2), data=blue, method="ML"),error = function(e) { skip_to_next <<- TRUE})
  #bind relevant data columns together
  output_id<-cbind(speciesCode, ecosystem, skip_to_next)
  #append new results to output dataset
  beeCSF_output<-rbind(beeCSF_output,output_id)
  
  skip_to_next <- FALSE
  
  tryCatch(lme(formula(paste(speciesCode, "~spring6SPEI +I(spring6SPEI^2) +I(spring6SPEI^3)")), random=~1|transect,correlation=corAR1(form=~year|transect), data=blue, method="ML"),error = function(e) { skip_to_next <<- TRUE})
  #bind relevant data columns together
  output_id<-cbind(speciesCode, ecosystem, skip_to_next)
  #append new results to output dataset
  beeCSF_output<-rbind(beeCSF_output,output_id)
  
  skip_to_next <- FALSE
  
  tryCatch(lme(formula(paste(speciesCode, "~spring6SPEI +I(spring6SPEI^2) +I(spring6SPEI^3)")), random=~1|transect, correlation=corARMA(form=~year|transect,p=2), data=blue, method="ML"),error = function(e) { skip_to_next <<- TRUE})
  #bind relevant data columns together
  output_id<-cbind(speciesCode, ecosystem, skip_to_next)
  #append new results to output dataset
  beeCSF_output<-rbind(beeCSF_output,output_id)
  
  skip_to_next <- FALSE
  
  tryCatch(lme(formula(paste(speciesCode, "~monsoon6SPEI")), random=~1|transect, data=blue, method="ML"),error = function(e) { skip_to_next <<- TRUE})
  #bind relevant data columns together
  output_id<-cbind(speciesCode, ecosystem, skip_to_next)
  #append new results to output dataset
  beeCSF_output<-rbind(beeCSF_output,output_id)
  
  skip_to_next <- FALSE
  
  tryCatch(lme(formula(paste(speciesCode, "~monsoon6SPEI +I(monsoon6SPEI^2)")), random=~1|transect, data=blue, method="ML"),error = function(e) { skip_to_next <<- TRUE})
  #bind relevant data columns together
  output_id<-cbind(speciesCode, ecosystem, skip_to_next)
  #append new results to output dataset
  beeCSF_output<-rbind(beeCSF_output,output_id)
  
  skip_to_next <- FALSE
  
  
  tryCatch(lme(formula(paste(speciesCode, "~monsoon6SPEI +I(monsoon6SPEI^2) +I(monsoon6SPEI^3)")), random=~1|transect, data=blue, method="ML"),error = function(e) { skip_to_next <<- TRUE})
  #bind relevant data columns together
  output_id<-cbind(speciesCode, ecosystem, skip_to_next)
  #append new results to output dataset
  beeCSF_output<-rbind(beeCSF_output,output_id)
  
  skip_to_next <- FALSE
  
  tryCatch(lme(formula(paste(speciesCode, "~monsoon6SPEI")), correlation=corAR1(form=~year|transect), random=~1|transect, data=blue, method="ML"),error = function(e) { skip_to_next <<- TRUE})
  #bind relevant data columns together
  output_id<-cbind(speciesCode, ecosystem, skip_to_next)
  #append new results to output dataset
  beeCSF_output<-rbind(beeCSF_output,output_id)
  
  skip_to_next <- FALSE
  
  tryCatch(lme(formula(paste(speciesCode, "~monsoon6SPEI")), correlation=corARMA(form=~year|transect,p=2), random=~1|transect, data=blue, method="ML"),error = function(e) { skip_to_next <<- TRUE})
  #bind relevant data columns together
  output_id<-cbind(speciesCode, ecosystem, skip_to_next)
  #append new results to output dataset
  beeCSF_output<-rbind(beeCSF_output,output_id)
  
  skip_to_next <- FALSE
  
  tryCatch(lme(formula(paste(speciesCode, "~monsoon6SPEI + I(monsoon6SPEI^2)")), correlation=corAR1(form=~year|transect), random=~1|transect, data=blue, method="ML"),error = function(e) { skip_to_next <<- TRUE})
  #bind relevant data columns together
  output_id<-cbind(speciesCode, ecosystem, skip_to_next)
  #append new results to output dataset
  beeCSF_output<-rbind(beeCSF_output,output_id)
  
  skip_to_next <- FALSE
  
  tryCatch(lme(formula(paste(speciesCode, "~monsoon6SPEI + I(monsoon6SPEI^2)")), correlation=corARMA(form=~year|transect,p=2), random=~1|transect, data=blue, method="ML"),error = function(e) { skip_to_next <<- TRUE})
  #bind relevant data columns together
  output_id<-cbind(speciesCode, ecosystem, skip_to_next)
  #append new results to output dataset
  beeCSF_output<-rbind(beeCSF_output,output_id)
  
  skip_to_next <- FALSE
  
  tryCatch(lme(formula(paste(speciesCode, "~monsoon6SPEI + I(monsoon6SPEI^2) + I(monsoon6SPEI^3)")), correlation=corAR1(form=~year|transect), random=~1|transect, data=blue, method="ML"),error = function(e) { skip_to_next <<- TRUE})
  #bind relevant data columns together
  output_id<-cbind(speciesCode, ecosystem, skip_to_next)
  #append new results to output dataset
  beeCSF_output<-rbind(beeCSF_output,output_id)
  
  skip_to_next <- FALSE
  
  tryCatch(lme(formula(paste(speciesCode, "~monsoon6SPEI + I(monsoon6SPEI^2) + I(monsoon6SPEI^3)")), correlation=corARMA(form=~year|transect,p=2), random=~1|transect, data=blue, method="ML"),error = function(e) { skip_to_next <<- TRUE})
  #bind relevant data columns together
  output_id<-cbind(speciesCode, ecosystem, skip_to_next)
  #append new results to output dataset
  beeCSF_output<-rbind(beeCSF_output,output_id)
  
  skip_to_next <- FALSE
  
  
  tryCatch(lme(formula(paste(speciesCode, "~spring6SPEI_prioryear")), random=~1|transect, data=blue, method="ML"),error = function(e) { skip_to_next <<- TRUE})
  #bind relevant data columns together
  output_id<-cbind(speciesCode, ecosystem, skip_to_next)
  #append new results to output dataset
  beeCSF_output<-rbind(beeCSF_output,output_id)
  
  skip_to_next <- FALSE
  
  
  tryCatch(lme(formula(paste(speciesCode, "~spring6SPEI_prioryear +I(spring6SPEI_prioryear^2)")), random=~1|transect, data=blue, method="ML"),error = function(e) { skip_to_next <<- TRUE})
  #bind relevant data columns together
  output_id<-cbind(speciesCode, ecosystem, skip_to_next)
  #append new results to output dataset
  beeCSF_output<-rbind(beeCSF_output,output_id)
  
  skip_to_next <- FALSE
  
  tryCatch(lme(formula(paste(speciesCode, "~spring6SPEI_prioryear +I(spring6SPEI_prioryear^2) + I(spring6SPEI_prioryear^3)")), random=~1|transect, data=blue, method="ML") ,error = function(e) { skip_to_next <<- TRUE})
  #bind relevant data columns together
  output_id<-cbind(speciesCode, ecosystem, skip_to_next)
  #append new results to output dataset
  beeCSF_output<-rbind(beeCSF_output,output_id)
  
  skip_to_next <- FALSE
  
  tryCatch(lme(formula(paste(speciesCode, "~spring6SPEI_prioryear")), correlation=corAR1(form=~year|transect), random=~1|transect, data=blue, method="ML"),error = function(e) { skip_to_next <<- TRUE})
  #bind relevant data columns together
  output_id<-cbind(speciesCode, ecosystem, skip_to_next)
  #append new results to output dataset
  beeCSF_output<-rbind(beeCSF_output,output_id)
  
  skip_to_next <- FALSE
  
  tryCatch(lme(formula(paste(speciesCode, "~spring6SPEI_prioryear")), correlation=corARMA(form=~year|transect,p=2), random=~1|transect, data=blue, method="ML"),error = function(e) { skip_to_next <<- TRUE})
  #bind relevant data columns together
  output_id<-cbind(speciesCode, ecosystem, skip_to_next)
  #append new results to output dataset
  beeCSF_output<-rbind(beeCSF_output,output_id)
  
  skip_to_next <- FALSE
  
  tryCatch(lme(formula(paste(speciesCode, "~spring6SPEI_prioryear +I(spring6SPEI_prioryear^2)")), random=~1|transect, correlation=corAR1(form=~year|transect), data=blue, method="ML"),error = function(e) { skip_to_next <<- TRUE})
  #bind relevant data columns together
  output_id<-cbind(speciesCode, ecosystem, skip_to_next)
  #append new results to output dataset
  beeCSF_output<-rbind(beeCSF_output,output_id)
  
  skip_to_next <- FALSE
  
  tryCatch(lme(formula(paste(speciesCode, "~spring6SPEI_prioryear +I(spring6SPEI_prioryear^2)")), random=~1|transect, correlation=corARMA(form=~year|transect,p=2), data=blue, method="ML"),error = function(e) { skip_to_next <<- TRUE})
  #bind relevant data columns together
  output_id<-cbind(speciesCode, ecosystem, skip_to_next)
  #append new results to output dataset
  beeCSF_output<-rbind(beeCSF_output,output_id)
  
  skip_to_next <- FALSE
  
  tryCatch(lme(formula(paste(speciesCode, "~spring6SPEI_prioryear +I(spring6SPEI_prioryear^2) +I(spring6SPEI_prioryear^3)")), random=~1|transect,correlation=corAR1(form=~year|transect), data=blue, method="ML"),error = function(e) { skip_to_next <<- TRUE})
  #bind relevant data columns together
  output_id<-cbind(speciesCode, ecosystem, skip_to_next)
  #append new results to output dataset
  beeCSF_output<-rbind(beeCSF_output,output_id)
  
  skip_to_next <- FALSE
  
  tryCatch(lme(formula(paste(speciesCode, "~spring6SPEI_prioryear +I(spring6SPEI_prioryear^2) +I(spring6SPEI_prioryear^3)")), random=~1|transect, correlation=corARMA(form=~year|transect,p=2), data=blue, method="ML"),error = function(e) { skip_to_next <<- TRUE})
  #bind relevant data columns together
  output_id<-cbind(speciesCode, ecosystem, skip_to_next)
  #append new results to output dataset
  beeCSF_output<-rbind(beeCSF_output,output_id)
  
  skip_to_next <- FALSE
  
  tryCatch(lme(formula(paste(speciesCode, "~monsoon6SPEI_prioryear")), random=~1|transect, data=blue, method="ML"),error = function(e) { skip_to_next <<- TRUE})
  #bind relevant data columns together
  output_id<-cbind(speciesCode, ecosystem, skip_to_next)
  #append new results to output dataset
  beeCSF_output<-rbind(beeCSF_output,output_id)
  
  skip_to_next <- FALSE
  
  tryCatch(lme(formula(paste(speciesCode, "~monsoon6SPEI_prioryear +I(monsoon6SPEI_prioryear^2)")), random=~1|transect, data=blue, method="ML"),error = function(e) { skip_to_next <<- TRUE})
  #bind relevant data columns together
  output_id<-cbind(speciesCode, ecosystem, skip_to_next)
  #append new results to output dataset
  beeCSF_output<-rbind(beeCSF_output,output_id)
  
  skip_to_next <- FALSE
  
  tryCatch(lme(formula(paste(speciesCode, "~monsoon6SPEI_prioryear +I(monsoon6SPEI_prioryear^2) +I(monsoon6SPEI_prioryear^3)")), random=~1|transect, data=blue, method="ML"),error = function(e) { skip_to_next <<- TRUE})
  #bind relevant data columns together
  output_id<-cbind(speciesCode, ecosystem, skip_to_next)
  #append new results to output dataset
  beeCSF_output<-rbind(beeCSF_output,output_id)
  
  skip_to_next <- FALSE
  
  tryCatch(lme(formula(paste(speciesCode, "~monsoon6SPEI_prioryear")), correlation=corAR1(form=~year|transect), random=~1|transect, data=blue, method="ML"),error = function(e) { skip_to_next <<- TRUE})
  #bind relevant data columns together
  output_id<-cbind(speciesCode, ecosystem, skip_to_next)
  #append new results to output dataset
  beeCSF_output<-rbind(beeCSF_output,output_id)
  
  skip_to_next <- FALSE
  
  tryCatch(lme(formula(paste(speciesCode, "~monsoon6SPEI_prioryear")), correlation=corARMA(form=~year|transect,p=2), random=~1|transect, data=blue, method="ML"),error = function(e) { skip_to_next <<- TRUE})
  #bind relevant data columns together
  output_id<-cbind(speciesCode, ecosystem, skip_to_next)
  #append new results to output dataset
  beeCSF_output<-rbind(beeCSF_output,output_id)
  
  skip_to_next <- FALSE
  
  tryCatch(lme(formula(paste(speciesCode, "~monsoon6SPEI_prioryear + I(monsoon6SPEI_prioryear^2)")), correlation=corAR1(form=~year|transect), random=~1|transect, data=blue, method="ML"),error = function(e) { skip_to_next <<- TRUE})
  #bind relevant data columns together
  output_id<-cbind(speciesCode, ecosystem, skip_to_next)
  #append new results to output dataset
  beeCSF_output<-rbind(beeCSF_output,output_id)
  
  skip_to_next <- FALSE
  
  tryCatch(lme(formula(paste(speciesCode, "~monsoon6SPEI_prioryear + I(monsoon6SPEI_prioryear^2)")), correlation=corARMA(form=~year|transect,p=2), random=~1|transect, data=blue, method="ML"),error = function(e) { skip_to_next <<- TRUE})
  #bind relevant data columns together
  output_id<-cbind(speciesCode, ecosystem, skip_to_next)
  #append new results to output dataset
  beeCSF_output<-rbind(beeCSF_output,output_id)
  
  skip_to_next <- FALSE
  
  tryCatch(lme(formula(paste(speciesCode, "~monsoon6SPEI_prioryear + I(monsoon6SPEI_prioryear^2) + I(monsoon6SPEI_prioryear^3)")), correlation=corAR1(form=~year|transect), random=~1|transect, data=blue, method="ML"),error = function(e) { skip_to_next <<- TRUE})
  #bind relevant data columns together
  output_id<-cbind(speciesCode, ecosystem, skip_to_next)
  #append new results to output dataset
  beeCSF_output<-rbind(beeCSF_output,output_id)
  
  skip_to_next <- FALSE
  
  tryCatch(lme(formula(paste(speciesCode, "~monsoon6SPEI_prioryear + I(monsoon6SPEI_prioryear^2) + I(monsoon6SPEI_prioryear^3)")), correlation=corARMA(form=~year|transect,p=2), random=~1|transect, data=blue, method="ML"),error = function(e) { skip_to_next <<- TRUE})
  #bind relevant data columns together
  output_id<-cbind(speciesCode, ecosystem, skip_to_next)
  #append new results to output dataset
  beeCSF_output<-rbind(beeCSF_output,output_id)
  
}

# Create output data frame
beeCSF_output <- data.frame(beeCSF_output)

# Write .csv file of output
# write.csv(beeCSF_output, "CSFs_which_models_run_no2016or2017_2023-12-01.csv",row.names = FALSE)

