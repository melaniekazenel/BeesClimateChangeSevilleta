################################################################################### 
# Aridity Index (SPEI) Calculations
# Historic Climate Data and Predicted Future Climate Data from the CanESM2 GCM

# Heat and desiccation tolerances predict bee abundance under climate change
# Melanie R. Kazenel, Karen W. Wright, Terry Griswold, Kenneth D. Whitney, and Jennifer A. Rudgers

# Date: 2023-08-29
# Corresponding author's email: melanie.kazenel@gmail.com
################################################################################### 

# Load required packages
library(SPEI)
library(dplyr)
library(tidyr)
library(zoo)


##### Initial Calculations: Historic Climate Data #####

# Read in historic monthly meteorological data, and subset to just include focal weather stations
met_monthly <- read.csv("SevMet_monthly_12Mar21.csv")
met_monthly <- subset(met_monthly, Station_Name == "Five Points" | Station_Name == "Blue Grama")

# Calculate potential evapotranspiration using the Thornthwaite method
met_monthly$PET <- thornthwaite(Tave = met_monthly$airtemp_mean, lat = 34.33532, na.rm = FALSE)

# Create column of climatic water balance values
met_monthly$CWB <- met_monthly$precip - met_monthly$PET

# Create a new data frame containing only columns of interest, and rename one column to facilitate later steps
met_data<-met_monthly[,c("Station_Name","year","month","CWB")]
names(met_data)[1]<-"station_or_scenario"


##### Initial Calculations: CanESM2 Projected Future Climate Data #####

# Read in projected future climate data from ClimateNA
rcp2.6<-read.csv("climatedata_CanESM2_RCP2.6_2020-12-22.csv")
rcp4.5<-read.csv("climatedata_CanESM2_RCP4.5_2020-12-22.csv")
rcp8.5<-read.csv("climatedata_CanESM2_RCP8.5_2020-12-22.csv")

# Add a scenario column to each data frame, and combine the data frames
rcp2.6$scenario <- "rcp2.6"
rcp4.5$scenario <- "rcp4.5"
rcp8.5$scenario <- "rcp8.5"
all_data<-bind_rows(rcp2.6,rcp4.5,rcp8.5)

# Create a new data frame containing only the desired columns
all_data_2 <- data_frame(Year=all_data$Year,all_data[,4:6],all_data[,31:54],scenario=all_data$scenario)

# Create separate data frames of temperature and precipitation data
temp <- data_frame(all_data_2[,1:16],all_data_2[,29])
precip <- data_frame(all_data_2[,1:4],all_data_2[,17:28],all_data_2[,29])

# Convert the data from wide to long format
temp_long <- temp %>% pivot_longer(cols=5:16,names_to="temp_month")
precip_long <- precip %>% pivot_longer(cols=5:16,names_to="precip_month")

# Create month columns
temp_long$month <- substr(temp_long$temp_month,5,6)
precip_long$month <- substr(precip_long$precip_month,4,5)

# Add columns indicating the variable of interest
temp_long$variable <- "Tave"
precip_long$variable <- "PPT"

# Remove unneeded columns
temp_long$temp_month <- NULL
precip_long$precip_month <- NULL

# Combine long-form temperature and precipitation data frames
temp_precip_long <- bind_rows(temp_long,precip_long)

# Convert the combined data frame back to wide format
temp_precip_wide <- temp_precip_long %>% pivot_wider(id_cols=c(Year,Latitude,Longitude,Elevation,month,scenario),names_from=variable,values_from=value)

# Calculate potential evapotranspiration using the Thornthwaite method
temp_precip_wide$PET <- thornthwaite(Tave = temp_precip_wide$Tave, lat = 34.33532, na.rm = FALSE)

# Create column of climatic water balance values
temp_precip_wide$CWB <- temp_precip_wide$PPT - temp_precip_wide$PET

# Create a new data frame of relevant columns
combined_data <- data_frame(station_or_scenario=temp_precip_wide$scenario,year=temp_precip_wide$Year,month=temp_precip_wide$month,CWB=temp_precip_wide$CWB)

# Treat select columns as numeric data
combined_data$month<-as.numeric(combined_data$month)
combined_data$CWB<-as.numeric(combined_data$CWB)
met_data$CWB<-as.numeric(met_data$CWB)


##### Aridity Index (SPEI) Calculations: Combined Historic and Projected Future Climate Data #####

# Merge historic and projected future climate data frames
combined_data <- bind_rows(met_data,combined_data)

# Subset to just include years of interest
combined_data <- subset(combined_data,year>1999)

# Create separate data frames for each site x future scenario combination, including the historic data in each data frame

# Five Points weather station, RCP 2.6
rcp2.6<-subset(combined_data,station_or_scenario=="rcp2.6" & year>2020) #future data
fp<-subset(combined_data,station_or_scenario=="Five Points") # historic data
fp2.6<-bind_rows(fp,rcp2.6) # combine the two
fp2.6$dataset<-"fp2.6" # add column indicating dataset

# Five Points weather station, RCP 4.5
rcp4.5<-subset(combined_data,station_or_scenario=="rcp4.5" & year>2020) # future data
fp<-subset(combined_data,station_or_scenario=="Five Points") # historic data
fp4.5<-bind_rows(fp,rcp4.5) # combine the two
fp4.5$dataset<-"fp4.5" # add column indicating dataset

# Five Points weather station, RCP 8.5
rcp8.5<-subset(combined_data,station_or_scenario=="rcp8.5" & year>2020) # future data
fp<-subset(combined_data,station_or_scenario=="Five Points") # historic data
fp8.5<-bind_rows(fp,rcp8.5) # combine the two
fp8.5$dataset<-"fp8.5" # add column indicating dataset

# Blue Grama weather station, RCP 2.6
rcp2.6<-subset(combined_data,station_or_scenario=="rcp2.6" & year>2020) # future data
bg<-subset(combined_data,station_or_scenario=="Blue Grama") # historic data
b2.6<-bind_rows(bg,rcp2.6) # combine the two
b2.6$dataset<-"b2.6" # add column indicating dataset

# Blue Grama weather station, RCP 4.5
rcp4.5<-subset(combined_data,station_or_scenario=="rcp4.5" & year>2020) # future data
bg<-subset(combined_data,station_or_scenario=="Blue Grama")  # historic data
b4.5<-bind_rows(bg,rcp4.5) # combine the two
b4.5$dataset<-"b4.5" # add column indicating dataset

# Blue Grama weather station, RCP 8.5
rcp8.5<-subset(combined_data,station_or_scenario=="rcp8.5" & year>2020) # future data
bg<-subset(combined_data,station_or_scenario=="Blue Grama")  # historic data
b8.5<-bind_rows(bg,rcp8.5) # combine the two
b8.5$dataset<-"b8.5" # add column indicating dataset

# Create data frame just for historic climate values
historic<-subset(combined_data,station_or_scenario=="Five Points" | station_or_scenario=="Blue Grama")
historic$dataset<-"historic" # add column indicating dataset
historic <- historic %>% arrange(station_or_scenario,year,month) # sort the data frame

# Combine all data frames created above
combined_cwb<-bind_rows(historic,b2.6,b4.5,b8.5,fp2.6,fp4.5,fp8.5)

# Calculate 6-month SPEI for each month of each year in the historic and future data, using 2002-2020 data from the Blue Grama weather station as a reference period.
SPEI_6mo <- spei(data=ts(combined_cwb$CWB,freq=12,start=c(2002,1)), scale=6, ref.start = c(2002,1), ref.end = c(2020,12)) 

# Add the calculated SPEI values to the data frame
combined_cwb$SPEI_6mo<-SPEI_6mo$fitted

# Add NAs for dummy-calculated values at the beginning of each site x scenario time series
combined_cwb$SPEI_6mo <- ifelse(combined_cwb$station_or_scenario=="Five Points" & combined_cwb$year==2000 & combined_cwb$month %in% 1:5 | combined_cwb$station_or_scenario=="Blue Grama" & combined_cwb$year==2002 & combined_cwb$month %in% 1:5, NA, combined_cwb$SPEI_6mo)


##### Subset and Format the Data #####

# Divide the data frame back into separate datasets for each weather station x scenario combination
fp2.6<-subset(combined_cwb,dataset=="fp2.6")
fp4.5<-subset(combined_cwb,dataset=="fp4.5")
fp8.5<-subset(combined_cwb,dataset=="fp8.5")
b2.6<-subset(combined_cwb,dataset=="b2.6")
b4.5<-subset(combined_cwb,dataset=="b4.5")
b8.5<-subset(combined_cwb,dataset=="b8.5")

# Format each data frame

# Five Points, RCP 2.6
# subset to just include spring and fall SPEI values
spei_subset<-subset(fp2.6,month==5|month==9)
# rename month to season, and recode the variable
names(spei_subset)[3]<-"season"
spei_subset$season<-as.factor(spei_subset$season)
spei_subset$season <- recode(spei_subset$season, "5" = "spring", "9" = "monsoon")
# convert data from long to wide format
spei_subset_wide <- spei_subset %>% pivot_wider(id_cols=c(year,station_or_scenario),names_from=season,values_from=SPEI_6mo)
# rename SPEI columns
names(spei_subset_wide)[3:4]<-c("spring6SPEI","monsoon6SPEI")
# add columns representing the previous year's SPEI values
spei_subset_wide_nextyear <- spei_subset_wide
spei_subset_wide_nextyear$next_year<- spei_subset_wide$year + 1
names(spei_subset_wide_nextyear)[3:4]<-c("spring6SPEI_prioryear","monsoon6SPEI_prioryear")
spei_subset_wide_nextyear$year<-NULL
names(spei_subset_wide_nextyear)[4]<-"year"
fp2.6_spei_subset_wide<-left_join(spei_subset_wide,spei_subset_wide_nextyear[,2:4],by=c("year"))

# Five Points, RCP 4.5
# subset to just include spring and fall SPEI values
spei_subset<-subset(fp4.5,month==5|month==9)
# rename month to season, and recode the variable
names(spei_subset)[3]<-"season"
spei_subset$season<-as.factor(spei_subset$season)
spei_subset$season <- recode(spei_subset$season, "5" = "spring", "9" = "monsoon")
# convert data from long to wide format
spei_subset_wide <- spei_subset %>% pivot_wider(id_cols=c(year,station_or_scenario),names_from=season,values_from=SPEI_6mo)
# rename SPEI columns
names(spei_subset_wide)[3:4]<-c("spring6SPEI","monsoon6SPEI")
# add columns representing the previous year's SPEI values
spei_subset_wide_nextyear <- spei_subset_wide 
spei_subset_wide_nextyear$next_year<- spei_subset_wide$year + 1
names(spei_subset_wide_nextyear)[3:4]<-c("spring6SPEI_prioryear","monsoon6SPEI_prioryear")
spei_subset_wide_nextyear$year<-NULL
names(spei_subset_wide_nextyear)[4]<-"year"
fp4.5_spei_subset_wide<-left_join(spei_subset_wide,spei_subset_wide_nextyear[,2:4],by=c("year"))

# Five Points, RCP 8.5
# subset to just include spring and fall SPEI values
spei_subset<-subset(fp8.5,month==5|month==9)
# rename month to season, and recode the variable
names(spei_subset)[3]<-"season"
spei_subset$season<-as.factor(spei_subset$season)
spei_subset$season <- recode(spei_subset$season, "5" = "spring", "9" = "monsoon")
# convert data from long to wide format
spei_subset_wide <- spei_subset %>% pivot_wider(id_cols=c(year,station_or_scenario),names_from=season,values_from=SPEI_6mo)
# rename SPEI columns
names(spei_subset_wide)[3:4]<-c("spring6SPEI","monsoon6SPEI")
# add columns representing the previous year's SPEI values
spei_subset_wide_nextyear <- spei_subset_wide
spei_subset_wide_nextyear$next_year<- spei_subset_wide$year + 1
names(spei_subset_wide_nextyear)[3:4]<-c("spring6SPEI_prioryear","monsoon6SPEI_prioryear")
spei_subset_wide_nextyear$year<-NULL
names(spei_subset_wide_nextyear)[4]<-"year"
fp8.5_spei_subset_wide<-left_join(spei_subset_wide,spei_subset_wide_nextyear[,2:4],by=c("year"))

# Blue Grama, RCP 2.6
# subset to just include spring and fall SPEI values
spei_subset<-subset(b2.6,month==5|month==9)
# rename month to season, and recode the variable
names(spei_subset)[3]<-"season"
spei_subset$season<-as.factor(spei_subset$season)
spei_subset$season <- recode(spei_subset$season, "5" = "spring", "9" = "monsoon")
# convert data from long to wide format
spei_subset_wide <- spei_subset %>% pivot_wider(id_cols=c(year,station_or_scenario),names_from=season,values_from=SPEI_6mo)
# rename SPEI columns
names(spei_subset_wide)[3:4]<-c("spring6SPEI","monsoon6SPEI")
# add columns representing the previous year's SPEI values
spei_subset_wide_nextyear <- spei_subset_wide
spei_subset_wide_nextyear$next_year<- spei_subset_wide$year + 1
names(spei_subset_wide_nextyear)[3:4]<-c("spring6SPEI_prioryear","monsoon6SPEI_prioryear")
spei_subset_wide_nextyear$year<-NULL
names(spei_subset_wide_nextyear)[4]<-"year"
b2.6_spei_subset_wide<-left_join(spei_subset_wide,spei_subset_wide_nextyear[,2:4],by=c("year"))

# Blue Grama, RCP 4.5
# subset to just include spring and fall SPEI values
spei_subset<-subset(b4.5,month==5|month==9)
# rename month to season, and recode the variable
names(spei_subset)[3]<-"season"
spei_subset$season<-as.factor(spei_subset$season)
spei_subset$season <- recode(spei_subset$season, "5" = "spring", "9" = "monsoon")
# convert data from long to wide format
spei_subset_wide <- spei_subset %>% pivot_wider(id_cols=c(year,station_or_scenario),names_from=season,values_from=SPEI_6mo)
# rename SPEI columns
names(spei_subset_wide)[3:4]<-c("spring6SPEI","monsoon6SPEI")
# add columns representing the previous year's SPEI values
spei_subset_wide_nextyear <- spei_subset_wide
spei_subset_wide_nextyear$next_year<- spei_subset_wide$year + 1
names(spei_subset_wide_nextyear)[3:4]<-c("spring6SPEI_prioryear","monsoon6SPEI_prioryear")
spei_subset_wide_nextyear$year<-NULL
names(spei_subset_wide_nextyear)[4]<-"year"
b4.5_spei_subset_wide<-left_join(spei_subset_wide,spei_subset_wide_nextyear[,2:4],by=c("year"))

# Blue Grama, RCP 8.5
# subset to just include spring and fall SPEI values
spei_subset<-subset(b8.5,month==5|month==9)
# rename month to season, and recode the variable
names(spei_subset)[3]<-"season"
spei_subset$season<-as.factor(spei_subset$season)
spei_subset$season <- recode(spei_subset$season, "5" = "spring", "9" = "monsoon")
# convert data from long to wide format
spei_subset_wide <- spei_subset %>% pivot_wider(id_cols=c(year,station_or_scenario),names_from=season,values_from=SPEI_6mo)
# rename SPEI columns
names(spei_subset_wide)[3:4]<-c("spring6SPEI","monsoon6SPEI")
# add columns representing the previous year's SPEI values
spei_subset_wide_nextyear <- spei_subset_wide
spei_subset_wide_nextyear$next_year<- spei_subset_wide$year + 1
names(spei_subset_wide_nextyear)[3:4]<-c("spring6SPEI_prioryear","monsoon6SPEI_prioryear")
spei_subset_wide_nextyear$year<-NULL
names(spei_subset_wide_nextyear)[4]<-"year"
b8.5_spei_subset_wide<-left_join(spei_subset_wide,spei_subset_wide_nextyear[,2:4],by=c("year"))


# Format the data frames above in preparation for merging them
fp2.6_spei_subset_wide$station<-"Five Points"
names(fp2.6_spei_subset_wide)[2]<-"source"
fp2.6_spei_subset_wide$source <- recode(fp2.6_spei_subset_wide$source, "Five Points" = "SEV-met")

fp4.5_spei_subset_wide$station<-"Five Points"
names(fp4.5_spei_subset_wide)[2]<-"source"
fp4.5_spei_subset_wide$source <- recode(fp4.5_spei_subset_wide$source, "Five Points" = "SEV-met")

fp8.5_spei_subset_wide$station<-"Five Points"
names(fp8.5_spei_subset_wide)[2]<-"source"
fp8.5_spei_subset_wide$source <- recode(fp8.5_spei_subset_wide$source, "Five Points" = "SEV-met")

b2.6_spei_subset_wide$station<-"Blue Grama"
names(b2.6_spei_subset_wide)[2]<-"source"
b2.6_spei_subset_wide$source <- recode(b2.6_spei_subset_wide$source, "Blue Grama" = "SEV-met")

b4.5_spei_subset_wide$station<-"Blue Grama"
names(b4.5_spei_subset_wide)[2]<-"source"
b4.5_spei_subset_wide$source <- recode(b4.5_spei_subset_wide$source, "Blue Grama" = "SEV-met")

b8.5_spei_subset_wide$station<-"Blue Grama"
names(b8.5_spei_subset_wide)[2]<-"source"
b8.5_spei_subset_wide$source <- recode(b8.5_spei_subset_wide$source, "Blue Grama" = "SEV-met")

# Remove unneeded rows
fp4.5_spei_subset_wide<-subset(fp4.5_spei_subset_wide, source!="SEV-met")
fp8.5_spei_subset_wide<-subset(fp8.5_spei_subset_wide, source!="SEV-met")

b4.5_spei_subset_wide<-subset(b4.5_spei_subset_wide, source!="SEV-met")
b8.5_spei_subset_wide<-subset(b8.5_spei_subset_wide, source!="SEV-met")

# Combine the data frames
combined_data_subset <- bind_rows(fp2.6_spei_subset_wide,fp4.5_spei_subset_wide,fp8.5_spei_subset_wide,b2.6_spei_subset_wide,b4.5_spei_subset_wide,b8.5_spei_subset_wide)


# Write a .csv file
# write.csv(combined_data_subset,"spei_historic_and_future_byscenario_CanESM2_2023-08-29.csv",row.names=FALSE)

