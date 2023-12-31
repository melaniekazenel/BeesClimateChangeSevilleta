################################################################################### 
# Code for Select Extended Data Figures

# Heat and desiccation tolerances predict bee abundance under climate change
# Melanie R. Kazenel, Karen W. Wright, Terry Griswold, Kenneth D. Whitney, and Jennifer A. Rudgers

# Date: 2023-12-01
# Corresponding author's email: melanie.kazenel@gmail.com
################################################################################### 


# Load required packages
library(car)
library(dplyr)
library(lubridate)
library(patchwork)
library(ggplot2)
library(plantecophys)
library(zoo)
library(tidyr)
library(grid)
library(lme4)
library(piecewiseSEM)
library(viridis)
library(SPEI)
library(patchwork)


setwd("/Users/mkazenel/Documents/Ph.D./Research/Bees_Climate/Bee_CSFs/Recent_Code/Data Associated with Final Code")


##### EXTENDED DATA FIG. 1 #####
### 1a and 1b ######

# Read in and format Socorro, NM SPEI data from Rudgers et al. 2018 (Ecology)
spei<-read.csv("socorro_1900_speiSept.csv") # SPEI 
cv<-read.csv("socorro_1900_speiCVSept.csv") # CV in SPEI

# Plot inverse SPEI as a function of year
p1<-ggplot(data=spei, aes(x=Year, y=spei6mo*(-1))) +
  geom_point(size=3) +
  stat_smooth(method = "lm", formula = y ~ x) +
  xlab("Year") + 
  ylab("Aridity index") + 
  theme_classic() + theme(axis.text.x = element_text(color="black",size=14)) + theme(axis.text.y = element_text(color="black",size=14)) + theme(axis.title=element_text(size=16))+
  xlim(1900,2015)
p1

# Save the plot
# ggsave("socorro_spei_2023-08-29.jpg", p1,
#        width=5.5,height=4,units = c("in"),
#        dpi = 300)

# Plot CV in SPEI as a function of year
p2<-ggplot(data=cv, aes(x=Year, y=spei6mo_adj5Ycv)) +
  geom_point(size=3) +
  stat_smooth(method = "lm", formula = y ~ x + I(x^2), size = 1) +
  xlab("Year") + 
  ylab("Variance in aridity index (CV)") + 
  theme_classic() + theme(axis.text.x = element_text(color="black",size=14)) + theme(axis.text.y = element_text(color="black",size=14)) + theme(axis.title=element_text(size=16))+
  xlim(1900,2015)
p2

# Save the plot
# ggsave("socorro_spei_cv_2023-08-29.jpg", p2,
#        width=5.5,height=4,units = c("in"),
#        dpi = 300)


### 1c #####

# Read in and format climate data (CanESM2 GCM)
clim_new<-read.csv("spei_historic_and_future_byscenario_CanESM2_2023-12-01.csv")

names(clim_new)[2]<-"scenario"
clim_new$scenario<-as.factor(clim_new$scenario)
levels(clim_new$scenario)[levels(clim_new$scenario)=="SEV-met"] <- "historic"

# Create data frames that assign the year 2020 to the projected future dataset, and bind these to the original climate dataset
rcp2.6<-subset(clim_new,scenario=="historic" & year==2020)
rcp2.6$scenario<-"rcp2.6"

rcp4.5<-subset(clim_new,scenario=="historic" & year==2020)
rcp4.5$scenario<-"rcp4.5"

rcp8.5<-subset(clim_new,scenario=="historic" & year==2020)
rcp8.5$scenario<-"rcp8.5"

clim_new<-bind_rows(clim_new,rcp2.6,rcp4.5,rcp8.5)

# Plot inverse monsoon SPEI as a function of year, for each climate scenario
clim_new_future<-subset(clim_new, scenario!="historic")
clim_new_future$scenario<-as.factor(clim_new_future$scenario)
levels(clim_new_future$scenario) <- c("RCP 2.6","RCP 4.5","RCP 8.5")

p1<-ggplot(data=clim_new_future, aes(x=year,y=monsoon6SPEI*(-1))) + geom_point(size=3,color="black")+ ylab("Aridity index") + xlab("Year") + theme_classic() +theme(axis.text.x = element_text(size=12))+ theme(axis.text.y = element_text(size=12)) + theme(legend.text=element_text(size=15), legend.title=element_text(size=14))+ theme(axis.title.x = element_text(size=18))+ theme(axis.title.y = element_text(size=18)) + geom_smooth(method="lm",formula=y~x,color="blue") + facet_wrap(~scenario) + theme(panel.spacing.x = unit(7, "mm")) + theme(strip.text.x = element_text(size = 16)) + ggtitle("CanESM2") + theme(plot.title = element_text(size = 20)) + theme(strip.background = element_blank(),strip.text.x = element_blank()) + theme(axis.title.x=element_blank(), axis.text.x = element_blank())
p1

# Read in and format climate data (ACCESS 1.0 GCM)
clim_new<-read.csv("spei_historic_and_future_byscenario_ACCESS1-0_2023-12-01.csv")

names(clim_new)[2]<-"scenario"
clim_new$scenario<-as.factor(clim_new$scenario)
levels(clim_new$scenario)[levels(clim_new$scenario)=="SEV-met"] <- "historic"

# Create data frames that assign the year 2020 to the projected future dataset, and bind these to the original climate dataset
rcp2.6<-subset(clim_new,scenario=="historic" & year==2020)
rcp2.6$scenario<-"rcp2.6"

rcp4.5<-subset(clim_new,scenario=="historic" & year==2020)
rcp4.5$scenario<-"rcp4.5"

rcp8.5<-subset(clim_new,scenario=="historic" & year==2020)
rcp8.5$scenario<-"rcp8.5"

clim_new<-bind_rows(clim_new,rcp2.6,rcp4.5,rcp8.5)

# Plot inverse monsoon SPEI as a function of year, for each climate scenario
clim_new_future<-subset(clim_new, scenario!="historic")
clim_new_future$scenario<-as.factor(clim_new_future$scenario)
levels(clim_new_future$scenario) <- c("RCP 2.6","RCP 4.5","RCP 8.5")

p2<-ggplot(data=clim_new_future, aes(x=year,y=monsoon6SPEI*(-1))) + geom_point(size=3,color="black")+ ylab("Aridity index") + xlab("Year") + theme_classic() +theme(axis.text.x = element_text(size=12))+ theme(axis.text.y = element_text(size=12)) + theme(legend.text=element_text(size=15), legend.title=element_text(size=14))+ theme(axis.title.x = element_text(size=18))+ theme(axis.title.y = element_text(size=18)) + geom_smooth(method="lm",formula=y~x,color="blue") + facet_wrap(~scenario) + theme(panel.spacing.x = unit(7, "mm")) + theme(strip.text.x = element_text(size = 16)) + ggtitle("ACCESS 1.0") + theme(plot.title = element_text(size = 20)) + theme(axis.title.x=element_blank(), axis.text.x = element_blank()) + theme(strip.background = element_blank(),strip.text.x = element_blank())
p2

# Read in and format climate data (CCSM4 GCM)
clim_new<-read.csv("spei_historic_and_future_byscenario_CCSM4_2023-12-01.csv")

names(clim_new)[2]<-"scenario"
clim_new$scenario<-as.factor(clim_new$scenario)
levels(clim_new$scenario)[levels(clim_new$scenario)=="SEV-met"] <- "historic"

# Create data frames that assign the year 2020 to the projected future dataset, and bind these to the original climate dataset
rcp2.6<-subset(clim_new,scenario=="historic" & year==2020)
rcp2.6$scenario<-"rcp2.6"

rcp4.5<-subset(clim_new,scenario=="historic" & year==2020)
rcp4.5$scenario<-"rcp4.5"

rcp8.5<-subset(clim_new,scenario=="historic" & year==2020)
rcp8.5$scenario<-"rcp8.5"

clim_new<-bind_rows(clim_new,rcp2.6,rcp4.5,rcp8.5)

# Plot inverse monsoon SPEI as a function of year, for each climate scenario
clim_new_future<-subset(clim_new, scenario!="historic")
clim_new_future$scenario<-as.factor(clim_new_future$scenario)
levels(clim_new_future$scenario) <- c("RCP 2.6","RCP 4.5","RCP 8.5")

p3<-ggplot(data=clim_new_future, aes(x=year,y=monsoon6SPEI*(-1))) + geom_point(size=3,color="black")+ ylab("Aridity index") + xlab("Year") + theme_classic() +theme(axis.text.x = element_text(size=12))+ theme(axis.text.y = element_text(size=12)) + theme(legend.text=element_text(size=15), legend.title=element_text(size=14))+ theme(axis.title.x = element_text(size=18))+ theme(axis.title.y = element_text(size=18)) + geom_smooth(method="lm",formula=y~x,color="blue") + facet_wrap(~scenario) + theme(panel.spacing.x = unit(7, "mm")) + theme(strip.text.x = element_text(size = 16)) + ggtitle("CCSM 4.0") + theme(plot.title = element_text(size = 20)) + theme(axis.title.x=element_blank(), axis.text.x = element_blank()) + theme(strip.background = element_blank(),strip.text.x = element_blank())
p3

# Read in and format climate data (CNRM GCM)
clim_new<-read.csv("spei_historic_and_future_byscenario_CNRM-CM5_2023-12-01.csv")

names(clim_new)[2]<-"scenario"
clim_new$scenario<-as.factor(clim_new$scenario)
levels(clim_new$scenario)[levels(clim_new$scenario)=="SEV-met"] <- "historic"

# Create data frames that assign the year 2020 to the projected future dataset, and bind these to the original climate dataset
rcp2.6<-subset(clim_new,scenario=="historic" & year==2020)
rcp2.6$scenario<-"rcp2.6"

rcp4.5<-subset(clim_new,scenario=="historic" & year==2020)
rcp4.5$scenario<-"rcp4.5"

rcp8.5<-subset(clim_new,scenario=="historic" & year==2020)
rcp8.5$scenario<-"rcp8.5"

clim_new<-bind_rows(clim_new,rcp2.6,rcp4.5,rcp8.5)

# Plot inverse monsoon SPEI as a function of year, for each climate scenario
clim_new_future<-subset(clim_new, scenario!="historic")
clim_new_future$scenario<-as.factor(clim_new_future$scenario)
levels(clim_new_future$scenario) <- c("RCP 2.6","RCP 4.5","RCP 8.5")

p4<-ggplot(data=clim_new_future, aes(x=year,y=monsoon6SPEI*(-1))) + geom_point(size=3,color="black")+ ylab("Aridity index") + xlab("Year") + theme_classic() +theme(axis.text.x = element_text(size=12))+ theme(axis.text.y = element_text(size=12)) + theme(legend.text=element_text(size=15), legend.title=element_text(size=14))+ theme(axis.title.x = element_text(size=18))+ theme(axis.title.y = element_text(size=18)) + geom_smooth(method="lm",formula=y~x,color="blue") + facet_wrap(~scenario) + theme(panel.spacing.x = unit(7, "mm")) + theme(strip.text.x = element_text(size = 16)) + ggtitle("CNRM-CM5") + theme(plot.title = element_text(size = 20)) + theme(axis.title.x=element_blank(), axis.text.x = element_blank()) + theme(strip.background = element_blank(),strip.text.x = element_blank())
p4

# Read in and format climate data (CSIRO GCM)
clim_new<-read.csv("spei_historic_and_future_byscenario_CSIRO-Mk3-6-0_2023-12-01.csv")

names(clim_new)[2]<-"scenario"
clim_new$scenario<-as.factor(clim_new$scenario)
levels(clim_new$scenario)[levels(clim_new$scenario)=="SEV-met"] <- "historic"

# Create data frames that assign the year 2020 to the projected future dataset, and bind these to the original climate dataset
rcp2.6<-subset(clim_new,scenario=="historic" & year==2020)
rcp2.6$scenario<-"rcp2.6"

rcp4.5<-subset(clim_new,scenario=="historic" & year==2020)
rcp4.5$scenario<-"rcp4.5"

rcp8.5<-subset(clim_new,scenario=="historic" & year==2020)
rcp8.5$scenario<-"rcp8.5"

clim_new<-bind_rows(clim_new,rcp2.6,rcp4.5,rcp8.5)

# Plot inverse monsoon SPEI as a function of year, for each climate scenario
clim_new_future<-subset(clim_new, scenario!="historic")
clim_new_future$scenario<-as.factor(clim_new_future$scenario)
levels(clim_new_future$scenario) <- c("RCP 2.6","RCP 4.5","RCP 8.5")

p5<-ggplot(data=clim_new_future, aes(x=year,y=monsoon6SPEI*(-1))) + geom_point(size=3,color="black")+ ylab("Aridity index") + xlab("Year") + theme_classic() +theme(axis.text.x = element_text(size=12))+ theme(axis.text.y = element_text(size=12)) + theme(legend.text=element_text(size=15), legend.title=element_text(size=14))+ theme(axis.title.x = element_text(size=18))+ theme(axis.title.y = element_text(size=18)) + geom_smooth(method="lm",formula=y~x,color="blue") + facet_wrap(~scenario) + theme(panel.spacing.x = unit(7, "mm")) + theme(strip.text.x = element_text(size = 16)) + ggtitle("CSIRO-Mk3.6.0") + theme(plot.title = element_text(size = 20)) + theme(axis.title.x=element_blank(), axis.text.x = element_blank()) + theme(strip.background = element_blank(),strip.text.x = element_blank())
p5

# Read in and format climate data (INM GCM)
clim_new<-read.csv("spei_historic_and_future_byscenario_INM-CM4_2023-12-01.csv")

names(clim_new)[2]<-"scenario"
clim_new$scenario<-as.factor(clim_new$scenario)
levels(clim_new$scenario)[levels(clim_new$scenario)=="SEV-met"] <- "historic"

# Create data frames that assign the year 2020 to the projected future dataset, and bind these to the original climate dataset
rcp2.6<-subset(clim_new,scenario=="historic" & year==2020)
rcp2.6$scenario<-"rcp2.6"

rcp4.5<-subset(clim_new,scenario=="historic" & year==2020)
rcp4.5$scenario<-"rcp4.5"

rcp8.5<-subset(clim_new,scenario=="historic" & year==2020)
rcp8.5$scenario<-"rcp8.5"

clim_new<-bind_rows(clim_new,rcp2.6,rcp4.5,rcp8.5)

# Plot inverse monsoon SPEI as a function of year, for each climate scenario
clim_new_future<-subset(clim_new, scenario!="historic")
clim_new_future$scenario<-as.factor(clim_new_future$scenario)
levels(clim_new_future$scenario) <- c("RCP 2.6","RCP 4.5","RCP 8.5")

p6<-ggplot(data=clim_new_future, aes(x=year,y=monsoon6SPEI*(-1))) + geom_point(size=3,color="black")+ ylab("Aridity index") + xlab("Year") + theme_classic() +theme(axis.text.x = element_text(size=12))+ theme(axis.text.y = element_text(size=12)) + theme(legend.text=element_text(size=15), legend.title=element_text(size=14))+ theme(axis.title.x = element_text(size=18))+ theme(axis.title.y = element_text(size=18)) + geom_smooth(method="lm",formula=y~x,color="blue") + facet_wrap(~scenario) + theme(panel.spacing.x = unit(7, "mm")) + theme(strip.text.x = element_text(size = 16)) + ggtitle("INM") + theme(plot.title = element_text(size = 20)) + theme(strip.background = element_blank(),strip.text.x = element_blank())
p6

# Combine graphs into a multi-figure plot
p<-p1 + p2 + p3 +p4 + p5 + p6 + plot_layout(ncol = 1, guides = "collect")
p

# save plot
#ggsave("future_monsoon_spei_2023-12-01.jpg", p, width=8, height=12, units="in", dpi=300)


##### EXTENDED DATA FIG. 2 #####

# t-tests and summary statistics for table of results related to climate differences between the Chihuahuan Desert and Plains ecosystems in the month of maximum difference between ecosystems

# maximum temperature
clim<-read.csv("maximum_temperature_largest_difference.csv")
t.test(clim$Station_49, clim$Station_50, paired = TRUE, alternative = "two.sided")

clim_long<-pivot_longer(data=clim, cols=3:4, names_to="station", values_to="climate_value")
clim_long %>% group_by(station) %>% summarise(mean=mean(climate_value), se=sd(climate_value)/sqrt(n()))


# mean temperature
clim<-read.csv("mean_temperature_largest_difference.csv")
t.test(clim$Station_49, clim$Station_50, paired = TRUE, alternative = "two.sided")

clim_long<-pivot_longer(data=clim, cols=3:4, names_to="station", values_to="climate_value")
clim_long %>% group_by(station) %>% summarise(mean=mean(climate_value), se=sd(climate_value)/sqrt(n()))


# minimum temperature
clim<-read.csv("minimum_temperature_largest_difference.csv")
t.test(clim$Station_49, clim$Station_50, paired = TRUE, alternative = "two.sided")

clim_long<-pivot_longer(data=clim, cols=3:4, names_to="station", values_to="climate_value")
clim_long %>% group_by(station) %>% summarise(mean=mean(climate_value), se=sd(climate_value)/sqrt(n()))


# precipitation
clim<-read.csv("precipitation_largest_difference.csv")
t.test(clim$Station_49, clim$Station_50, paired = TRUE, alternative = "two.sided")

clim_long<-pivot_longer(data=clim, cols=3:4, names_to="station", values_to="climate_value")
clim_long %>% group_by(station) %>% summarise(mean=mean(climate_value), se=sd(climate_value)/sqrt(n()))


##### EXTENDED DATA FIG. 4a #####

# Read in historic monthly meteorological data, and filter to just include years and weather stations of interest
met_monthly <- read.csv("SevMet_monthly_12Mar21.csv")
met_monthly <- subset(met_monthly, Station_Name == "Five Points" | Station_Name == "Blue Grama")

# Calculate rolling 6-month max temperature and mean precipitation
met_monthly$airt_max6mo<-rollmax(x=met_monthly$airtemp_max,k=6,fill=NA,align="right")
met_monthly$precip_mean6mo<-rollmean(x=met_monthly$precip,k=6,fill=NA,align="right")

# Subset the data to just include spring and monsoon values
met_subset<-subset(met_monthly,month==5|month==9)

# Select the columns of interest
met_subset<-met_subset[,c(2,7,8,19,20)]

# Rename columns and recode values
names(met_subset)[c(1,3)]<-c("station","season")
met_subset$season<-as.factor(met_subset$season)
met_subset$season <- recode(met_subset$season, "5" = "spring", "9" = "monsoon")

# Transform data from long to wide format
met_subset_wide <- met_subset %>% pivot_wider(id_cols=c(year,station),names_from=season,values_from=c(airt_max6mo,precip_mean6mo))

# Read in projected future climate data from ClimateNA for a given GCM (two files per GCM)
rcp4.5<-read.csv("climatedata_ACCESS1-0_RCP45_2011-2100MP.csv")
rcp8.5<-read.csv("climatedata_ACCESS1-0_RCP85_2011-2100MP.csv")
# rcp4.5<-read.csv("climatedata_CCSM4_RCP45_2011-2100MP.csv")
# rcp8.5<-read.csv("climatedata_CCSM4_RCP85_2011-2100MP.csv")
# rcp4.5<-read.csv("climatedata_CNRM-CM5_RCP45_2011-2100MP.csv")
# rcp8.5<-read.csv("climatedata_CNRM-CM5_RCP85_2011-2100MP.csv")
# rcp4.5<-read.csv("climatedata_CSIRO-Mk3-6-0_RCP45_2011-2100MP.csv")
# rcp8.5<-read.csv("climatedata_CSIRO-Mk3-6-0_RCP85_2011-2100MP.csv")
# rcp4.5<-read.csv("climatedata_INM-CM4_RCP45_2011-2100MP.csv")
# rcp8.5<-read.csv("climatedata_INM-CM4_RCP85_2011-2100MP.csv")

# Add a scenario column to each data frame, and combine the data frames
rcp4.5$scenario <- "rcp4.5"
rcp8.5$scenario <- "rcp8.5"
all_data<-bind_rows(rcp4.5,rcp8.5)

# Create new data frame containing only the columns of interest
all_data_2 <- data_frame(Year=all_data$Year,all_data[,4:6],all_data[,7:18],all_data[,43:54],scenario=all_data$scenario)

# Create separate data frames for temperature and precipitation variables
temp <- data_frame(all_data_2[,1:16],all_data_2[,29])
precip <- data_frame(all_data_2[,1:4],all_data_2[,17:28],all_data_2[,29])

# Transform data from long to wide format
temp_long <- temp %>% pivot_longer(cols=5:16,names_to="temp_month")
precip_long <- precip %>% pivot_longer(cols=5:16,names_to="precip_month")

# Create month columns in each data frame
temp_long$month <- substr(temp_long$temp_month,5,6)
precip_long$month <- substr(precip_long$precip_month,4,5)

# Create variable columns in each data frame
temp_long$variable <- "Tmax"
precip_long$variable <- "PPT"

# Remove unneeded columns
temp_long$temp_month <- NULL
precip_long$precip_month <- NULL

# Combine the temperature and precipitation data frames
temp_precip_long <- bind_rows(temp_long,precip_long)

# Tranform the data from long to wide format
temp_precip_wide <- temp_precip_long %>% pivot_wider(id_cols=c(Year,Latitude,Longitude,Elevation,month,scenario),names_from=variable,values_from=value)

### For RCP 4.5: ###
# Create data frame of just the RCP 4.5 data
rcp4.5<-subset(temp_precip_wide,scenario=="rcp4.5")

# Calculate rolling 6-month max temperature and mean precipitation
rcp4.5$airt_max6mo<-rollmax(x=rcp4.5$Tmax,k=6,fill=NA,align="right")
rcp4.5$precip_mean6mo<-rollmean(x=rcp4.5$PPT,k=6,fill=NA,align="right")

# Subset the data to just include spring and monsoon values
rcp4.5_subset<-subset(rcp4.5,month=="05"|month=="09")

# Select the columns of interest
rcp4.5_subset<-rcp4.5_subset[,c(1,5,6,9,10)]

# Rename columns and recode values
names(rcp4.5_subset)[c(1,2)]<-c("year","season")
rcp4.5_subset$season<-as.factor(rcp4.5_subset$season)
rcp4.5_subset$season <- recode(rcp4.5_subset$season, "05" = "spring", "09" = "monsoon")

# Transform data from long to wide format
rcp4.5_subset_wide <- rcp4.5_subset %>% pivot_wider(id_cols=c(year,scenario),names_from=season,values_from=c(airt_max6mo,precip_mean6mo))

### For RCP 8.5: ###
# Create data frame of just the RCP 8.5 data
rcp8.5<-subset(temp_precip_wide,scenario=="rcp8.5")

# Calculate rolling 6-month max temperature and mean precipitation
rcp8.5$airt_max6mo<-rollmax(x=rcp8.5$Tmax,k=6,fill=NA,align="right")
rcp8.5$precip_mean6mo<-rollmean(x=rcp8.5$PPT,k=6,fill=NA,align="right")

# Subset the data to just include spring and monsoon values
rcp8.5_subset<-subset(rcp8.5,month=="05"|month=="09")

# Select the columns of interest
rcp8.5_subset<-rcp8.5_subset[,c(1,5,6,9,10)]

# Rename columns and recode values
names(rcp8.5_subset)[c(1,2)]<-c("year","season")
rcp8.5_subset$season<-as.factor(rcp8.5_subset$season)
rcp8.5_subset$season <- recode(rcp8.5_subset$season, "05" = "spring", "09" = "monsoon")

# Transform data from long to wide format
rcp8.5_subset_wide <- rcp8.5_subset %>% pivot_wider(id_cols=c(year,scenario),names_from=season,values_from=c(airt_max6mo,precip_mean6mo))

# Combine the data frames for each scenario, and subset to include just the years after 2020
future<-bind_rows(rcp4.5_subset_wide,rcp8.5_subset_wide)
future<-subset(future,year>2020)

# Rename column and create new station column
names(future)[2]<-"source"
future$station<-"both"

# Read in SPEI data for the focal GCM
spei<-read.csv("spei_historic_and_future_byscenario_ACCESS1-0_2023-12-01.csv")
#spei<-read.csv("spei_historic_and_future_byscenario_CCSM4_2023-12-01.csv")
#spei<-read.csv("spei_historic_and_future_byscenario_CNRM-CM5_2023-12-01.csv")
#spei<-read.csv("spei_historic_and_future_byscenario_CSIRO-Mk3-6-0_2023-12-01.csv")
#spei<-read.csv("spei_historic_and_future_byscenario_INM-CM4_2023-12-01.csv")

# Subset SPEI data to only include historic records, and join the historic temperature/precipitation and SPEI data frames
spei_past<-subset(spei,year<=2020)
past<-left_join(met_subset_wide,spei_past,by=c("year","station"))

# Subset the SPEI data to only include projected future values, for only one weather station, since the future values are the same for both stations (we projected future climate values for the midpoint between them)
spei_future<-subset(spei,year>2020 & station=="Five Points")

# Remove unneeded column, and join the future temperature/precipitation and SPEI data frames
spei_future$station<-NULL
future<-left_join(future,spei_future,by=c("year","source"))

# Select focal columns from the past climate data frame, and select focal years of data
colnames(past)
past2<-past[,c(1,2,7,3:6,8:10)]
colnames(past2)
past2<-subset(past2,year>=2002)

# Select focal columns from the future climate data frame
colnames(future)
future2<-future[,c(1,7,2,3:6,8:10)]
colnames(future2)

# Combine the focal past and future data frames
all<-bind_rows(past2,future2)

# write .csv file of data
#write.csv(all, "data_extdatafig4a_ACCESS1-0_2023-12-01.csv", row.names=FALSE)
#write.csv(all, "data_extdatafig4a_CCSM4_2023-12-01.csv", row.names=FALSE)
#write.csv(all, "data_extdatafig4a_CNRM-CM5_2023-12-01.csv", row.names=FALSE)
#write.csv(all, "data_extdatafig4a_CSIRO-Mk3-6-0_2023-12-01.csv", row.names=FALSE)
#write.csv(all, "data_extdatafig4a_INM_2023-12-01.csv", row.names=FALSE)

# Read in data from each GCM, and add a column indicating the GCM to each data frame
all1<-read.csv("data_extdatafig4a_CanESM2_2023-12-01.csv")
all1$gcm<-"CanESM2"
all2<-read.csv("data_extdatafig4a_ACCESS1-0_2023-12-01.csv")
all2$gcm<-"ACCESS1-0"
all3<-read.csv("data_extdatafig4a_CCSM4_2023-12-01.csv")
all3$gcm<-"CCSM4"
all4<-read.csv("data_extdatafig4a_CNRM-CM5_2023-12-01.csv")
all4$gcm<-"CNRM"
all5<-read.csv("data_extdatafig4a_CSIRO-Mk3-6-0_2023-12-01.csv")
all5$gcm<-"CSIRO"
all6<-read.csv("data_extdatafig4a_INM_2023-12-01.csv")
all6$gcm<-"INM"

# Combine data frames
all<-bind_rows(all1,all2,all3,all4,all5,all6)

# Create a new data frame and recode the source column in preparation for graphing
forgraph<-all
forgraph$source<-as.factor(forgraph$source)
levels(forgraph$source)<-c("RCP 2.6","RCP 4.5","RCP 8.5","Historic")

# Graph the relationship between inverse SPEI and maximum air temperature for the period leading up to the monsoon season
(p1<-ggplot(aes(x=airt_max6mo_monsoon,y=monsoon6SPEI*(-1)),data=forgraph) + 
    geom_point(aes(colour=source)) + 
    labs(color="Data source") +
    geom_smooth(method="lm", color="black",size=0.75) +
    xlab("Maximum air temperature (degrees C)") +
    ylab("Aridity index"))

# Add a red bar and asterisk to the graph indicating the the critical thermal maximum (CTMax) of the least thermally tolerant bee taxon in the dataset
gtext = textGrob("*", y = -.05, gp = gpar(col = "red"))
gline = linesGrob(y = c(-.02, .02),  gp = gpar(col = "red", lwd = 2)) 
p1 = p1 + annotation_custom(gtext, xmin=39.09, xmax=39.09, ymin=-Inf, ymax=Inf) +
  annotation_custom(gline, xmin=39.09, xmax=39.09, ymin=-Inf, ymax=Inf)
g = ggplotGrob(p1)
g$layout$clip[g$layout$name=="panel"] <- "off"
grid.draw(g)

# Save the plot
#ggsave("spei-by-max-air-temp_2023-12-29.jpg", g, width=5,height=4,units = c("in"), dpi = 300)

# Model the graphed relationship (linear regression)
lm(monsoon6SPEI*(-1) ~ airt_max6mo_monsoon, data=all) 
# equation of the line: spei = 0.3979*maxairtemp -14.2091 

# Use the equation above to calculate the SPEI value corresponding with the CTMax value of the least thermally tolerant bee taxon in the dataset
0.3979*39.09 -14.2091 # output: 1.344811

# ACCESS 1.0 GCM: For each projected future climate scenario, calculate the percentage of future years for which monsoon conditions will exceed the thermal maximum of the least thermally tolerant bee taxon 
future4.5<-filter(all,source=="rcp4.5" & gcm=="ACCESS1-0" & monsoon6SPEI*(-1)>=1.344811)
nrow(future4.5)/80 # 0.2625
future8.5<-filter(all,source=="rcp8.5" & monsoon6SPEI*(-1)>=1.344811)
nrow(future8.5)/80

# CanESM2 GCM: For each projected future climate scenario, calculate the percentage of future years for which monsoon conditions will exceed the thermal maximum of the least thermally tolerant bee taxon 
future4.5<-filter(all,source=="rcp4.5" & gcm=="CanESM2" & monsoon6SPEI*(-1)>=1.344811)
nrow(future4.5)/80 # 0.2375
future8.5<-filter(all,source=="rcp8.5" & monsoon6SPEI*(-1)>=1.344811)
nrow(future8.5)/80

# CCSM4 GCM: For each projected future climate scenario, calculate the percentage of future years for which monsoon conditions will exceed the thermal maximum of the least thermally tolerant bee taxon 
future4.5<-filter(all,source=="rcp4.5" & gcm=="CCSM4" & monsoon6SPEI*(-1)>=1.344811)
nrow(future4.5)/80 # 0.1375
future8.5<-filter(all,source=="rcp8.5" & monsoon6SPEI*(-1)>=1.344811)
nrow(future8.5)/80

# CNRM GCM: For each projected future climate scenario, calculate the percentage of future years for which monsoon conditions will exceed the thermal maximum of the least thermally tolerant bee taxon 
future4.5<-filter(all,source=="rcp4.5" & gcm=="CNRM" & monsoon6SPEI*(-1)>=1.344811)
nrow(future4.5)/80 # 0.075
future8.5<-filter(all,source=="rcp8.5" & monsoon6SPEI*(-1)>=1.344811)
nrow(future8.5)/80

# CSIRO GCM: For each projected future climate scenario, calculate the percentage of future years for which monsoon conditions will exceed the thermal maximum of the least thermally tolerant bee taxon 
future4.5<-filter(all,source=="rcp4.5" & gcm=="CSIRO" & monsoon6SPEI*(-1)>=1.344811)
nrow(future4.5)/80 # 0.2375
future8.5<-filter(all,source=="rcp8.5" & monsoon6SPEI*(-1)>=1.344811)
nrow(future8.5)/80

# INM GCM: For each projected future climate scenario, calculate the percentage of future years for which monsoon conditions will exceed the thermal maximum of the least thermally tolerant bee taxon 
future4.5<-filter(all,source=="rcp4.5" & gcm=="INM" & monsoon6SPEI*(-1)>=1.344811)
nrow(future4.5)/80 # 0.05
future8.5<-filter(all,source=="rcp8.5" & monsoon6SPEI*(-1)>=1.344811)
nrow(future8.5)/80

# Calculate the mean percentage of future years across GCMs
mean(c(0.2625,0.2375,0.1375,0.075,0.2375,0.05))

##### EXTENDED DATA FIG. 7 #####

# read in body mass data
data <- read.csv("SEVBeeBodyMassData_2023-08-29_forpub.csv")

# subset body mass data to just include focal species
data2 <- subset(data, 
                code=="MEOSMWAT"|
                  code=="APANTLES"|
                  code=="HAHALLIG"|
                  code=="APDIARIN"|
                  code=="APDIAAUS"|
                  code=="HAHALTRI"|
                  code=="HAAGAANG"|
                  code=="APANTAFF"|
                  code=="HALASHUD"|
                  code=="HALASDEL"|
                  code=="APANTMON"|
                  code=="HALASSIS"|
                  code=="MEASHMEL"|
                  code=="ANPERCAL"|
                  code=="HALASSEM"|
                  code=="HALASDIA")

# calculate mean mass for each species x year combination
summary<-data2 %>% group_by(code,genus,species,year) %>% summarise(mean_mass_mg=mean(mass_mg),se_mass_mg=sd(mass_mg)/sqrt(n()))

# create a genus_species column
summary$genus_species<-paste(summary$genus,summary$species,sep=" ")

# edit genus_species name for one species
summary$genus_species[summary$genus_species=="Lasioglossum A"] <- "Lasioglossum sp. A"


# plot body mass as a function of year
p<-ggplot(data=summary, aes(x=year, y=mean_mass_mg)) +
  geom_point(size=1.5) +
  facet_wrap(~genus_species, scales = "free_y") +
  geom_errorbar(aes(ymin=mean_mass_mg-se_mass_mg, ymax=mean_mass_mg+se_mass_mg), width=.1)+
  ylab("Mean body mass (mg)") +
  xlab("Year")+
  theme(axis.text.x = element_text(color="black",size=11)) +
  theme(axis.text.y = element_text(color="black",size=11)) +
  theme(axis.title=element_text(size=15)) +
  theme(panel.spacing.x = unit(2, "mm")) +
  theme(strip.text.x = element_text(size = 10, face = "italic"))
p

# save the plot
# ggsave("with_species_mass_change_2023-08-29.jpg", p,
#        width=11,height=6,units = c("in"),
#        dpi = 300)


# format data to enable running linear regressions
summary2<-data2 %>% group_by(code,year) %>% summarise(mean_mass_mg=mean(mass_mg))
summary2_wide<-pivot_wider(data=summary2,names_from = code, values_from = mean_mass_mg)

# linear regression models: body mass as a function of year, for each species
summary(lm(HAAGAANG~year,data=summary2_wide))
summary(lm(APANTAFF~year,data=summary2_wide))
summary(lm(APANTLES~year,data=summary2_wide))
summary(lm(APANTMON~year,data=summary2_wide))
summary(lm(MEASHMEL~year,data=summary2_wide))
summary(lm(APDIAAUS~year,data=summary2_wide))
summary(lm(APDIARIN~year,data=summary2_wide))
summary(lm(HAHALLIG~year,data=summary2_wide))
summary(lm(HAHALTRI~year,data=summary2_wide))
summary(lm(HALASDEL~year,data=summary2_wide))
summary(lm(HALASHUD~year,data=summary2_wide))
summary(lm(HALASSEM~year,data=summary2_wide))
summary(lm(HALASSIS~year,data=summary2_wide))
summary(lm(HALASDIA~year,data=summary2_wide))
summary(lm(MEOSMWAT~year,data=summary2_wide))
summary(lm(ANPERCAL~year,data=summary2_wide))


##### EXTENDED DATA FIG. 8 #####

# read in climate data, plant phenology data, and plant species list
clim<-read.csv("spei_historic_and_future_byscenario_CanESM2_2023-12-01.csv")
phen<-read.csv("sev137_plant_phenology.csv")
plantlist<-read.csv("sev051_plantspecieslist_20140818.txt")

# subset plant phenology dataset to just include relevant years, months, sites, and plant taxa
phen<-subset(phen,site=="B" | site=="C" | site=="G")
phen$year<-substr(phen$date,1,4)
phen<-subset(phen,year>2000)
phen<-subset(phen,month>2 & month<11)
names(phen)[4]<-"species_code"
phen<-left_join(phen,plantlist,by="species_code")
phen<-subset(phen,life_form!="G") # remove grasses

# calculate proportion of forb individuals flowering on each transect (web)
# spring
spring<-subset(phen,month<7)

spring_summary<-spring %>% group_by(year,site,web) %>% summarise(total_obs=n())

spring_flow<-subset(spring,reproduction=="FL" | reproduction=="FF" | reproduction=="BFR" | reproduction=="B")
flow_summary_spring<-spring_flow %>% group_by(year,site,web) %>% summarise(indivs_flowering=n())

spring_summary2<-left_join(spring_summary,flow_summary_spring,by=c("year","site","web"))
spring_summary2$indivs_flowering[is.na(spring_summary2$indivs_flowering)] <- 0
spring_summary2$prop_flowering_spring<-spring_summary2$indivs_flowering/spring_summary2$total_obs

summary(spring_summary2$prop_flowering_spring)

# monsoon
monsoon<-subset(phen,month>=7)

monsoon_summary<-monsoon %>% group_by(year,site,web) %>% summarise(total_obs=n())

monsoon_flow<-subset(monsoon,reproduction=="FL" | reproduction=="FF" | reproduction=="BFR" | reproduction=="B")
flow_summary_monsoon<-monsoon_flow %>% group_by(year,site,web) %>% summarise(indivs_flowering=n())

monsoon_summary2<-left_join(monsoon_summary,flow_summary_monsoon,by=c("year","site","web"))
monsoon_summary2$indivs_flowering[is.na(monsoon_summary2$indivs_flowering)] <- 0
monsoon_summary2$prop_flowering_monsoon<-monsoon_summary2$indivs_flowering/monsoon_summary2$total_obs

summary(monsoon_summary2$prop_flowering_monsoon)

# join data frames
summary_all<-left_join(spring_summary2[,c(1:3,6)],monsoon_summary2[,c(1:3,6)],by=c("year","site","web"))

# add transect column
summary_all$transect<-paste(summary_all$site,summary_all$web,sep="")
summary_all2<-summary_all[,c(1,2,6,4,5)]
names(summary_all2)[2]<-"ecosystem"

# write csv
#write.csv(summary_all2, "prop_flowering_by_season_2021-12-22.csv",row.names=FALSE)

# add weather station column
summary_all2$station<-ifelse(summary_all2$ecosystem=="B","Blue Grama", "Five Points")

# change year to numeric
summary_all2$year<-as.numeric(summary_all2$year)

# join climate and phenology data frames
summary_all3<-left_join(summary_all2,clim[,c(1,3,4,7)],by=c("year","station"))

# format data in preparation for graphing
plotdata<-summary_all3
plotdata$ecosystem<-as.factor(plotdata$ecosystem)
levels(plotdata$ecosystem) <- c("Plains grassland","Chihuahuan Desert shrubland","Chihuahuan Desert grassland")
plotdata$ecosystem <- factor(plotdata$ecosystem, levels = c("Plains grassland","Chihuahuan Desert grassland","Chihuahuan Desert shrubland"))

# graph: flowering as a function of inverse spring SPEI
p<-ggplot(data=plotdata, aes(x=spring6SPEI*(-1),y=prop_flowering_spring)) + geom_point(size=3,color="black")+ geom_smooth(method = "lm") + facet_wrap(~ecosystem) + xlab("Spring aridity") + ylab("Proportion of forb and shrub \nindividuals in flower")
p

# save plot
# ggsave("plant_phenology_spring_spei_2023-12-01.jpg", p,
#        width=7,height=3,units = c("in"),
#        dpi = 300)

# subset data by ecosystem to run statistics
B<-subset(summary_all3, ecosystem=="B")
C<-subset(summary_all3, ecosystem=="C")
G<-subset(summary_all3, ecosystem=="G")

# statistics: flowering as a function of inverse spring SPEI
m1<-lmer(prop_flowering_spring~spring6SPEI*(-1)+(1|transect),data=B,na.action=na.omit)
Anova(m1,type = 3)
rsquared(m1)

m1<-lmer(prop_flowering_spring~spring6SPEI*(-1)+(1|transect),data=G,na.action=na.omit)
Anova(m1,type = 3)
rsquared(m1)

m1<-lmer(prop_flowering_spring~spring6SPEI*(-1)+(1|transect),data=C,na.action=na.omit)
Anova(m1,type = 3)
rsquared(m1)

# graph: flowering as a function of inverse monsoon SPEI
p<-ggplot(data=plotdata, aes(x=monsoon6SPEI*(-1),y=prop_flowering_monsoon)) + geom_point(size=3,color="black")+ geom_smooth(method = "lm") + facet_wrap(~ecosystem) + xlab("Monsoon aridity") + ylab("Proportion of forb and shrub \nindividuals in flower")
p

# save graph
# ggsave("plant_phenology_monsoon_spei_2023-12-01.jpg", p,
#        width=7,height=3,units = c("in"),
#        dpi = 300)

# statistics: flowering as a function of inverse monsoon SPEI
m1<-lmer(prop_flowering_monsoon~monsoon6SPEI*(-1)+(1|transect),data=B,na.action=na.omit)
Anova(m1,type = 3)
rsquared(m1)

m1<-lmer(prop_flowering_monsoon~monsoon6SPEI*(-1)+(1|transect),data=G,na.action=na.omit)
Anova(m1,type = 3)
rsquared(m1)

m1<-lmer(prop_flowering_monsoon~monsoon6SPEI*(-1)+(1|transect),data=C,na.action=na.omit)
Anova(m1,type = 3)
rsquared(m1)


##### EXTENDED DATA FIG. 9 #####

# Read in historic monthly meteorological data, and subset to just include focal weather stations and years
met_monthly <- read.csv("SevMet_monthly_12Mar21.csv")
met_monthly <- filter(met_monthly, Station_Name == "Five Points" | Station_Name == "Blue Grama")
met_monthly <- filter(met_monthly, year >= 2002 & year <= 2019)

# Create month_year column, and format month column
met_monthly$month_year<-ymd(paste(met_monthly$year,met_monthly$month,sep="-","01"))
met_monthly$month<-as.factor(met_monthly$month)

# FOR MANUSCRIPT TEXT: Calculate mean monthly precipitation and temperature across years for each weather station
met_monthly %>% group_by(sta) %>% summarise(meanprecip=mean(precip),meantemp=mean(airtemp_mean))

# For each weather station, calculate monthly mean temperature and mean total precipitation, averaging across years
monthly_avg<-met_monthly%>%group_by(month,sta)%>%summarise(mean_temp=mean(airtemp_mean),mean_totalprecip=mean(precip))

# Graph mean temperature as a function of month for the Plains ecosystem, averaging across years
(bluegrama_temp <- monthly_avg %>%
    filter(sta == "50") %>% 
    ggplot() +
    geom_col(aes(x = month, y= mean_temp,fill=month),color="black",linewidth=0.1) +
    labs(title = "Plains ecosystem",
         x = "Month",
         y = "Temperature (C)") +
    scale_fill_viridis(discrete=TRUE) +
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
)

# Graph mean temperature as a function of month for the Chihuahuan Desert ecosystems, averaging across years
(fivepoints_temp <- monthly_avg %>%
    filter(sta == "49") %>% 
    ggplot() +
    geom_col(aes(x = month, y= mean_temp, fill=month),color="black",linewidth=0.1) +
    labs(title = "Chihuahuan Desert ecosystems",
         x = "Month",
         y = "Temperature (C)") +
    scale_fill_viridis(discrete=TRUE) +
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
)

# Graph mean total precipitation as a function of month for the Plains ecosystem, averaging across years
(bluegrama_precip <- monthly_avg %>%
    filter(sta == "50") %>% 
    ggplot() +
    geom_col(aes(x = month, y= mean_totalprecip,fill=month),color="black",linewidth=0.1) +
    labs(title = "Plains ecosystem",
         x = "Month",
         y = "Precipitation (mm)") +
    scale_fill_viridis(discrete=TRUE) +
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
)

# Graph mean total precipitation as a function of month for the Chihuahuan Desert ecosystems, averaging across years
(fivepoints_precip <- monthly_avg %>%
    filter(sta == "49") %>% 
    ggplot() +
    geom_col(aes(x = month, y= mean_totalprecip, fill=month),color="black",linewidth=0.1) +
    labs(title = "Chihuahuan Desert ecosystems",
         x = "Month",
         y = "Precipitation (mm)") +
    scale_fill_viridis(discrete=TRUE) +
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
)

# Graph change in mean monthly temperature from 2002–2019, for the Plains ecosystem
(bluegrama_temperature2 <- met_monthly %>% 
    filter(sta == "50") %>% 
    ggplot() +
    geom_col(aes(x = month_year, y= airtemp_mean,fill=month),color="black",linewidth=0.1) +
    labs(title = "Plains ecosystem",
         x = "Year",
         y = "Temperature (C)") +
    scale_fill_viridis(discrete=TRUE) +
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
)

# Graph change in mean monthly precipitation from 2002–2019, for the Plains ecosystem
(bluegrama_precipitation2 <- met_monthly %>% 
    filter(sta == "50") %>% 
    ggplot() +
    geom_col(aes(x = month_year, y= precip,fill=month),color="black",linewidth=0.1) +
    labs(title = "Plains ecosystem",
         x = "Year",
         y = "Precipitation (mm)") +
    scale_fill_viridis(discrete=TRUE) +
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
)

# Graph change in mean monthly temperature from 2002–2019, for the Chihuahuan Desert ecosystems
(fivepoints_temperature2 <- met_monthly %>% 
    filter(sta == "49") %>% 
    ggplot() +
    geom_col(aes(x = month_year, y= airtemp_mean,fill=month),color="black",linewidth=0.1) +
    labs(title = "Chihuahuan Desert ecosystems",
         x = "Year",
         y = "Temperature (C)") +
    scale_fill_viridis(discrete=TRUE) +
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
)

# Graph change in mean monthly precipitation from 2002–2019, for the Chihuahuan Desert ecosystems
(fivepoints_precipitation2 <- met_monthly %>% 
    filter(sta == "49") %>% 
    ggplot() +
    geom_col(aes(x = month_year, y= precip,fill=month),color="black",linewidth=0.1) +
    labs(title = "Chihuahuan Desert ecosystems",
         x = "Year",
         y = "Precipitation (mm)") +
    theme_minimal()+
    scale_fill_viridis(discrete=TRUE) +
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
)

# Combine graphs into a multi-figure plot
p1<-bluegrama_temp + fivepoints_temp + bluegrama_precip +fivepoints_precip + bluegrama_temperature2 + fivepoints_temperature2 + bluegrama_precipitation2 +fivepoints_precipitation2 + plot_layout(ncol = 2, guides = "collect")
p1

# Save the plot
# ggsave("metstationtempprecip_2023-08-29.jpg", p1, width=10,height=12,units = c("in"),dpi = 600)


##### EXTENDED DATA FIG. 10 #####

# Read in historic monthly meteorological data, and subset to just include focal weather stations and years
met_monthly <- read.csv("SevMet_monthly_12Mar21.csv") %>% 
  mutate(sta = as.factor(sta))
met_monthly <- filter(met_monthly, Station_Name == "Five Points" | Station_Name == "Blue Grama")
met_monthly <- filter(met_monthly, year >= 2002 & year <= 2019)

# Calculate vapor pressure deficit
met_monthly$vpd<-RHtoVPD(RH=met_monthly$RH_mean,TdegC=met_monthly$airtemp_mean)

# Rename variables
names(met_monthly)[c(9,12,18)] <- c("airt","rh","ppt")

# Create data frame of relevant columns for the Chihuahuan Desert ecosystems
met49 <- met_monthly %>% 
  filter(sta == "49") %>% 
  select(sta, year, month, airt, ppt, rh, vpd)

# Create data frame of relevant columns for the Plains ecosystem
met50 <- met_monthly %>% 
  filter(sta == "50") %>% 
  select(sta, year, month, airt, ppt, rh, vpd)

# Write function to calculate potential evapotranspiration (PET) and climatic water balance for each month x year combination, using data from a given weather station
thornthwaite_by_sta <- function(data, station, latitude) {
  calc_for_sta <- data %>% 
    filter(sta == station) %>% 
    mutate(pet = thornthwaite(Tave  = airt, # Thornthwaite method used in PET calculation
                              lat   = latitude),
           balance = ppt - pet)
  
  return(calc_for_sta)
}

# Run the function for both weather stations
met49_thornth <- thornthwaite_by_sta(met49, "49", 34.335)
met50_thornth <- thornthwaite_by_sta(met50, "50", 34.335)

# Write function to calculate SPEI for each month x year combination, using data from a given weather station
spei_by_sta <- function(data) {
  spei_calc <- data %>% 
    mutate(spei1  = as.vector(spei(scale =  1, na.rm = TRUE, data = .$balance)$fitted)) 
  
  return(spei_calc)
}

# Run the function for both weather stations
met49_spei <- spei_by_sta(met49_thornth)
met50_spei <- spei_by_sta(met50_thornth)

# Sort the data by date
met49_spei <- met49_spei %>% arrange(year,month)
met50_spei <- met50_spei %>% arrange(year,month)

# For the Chihuahuan Desert ecosystems, plot monthly inverse SPEI as a function of temperature
(p1<-ggplot(aes(x=airt,y=spei1*(-1)),data=met49_spei) + 
    geom_point() + facet_wrap(~month) + 
    geom_smooth(method="lm", color="black",size=0.75) +
    xlab("Air temperature (degrees C)") +
    ylab("Aridity index") +
    labs(title = "Chihuahuan Desert ecosystems") +
    theme_bw())

# For the Plains ecosystem, plot monthly inverse SPEI as a function of temperature
(p2<-ggplot(aes(x=airt,y=spei1*(-1)),data=met50_spei) + 
    geom_point() + facet_wrap(~month) + 
    geom_smooth(method="lm", color="black",size=0.75) +
    xlab("Air temperature (degrees C)") +
    ylab("Aridity index") +
    labs(title = "Plains ecosystem") +
    theme_bw())

# For the Chihuahuan Desert ecosystems, plot monthly inverse SPEI as a function of precipitation
(p3<-ggplot(aes(x=ppt,y=spei1*(-1)),data=met49_spei) + 
    geom_point() + facet_wrap(~month) + 
    geom_smooth(method="lm", color="black",size=0.75) +
    xlab("Precipitation (mm)") +
    ylab("Aridity index") +
    theme_bw())

# For the Plains ecosystem, plot monthly inverse SPEI as a function of precipitation
(p4<-ggplot(aes(x=ppt,y=spei1*(-1)),data=met50_spei) + 
    geom_point() + facet_wrap(~month) + 
    geom_smooth(method="lm", color="black",size=0.75) +
    xlab("Precipitation (mm)") +
    ylab("Aridity index") +
    theme_bw())

# For the Chihuahuan Desert ecosystems, plot monthly inverse SPEI as a function of relative humidity
(p5<-ggplot(aes(x=rh,y=spei1*(-1)),data=met49_spei) + 
    geom_point() + facet_wrap(~month) + 
    geom_smooth(method="lm", color="black",size=0.75) +
    xlab("Relative humidity (%)") +
    ylab("Aridity index") +
    theme_bw())

# For the Plains ecosystem, plot monthly inverse SPEI as a function of relative humidity
(p6<-ggplot(aes(x=rh,y=spei1*(-1)),data=met50_spei) + 
    geom_point() + facet_wrap(~month) + 
    geom_smooth(method="lm", color="black",size=0.75) +
    xlab("Relative humidity (%)") +
    ylab("Aridity index") +
    theme_bw())

# For the Chihuahuan Desert ecosystems, plot monthly inverse SPEI as a function of vapor pressure deficit
(p7<-ggplot(aes(x=vpd,y=spei1*(-1)),data=met49_spei) + 
    geom_point() + facet_wrap(~month) + 
    geom_smooth(method="lm", color="black",size=0.75) +
    ylab("Aridity index") +
    xlab("Vapor pressure deficit (kPa)") +
    theme_bw())

# For the Plains ecosystem, plot monthly inverse SPEI as a function of vapor pressure deficit
(p8<-ggplot(aes(x=vpd,y=spei1*(-1)),data=met50_spei) + 
    geom_point() + facet_wrap(~month) + 
    geom_smooth(method="lm", color="black",size=0.75) +
    ylab("Aridity index") +
    xlab("Vapor pressure deficit (kPa)") +
    theme_bw())

# Combine graphs into a multi-figure plot
full<-p1 + p2 + p3 + p4 + p5 + p6 + p7 + p8 + plot_layout(ncol = 2)
full

# Save the plot
# ggsave("aridity-by-other-climate-variables_2023-08-29.jpg", full,
#       width=10,height=20,units = c("in"),dpi = 600)

