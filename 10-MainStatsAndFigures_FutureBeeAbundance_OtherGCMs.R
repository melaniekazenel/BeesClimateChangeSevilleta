################################################################################### 
# Statistics related to predicted future bee community change, using data from GCMs other than CanESM2 

# Heat and desiccation tolerances predict bee abundance under climate change
# Melanie R. Kazenel, Karen W. Wright, Terry Griswold, Kenneth D. Whitney, and Jennifer A. Rudgers

# Date: 2023-08-29
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



##### PREDICTED FUTURE CHANGE IN INDIVIDUAL SPECIES ABUNDANCES #####
## Data manipulation #####

## Format predicted future bee abundance data ##

# read in data frame of projected future bee data, completing this step and running the code below separately for each GCM
abund<-read.csv("predicted_abundances_ACCESS1-0-2023-08-29.csv")
# abund<-read.csv("predicted_abundances_CCSM4-2023-08-29.csv")
# abund<-read.csv("predicted_abundances_CNRM-CM5-2023-08-29.csv")
# abund<-read.csv("predicted_abundances_CSIRO-Mk3-6-0-2023-08-29.csv")
# abund<-read.csv("predicted_abundances_INM-CM4-2023-08-29.csv")

# remove NA and Inf values
abund<-subset(abund,predicted_max_abundance_per_transect>-Inf & predicted_max_abundance_per_transect<Inf)
abund$spei_value<-as.numeric(abund$spei_value)
abund<-na.omit(abund)

# remove outlier bee abundance values by excluding records that are more than five times the maximum number of bees of a single species recorded for a given transect x year combination 
# max. number recorded = 822 individuals; 822*5 = 4110
abund_updated<-subset(abund,predicted_max_abundance_per_transect<=4110)

# create list of species for which we have predicted future abundance data (i.e., species for which CSFs ran)
species_future<-unique(abund_updated$code)


## Format historic bee abundance data ##

# read in data
abund_hist_wide<-read.csv("bee_wide_year_2002-2019_no2016or2017_maxabund_2023-08-29.csv")

# add "scenario" column
abund_hist_wide <- abund_hist_wide %>%  mutate(scenario="historic", .after = monsoon6SPEI_prioryear)

# pivot to long form
abund_historic_melt<-pivot_longer(data=abund_hist_wide, cols = 11:349, names_to = "code", values_to = "abund")

# subset to just include species in predicted future abundance dataset (i.e., species for which CSFs ran)
abund_historic_2_melt <- abund_historic_melt %>% filter(code %in% species_future)


## Combine historic and predicted future data frames ##

# select columns and format data frames in preparation for combining them
abund_historic_final<-abund_historic_2_melt[,c(2,3,10,1,11,12)]
abund_future<-abund_updated[,c(2,3,5:6,1,4)]
names(abund_future)[6]<-"abund"

# combine data frames
abund_combined<-bind_rows(abund_future,abund_historic_final)

# add station column
abund_combined$station<-abund_combined$ecosystem

# rename levels of the "station" variable
abund_combined$station<-as.factor(abund_combined$station)
levels(abund_combined$station)[levels(abund_combined$station)=="B"] <- "Blue Grama"
levels(abund_combined$station)[levels(abund_combined$station)=="C"] <- "Five Points"
levels(abund_combined$station)[levels(abund_combined$station)=="G"] <- "Five Points"


## Add climate data ##

# read in SPEI data, completing this step and running the code below separately for each GCM
clim_new<-read.csv("spei_historic_and_future_byscenario_ACCESS1-0_2023-08-29.csv")
# clim_new<-read.csv("spei_historic_and_future_byscenario_CCSM4_2023-08-29.csv")
# clim_new<-read.csv("spei_historic_and_future_byscenario_CNRM-CM5_2023-08-29.csv")
# clim_new<-read.csv("spei_historic_and_future_byscenario_CSIRO-Mk3-6-0_2023-08-29.csv")
# clim_new<-read.csv("spei_historic_and_future_byscenario_INM-CM4_2023-08-29.csv")

# create data frame for historic data, and add climate data to it
colnames(abund_combined)
colnames(clim_new)
historic<-subset(abund_combined,year<=2020)
clim_historic<-subset(clim_new,source=="SEV-met")
historic2<-left_join(historic,clim_historic,by=c("year","station"))
historic2$source<-NULL

# create data frame for future bee data, and add climate data to it
future<-subset(abund_combined,year>2020)
clim_future<-subset(clim_new,source!="SEV-met")
names(clim_future)[2]<-"scenario"
future2<-left_join(future,clim_future,by=c("year","scenario","station"))

# combine the data frames
colnames(historic2)
colnames(future2)
abund_combined_2<-bind_rows(historic2,future2)

# sort the resulting data frame
abund_combined_2<-abund_combined_2[order(abund_combined_2$code, abund_combined_2$ecosystem,abund_combined_2$transect,abund_combined_2$scenario,abund_combined_2$year),]


## Loop: Which species are predicted to increase, decrease, or remain stable in abundance over time? #####

# create a vector of species codes
species<-levels(as.factor(abund_combined_2$code))

# create a vector of scenario names
scenarios <- c("rcp4.5","rcp8.5")

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
#write.csv(specieslist2,"slopes_predicted_future_change_ACCESS1-0_2023-08-29.csv",row.names = FALSE)
#write.csv(specieslist2,"slopes_predicted_future_change_CCSM4_2023-08-29.csv",row.names = FALSE)
#write.csv(specieslist2,"slopes_predicted_future_change_CNRM-CM5_2023-08-29.csv",row.names = FALSE)
#write.csv(specieslist2,"slopes_predicted_future_change_CSIRO-Mk3-6-0_2023-08-29.csv",row.names = FALSE)
#write.csv(specieslist2,"slopes_predicted_future_change_INM-CM4_2023-08-29.csv",row.names = FALSE)


## Summary statistics: Which species predicted to increase, decrease, or remain stable in abundance over time? #####

# Read in slopes data from calculations above
slopes <- read.csv("slopes_predicted_future_change_ACCESS1-0_2023-08-29.csv")
#slopes <- read.csv("slopes_predicted_future_change_CCSM4_2023-08-29.csv")
#slopes <- read.csv("slopes_predicted_future_change_CNRM-CM5_2023-08-29.csv")
#slopes <- read.csv("slopes_predicted_future_change_CSIRO-Mk3-6-0_2023-08-29.csv")
#slopes <- read.csv("slopes_predicted_future_change_INM-CM4_2023-08-29.csv")

# RCP 4.5
# create separate data frames for species predicted to increase, decrease, or not change in abundance
decrease_4 <- subset(slopes, scenario=="rcp4.5" & slope_year<0 & p_value<=0.05)
increase_4<- subset(slopes, scenario=="rcp4.5" & slope_year>0 & p_value<=0.05)
no_change_4 <- subset(slopes, scenario=="rcp4.5" & p_value>0.05)

# create single data frame with counts of species predicted to increase, decrease, or not change in abundance
df<-data.frame(rbind(c("decrease",nrow(decrease_4)),c("increase",nrow(increase_4)),c("no_change",nrow(no_change_4))))
names(df)<-c("direction","count")
df$scenario<-"rcp4.5"

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
df$proportion=as.numeric(df$count)/245

# Format for table in manuscript
df$count<-NULL
df_wide<-pivot_wider(df,names_from="direction",values_from="proportion") 
# Uncomment correct line below to specify GCM
#df_wide<-df_wide %>% mutate(GCM="ACCESS1-0", .before="scenario") 
#df_wide<-df_wide %>% mutate(GCM="CCSM4", .before="scenario") 
#df_wide<-df_wide %>% mutate(GCM="CNRM", .before="scenario") 
#df_wide<-df_wide %>% mutate(GCM="CSIRO", .before="scenario") 
#df_wide<-df_wide %>% mutate(GCM="INM", .before="scenario") 

# Uncomment correct line below to write .csv file
#write.csv(df_wide,"table_predicted_abundance_change_ACCESS1-0.csv",row.names=FALSE)
#write.csv(df_wide,"table_predicted_abundance_change_CCSM4.csv",row.names=FALSE)
#write.csv(df_wide,"table_predicted_abundance_change_CNRM-CM5.csv",row.names=FALSE)
#write.csv(df_wide,"table_predicted_abundance_change_CSIRO-Mk3-6-0.csv",row.names=FALSE)
#write.csv(df_wide,"table_predicted_abundance_change_INM-CM4.csv",row.names=FALSE)


##### PREDICTED FUTURE CHANGE IN TOTAL ABUNDANCE ACROSS SPECIES #####
## How has total bee abundance changed over time? #####

# calculate total abundance across bee species for each ecosystem x year combination in the historic and projected future data
abund_summary <- abund_combined_2 %>% group_by(ecosystem,scenario,year,monsoon6SPEI) %>% summarise(total_abund=sum(abund))

# remove very high outlier values
abund_summary<-subset(abund_summary,total_abund<10000)

# Plot predicted change in total abundance over time for RCP 4.5, with points colored by SPEI

# format data
rcp4.5<-subset(abund_summary, scenario=="historic" | scenario=="rcp4.5")
rcp4.5$ecosystem<-as.factor(rcp4.5$ecosystem)
levels(rcp4.5$ecosystem) <- c("Plains grassland","Chihuahuan Desert shrubland","Chihuahuan Desert grassland")
rcp4.5$scenario<-as.factor(rcp4.5$scenario)
levels(rcp4.5$scenario)<-c("long-term historic","predicted future")

# plot the trend
plot<-ggplot(data=rcp4.5, aes(x=year,y=total_abund)) + 
  geom_point(size=3,aes(fill=monsoon6SPEI*(-1),shape=scenario))+ 
  scale_shape_manual(values=c(22,21)) +
  xlab("Year") + ylab("Total abundance") + 
  theme_classic() + theme(axis.text.x = element_text(size=16))+ theme(axis.text.y = element_text(size=16)) + theme(legend.text=element_text(size=12), legend.title=element_text(size=15))+
  theme(axis.title.x = element_text(size=20))+ theme(axis.title.y = element_text(size=20)) +
  geom_smooth(method="lm",color="black") + 
  scale_fill_gradient2(high = "red3", low = "royalblue",mid = "white",midpoint = 0) + labs(fill="Aridity index", shape="Data type") 
plot

# Plot predicted change in abundance over time for RCP 8.5, with points colored by SPEI

# format data
rcp8.5<-subset(abund_summary, scenario=="historic" | scenario=="rcp8.5")
rcp8.5$ecosystem<-as.factor(rcp8.5$ecosystem)
levels(rcp8.5$ecosystem) <- c("Plains grassland","Chihuahuan Desert shrubland","Chihuahuan Desert grassland")
rcp8.5$scenario<-as.factor(rcp8.5$scenario)
levels(rcp8.5$scenario)<-c("long-term historic","predicted future")

# plot the trend
plot<-ggplot(data=rcp8.5, aes(x=year,y=total_abund)) + 
  geom_point(size=3,aes(fill=monsoon6SPEI*(-1),shape=scenario))+ 
  scale_shape_manual(values=c(22,21)) +
  xlab("Year") + ylab("Total abundance") + 
  theme_classic() + theme(axis.text.x = element_text(size=16))+ theme(axis.text.y = element_text(size=16)) + theme(legend.text=element_text(size=12), legend.title=element_text(size=15))+
  theme(axis.title.x = element_text(size=20))+ theme(axis.title.y = element_text(size=20)) +
  geom_smooth(method="lm",color="black") + 
  scale_fill_gradient2(high = "red3", low = "royalblue",mid = "white",midpoint = 0) + labs(fill="Aridity index", shape="Data type") 
plot


# Linear regressions corresponding with trends above (how has total abundance changed over time under each climate scenario?)

# RCP 4.5
rcp4.5<-subset(abund_summary, scenario=="historic" | scenario=="rcp4.5")
anova(lm(total_abund~year*ecosystem, data=rcp4.5))
summary(lm(total_abund~year*ecosystem, data=rcp4.5))

# RCP 8.5
rcp8.5<-subset(abund_summary, scenario=="historic" | scenario=="rcp8.5")
anova(lm(total_abund~year*ecosystem, data=rcp8.5))
summary(lm(total_abund~year*ecosystem, data=rcp8.5))


##### CHANGE IN COMMUNITY-WEIGHTED MEAN BODY MASS WITH ARIDITY AND OVER TIME #####
## Format data #####

# read in mass data for each bee species
mass<-read.csv("SEVBeeBodyMassData_2023-08-29_forpub.csv")

# calculate average mass for each species
avg_mass<-mass %>% group_by(code) %>% summarise(mean_mass_mg=mean(mass_mg),se_mass_mg=sd(mass_mg)/sqrt(n()))

# Calculate community-weighted mean body mass:

# merge historic and future bee abundance data with mass data
abund_focal_merged <- left_join(abund_combined_2, avg_mass, by = c("code"))

# remove rows with NAs (species for which we lack body mass data)
abund_focal_merged<-subset(abund_focal_merged,!is.na(mean_mass_mg))

# calculate total bee abundance for each ecosystem x transect x year x scenario combination
summary <- abund_focal_merged %>% group_by(ecosystem, transect, year, scenario) %>% summarise(site_year_total_abund = sum(abund))

# join the two data frames from above
cwdata <- left_join(abund_focal_merged, summary, by = c("ecosystem","transect", "year","scenario"))

# add a column of mass-weighted abundance 
cwdata$prop_abund <- cwdata$abund/cwdata$site_year_total_abund
cwdata$cw_mass <- cwdata$prop_abund*cwdata$mean_mass_mg

# calculate community-weighted mean (CWM) body mass
cwdata_final <- cwdata %>% group_by(ecosystem, transect, year, scenario) %>% summarise(cwm_mass = sum(cw_mass, na.rm = TRUE))

# Add climate data of interest:
# get climate data
climate<-abund_combined_2 %>% group_by(ecosystem, transect, year, scenario, monsoon6SPEI) %>% summarise(n=n())
# add to community-weighted mean body mass dataset
cwdata_final <- left_join(cwdata_final,climate, by=c("ecosystem","transect","year","scenario"))
# create inverse monsoon SPEI column
cwdata_final$monsoon6SPEI_positivized<-cwdata_final$monsoon6SPEI*(-1)

## Relationship between CWM mass and year in combined historic and predicted future datasets #####

# RCP 4.5

# create data frame of historic data and RCP 4.5 data
rcp4.5<-subset(cwdata_final,scenario == "historic" | scenario == "rcp4.5")

# create separate data frames for each ecosystem, and relevel factor
rcp4.5$scenario<-as.factor(rcp4.5$scenario)
levels(rcp4.5$scenario)<-c("long-term historic","predicted future")
B <- subset(rcp4.5,ecosystem=="B")
C <- subset(rcp4.5,ecosystem=="C")
G <- subset(rcp4.5,ecosystem=="G")

# Plains grassland: plot relationship between CWM body mass and year
p<-ggplot(data=B, aes(x=year,y=cwm_mass)) + 
  geom_point(size=3,aes(fill=monsoon6SPEI*(-1),shape=scenario))+ 
  xlab("Year") + ylab("Community-weighted \nmean body mass (mg)") + 
  theme_classic() + theme(axis.text.x = element_text(size=12))+ theme(axis.text.y = element_text(size=12)) + theme(legend.text=element_text(size=15), legend.title=element_text(size=14))+
  theme(axis.title.x = element_text(size=18))+ theme(axis.title.y = element_text(size=18)) +
  stat_smooth(method = "lm", formula = y ~ x, color="black") + 
  scale_shape_manual(values=c(22,21)) +
  scale_fill_gradient2(high = "red3", low = "royalblue",mid = "white",midpoint = 0) + labs(fill="Aridity index", shape="Data type") +
  ylim(0,40)
p


# Chihuahuan Desert shrubland: plot relationship between CWM body mass and year
p<-ggplot(data=C, aes(x=year,y=cwm_mass)) + 
  geom_point(size=3,aes(fill=monsoon6SPEI*(-1),shape=scenario))+ 
  xlab("Year") + ylab("Community-weighted \nmean body mass (mg)") + 
  theme_classic() + theme(axis.text.x = element_text(size=12))+ theme(axis.text.y = element_text(size=12)) + theme(legend.text=element_text(size=15), legend.title=element_text(size=14))+
  theme(axis.title.x = element_text(size=18))+ theme(axis.title.y = element_text(size=18)) +
  stat_smooth(method = "lm", formula = y ~ x, color="black") + 
  scale_shape_manual(values=c(22,21)) +
  scale_fill_gradient2(high = "red3", low = "royalblue",mid = "white",midpoint = 0) + labs(fill="Aridity index", shape="Data type")+
  ylim(0,30)
p


# Chihuahuan Desert grassland: plot relationship between CWM body mass and year
p<-ggplot(data=G, aes(x=year,y=cwm_mass)) + 
  geom_point(size=3,aes(fill=monsoon6SPEI*(-1),shape=scenario))+ 
  xlab("Year") + ylab("Community-weighted \nmean body mass (mg)") + 
  theme_classic() + theme(axis.text.x = element_text(size=12))+ theme(axis.text.y = element_text(size=12)) + theme(legend.text=element_text(size=15), legend.title=element_text(size=14))+
  theme(axis.title.x = element_text(size=18))+ theme(axis.title.y = element_text(size=18)) +
  stat_smooth(method = "lm", formula = y ~ x, color="black") +  
  scale_shape_manual(values=c(22,21)) +
  scale_fill_gradient2(high = "red3", low = "royalblue",mid = "white",midpoint = 0) + labs(fill="Aridity index", shape="Data type")+
  ylim(0,30)

p


# Statistical analyses corresponding with graphs above

# Plains grassland: subset data, run mixed effects model, and obtain statistical results
rcp4.5_B<-subset(cwdata_final, scenario=="rcp4.5" & ecosystem=="B")
m1<-lmer(cwm_mass~year+(1|transect),data=rcp4.5_B,na.action=na.omit)
Anova(m1,type = 3)
summary(m1)
rsquared(m1)

# Chihuahuan Desert grassland: subset data, run mixed effects model, and obtain statistical results
rcp4.5_G<-subset(cwdata_final, scenario=="rcp4.5" & ecosystem=="G")
m1<-lmer(cwm_mass~year+(1|transect),data=rcp4.5_G,na.action=na.omit)
Anova(m1,type = 3)
summary(m1) 
rsquared(m1)

# Chihuahuan Desert shrubland: subset data, run mixed effects model, and obtain statistical results
rcp4.5_C<-subset(cwdata_final, scenario=="rcp4.5" & ecosystem=="C")
m1<-lmer(cwm_mass~year+(1|transect),data=rcp4.5_C,na.action=na.omit)
Anova(m1,type = 3)
summary(m1) 
rsquared(m1)


# RCP 8.5

# create data frame of historic and projected future data
rcp8.5<-subset(cwdata_final,scenario == "historic" | scenario == "rcp8.5")

# create separate data frames for each ecosystem, and relevel factor
rcp8.5$scenario<-as.factor(rcp8.5$scenario)
levels(rcp8.5$scenario)<-c("long-term historic","predicted future")
B <- subset(rcp8.5,ecosystem=="B")
C <- subset(rcp8.5,ecosystem=="C")
G <- subset(rcp8.5,ecosystem=="G")

# Plains grassland: plot relationship between CWM body mass and year
p<-ggplot(data=B, aes(x=year,y=cwm_mass)) + 
  geom_point(size=3,aes(fill=monsoon6SPEI*(-1),shape=scenario))+ 
  xlab("Year") + ylab("Community-weighted \nmean body mass (mg)") + 
  theme_classic() +theme(axis.text.x = element_text(size=12))+ theme(axis.text.y = element_text(size=12)) + theme(legend.text=element_text(size=15), legend.title=element_text(size=14))+ theme(axis.title.x = element_text(size=18))+ theme(axis.title.y = element_text(size=18)) +
  stat_smooth(method = "lm", formula = y ~ x, color="black") + 
  scale_shape_manual(values=c(22,21)) +
  scale_fill_gradient2(high = "red3", low = "royalblue",mid = "white",midpoint = 0) + labs(fill="Aridity index", shape="Data type") +
  ylim(0,40)
p


# Chihuahuan Desert grassland: plot relationship between CWM body mass and year
p<-ggplot(data=G, aes(x=year,y=cwm_mass)) + 
  geom_point(size=3,aes(fill=monsoon6SPEI*(-1),shape=scenario))+ 
  xlab("Year") + ylab("Community-weighted \nmean body mass (mg)") + 
  theme_classic() +theme(axis.text.x = element_text(size=12))+ theme(axis.text.y = element_text(size=12)) + theme(legend.text=element_text(size=15), legend.title=element_text(size=14))+ theme(axis.title.x = element_text(size=18))+ theme(axis.title.y = element_text(size=18)) +
  stat_smooth(method = "lm", formula = y ~ x, color="black") +  
  scale_shape_manual(values=c(22,21)) +
  scale_fill_gradient2(high = "red3", low = "royalblue",mid = "white",midpoint = 0) + labs(fill="Aridity index", shape="Data type")+
  ylim(0,40)
p


# Chihuahuan Desert shrubland: plot relationship between CWM body mass and year
p<-ggplot(data=C, aes(x=year,y=cwm_mass)) + 
  geom_point(size=3,aes(fill=monsoon6SPEI*(-1),shape=scenario))+ 
  xlab("Year") + ylab("Community-weighted \nmean body mass (mg)") + 
  theme_classic() +theme(axis.text.x = element_text(size=12))+ theme(axis.text.y = element_text(size=12)) + theme(legend.text=element_text(size=15), legend.title=element_text(size=14))+ theme(axis.title.x = element_text(size=18))+ theme(axis.title.y = element_text(size=18)) +
  stat_smooth(method = "lm", formula = y ~ x, color="black") + 
  scale_shape_manual(values=c(22,21)) +
  scale_fill_gradient2(high = "red3", low = "royalblue",mid = "white",midpoint = 0) + labs(fill="Aridity index", shape="Data type")+
  ylim(0,40)
p


# Statistical analyses corresponding with graphs above

# Plains grassland: subset data, run mixed effects model, and obtain statistical results
rcp8.5_B<-subset(cwdata_final, scenario=="rcp8.5" & ecosystem=="B")
m1<-lmer(cwm_mass~year+(1|transect),data=rcp8.5_B,na.action=na.omit)
Anova(m1,type = 3)
summary(m1) 
rsquared(m1)

# Chihuahuan Desert grassland: subset data, run mixed effects model, and obtain statistical results
rcp8.5_G<-subset(cwdata_final, scenario=="rcp8.5" & ecosystem=="G")
m1<-lmer(cwm_mass~year+(1|transect),data=rcp8.5_G,na.action=na.omit)
Anova(m1,type = 3)
summary(m1) 
rsquared(m1)

# Chihuahuan Desert shrubland: subset data, run mixed effects model, and obtain statistical results
rcp8.5_C<-subset(cwdata_final, scenario=="rcp8.5" & ecosystem=="C")
m1<-lmer(cwm_mass~year+(1|transect),data=rcp8.5_C,na.action=na.omit)
Anova(m1,type = 3)
summary(m1) 
rsquared(m1)


## SUMMARY OF RESULTS #####

# Read in file that output from above was transferred into
data<-read.csv("FutureBeeAbundanceStats_ComparisonsAmongGCMs.csv")

# Create summaries

# proportion of species predicted to increase/decrease/not change in abundance
summary<-data %>% group_by(scenario) %>% summarise(prop_decrease=mean(prop_decrease),
                                                   prop_increase=mean(prop_increase),
                                                   prop_nochange=mean(prop_nochange))

incrdecr<-data %>% filter(scenario=="RCP 4.5") %>% summarise(max_decrease=max(prop_decrease),
                                                            min_decrease=min(prop_decrease),
                                                            max_increase=max(prop_increase),
                                                            min_increase=min(prop_increase),
                                                            max_nochange=max(prop_nochange),
                                                            min_nochange=min(prop_nochange))

# total change in abundance over time
abundchange<-data %>% mutate(count=1) %>% group_by(scenario, abund_trend) %>% summarise(count=sum(count)) %>% pivot_wider(names_from="abund_trend",values_from="count")

# community-weighted mean body mass change
mass_b<-data %>% mutate(count=1) %>% group_by(scenario, body_mass_b) %>% summarise(count=sum(count)) %>% pivot_wider(names_from="body_mass_b",values_from="count")

mass_g<-data %>% mutate(count=1) %>% group_by(scenario, body_mass_g) %>% summarise(count=sum(count)) %>% pivot_wider(names_from="body_mass_g",values_from="count")

mass_c<-data %>% mutate(count=1) %>% group_by(scenario, body_mass_c) %>% summarise(count=sum(count)) %>% pivot_wider(names_from="body_mass_c",values_from="count")
