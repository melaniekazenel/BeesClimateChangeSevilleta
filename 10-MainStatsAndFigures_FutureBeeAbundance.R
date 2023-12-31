################################################################################### 
# Statistics related to predicted bee community change under projected future climate scenarios

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

##### PREDICTED FUTURE CHANGE IN INDIVIDUAL SPECIES ABUNDANCES: CANESM2 GCM #####
## Data manipulation #####

## Format predicted future bee abundance data ##

# read in data
abund<-read.csv("predicted_abundances_CanESM2-2023-12-01.csv")

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

# read in SPEI data
clim_new<-read.csv("spei_historic_and_future_byscenario_CanESM2_2023-12-01.csv")

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
#write.csv(specieslist2,"slopes_predicted_future_change_CanESM2_2023-12-01.csv",row.names = FALSE)


## Summary statistics: Which species predicted to increase, decrease, or remain stable in abundance over time? #####

# Read in slope data (calculated above)
slopes <- read.csv("slopes_predicted_future_change_CanESM2_2023-12-01.csv")

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
df_wide<-df_wide %>% mutate(GCM="CanESM2", .before="scenario") 
#write.csv(df_wide,"table_predicted_abundance_change_CanESM2_2023-12-01.csv",row.names=FALSE)


##### PREDICTED FUTURE CHANGE IN TOTAL ABUNDANCE ACROSS SPECIES #####
## How will total bee abundance change over time? #####

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
plot_abundchange<-ggplot(data=rcp4.5, aes(x=year,y=total_abund)) + 
  geom_point(size=3,aes(fill=monsoon6SPEI*(-1),shape=scenario))+ 
  scale_shape_manual(values=c(22,21)) +
  xlab("Year") + ylab("Total abundance") + 
  theme_classic() + theme(axis.text.x = element_text(size=16))+ theme(axis.text.y = element_text(size=16)) + theme(legend.text=element_text(size=13), legend.title=element_text(size=15))+
  theme(axis.title.x = element_text(size=20))+ theme(axis.title.y = element_text(size=20)) +
  geom_smooth(method="lm",color="black") + 
  scale_fill_gradient2(high = "red3", low = "royalblue",mid = "white",midpoint = 0) + labs(fill="Aridity index", shape="Data type") +
  annotate("text", x = 2011, y = 6000, label = expression(atop(italic(P)*" = 0.6324",italic(r^2)=="0.10")),size=4.5) +
  theme(plot.margin = unit(c(0,0,0,20), "pt")) +
  scale_y_continuous(breaks = c(0,1000,2000,3000,4000,5000,6000))
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
  geom_point(size=3,aes(fill=monsoon6SPEI*(-1),shape=scenario))+ 
  scale_shape_manual(values=c(22,21)) +
  xlab("Year") + ylab("Total abundance") + 
  theme_bw() + theme(axis.text.x = element_text(size=16))+ theme(axis.text.y = element_text(size=16)) + theme(legend.text=element_text(size=12), legend.title=element_text(size=15))+
  theme(axis.title.x = element_text(size=20))+ theme(axis.title.y = element_text(size=20)) +
  geom_smooth(method="lm",color="black") + 
  scale_fill_gradient2(high = "red3", low = "royalblue",mid = "white",midpoint = 0) + labs(fill="Aridity index", shape="Data type") + ylim(0,10000) + facet_wrap(~scenario2) + theme(panel.spacing.x = unit(10, "mm")) + theme(strip.text.x = element_text(size = 18)) + theme(axis.line = element_line(color='black'), plot.background = element_blank(),panel.grid.minor = element_blank(),panel.grid.major = element_blank())
plot

# save plot
# ggsave("change_total_abundance_rcp2.6_and_8.5_2023-08-29.jpg", plot,
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


## Relationship between CWM mass and SPEI or year in combined historic and predicted future datasets #####

# create a vector of ecosystem names
ecosystems <- c("B","G","C")

# create a vector of scenario names
scenarios <- c("rcp2.6","rcp4.5","rcp8.5")

# create a new matrix to hold the output of the loop below
output<-matrix(nrow=1,ncol=19)

# Loop: CWM mass ~ aridity
for (i in 1:length(ecosystems)){
  for (j in 1:length(scenarios)){
    
    variable <- "aridity"
    
    # create ecosystem and scenario ID objects
    eco_id<-ecosystems[i]
    scenario_id<-scenarios[j]
    
    # subset the data
    model_data<-subset(cwdata_final,ecosystem==eco_id)
    model_data<-subset(model_data,scenario==scenario_id | scenario=="historic")
    
    # mixed effects models
    m1<-lmer(cwm_mass~monsoon6SPEI_positivized+(1|transect)+(1|year),data=model_data,na.action=na.omit)
    m2<-lmer(cwm_mass~monsoon6SPEI_positivized+I(monsoon6SPEI_positivized^2)+(1|transect)+(1|year),data=model_data,na.action=na.omit)
    
    # get statistical output
    AICc_lin <-AICc(m1)
    rsquared_lin <- rsquared(m1)$Conditional
    xsquared_lin <- Anova(m1,type = 3)[2,1]
    p_lin <- Anova(m1,type = 3)[2,3]
    lin_param <- summary(m1)$coefficients[2,1]
    lin_param_se <- summary(m1)$coefficients[2,2]
    
    AICc_quad <-AICc(m2)
    rsquared_quad <- rsquared(m2)$Conditional
    xsquared_quad_lin <- Anova(m2,type = 3)[2,1]
    p_quad_lin <- Anova(m2,type = 3)[2,3]
    xsquared_quad <-  Anova(m2,type = 3)[3,1]
    p_quad <- Anova(m2,type = 3)[3,3]
    quad_lin_param <- summary(m2)$coefficients[2,1]
    quad_lin_param_se <- summary(m2)$coefficients[2,2]
    quad_param <- summary(m2)$coefficients[3,1]
    quad_param_se <- summary(m2)$coefficients[3,2]
    
    # bind ID and statistical output values
    output2<-cbind(variable,scenario_id,eco_id,
                   AICc_lin,rsquared_lin,xsquared_lin,p_lin,lin_param,lin_param_se,
                   AICc_quad,rsquared_quad,xsquared_quad_lin,p_quad_lin,xsquared_quad,p_quad, 
                   quad_lin_param,quad_lin_param_se,quad_param,quad_param_se)
    
    # bind species output to main output data frame
    output <-rbind(output,output2)
    
  }
  
}

# create and format data frame containing output
output<-output[-1,]

# Loop: CWM mass ~ year
for (i in 1:length(ecosystems)){
  for (j in 1:length(scenarios)){
    
    variable <- "year"
    
    # create ecosystem and scenario ID objects
    eco_id<-ecosystems[i]
    scenario_id<-scenarios[j]
    
    # subset the data
    model_data<-subset(cwdata_final,ecosystem==eco_id)
    model_data<-subset(model_data,scenario==scenario_id | scenario=="historic")
    
    # mixed effects models
    m1<-lmer(cwm_mass~year+(1|transect),data=model_data,na.action=na.omit)
    
    # get statistical output
    AICc_lin <-AICc(m1)
    rsquared_lin <- rsquared(m1)$Conditional
    xsquared_lin <- Anova(m1,type = 3)[2,1]
    p_lin <- Anova(m1,type = 3)[2,3]
    lin_param <- summary(m1)$coefficients[2,1]
    lin_param_se <- summary(m1)$coefficients[2,2]
    
    AICc_quad <-NA
    rsquared_quad <- NA
    xsquared_quad_lin <- NA
    p_quad_lin <- NA
    xsquared_quad <-  NA
    p_quad <- NA
    quad_lin_param <- NA
    quad_lin_param_se <- NA
    quad_param <- NA
    quad_param_se <- NA
    
    # bind ID and statistical output values
    output2<-cbind(variable,scenario_id,eco_id,
                   AICc_lin,rsquared_lin,xsquared_lin,p_lin,lin_param,lin_param_se,
                   AICc_quad,rsquared_quad,xsquared_quad_lin,p_quad_lin,xsquared_quad,p_quad, 
                   quad_lin_param,quad_lin_param_se,quad_param,quad_param_se)
    
    # bind species output to main output data frame
    output <-rbind(output,output2)
    
  }
  
}

# create and format data frame containing output
output<-as.data.frame(output)

# write .csv file of results
#write.csv(output,"cwbodymassresults_CanESM2_2023-12-29.csv",row.names=FALSE)



##### PREDICTED FUTURE CHANGE IN INDIVIDUAL SPECIES ABUNDANCES: OTHER GCMS #####
## Data manipulation #####

## Format predicted future bee abundance data ##

# read in data frame of projected future bee data, completing this step and running the code below separately for each GCM
abund<-read.csv("predicted_abundances_ACCESS1-0-2023-12-01.csv")
#abund<-read.csv("predicted_abundances_CCSM4-2023-12-01.csv")
#abund<-read.csv("predicted_abundances_CNRM-CM5-2023-12-01.csv")
#abund<-read.csv("predicted_abundances_CSIRO-Mk3-6-0-2023-12-01.csv")
#abund<-read.csv("predicted_abundances_INM-CM4-2023-12-01.csv")

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
clim_new<-read.csv("spei_historic_and_future_byscenario_ACCESS1-0_2023-12-01.csv")
#clim_new<-read.csv("spei_historic_and_future_byscenario_CCSM4_2023-12-01.csv")
#clim_new<-read.csv("spei_historic_and_future_byscenario_CNRM-CM5_2023-12-01.csv")
#clim_new<-read.csv("spei_historic_and_future_byscenario_CSIRO-Mk3-6-0_2023-12-01.csv")
#clim_new<-read.csv("spei_historic_and_future_byscenario_INM-CM4_2023-12-01.csv")

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
#write.csv(specieslist2,"slopes_predicted_future_change_ACCESS1-0_2023-12-01.csv",row.names = FALSE)
#write.csv(specieslist2,"slopes_predicted_future_change_CCSM4_2023-12-01.csv",row.names = FALSE)
#write.csv(specieslist2,"slopes_predicted_future_change_CNRM-CM5_2023-12-01.csv",row.names = FALSE)
#write.csv(specieslist2,"slopes_predicted_future_change_CSIRO-Mk3-6-0_2023-12-01.csv",row.names = FALSE)
#write.csv(specieslist2,"slopes_predicted_future_change_INM-CM4_2023-12-01.csv",row.names = FALSE)


## Summary statistics: Which species predicted to increase, decrease, or remain stable in abundance over time? #####

# Read in slopes data from calculations above
slopes <- read.csv("slopes_predicted_future_change_ACCESS1-0_2023-12-01.csv")
#slopes <- read.csv("slopes_predicted_future_change_CCSM4_2023-12-01.csv")
#slopes <- read.csv("slopes_predicted_future_change_CNRM-CM5_2023-12-01.csv")
#slopes <- read.csv("slopes_predicted_future_change_CSIRO-Mk3-6-0_2023-12-01.csv")
#slopes <- read.csv("slopes_predicted_future_change_INM-CM4_2023-12-01.csv")

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
df$proportion=as.numeric(df$count)/243

# Format results for table in manuscript
df$count<-NULL
df_wide<-pivot_wider(df,names_from="direction",values_from="proportion") 
# Uncomment correct line below to specify GCM
df_wide<-df_wide %>% mutate(GCM="ACCESS1-0", .before="scenario") 
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

## Relationship between CWM mass and SPEI or year in combined historic and predicted future datasets #####

# create a vector of ecosystems
ecosystems <- c("B","G","C")

# create a vector of scenario names
scenarios <- c("rcp4.5","rcp8.5")

# create a new matrix to hold the output of the loop below
output<-matrix(nrow=1,ncol=19)

# Loop: CWM mass ~ aridity
for (i in 1:length(ecosystems)){
  for (j in 1:length(scenarios)){
    
    variable <- "aridity"
    
    # create ecosystem and scenario ID objects
    eco_id<-ecosystems[i]
    scenario_id<-scenarios[j]
    
    # subset the data
    model_data<-subset(cwdata_final,ecosystem==eco_id)
    model_data<-subset(model_data,scenario==scenario_id | scenario=="historic")
    
    # mixed effects models
    m1<-lmer(cwm_mass~monsoon6SPEI_positivized+(1|transect)+(1|year),data=model_data,na.action=na.omit)
    m2<-lmer(cwm_mass~monsoon6SPEI_positivized+I(monsoon6SPEI_positivized^2)+(1|transect)+(1|year),data=model_data,na.action=na.omit)
    
    # get statistical output
    AICc_lin <-AICc(m1)
    rsquared_lin <- rsquared(m1)$Conditional
    xsquared_lin <- Anova(m1,type = 3)[2,1]
    p_lin <- Anova(m1,type = 3)[2,3]
    lin_param <- summary(m1)$coefficients[2,1]
    lin_param_se <- summary(m1)$coefficients[2,2]
    
    AICc_quad <-AICc(m2)
    rsquared_quad <- rsquared(m2)$Conditional
    xsquared_quad_lin <- Anova(m2,type = 3)[2,1]
    p_quad_lin <- Anova(m2,type = 3)[2,3]
    xsquared_quad <-  Anova(m2,type = 3)[3,1]
    p_quad <- Anova(m2,type = 3)[3,3]
    quad_lin_param <- summary(m2)$coefficients[2,1]
    quad_lin_param_se <- summary(m2)$coefficients[2,2]
    quad_param <- summary(m2)$coefficients[3,1]
    quad_param_se <- summary(m2)$coefficients[3,2]
    
    # bind ID and statistical output values
    output2<-cbind(variable,scenario_id,eco_id,
                   AICc_lin,rsquared_lin,xsquared_lin,p_lin,lin_param,lin_param_se,
                   AICc_quad,rsquared_quad,xsquared_quad_lin,p_quad_lin,xsquared_quad,p_quad, 
                   quad_lin_param,quad_lin_param_se,quad_param,quad_param_se)
    
    # bind species output to main output data frame
    output <-rbind(output,output2)
    
  }
  
}

# create and format data frame containing output
output<-output[-1,]

# Loop: CWM mass ~ year
for (i in 1:length(ecosystems)){
  for (j in 1:length(scenarios)){
    
    variable <- "year"
    
    # create ecosystem and scenario ID objects
    eco_id<-ecosystems[i]
    scenario_id<-scenarios[j]
    
    # subset the data
    model_data<-subset(cwdata_final,ecosystem==eco_id)
    model_data<-subset(model_data,scenario==scenario_id | scenario=="historic")
    
    # mixed effects models
    m1<-lmer(cwm_mass~year+(1|transect),data=model_data,na.action=na.omit)
    
    # get statistical output
    AICc_lin <-AICc(m1)
    rsquared_lin <- rsquared(m1)$Conditional
    xsquared_lin <- Anova(m1,type = 3)[2,1]
    p_lin <- Anova(m1,type = 3)[2,3]
    lin_param <- summary(m1)$coefficients[2,1]
    lin_param_se <- summary(m1)$coefficients[2,2]
    
    AICc_quad <-NA
    rsquared_quad <- NA
    xsquared_quad_lin <- NA
    p_quad_lin <- NA
    xsquared_quad <-  NA
    p_quad <- NA
    quad_lin_param <- NA
    quad_lin_param_se <- NA
    quad_param <- NA
    quad_param_se <- NA
    
    # bind ID and statistical output values
    output2<-cbind(variable,scenario_id,eco_id,
                   AICc_lin,rsquared_lin,xsquared_lin,p_lin,lin_param,lin_param_se,
                   AICc_quad,rsquared_quad,xsquared_quad_lin,p_quad_lin,xsquared_quad,p_quad, 
                   quad_lin_param,quad_lin_param_se,quad_param,quad_param_se)
    
    # bind species output to main output data frame
    output <-rbind(output,output2)
    
  }
  
}

# create and format data frame containing output
output<-as.data.frame(output)

# write output as .csv file
#write.csv(output,"cwbodymassresults_ACCESS1-0_2023-12-29.csv",row.names=FALSE)
#write.csv(output,"cwbodymassresults_CCSM_2023-12-29.csv",row.names=FALSE)
#write.csv(output,"cwbodymassresults_CNRM_2023-12-29.csv",row.names=FALSE)
#write.csv(output,"cwbodymassresults_CSIRO_2023-12-29.csv",row.names=FALSE)
#write.csv(output,"cwbodymassresults_INM_2023-12-29.csv",row.names=FALSE)


##### SUMMARY OF RESULTS #####
## Summary of results from the six GCMs #####
# Read in file that output from above was transferred into
data<-read.csv("FutureBeeAbundanceStats_ComparisonsAmongGCMs_2023-12-01.csv")

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

## Figure: Parameter estimates from models of CWM body mass ~ aridity or year, for each GCM - RCP 4.5 #####

# read in results from aridity models and recode variables
mass_spei<-read.csv("paramest_cwmbodymass_aridity_diffgcms_rcp4-5_2023-12-01.csv")

mass_spei$ecosystem<-as.factor(mass_spei$ecosystem)
mass_spei$ecosystem <- factor(mass_spei$ecosystem, levels=c("B","G","C"))

levels(mass_spei$ecosystem)[levels(mass_spei$ecosystem)=="B"] <- "Plains grassland"
levels(mass_spei$ecosystem)[levels(mass_spei$ecosystem)=="C"] <- "Desert shrubland"
levels(mass_spei$ecosystem)[levels(mass_spei$ecosystem)=="G"] <- "Desert grassland"

mass_spei$gcm<-as.factor(mass_spei$gcm)
levels(mass_spei$gcm)
levels(mass_spei$gcm) <- c("ACCESS 1.0", "CanESM2", "CCSM 4.0", "CNRM-CM5", "CSIRO-Mk3.6.0", "INM-CM4")
mass_spei$gcm <- factor(mass_spei$gcm, levels=c( "INM-CM4", "CSIRO-Mk3.6.0", "CNRM-CM5", "CCSM 4.0","CanESM2", "ACCESS 1.0" ))

# plot: linear parameter estimate
lin <- ggplot(data=mass_spei, aes(x=lin_param, y=gcm)) +
  geom_errorbarh(height = 0.2, aes(xmax = lin_param + lin_param_se, xmin = lin_param - lin_param_se))+
  geom_vline(xintercept=0, linetype="dashed", color = "grey20") +
  geom_point(size=4) +
  theme_bw() + theme(axis.text.x = element_text(size=15, color="black"))+ theme(axis.text.y = element_text(size=16, color="black"))+
  theme(axis.title.x = element_text(size=17))+ 
  xlab("Linear parameter \nestimate (aridity)") +
  theme(axis.title.y = element_blank()) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  facet_wrap(~ecosystem, ncol=1) + 
  theme(strip.text.x = element_text(size = 15))
lin

# plot: quadratic parameter estimate
quad <- ggplot(data=mass_spei, aes(x=quad_param, y=gcm)) +
  geom_errorbarh(height = 0.2, aes(xmax = quad_param + quad_param_se, xmin = quad_param - quad_param_se))+
  geom_vline(xintercept=0, linetype="dashed", color = "grey20") +
  geom_point(size=4) +
  theme_bw() + theme(axis.text.x = element_text(size=15, color="black"))+ theme(axis.text.y = element_text(size=16, color="black"))+
  theme(axis.title.x = element_text(size=17))+ 
  xlab("Quadratic parameter \nestimate (aridity)") +
  theme(axis.title.y = element_blank()) +
  theme(axis.text.y = element_blank()) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  facet_wrap(~ecosystem, ncol=1) + 
  theme(strip.text.x = element_text(size = 15))
quad

# read in results from year models and recode variables
mass_year<-read.csv("paramest_cwmbodymass_year_diffgcms_rcp4-5_2023-12-01.csv")

# rename levels of the "ecosystem" variable
mass_year$ecosystem<-as.factor(mass_year$ecosystem)
mass_year$ecosystem <- factor(mass_year$ecosystem, levels=c("B","G","C"))

levels(mass_year$ecosystem)[levels(mass_year$ecosystem)=="B"] <- "Plains grassland"
levels(mass_year$ecosystem)[levels(mass_year$ecosystem)=="C"] <- "Desert shrubland"
levels(mass_year$ecosystem)[levels(mass_year$ecosystem)=="G"] <- "Desert grassland"

mass_year$gcm<-as.factor(mass_year$gcm)
levels(mass_year$gcm)
levels(mass_year$gcm) <- c("ACCESS 1.0", "CanESM2", "CCSM 4.0", "CNRM-CM5", "CSIRO-Mk3.6.0", "INM-CM4")
mass_year$gcm <- factor(mass_year$gcm, levels=c( "INM-CM4", "CSIRO-Mk3.6.0", "CNRM-CM5", "CCSM 4.0","CanESM2", "ACCESS 1.0" ))

# plot: linear parameter estimate
linyear <- ggplot(data=mass_year, aes(x=lin_param, y=gcm)) +
  geom_errorbarh(height = 0.2, aes(xmax = lin_param + lin_param_se, xmin = lin_param - lin_param_se))+
  geom_vline(xintercept=0, linetype="dashed", color = "grey20") +
  geom_point(size=4) +
  theme_bw() + theme(axis.text.x = element_text(size=15, color="black"))+ theme(axis.text.y = element_text(size=16, color="black"))+
  theme(axis.title.x = element_text(size=17))+ 
  xlab("Linear parameter \nestimate (year)") +
  theme(axis.title.y = element_blank()) +
  theme(axis.text.y = element_blank()) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  facet_wrap(~ecosystem, ncol=1) + 
  theme(strip.text.x = element_text(size = 15))
linyear

p1<-lin + quad + linyear + plot_layout(ncol = 3)
p1

#ggsave("fig_cwm_mass_aridityandyear_paramestimates_2023-12-29.pdf", p1, width=13,height=10,units = c("in"),dpi = 600)


## Figure: Parameter estimates from models of CWM body mass ~ aridity or year, for each GCM - RCP 2.6 #####

# read in results from aridity models and recode variables
mass_spei<-read.csv("paramest_cwmbodymass_aridity_diffgcms_rcp2-6_2023-12-01.csv")

mass_spei$ecosystem<-as.factor(mass_spei$ecosystem)
mass_spei$ecosystem <- factor(mass_spei$ecosystem, levels=c("B","G","C"))

levels(mass_spei$ecosystem)[levels(mass_spei$ecosystem)=="B"] <- "Plains grassland"
levels(mass_spei$ecosystem)[levels(mass_spei$ecosystem)=="C"] <- "Desert shrubland"
levels(mass_spei$ecosystem)[levels(mass_spei$ecosystem)=="G"] <- "Desert grassland"

# plot: linear parameter estimate
lin <- ggplot(data=mass_spei, aes(x=lin_param, y=gcm)) +
  geom_errorbarh(height = 0.2, aes(xmax = lin_param + lin_param_se, xmin = lin_param - lin_param_se))+
  geom_vline(xintercept=0, linetype="dashed", color = "grey20") +
  geom_point(size=4) +
  theme_bw() + theme(axis.text.x = element_text(size=15, color="black"))+ theme(axis.text.y = element_text(size=16, color="black"))+
  theme(axis.title.x = element_text(size=17))+ 
  xlab("Linear parameter \nestimate (aridity)") +
  theme(axis.title.y = element_blank()) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  facet_wrap(~ecosystem, ncol=1) + 
  theme(strip.text.x = element_text(size = 15))
lin

# plot: quadratic parameter estimate
quad <- ggplot(data=mass_spei, aes(x=quad_param, y=gcm)) +
  geom_errorbarh(height = 0.2, aes(xmax = quad_param + quad_param_se, xmin = quad_param - quad_param_se))+
  geom_vline(xintercept=0, linetype="dashed", color = "grey20") +
  geom_point(size=4) +
  theme_bw() + theme(axis.text.x = element_text(size=15, color="black"))+ theme(axis.text.y = element_text(size=16, color="black"))+
  theme(axis.title.x = element_text(size=17))+ 
  xlab("Quadratic parameter \nestimate (aridity)") +
  theme(axis.title.y = element_blank()) +
  theme(axis.text.y = element_blank()) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  facet_wrap(~ecosystem, ncol=1) + 
  theme(strip.text.x = element_text(size = 15))
quad

# read in results from year models and recode variables
mass_year<-read.csv("paramest_cwmbodymass_year_diffgcms_rcp2-6_2023-12-01.csv")

# rename levels of the "ecosystem" variable
mass_year$ecosystem<-as.factor(mass_year$ecosystem)
mass_year$ecosystem <- factor(mass_year$ecosystem, levels=c("B","G","C"))

levels(mass_year$ecosystem)[levels(mass_year$ecosystem)=="B"] <- "Plains grassland"
levels(mass_year$ecosystem)[levels(mass_year$ecosystem)=="C"] <- "Desert shrubland"
levels(mass_year$ecosystem)[levels(mass_year$ecosystem)=="G"] <- "Desert grassland"

# plot: linear parameter estimate
linyear <- ggplot(data=mass_year, aes(x=lin_param, y=gcm)) +
  geom_errorbarh(height = 0.2, aes(xmax = lin_param + lin_param_se, xmin = lin_param - lin_param_se))+
  geom_vline(xintercept=0, linetype="dashed", color = "grey20") +
  geom_point(size=4) +
  theme_bw() + theme(axis.text.x = element_text(size=15, color="black"))+ theme(axis.text.y = element_text(size=16, color="black"))+
  theme(axis.title.x = element_text(size=17))+ 
  xlab("Linear parameter \nestimate (year)") +
  theme(axis.title.y = element_blank()) +
  theme(axis.text.y = element_blank()) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  facet_wrap(~ecosystem, ncol=1) + 
  theme(strip.text.x = element_text(size = 15))
linyear

p1<-lin + quad + linyear + plot_layout(ncol = 3)
p1

#ggsave("fig_cwm_mass_aridityandyear_paramestimates_rcp2-6_2023-12-29.pdf", p1, width=13,height=4,units = c("in"),dpi = 600)

## Figure: Parameter estimates from models of CWM body mass ~ aridity or year, for each GCM - RCP 8.5 #####

# read in results from aridity models and recode variables
mass_spei<-read.csv("paramest_cwmbodymass_aridity_diffgcms_rcp8-5_2023-12-01.csv")

mass_spei$ecosystem<-as.factor(mass_spei$ecosystem)
mass_spei$ecosystem <- factor(mass_spei$ecosystem, levels=c("B","G","C"))

levels(mass_spei$ecosystem)[levels(mass_spei$ecosystem)=="B"] <- "Plains grassland"
levels(mass_spei$ecosystem)[levels(mass_spei$ecosystem)=="C"] <- "Desert shrubland"
levels(mass_spei$ecosystem)[levels(mass_spei$ecosystem)=="G"] <- "Desert grassland"

mass_spei$gcm<-as.factor(mass_spei$gcm)
levels(mass_spei$gcm)
levels(mass_spei$gcm) <- c("ACCESS 1.0", "CanESM2", "CCSM 4.0", "CNRM-CM5", "CSIRO-Mk3.6.0", "INM-CM4")
mass_spei$gcm <- factor(mass_spei$gcm, levels=c( "INM-CM4", "CSIRO-Mk3.6.0", "CNRM-CM5", "CCSM 4.0","CanESM2", "ACCESS 1.0" ))

# plot: linear parameter estimate
lin <- ggplot(data=mass_spei, aes(x=lin_param, y=gcm)) +
  geom_errorbarh(height = 0.2, aes(xmax = lin_param + lin_param_se, xmin = lin_param - lin_param_se))+
  geom_vline(xintercept=0, linetype="dashed", color = "grey20") +
  geom_point(size=4) +
  theme_bw() + theme(axis.text.x = element_text(size=15, color="black"))+ theme(axis.text.y = element_text(size=16, color="black"))+
  theme(axis.title.x = element_text(size=17))+ 
  xlab("Linear parameter \nestimate (aridity)") +
  theme(axis.title.y = element_blank()) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  facet_wrap(~ecosystem, ncol=1) + 
  theme(strip.text.x = element_text(size = 15))
lin

# plot: quadratic parameter estimate
quad <- ggplot(data=mass_spei, aes(x=quad_param, y=gcm)) +
  geom_errorbarh(height = 0.2, aes(xmax = quad_param + quad_param_se, xmin = quad_param - quad_param_se))+
  geom_vline(xintercept=0, linetype="dashed", color = "grey20") +
  geom_point(size=4) +
  theme_bw() + theme(axis.text.x = element_text(size=15, color="black"))+ theme(axis.text.y = element_text(size=16, color="black"))+
  theme(axis.title.x = element_text(size=17))+ 
  xlab("Quadratic parameter \nestimate (aridity)") +
  theme(axis.title.y = element_blank()) +
  theme(axis.text.y = element_blank()) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  facet_wrap(~ecosystem, ncol=1) + 
  theme(strip.text.x = element_text(size = 15))
quad

# read in results from year models and recode variables
mass_year<-read.csv("paramest_cwmbodymass_year_diffgcms_rcp8-5_2023-12-01.csv")

# rename levels of the "ecosystem" variable
mass_year$ecosystem<-as.factor(mass_year$ecosystem)
mass_year$ecosystem <- factor(mass_year$ecosystem, levels=c("B","G","C"))

levels(mass_year$ecosystem)[levels(mass_year$ecosystem)=="B"] <- "Plains grassland"
levels(mass_year$ecosystem)[levels(mass_year$ecosystem)=="C"] <- "Desert shrubland"
levels(mass_year$ecosystem)[levels(mass_year$ecosystem)=="G"] <- "Desert grassland"

mass_year$gcm<-as.factor(mass_year$gcm)
levels(mass_year$gcm)
levels(mass_year$gcm) <- c("ACCESS 1.0", "CanESM2", "CCSM 4.0", "CNRM-CM5", "CSIRO-Mk3.6.0", "INM-CM4")
mass_year$gcm <- factor(mass_year$gcm, levels=c( "INM-CM4", "CSIRO-Mk3.6.0", "CNRM-CM5", "CCSM 4.0","CanESM2", "ACCESS 1.0" ))

# plot: linear parameter estimate
linyear <- ggplot(data=mass_year, aes(x=lin_param, y=gcm)) +
  geom_errorbarh(height = 0.2, aes(xmax = lin_param + lin_param_se, xmin = lin_param - lin_param_se))+
  geom_vline(xintercept=0, linetype="dashed", color = "grey20") +
  geom_point(size=4) +
  theme_bw() + theme(axis.text.x = element_text(size=15, color="black"))+ theme(axis.text.y = element_text(size=16, color="black"))+
  theme(axis.title.x = element_text(size=17))+ 
  xlab("Linear parameter \nestimate (year)") +
  theme(axis.title.y = element_blank()) +
  theme(axis.text.y = element_blank()) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  facet_wrap(~ecosystem, ncol=1) + 
  theme(strip.text.x = element_text(size = 15))
linyear

p1<-lin + quad + linyear + plot_layout(ncol = 3)
p1

#ggsave("fig_cwm_mass_aridityandyear_paramestimates_rcp8-5_2023-12-29.pdf", p1, width=13,height=10,units = c("in"),dpi = 600)
