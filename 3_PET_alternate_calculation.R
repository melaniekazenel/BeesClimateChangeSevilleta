################################################################################### 
# Year-to-year variation in the aridity index SPEI calculated using two different PET estimation methods (Thornthwaite and Penman) 

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
library(lubridate)


# Read in data from Ameriflux towers
grassland<-read.csv("US-Seg_daily_aflx.csv") # Plains ecosystem
shrubland<-read.csv("US-Ses_daily_aflx.csv") # Chihuahuan Desert ecosystems


##### SPEI Calculations: Plains Ecosystem #####

# Convert date to standard format
names(grassland)[1]<-"date"
grassland$date<-mdy(grassland$date)

# Create year and month columns
grassland$year<-year(grassland$date)
grassland$month<-month(grassland$date)

# For variables of interest, calculate mean values for each month x year combination
grass_month<-grassland %>% group_by(month,year) %>% summarise(P_sum=sum(P_int, na.rm=TRUE),PET1pen_mean=mean(PET1_int, na.rm=TRUE),PET2pen_mean=mean(PET2_int, na.rm=TRUE),Tair_mean=mean(TA_avg, na.rm=TRUE))
grass_month<-na.omit(grass_month)

# Calculate PET using the Thornthwaite method
grass_month$PET_thorn <- as.numeric(thornthwaite(Tave = grass_month$Tair_mean, lat = 34.35, na.rm = FALSE))

# Calculate climatic water balance (CWB) using both Penman and Thornthwaite PET values, and sort the data frame. Two different Penman values were used (see Supplementary Information for details)
grass_month <- grass_month %>% mutate(CWB1_pen=P_sum-PET1pen_mean, CWB2_pen=P_sum-PET2pen_mean, CWB_thorn=P_sum-PET_thorn)
grass_month<-arrange(grass_month, year, month)

# Subset the data to just include focal years
grass_month2<-filter(grass_month, year>2007)

# Calculate 6-month SPEI using the three different CWB values calculated above
SPEI_6mo1 <- spei(data=ts(grass_month2$CWB1_pen,freq=12,start=c(2008,1)), scale=6)
grass_month2$SPEI1_pen<-as.numeric(SPEI_6mo1$fitted)

SPEI_6mo2 <- spei(data=ts(grass_month2$CWB2_pen,freq=12,start=c(2008,1)), scale=6)
grass_month2$SPEI2_pen<-as.numeric(SPEI_6mo2$fitted)

SPEI_6mo3 <- spei(data=ts(grass_month2$CWB_thorn,freq=12,start=c(2008,1)), scale=6)
grass_month2$SPEI_thorn<-as.numeric(SPEI_6mo3$fitted)

# Add a column indicating station (Ameriflux tower) identity
grass_month2$station<-"grassland"


##### SPEI Calcuations: Chihuahuan Desert Ecosystems #####

# Convert date to standard format
names(shrubland)[1]<-"date"
shrubland$date<-mdy(shrubland$date)

# Create year and month columns
shrubland$year<-year(shrubland$date)
shrubland$month<-month(shrubland$date)

# For variables of interest, calculate mean values for each month x year combination
shrub_month<-shrubland %>% group_by(month,year) %>% summarise(P_sum=sum(P_int, na.rm=TRUE),PET1pen_mean=mean(PET1_int, na.rm=TRUE),PET2pen_mean=mean(PET2_int, na.rm=TRUE),Tair_mean=mean(TA_avg, na.rm=TRUE))
shrub_month<-na.omit(shrub_month)

# Calculate PET using the Thornthwaite method
shrub_month$PET_thorn <- as.numeric(thornthwaite(Tave = shrub_month$Tair_mean, lat = 34.35, na.rm = FALSE))

# Calculate climatic water balance (CWB) using both Penman and Thornthwaite PET values, and sort the data frame. Two different Penman values were used (see Supplementary Information for details)
shrub_month <- shrub_month %>% mutate(CWB1_pen=P_sum-PET1pen_mean,CWB2_pen=P_sum-PET2pen_mean, CWB_thorn=P_sum-PET_thorn)
shrub_month<-arrange(shrub_month, year, month)

# Subset the data to just include focal years
shrub_month2<-filter(shrub_month, year>2007)

# Calculate 6-month SPEI using the three different CWB values calculated above
SPEI_6mo1 <- spei(data=ts(shrub_month2$CWB1_pen,freq=12,start=c(2008,1)), scale=6)
shrub_month2$SPEI1_pen<-as.numeric(SPEI_6mo1$fitted)

SPEI_6mo2 <- spei(data=ts(shrub_month2$CWB2_pen,freq=12,start=c(2008,1)), scale=6)
shrub_month2$SPEI2_pen<-as.numeric(SPEI_6mo2$fitted)

SPEI_6mo3 <- spei(data=ts(shrub_month2$CWB_thorn,freq=12,start=c(2008,1)), scale=6)
shrub_month2$SPEI_thorn<-as.numeric(SPEI_6mo3$fitted)

# Add a column indicating station (Ameriflux tower) identity
shrub_month2$station<-"shrubland"


###### Combine the data and examine differences in SPEI between the two PET estimation methods #####

# Combine the data frames from the two Ameriflux towers, and omit rows containing NA values
spei<-bind_rows(grass_month2,shrub_month2)
spei<-na.omit(spei)

# Calculate inverse SPEI for the sake of analyses
spei<- spei %>% mutate(SPEI1_pen=-1*SPEI1_pen, SPEI2_pen=-1*SPEI2_pen, SPEI_thorn=-1*SPEI_thorn)

# Create long form version of data
spei2<-pivot_longer(spei, cols=11:13, names_to = "method", values_to = "spei")

# Calculate a mean Penman-method PET value
spei3<-spei %>% mutate(SPEI_pen_mean=(SPEI1_pen+SPEI2_pen)/2)

# For the spring and monsoon season of each year, calculate the difference between the SPEI values derived from the Thornthwaite and Penman PET values
summary<-spei3 %>% filter(month==5 | month == 9) %>% group_by(year,month,station) %>% summarise(SPEI_diff=SPEI_thorn-SPEI_pen_mean)

# Calculate, across years, the maximum difference between the Penman and Thornthwaite-derived SPEI values for the spring and monsoon seasons in each ecosystem
summary_speidiff<-summary %>% group_by(month, station) %>% summarise(max_diff=max(abs(SPEI_diff)))

# Create a column indicating whether the Thornthwaite-derived SPEI value was higher or lower than the Penman-derived value
summary$diff_direction<-ifelse(summary$SPEI_diff>0,"higher","lower")

# Summarize the results
summary2<- summary %>% group_by(month,station, diff_direction) %>% summarise(count=n())

# Calculate proportion of estimates that were higher vs. lower
summary2$year_count<-ifelse(summary2$station=="grassland",14,15)
summary2$proportion<-summary2$count/summary2$year_count

# Get a summary of the mininum and maximum Penman and Thornthwaite-estimated SPEI values for the spring and monsoon season in each ecosystem
summary <- spei2 %>% filter(month==5 | month == 9) %>% group_by(month,station,method) %>% summarise(max_spei=max(spei, na.rm=TRUE), min_spei=min(spei, na.rm=TRUE), range_spei=max(spei, na.rm=TRUE)-min(spei, na.rm=TRUE))

# Format the data in preparation for graphing trends
forgraph<-spei %>% mutate(SPEI_pen_mean=(SPEI1_pen+SPEI2_pen)/2)
forgraph<-pivot_longer(forgraph, cols=c(13,15), names_to = "method", values_to = "spei")
forgraph <- forgraph %>% filter(month==9 | month == 5)

forgraph$method<-as.factor(forgraph$method)
levels(forgraph$method)
levels(forgraph$method) <- c("Penman","Thornthwaite")

forgraph$month<-as.factor(forgraph$month)
levels(forgraph$month)
levels(forgraph$month) <- c("Spring","Monsoon")

forgraph$station<-as.factor(forgraph$station)
levels(forgraph$station)
levels(forgraph$station) <- c("Plains ecosystem","Chihuahuan Desert ecosystems")

# Penman and Thornthwaite aridity index (SPEI) values as a function of year, for the spring and monsoon season in each ecosytem 
p <- forgraph %>%
  ggplot(aes(x=year,y=spei, color=method)) + geom_point(size=2) + theme_bw() + xlab("Year") + ylab("Aridity index") + labs(color="PET method") + geom_line() + facet_grid(month~station)
p

# save graph
# ggsave("petgraph.jpg",p,width=7,height=5.5,units = c("in"),dpi = 300)
