# this file contains the models for RQ3
# Maartje Oostdijk, Laura G. Elsler

### set up
# clear workspace
rm(list = ls())
graphics.off()

# Load libraries
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(glmmTMB)
library(sandwich)
library(lmtest)

# Directory
datadir <- "~/Documents/SESYNC/Files/FISHMAR-data/rq3/final_push/analysis/data"

# Read data
data<- read.csv(file.path(datadir, "20200830_stage2.csv"), as.is=T) %>% rename(coopid=X...coopid) # Annual data 2-stage regression
data_trade <- read.csv(file.path(datadir, "CT_fish_trade92.csv"), as.is=T)%>% #international trade data
  filter(group_name=="rocklobster") %>% mutate(q=q*1000,v=v*1000) # tons -> kg
fmsy <- read.csv(file.path(datadir, "results_timeseries_MSY_rq3.csv"), as.is=T) %>% rename(coopid=Coop.Id) # FMSY timeseries
inflation <- read.csv(file.path(datadir, "inflation.csv"), as.is=T) # Extract functionality from here
data_functionality <- read.csv(file.path(datadir, "data_monthly_lobsters.csv"), as.is=T) # Extract functionality from here
#gammas = read.csv('~/Documents/SESYNC/Files/FISHMAR-data/rq3/final_push/data/gammas2.csv') %>% rename(coopid = Coop_Id)      # @MO: still needed?

# Save for statistics section
statistics1 = data_trade
statistics2 = data_functionality

###################################################################################################################
###########################################   DATA PREPARATION   ##################################################
###################################################################################################################


# Clean international trade data and inflation correct price
data_trade= data_trade%>%
  rename(year=t)%>%
  mutate(price = v/q)%>%
  filter(!is.na(price))%>%
  left_join(inflation)%>%
  mutate(price = price*inflation_corrector_2016)%>%
  group_by(year)%>%
  summarise(av_int_price = mean(price))%>%
  filter(!is.na(av_int_price))

# Create list of cooperatives and functionality
functionality <- data_functionality %>% 
  distinct(coopid,Functionality)

functionality$FBins <- cut(functionality$Functionality, breaks = c(5,6,7,8,9,20), 
                  labels = c("5 < F < 6","6 < F < 7","7 < F < 8","8 < F < 9","F > 9"))

# Join with fmsy data
fmsy <- fmsy %>% filter(Species=='Lobsters') %>% rename(Bbmsy=B.Bmsy,Ffmsy=F.Fmsy,year=yr) %>% select(coopid,year,Bbmsy,Ffmsy)
data = data%>%
  left_join(fmsy)

# Left join functionality to Mexican data
names(data)[1] <- "coopid"
data<- data%>% 
  left_join(functionality, by = c('coopid')) %>% 
  #left_join(gammas)%>%                                                 # @MO: still needed? 
  left_join(data_trade)%>%
  mutate(overfished = as.factor(ifelse(Bbmsy <0.8, 1, 0)),
         high_overfished = as.factor(ifelse(Bbmsy <0.5, 1, 0)))




#######################################   TEST BIOMASS DATA   #####################################################

hist(data$Bbmsy, breaks=200) # histogram of BBmsy
shapiro.test(data$Bbmsy) # normally distributed? 
# result: W = 0.97435, p-value = 5.609e-05, very small values reject the null hypothesis 
qqnorm(data$Bbmsy) 
# result: heavy tails, not normal

###################################################################################################################
###########################################   REGRESSION MODELS   #################################################
###################################################################################################################

################################################   TABLE 1    #####################################################

# Model bmsy
m1= lm(log((Bbmsy)) ~ scale(Functionality) * scale(av_int_price) + as.factor(year), data)
coeftest(m1, vcov = NeweyWest(m1, lag = 3, prewhite = FALSE))
coeftest(m1, vcov = kernHAC(m1, kernel = "Bartlett", bw = 3, prewhite = FALSE, adjust = FALSE)) 

# Model total catch lobsters
m2 <- glm(total_catch_lobsters~ scale(Functionality) * scale(av_int_price) + as.factor(year) , data, family=Gamma(link="log"))
summary(m2)
coeftest(m2, vcov = NeweyWest(m2, lag = 3, prewhite = FALSE))

# Model total value lobsters
m3<- glm(total_value_lobsters~ scale(Functionality) * scale(av_int_price) + as.factor(year), data, family=Gamma(link="log"))
summary(m3)
coeftest(m3, vcov = NeweyWest(m3, lag = 3, prewhite = FALSE))

################################################   TABLE 2    #####################################################

# Model local price
m4= lm(log((ic_price_lobsters)) ~ scale(av_int_price) * scale(Functionality) + as.factor(year), data)
coeftest(m4, vcov = NeweyWest(m4, lag = 3, prewhite = FALSE))

# Model local price
m5= lm(log((ic_price_lobsters)) ~ scale(av_int_price) + as.factor(year), data)
coeftest(m5, vcov = NeweyWest(m5, lag = 3, prewhite = FALSE))

rm(m1,m2,m3,m4,m5)

###################################################################################################################
#############################################   PLOTTING DATA   ###################################################
###################################################################################################################

# Color scale
colorF <- c('#9b7e2e','#7f7d32','#5d7a38','#3b784f','#1d7662')

################################################   FIGURE 1    ####################################################
# Conceptual figure called: Co-op_Infographic_DRAFT-V3 (1).jpg

################################################   FIGURE 2    ####################################################
# Cooperative functionality and interannual standard deviation of local price (MXN) per cooperative. 

#    monthly data has no inflation corrected price.                                                        @LE: go do it!
a = statistics2%>% mutate(price_lobsters=value_lobsters/catch_lobsters) %>%
  group_by(coopid, year)%>%
  filter(!is.na(Functionality) & !is.na(price_lobsters))%>%
  summarise(sdprice = sd(price_lobsters), av_price=mean(price_lobsters), av_catch=mean(catch_lobsters))%>% 
  left_join(functionality) 

f2 <- ggplot(data = a, aes(x=Functionality, y=sdprice, color = FBins)) +
  geom_point(size=2) +
  theme_classic() +
  ylab("Standard deviation local price") +
  xlab("Functionality") +
  scale_color_manual(values=colorF) +
  theme(text=element_text(size=18), #change font size of all text
        legend.title = element_blank())
# f2
# ggsave(plot = f2, width = 7, height = 5, dpi = 600, filename ='~/Documents/SESYNC/Files/FISHMAR-data/rq3/final_push/analysis/figures/2_sdprice_funct.png')



##############################################    FIGURE 3A    ##################################################
# Average monthly local price of lobster per cooperative from 2000 to 2016. 
f3a <- ggplot(data = a, aes(x=Functionality, y=av_price, color = FBins)) +
  geom_point(size=2) +
  theme_classic() +
  ylab("Average monthly price per year (MXN)") +
  xlab("Functionality") +
  ggtitle("A.") +
  scale_color_manual(values=colorF) +
  theme(text=element_text(size=18), #change font size of all text
        plot.title = element_text(hjust = 0),
        legend.position = "none")


# f3a
# ggsave(plot = f3a, width = 7, height = 5, dpi = 600, filename ='~/Documents/SESYNC/Files/FISHMAR-data/rq3/final_push/analysis/figures/3A_price_funct.png')

##############################################    FIGURE 3B    ##################################################
# Average monthly catch of lobster per cooperative from 2000 to 2016. 
f3b <- ggplot(data = a, aes(x=Functionality, y=av_catch, color = FBins)) +
  geom_point(size=2) +
  theme_classic() +
  ylab("Average monthly catch per year (kg)") +
  xlab("Functionality") +
  ggtitle("B.") +
  scale_color_manual(values=colorF) +
  theme(text=element_text(size=18), #change font size of all text
        plot.title = element_text(hjust = 0),
        legend.title = element_blank(),
        legend.position = "none")
# f3b
# ggsave(plot = f3b, width = 7, height = 5, dpi = 600, filename ='~/Documents/SESYNC/Files/FISHMAR-data/rq3/final_push/analysis/figures/3B_catch_funct.png')



##############################################   SI FIGURE 1    ###################################################
# Average annual lobster price per cooperative functionality level from 2000 to 2016. 
b <- data %>% filter(!is.na(Bbmsy)) %>% 
  group_by(FBins, year) %>%
  summarise(av_price = mean(ic_price_lobsters),av_bmsy=mean(Bbmsy))
b <- b %>% rename(Functionality=FBins)

si1 <- ggplot(b, mapping = aes(year, av_price, color = Functionality)) +
  geom_line(size = 2) + ylab('Average annual price per kg (MXN)') + xlab('Year') + 
  scale_color_manual(values=colorF) +
  theme_bw() + theme(strip.background = element_blank(), strip.text.x = element_blank()) +
  theme(legend.position="top",
        text=element_text(size=18), #change font size of all text
        legend.text=element_text(size=15), #change font size of legend text
        legend.title=element_text(size=15)) #change font size of legend title  
# si1
# ggsave(plot = si1, width = 10, height = 7, dpi = 600, filename ='~/Documents/SESYNC/Files/FISHMAR-data/rq3/final_push/analysis/figures/SI1_price_ts.png')

##############################################   SI FIGURE 2    ###################################################
# Average annual lobster population levels per cooperative functionality level from 2000 to 2016. 

si2 <- ggplot(b, mapping = aes(year, av_bmsy, color = Functionality)) +
  geom_line(size = 2) + ylab('B/Bmsy') + xlab('Year') + 
  scale_color_manual(values=colorF) +
  theme_bw() + theme(strip.background = element_blank(), strip.text.x = element_blank()) +
  theme(legend.position="top",
        text=element_text(size=18), #change font size of all text
        legend.text=element_text(size=15), #change font size of legend text
        legend.title=element_text(size=15)) #change font size of legend title  
# si2
# ggsave(plot = si2, width = 10, height = 7, dpi = 600, filename ='~/Documents/SESYNC/Files/FISHMAR-data/rq3/final_push/analysis/figures/SI2_bmsy_ts.png')

rm(a,b,f2,si1,si2,f3a,f3b)

###################################################################################################################
#############################################   STATISTICS FOR TEXT   #############################################
###################################################################################################################

##############################################   INTRODUCTION    ##################################################
# we selected lobster-trading cooperatives 
unique(functionality$coopid)

# Lobster is a high-value species (avg. price USD XX)   
mean(data_trade$av_int_price)
#! we could use the weighted average - maybe better                                                                 @LE: go do it!

################################################   METHODS    ######################################################
# we selected cooperatives that reported harvesting lobster (N=40) 
unique(functionality$coopid)

#  reduces volatility of individual transactions (std. before averaging = # USD, std. after averaging = # USD).   
# international price averaged
sd(data_trade$av_int_price) # averaged price

# international price original
statistics1%>% mutate(price=v/q) %>%
  filter(!is.na(price))%>%
  summarise(sdprice = sd(price)) 

# average international traded volume = #; average Mexican caught volume = #
# Mean total catch volume/annum of Mexican cooperatives 
a <- aggregate( total_catch_lobsters ~ year, data, sum ) # calculates catch per year
mean(a$total_catch_lobsters)
# Mean total international trade volume of rock lobster per year 
b <- aggregate( q ~ t, statistics1, sum ) # calculates trade per year
mean(b$q) 


###############################################   RESULTS 3.1    ######################################################
# The average functionality of lobster cooperatives was high compared to cooperatives targeting other species (lobster cooperatives averaged 7.79; avg. across all species was 7.19)
x <- data_functionality %>% select(coopid,Functionality) %>% distinct()
mean(x$Functionality)

# The lowest functioning lobster cooperative had a functionality level of 4.71 and the highest functioning cooperative of 9.61
min(x$Functionality)
max(x$Functionality)

# Average B/BMSY of cooperatives was # (std. #); increasing functionality by one point was associated with a # increase in B/BMSY.
mean(data$Bbmsy,na.rm=TRUE)
sd(data$Bbmsy,na.rm=TRUE)

# An increase in cooperative functionality is associated with an increase in catch volume of # kg (R squared=#).
c <- aggregate( total_catch_lobsters ~ coopid, data, mean ) # calculates catch per functionality
# either use lowest and highest functioning coops' catch for illustration - manually
# or use this linear model
c <- c %>% left_join(functionality)
functionality_catch <- lm(total_catch_lobsters ~ Functionality, data = c)
summary(functionality_catch)

###############################################   RESULTS 3.3    ######################################################

# Across the years the local prices cooperatives received were equivalent to #% of the international trade price. 
# calculate international lobster price per year
d <- statistics1%>% mutate(price=v/q) %>% rename(year=t)
d <- aggregate( price ~ year, d, mean ) 
# convert international lobster price into MXN
conversion <- read.csv(file.path(datadir, "Consulta_20220119-015051678.csv"), as.is=T) %>% # Original file 'BDM_USD_90_20.xlsx' timeseries of MXN to USD from the Banco de Mexico
  rename(rate=X...SF63528) %>% 
  group_by(year) %>%  
  summarise(conversion = mean(rate)) %>%
  filter(!is.na(conversion))
# convert
d <- d %>% left_join(conversion) %>%
  mutate(price_mxn=price*conversion) 
# calculates local mean price per year
e <- aggregate( ic_price_lobsters ~ year, data, mean ) 
# compare local and international price
d <- d %>% left_join(e) %>% mutate(price_rate=ic_price_lobsters/price_mxn) %>% filter(!is.na(price_rate))
mean(d$price_rate)


# The average price cooperatives received for lobsters was # MXN/kg and the international trade price was # MXN/kg. 
mean(d$ic_price_lobsters)
mean(d$price_mxn)

###############################################   RESULTS 3.4    ######################################################

# The yearly average local prices of all cooperatives increased between 2000-2016, from # to # MXN. 
f <- aggregate( ic_price_lobsters ~ year, data, mean ) 


##############################################   SUPPLEMENTARY    #####################################################
#   SI Table 1. Summary statistics of the lobster cooperatives unless otherwise specified.
# counts 
unique(functionality$coopid)
length(statistics2$catch_lobsters)
length(data$Bbmsy)
length(data$Ffmsy)

# standard deviation
g <- data %>% 
  filter(!is.na(Bbmsy),!is.na(Bbmsy)) %>% 
  summarise(sd_funct = sd(Functionality), sd_catch=sd(total_catch_lobsters), sd_bbmsy=sd(Bbmsy), sd_ffmsy=sd(Bbmsy))

# avg, min, max
f <- summary(data)
f

rm(a,b,c,d,e,f,g)



###################################################################################################################
######################## USEFUL ADDITIONAL ANALYSIS NOT INCLUDED IN THE ARTICLE   #################################
###################################################################################################################

# Linear trend + confidence interval for functionality and bbmsy
s1 <- ggplot(data, aes(x=Functionality, y=Bbmsy)) +
  geom_point() +geom_smooth(color="black", method=lm)+theme_classic() 
s1


###########################################     REMOVE

test.1 = read.csv('~/Documents/SESYNC/Files/FISHMAR-data/mexico/processed/RQ3_monthly_dataset_lobster.csv', as.is=T)
test.2 = read.csv('~/Documents/SESYNC/Files/FISHMAR-data/mexico/processed/laura/RQ3_monthly_dataset_lobster_lge_ag.csv', as.is=T)
test.3 = read.csv('~/Documents/SESYNC/Files/FISHMAR-data/mexico/processed/data_monthly.csv', as.is=T) # full data of data_monthly_lobsters.csv

# total global aggregates, comtrade
a <- aggregate( v ~ t, data_trade, sum ) # calculates value per year
b <- aggregate( q ~ t, data_trade, sum ) # calculates catch per year

# average cooperative functionality
c <- data_functionality %>% select(coopid,Functionality) %>% distinct()
mean(c$Functionality)

# are the annual catch data from original to stage2 data
keeps <- c('Coop_Id','functionality.x','year','month','catch_lobsters','value_lobsters')
test.3 <- test.3[keeps]
f <- data %>% select(coopid,year,total_catch_lobsters,total_value_lobsters)
d <- test.3 %>% group_by(Coop_Id,year) %>% filter(!is.na(catch_lobsters),!is.na(value_lobsters)) %>% # remove if column is NA
  summarise_at(vars(catch_lobsters),list(newcol = sum)) %>% 
  left_join(f) 
e <- test.3 %>% group_by(Coop_Id,year) %>% filter(!is.na(catch_lobsters),!is.na(value_lobsters)) %>% # remove if column is NA
  summarise_at(vars(value_lobsters),list(newcol = sum)) %>% rename(coopid=Coop_Id) %>% 
  left_join(f) # all correct!

# exchange rate new file: Consulta_20220119-015051678.csv
# inflation corrected? 

rm(a,b,c,d,e,f,test.1,test.2,test.3)

