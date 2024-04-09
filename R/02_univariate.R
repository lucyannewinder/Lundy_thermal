rm(list=ls())

# Proportion of variance models + plot code
# Each age modeled separately 

library(dplyr)
library(MCMCglmm)
library(MasterBayes)
library(tidyr)
library(survival)
library(lubridate)
library(ggplot2)


## --------
## Read in data
## --------
Tb.ped <- read.csv("Data/Processed/ped_for_analysis.csv")
Tb.dat.na <- read.csv("Data/Processed/data_for_analysis.csv")


## ---------
## DAY 2 model
## ---------

# subset data
Age2.tb<-subset(Tb.dat.na,Age==2)
Age2.tb<-droplevels(Age2.tb)

# add brood size
# does every chick have a record? is the number of phenotyped chicks a good measure of brood size?
Age2.tb <- Age2.tb %>%
  group_by(Brood.natal) %>%
  mutate(Natal_brood_size = n())

## what is this exclusion based on? There is no mention of this in the text
Age2.tb <- as.data.frame(subset(Age2.tb, BirdID != "10179"))

# mean center covariates used in models 
Age2.tb$Air.tempC <- scale(Age2.tb$Air.temp, scale=FALSE)
Age2.tb$Natal_brood_sizeC <- scale(Age2.tb$Natal_brood_size, scale=FALSE)


## get inverse A matrix for day 2 chicks
Tb.ped.day2<- prunePed(Tb.ped,unique(Age2.tb$BirdID), make.base = TRUE)
day2.Ainv<-inverseA(Tb.ped.day2)$Ainv

# set priors
prior_day2 <- list(
  R=list(V=1, nu=0.002), 
  G=list(
    BirdID=list(V=1, nu=1, alpha.mu=0, alpha.V=1000),
    Brood.natal=list(V=1, nu=1, alpha.mu=0, alpha.V=1000)
  )
)
a=15

mod_uni_day2 <- MCMCglmm(max.temp~Air.tempC+
                            as.factor(Year)+
                            Natal_brood_sizeC,
                     random=~BirdID+
                       Brood.natal,
                     data=Age2.tb,
                     nitt=13000*a, thin=10*a, burnin=3000*a, 
                     ginverse=list(BirdID=day2.Ainv), prior=prior_day2, verbose=FALSE)
summary(mod_uni_day2)
#plot(mod_uni_day2)


### -----
## DAY 5 model
### -----

# subset data
Age5.tb<-subset(Tb.dat.na,Age==5)
Age5.tb<-droplevels(Age5.tb)

# add brood size
Age5.tb <- Age5.tb %>%
  group_by(Brood.rearing) %>%
  mutate(Rearing_brood_size = n())

### ?????
values_to_remove<-c("10179" ,"10821", "10827", "10829" ,"10911")
Age5.tb <- subset(Age5.tb, !BirdID %in% values_to_remove)
Age5.tb<-as.data.frame(Age5.tb)

# mean center covariates used in models 
Age5.tb$Air.tempC <- scale(Age5.tb$Air.temp, scale=FALSE)
Age5.tb$Rearing_brood_sizeC <- scale(Age5.tb$Rearing_brood_size, scale=FALSE)

## get inverse A matrix for day 5 chicks
Tb.ped.day5<- prunePed(Tb.ped,unique(Age5.tb$BirdID), make.base = TRUE)
day5.Ainv<-inverseA(Tb.ped.day5)$Ainv

## Priors for day 5 plus models
prior_day5plus<-list(
  R=list(V=1, nu=0.002), 
  G=list(
    BirdID=list(V=1, nu=1, alpha.mu=0, alpha.V=1000),
    Brood.natal=list(V=1, nu=1, alpha.mu=0, alpha.V=1000),
    Brood.rearing=list(V=1, nu=1, alpha.mu=0, alpha.V=1000)
  )
)

mod_uni_day5 <- MCMCglmm(max.temp~Air.tempC+
                            as.factor(Year)+
                            Rearing_brood_sizeC,
                          random=~BirdID+ 
                            Brood.rearing+
                            Brood.natal,
                          data=Age5.tb,
                          nitt=13000*a, thin=10*a, burnin=3000*a, 
                          ginverse=list(BirdID=day5.Ainv), prior=prior_day5plus, verbose=FALSE)
summary(mod_uni_day5)
#plot(mod_uni_day5)


#### ---------
## DAY 10 model
#### ---------

# subset data
Age10.tb<-subset(Tb.dat.na,Age==10)
Age10.tb<-droplevels(Age10.tb)

# add brood size
Age10.tb <- Age10.tb %>%
  group_by(Brood.rearing) %>%
  mutate(Rearing_brood_size = n())

### ?????
values_to_remove<-c("9989" , "9993" , "9998"  ,"10011" ,"10155" ,"10313", "10314" ,"10315" ,"10806" ,"10857" ,"10869", "10887", "10894")
Age10.tb <- subset(Age10.tb, !BirdID %in% values_to_remove)
Age10.tb<-as.data.frame(Age10.tb)

# mean center covariates used in models 
Age10.tb$Air.tempC <- scale(Age10.tb$Air.temp, scale=FALSE)
Age10.tb$Rearing_brood_sizeC <- scale(Age10.tb$Rearing_brood_size, scale=FALSE)

## get inverse A matrix for day 10 chicks
Tb.ped.day10<- prunePed(Tb.ped,unique(Age10.tb$BirdID), make.base = TRUE)
day10.Ainv<-inverseA(Tb.ped.day10)$Ainv


mod_uni_day10 <- MCMCglmm(max.temp~Air.tempC+
                            as.factor(Year)+
                            Rearing_brood_sizeC,
                          random=~BirdID+ 
                            Brood.rearing+
                            Brood.natal,
                          data=Age10.tb,
                          nitt=13000*a, thin=10*a, burnin=3000*a, 
                          ginverse=list(BirdID=day10.Ainv), prior=prior_day5plus, verbose=FALSE)
summary(mod_uni_day10)
#plot(mod_uni_day10)


### ---------
## DAY 12 model
### ---------

# subset data
Age12.tb<-subset(Tb.dat.na,Age==12)
Age12.tb<-droplevels(Age12.tb)

# add brood size
Age12.tb <- Age12.tb %>%
  group_by(Brood.rearing) %>%
  mutate(Rearing_brood_size = n())

### ?????
values_to_remove<-c("10011", "10155", "10313" ,"10314", "10315", "10806", "10857", "10869", "10887", "10894")
Age12.tb <- subset(Age12.tb, !BirdID %in% values_to_remove)
Age12.tb<-as.data.frame(Age12.tb)

# mean center covariates used in models 
Age12.tb$Air.tempC <- scale(Age12.tb$Air.temp, scale=FALSE)
Age12.tb$Rearing_brood_sizeC <- scale(Age12.tb$Rearing_brood_size, scale=FALSE)

## get inverse A matrix for day 10 chicks
Tb.ped.day12<- prunePed(Tb.ped,unique(Age12.tb$BirdID), make.base = TRUE)
day12.Ainv<-inverseA(Tb.ped.day12)$Ainv


mod_uni_day12 <- MCMCglmm(max.temp~Air.tempC+
                            as.factor(Year)+
                            Rearing_brood_sizeC,
                          random=~BirdID+ 
                            Brood.rearing+
                            Brood.natal,
                          data=Age12.tb,
                          nitt=13000*a, thin=10*a, burnin=3000*a, 
                          ginverse=list(BirdID=day12.Ainv), prior=prior_day5plus, verbose=FALSE)
summary(mod_uni_day12)
#plot(mod_uni_day12)



### ---------
### Save all model outputs
### ---------


save(
  mod_uni_day2,
  mod_uni_day5,
  mod_uni_day10,
  mod_uni_day12,
  file = "Data/Output/uni_models.Rdata")


