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
library(pedtricks)
library(beepr)


## --------
## Read in data
## --------
Tb.ped <- read.csv("Data/Processed/ped_for_analysis.csv")
Tb.dat.na <- read.csv("Data/Processed/data_for_analysis.csv")

Tb.dat.na$BirdID<-as.factor(Tb.dat.na$BirdID)
Tb.dat.na$doy<-as.factor(Tb.dat.na$doy)


## ---------
## DAY 2 model
## ---------

# subset data
Age2.tb<-subset(Tb.dat.na,Age==2)
Age2.tb<-droplevels(Age2.tb)

# mean center covariates used in models 
Age2.tb$Air.tempC <- scale(Age2.tb$Air.temp, scale=FALSE)
Age2.tb$Currentage_brood_sizeC <- scale(Age2.tb$Currentage_brood_size, scale=FALSE)


## get inverse A matrix for day 2 chicks
Tb.ped.day2<- prunePed(Tb.ped,unique(Age2.tb$BirdID), make.base = TRUE)
day2.Ainv<-inverseA(Tb.ped.day2)$Ainv

# # pedigree stats
# ped_stats(Tb.ped.day2,retain="informative")
# ped_table<-summary(ped_stats(Tb.ped.day2,retain="informative"))
# write.csv(ped_table,"ped_summary2.csv")

# set priors
prior_day2 <- list(
  R=list(V=1, nu=0.002), 
  G=list(
    BirdID=list(V=1, nu=1, alpha.mu=0, alpha.V=1000),
    Brood.natal=list(V=1, nu=1, alpha.mu=0, alpha.V=1000),
    doy=list(V=1, nu=1, alpha.mu=0, alpha.V=1000)
  )
)
a=15

mod_uni_day2 <- MCMCglmm(max.temp~Air.tempC+ #max.temp+273.15~Air.tempC+  #for conversion to kelvin to calculate evolvability 
                            as.factor(Year)+
                           Currentage_brood_sizeC,
                     random=~BirdID+
                       Brood.natal+
                       doy,
                     data=Age2.tb,
                     nitt=13000*a, thin=10*a, burnin=3000*a, 
                     ginverse=list(BirdID=day2.Ainv), prior=prior_day2, verbose=FALSE)
summary(mod_uni_day2)
apply(mod_uni_day2$Sol, 2, median)
beep(1)
#plot(mod_uni_day2)


### -----
## DAY 5 model
### -----

# subset data
Age5.tb<-subset(Tb.dat.na,Age==5)
Age5.tb<-droplevels(Age5.tb)

# mean center covariates used in models 
Age5.tb$Air.tempC <- scale(Age5.tb$Air.temp, scale=FALSE)
Age5.tb$Currentage_brood_sizeC <- scale(Age5.tb$Currentage_brood_size, scale=FALSE)

## get inverse A matrix for day 5 chicks
Tb.ped.day5<- prunePed(Tb.ped,unique(Age5.tb$BirdID), make.base = TRUE)
day5.Ainv<-inverseA(Tb.ped.day5)$Ainv

# # pedigree stats
# ped_stats(Tb.ped.day5,retain="informative")
# ped_table<-summary(ped_stats(Tb.ped.day5,retain="informative"))
# write.csv(ped_table,"ped_summary5.csv")

## Priors for day 5 plus models
prior_day5plus<-list(
  R=list(V=1, nu=0.002), 
  G=list(
    BirdID=list(V=1, nu=1, alpha.mu=0, alpha.V=1000),
    Brood.natal=list(V=1, nu=1, alpha.mu=0, alpha.V=1000),
    Brood.rearing=list(V=1, nu=1, alpha.mu=0, alpha.V=1000),
    doy=list(V=1, nu=1, alpha.mu=0, alpha.V=1000)
  )
)

mod_uni_day5 <- MCMCglmm(max.temp~Air.tempC+ #max.temp+273.15~Air.tempC+  #for conversion to kelvin to calculate evolvability  
                            as.factor(Year)+
                           Currentage_brood_sizeC,
                          random=~BirdID+ 
                            Brood.rearing+
                            Brood.natal+
                           doy,
                          data=Age5.tb,
                          nitt=13000*a, thin=10*a, burnin=3000*a, 
                          ginverse=list(BirdID=day5.Ainv), prior=prior_day5plus, verbose=FALSE)
summary(mod_uni_day5)
apply(mod_uni_day5$Sol, 2, median)
beep(1)
#plot(mod_uni_day5)


#### ---------
## DAY 10 model
#### ---------

# subset data
Age10.tb<-subset(Tb.dat.na,Age==10)
Age10.tb<-droplevels(Age10.tb)

# mean center covariates used in models 
Age10.tb$Air.tempC <- scale(Age10.tb$Air.temp, scale=FALSE)
Age10.tb$Currentage_brood_sizeC <- scale(Age10.tb$Currentage_brood_size, scale=FALSE)

## get inverse A matrix for day 10 chicks
Tb.ped.day10<- prunePed(Tb.ped,unique(Age10.tb$BirdID), make.base = TRUE)
day10.Ainv<-inverseA(Tb.ped.day10)$Ainv

# # pedigree stats
# ped_stats(Tb.ped.day10,retain="informative")
# ped_table<-summary(ped_stats(Tb.ped.day10,retain="informative"))
# write.csv(ped_table,"ped_summary10.csv")


mod_uni_day10 <- MCMCglmm(max.temp~Air.tempC+ #max.temp+273.15~Air.tempC+  #for conversion to kelvin to calculate evolvability.  
                            as.factor(Year)+
                            Currentage_brood_sizeC,
                          random=~BirdID+ 
                            Brood.rearing+
                            Brood.natal+
                            doy,
                          data=Age10.tb,
                          nitt=13000*a, thin=10*a, burnin=3000*a, 
                          ginverse=list(BirdID=day10.Ainv), prior=prior_day5plus, verbose=FALSE)
summary(mod_uni_day10)
apply(mod_uni_day10$Sol, 2, median)
beep(1)
#plot(mod_uni_day10)


### ---------
## DAY 12 model
### ---------

# subset data
Age12.tb<-subset(Tb.dat.na,Age==12)
Age12.tb<-droplevels(Age12.tb)

# mean center covariates used in models 
Age12.tb$Air.tempC <- scale(Age12.tb$Air.temp, scale=FALSE)
Age12.tb$Currentage_brood_sizeC <- scale(Age12.tb$Currentage_brood_size, scale=FALSE)

## get inverse A matrix for day 10 chicks
Tb.ped.day12<- prunePed(Tb.ped,unique(Age12.tb$BirdID), make.base = TRUE)
day12.Ainv<-inverseA(Tb.ped.day12)$Ainv

# # pedigree stats
# ped_stats(Tb.ped.day12,retain="informative")
# ped_table<-summary(ped_stats(Tb.ped.day12,retain="informative"))
# write.csv(ped_table,"ped_summary12.csv")


mod_uni_day12 <- MCMCglmm(max.temp~Air.tempC+ #max.temp+273.15~Air.tempC+  #for conversion to kelvin to calculate evolvability.  
                            as.factor(Year)+
                            Currentage_brood_sizeC,
                          random=~BirdID+ 
                            Brood.rearing+
                            Brood.natal+
                            doy,
                          data=Age12.tb,
                          nitt=13000*a, thin=10*a, burnin=3000*a, 
                          ginverse=list(BirdID=day12.Ainv), prior=prior_day5plus, verbose=FALSE)
summary(mod_uni_day12)
apply(mod_uni_day12$Sol, 2, median)
beep(1)
#plot(mod_uni_day12)



#evolvability - change above models so max.temp+273.15 to convert to K.
VA_2<-(mod_uni_day2$VCV[,"BirdID"])
z_2<-mean(Age2.tb$max.temp + 273.15, na.rm = TRUE)
e_2<-VA_2/(z_2^2)
ci_2<-quantile(e_2, c(0.025, 0.975))
mean(e_2)
ci_2

VA_5<-(mod_uni_day5$VCV[,"BirdID"])
z_5<-mean(Age5.tb$max.temp + 273.15, na.rm = TRUE)
e_5<-VA_5/(z_5^2)
ci_5<-quantile(e_5, c(0.025, 0.975))
mean(e_5)
ci_5

VA_10<-(mod_uni_day10$VCV[,"BirdID"])
z_10<-mean(Age10.tb$max.temp + 273.15, na.rm = TRUE)
e_10<-VA_10/(z_10^2)
ci_10<-quantile(e_10, c(0.025, 0.975))
mean(e_10)
ci_10

VA_12<-(mod_uni_day12$VCV[,"BirdID"])
z_12<-mean(Age12.tb$max.temp + 273.15, na.rm = TRUE)
e_12<-VA_12/(z_12^2)
ci_12<-quantile(e_12, c(0.025, 0.975))
mean(e_12)
ci_12

### ---------
### Save all model outputs
### ---------


save(
  mod_uni_day2,
  mod_uni_day5,
  mod_uni_day10,
  mod_uni_day12,
  file = "Data/Output/uni_models.Rdata")


# cor(modM1_Tbher_2$VCV)
# cor(modM1_Tbher_5$VCV)
# cor(modM1_Tbher_10$VCV)
# cor(modM1_Tbher_12$VCV)

