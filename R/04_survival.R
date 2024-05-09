#---- survival analysis new ----
# selective disappearance for body temp + mass at each age

library(dplyr)
library(MCMCglmm)
library(MasterBayes)
library(tidyr)
library(survival)
library(lubridate)
library(lme4)

Tb.dat.na <- read.csv("Data/Processed/data_for_analysis.csv")

# subset into ages
Age2.tbsurv<-subset(Tb.dat.na,Age==2)
Age2.tbsurv<-droplevels(Age2.tbsurv)
Age5.tbsurv<-subset(Tb.dat.na,Age==5)
Age5.tbsurv<-droplevels(Age5.tbsurv)
Age10.tbsurv<-subset(Tb.dat.na,Age==10)
Age10.tbsurv<-droplevels(Age10.tbsurv)
Age12.tbsurv<-subset(Tb.dat.na,Age==12)
Age12.tbsurv<-droplevels(Age12.tbsurv)

# change survival so 1=alive next year 0=dead (currently other way around)
Age2.tbsurv$survnew<-ifelse(Age2.tbsurv$surv==1,0,1)
Age5.tbsurv$survnew<-ifelse(Age5.tbsurv$surv==1,0,1)
Age10.tbsurv$survnew<-ifelse(Age10.tbsurv$surv==1,0,1)
Age12.tbsurv$survnew<-ifelse(Age12.tbsurv$surv==1,0,1)


mod_age2_tbsurv<-glmer(survnew~max.temp+
                         Mass+
                         as.factor(Year)+
                         (1|Brood.natal),
                       data=Age2.tbsurv,
                       family = "binomial")
summary(mod_age2_tbsurv)

mod_age5_tbsurv<-glmer(survnew~max.temp+
                         Mass+
                         as.factor(Year)+
                         (1|Brood.natal)+
                         (1|Brood.rearing),
                       data=Age5.tbsurv,
                       family = "binomial")
summary(mod_age5_tbsurv)

mod_age10_tbsurv<-glmer(survnew~max.temp+
                          Mass+
                          as.factor(Year)+
                          (1|Brood.natal)+
                          (1|Brood.rearing),
                        data=Age10.tbsurv,
                        family = "binomial")
summary(mod_age10_tbsurv)

mod_age12_tbsurv<-glmer(survnew~max.temp+
                          Mass+
                          as.factor(Year)+
                          (1|Brood.natal)+
                          (1|Brood.rearing),
                        data=Age12.tbsurv,
                        family = "binomial")
summary(mod_age12_tbsurv)

