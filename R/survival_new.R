rm(list=ls())

#---- survival analysis new ----
# selective disappearance for body temp + mass at each age

library(dplyr)
library(lme4)
library(tidyr)
library(survival)
library(lubridate)



## --------
## Read in data
## --------

Tb.dat.na <- read.csv("Data/Processed/data_for_analysis.csv")

Tb.dat.na$Age<-as.factor(Tb.dat.na$Age)

# mean center covariates used in models 
Tb.dat.na$MassC <- scale(Tb.dat.na$Mass, scale=FALSE)
Tb.dat.na$max.tempC <- scale(Tb.dat.na$max.temp, scale=FALSE)
Tb.dat.na$Air.tempC <- scale(Tb.dat.na$Air.temp, scale=FALSE)

head(Tb.dat.na)

#Cuts for survival analysis

Tb.dat.na$Age.num<-as.numeric(Tb.dat.na$Age)
Tb.dat.cut<-survSplit(Surv(Age.num,surv)~.,data=Tb.dat.na,cut=c(2,5,10,12),episode="interval")
Tb.dat.cut$survnew<-ifelse(Tb.dat.cut$surv==1,0,1)
head(Tb.dat.cut)
nrow(Tb.dat.na)
nrow(Tb.dat.cut)
subset(Tb.dat.na,BirdID=="9989")
subset(Tb.dat.cut,BirdID=="9989")


Tb.dat.na$survival <- ifelse(Tb.dat.na$Age == Tb.dat.na$maxage,0,1)


Age2.tbsurv<-subset(Tb.dat.na,Age==2)
Age2.tbsurv<-droplevels(Age2.tbsurv)

Age5.tbsurv<-subset(Tb.dat.na,Age==5)
Age5.tbsurv<-droplevels(Age5.tbsurv)

Age10.tbsurv<-subset(Tb.dat.na,Age==10)
Age10.tbsurv<-droplevels(Age10.tbsurv)

Age12.tbsurv<-subset(Tb.dat.na,Age==12)
Age12.tbsurv<-droplevels(Age12.tbsurv)




mod_age2_tbsurv<-glmer(surv~max.temp+
                         Mass+
                         as.factor(Year)+
                         (1|Brood.natal),
                       data=Age2.tbsurv,
                       family = "binomial")
summary(mod_age2_tbsurv)

mod_age5_tbsurv<-glmer(surv~max.temp+
                         Mass+
                         as.factor(Year)+
                         (1|Brood.natal)+
                         (1|Brood.rearing),
                       data=Age5.tbsurv,
                       family = "binomial")
summary(mod_age5_tbsurv)

mod_age10_tbsurv<-glmer(surv~max.temp+
                          Mass+
                          as.factor(Year)+
                          (1|Brood.natal)+
                          (1|Brood.rearing),
                        data=Age10.tbsurv,
                        family = "binomial")
summary(mod_age10_tbsurv)

mod_age12_tbsurv<-glmer(surv~max.temp+
                          Mass+
                          # as.factor(Year)+
                          (1|Brood.natal)+
                          (1|Brood.rearing),
                        data=Age12.tbsurv,
                        family = "binomial")
summary(mod_age12_tbsurv)
