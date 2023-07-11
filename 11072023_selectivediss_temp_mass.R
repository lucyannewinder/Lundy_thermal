#---- survival analysis new ----
# selective dissappearance for body temp + mass at each age

library(dplyr)
library(MCMCglmm)
library(MasterBayes)
library(tidyr)
library(survival)
library(lubridate)




Tb.dat<-read.csv("Tb.dat.csv")
str(Tb.dat)
Tb.dat$Age<-as.factor(Tb.dat$Age)

age2.nat.cf<-subset(Tb.dat,Age==2) %>%
  dplyr::select(BirdID,
                Brood.natal=Brood.code,
                Brood.rearing=CF.brood.no)

age2.nat.cf<-age2.nat.cf %>%
  mutate(Brood.rearing=coalesce(Brood.rearing,Brood.natal))

Tb.dat<-merge(Tb.dat,age2.nat.cf,by="BirdID")

Tb.dat$BirdID<-as.character(Tb.dat$BirdID)
Tb.dat$ind.id<-Tb.dat$BirdID

Tb.dat.na<-Tb.dat[!is.na(Tb.dat$max.temp), ]
Tb.dat.na$BirdID<-as.factor(Tb.dat.na$BirdID)
Tb.dat.na$ind.id<-as.factor(Tb.dat.na$ind.id)

#Cuts for survival analysis
Tb.dat.cut<-Tb.dat.na
Tb.dat.cut$Age.num<-as.numeric(Tb.dat.cut$Age)
Tb.dat.cut<-survSplit(Surv(Age.num,surv)~.,data=Tb.dat.cut,cut=c(2,5,10,12),episode="interval")
Tb.dat.cut$agestart<-ifelse(Tb.dat.cut$Age==2,0,
                            ifelse(Tb.dat.cut$Age==5,2,
                                   ifelse(Tb.dat.cut$Age==10,5,
                                          ifelse(Tb.dat.cut$Age==12,10,NA))))
#write.csv(Tb.dat.cut,"Tb.dat.cut.csv")

Tb.dat.cut$BirdID<-as.character(Tb.dat.cut$BirdID)
Tb.dat.cut$Age<-as.numeric(Tb.dat.cut$Age)



### ***!!!!! NOTE - ages are 1=2, 2=5, 3=10, 4=12 !!!!****
Age2.tbsurv<-subset(Tb.dat.cut,Age==1)
Age2.tbsurv<-droplevels(Age2.tbsurv)
Age2.tbsurv$survnew<-ifelse(Age2.tbsurv$surv==1,0,1)

Age5.tbsurv<-subset(Tb.dat.cut,Age==2)
Age5.tbsurv<-droplevels(Age5.tbsurv)
Age5.tbsurv$survnew<-ifelse(Age5.tbsurv$surv==1,0,1)

Age10.tbsurv<-subset(Tb.dat.cut,Age==3)
Age10.tbsurv<-droplevels(Age10.tbsurv)
Age10.tbsurv$survnew<-ifelse(Age10.tbsurv$surv==1,0,1)

Age12.tbsurv<-subset(Tb.dat.cut,Age==4)
Age12.tbsurv<-droplevels(Age12.tbsurv)
Age12.tbsurv$survnew<-ifelse(Age12.tbsurv$surv==1,0,1)

library(lme4)

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
