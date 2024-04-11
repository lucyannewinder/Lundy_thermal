rm(list=ls())

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

Tb.dat1<-read.csv("Data/Raw/Tb.dat.csv")
str(Tb.dat1)
Tb.dat1$Age<-as.factor(Tb.dat1$Age)

# rearing brood is wrong in dataset for day 2 as it shows xfoster (which is at age 2)
age2.nat.cf<-subset(Tb.dat1,Age==2) %>%
  dplyr::select(BirdID,
                Brood.natal=Brood.code,
                Brood.rearing=CF.brood.no)

age2.nat.cf<-age2.nat.cf %>%
  mutate(Brood.rearing=coalesce(Brood.rearing,Brood.natal))

Tb.dat<-merge(Tb.dat1,age2.nat.cf,by="BirdID")

nrow(Tb.dat1)
nrow(Tb.dat)
## why is there a different number of rows here? 
## L - this is removing birds that don't have a day 2 record. Basically these are individuals who hatched later than the rest of the brood.

Tb.dat$BirdID<-as.character(Tb.dat$BirdID)

Tb.dat.na<-Tb.dat[!is.na(Tb.dat$max.temp), ]
nrow(Tb.dat.na)

# just get columns we use
Tb.dat.na<-Tb.dat.na[,c("BirdID","Brood.natal","Brood.rearing","Age","Year","Mass","max.temp","Air.temp","surv")]


## --------
## Read in and fix Pedigree
## --------

Tb.ped<-read.csv("Data/Raw/pedto19_new.csv")[,1:3]
colnames(Tb.ped)[1] <- "BirdID"

# add in missing parents
missing.birds<-unique(Tb.dat.na$BirdID[which(!Tb.dat.na$BirdID%in%Tb.ped[,1])])
Tb.ped<-orderPed(insertPed(Tb.ped,missing.birds))


## --------
## Write data to csvs
## --------

write.csv(Tb.ped,file="Data/Processed/ped_for_analysis.csv",row.names=FALSE)
write.csv(Tb.dat.na,file="Data/Processed/data_for_analysis.csv",row.names=FALSE)

