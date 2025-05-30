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

# change rearing brood for age 2
# rearing brood is wrong in dataset for day 2 as it shows xfoster (which is at age 2)
age2.nat.cf<-subset(Tb.dat1,Age==2) %>%
  dplyr::select(BirdID,
                Brood.natal=Brood.code,
                Brood.rearing=CF.brood.no,
                Date_chr=as.character("Date"))

age2.nat.cf<-age2.nat.cf %>%
  mutate(Brood.rearing=coalesce(Brood.rearing,Brood.natal),
         Date_parsed = as.Date(Date_chr, format = "%d/%m/%Y"),
         doy=as.numeric(format(Date_parsed,"%j")))

Tb.dat<-merge(Tb.dat1,age2.nat.cf,by="BirdID",all = TRUE)

# make column to filter out birds with missing brood information later. These birds are those who hatch significantly later than their siblings.
Tb.dat$keep_row<-!is.na(Tb.dat$Brood.natal)
# assign brood information for above birds - this will be used to calculate brood number.
Tb.dat$Brood.natal<-ifelse(is.na(Tb.dat$Brood.natal),Tb.dat$Brood.code,Tb.dat$Brood.natal)
Tb.dat$Brood.rearing<-ifelse(is.na(Tb.dat$Brood.rearing),Tb.dat$Brood.code,Tb.dat$Brood.rearing)

# make brood size variable
# age 2 needs to be calculated from nqtal brood whereas other ages from rearing brood
age2_Tb.dat<-subset(Tb.dat,Age==2)
age2_Tb.dat<- age2_Tb.dat %>%
  group_by(Age, Brood.natal) %>%
  mutate(Currentage_brood_size = n())

age5up_Tb.dat<-subset(Tb.dat,Age==5 | Age==10 | Age==12)
age5up_Tb.dat<- age5up_Tb.dat %>%
  group_by(Age, Brood.rearing) %>%
  mutate(Currentage_brood_size = n())

Tb.dat<-rbind(age2_Tb.dat,age5up_Tb.dat)

# remove birds using filter for missing brood information (above)
Tb.dat <- subset(Tb.dat,keep_row=="TRUE")

Tb.dat$BirdID<-as.character(Tb.dat$BirdID)

Tb.dat.na<-Tb.dat[!is.na(Tb.dat$max.temp), ]
nrow(Tb.dat.na)

# just get columns we use
Tb.dat.na<-Tb.dat.na[,c("BirdID","Brood.natal","Brood.rearing","Age","Year","Mass","max.temp","Air.temp","surv","Currentage_brood_size","doy")]


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
