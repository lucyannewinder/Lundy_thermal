#temp by age - brood rearing for ages 5-12 only

library(dplyr)
library(MCMCglmm)
library(MasterBayes)
library(tidyr)
library(survival)
library(lubridate)

#data set up
Tb.dat<-read.csv("Tb.dat.csv")
Tb.dat$Age<-as.factor(Tb.dat$Age)


# rearing brood is wrong in dataset for day 2 as it shows xfoster (which is at age 2)
age2.nat.cf<-subset(Tb.dat,Age==2) %>%
  dplyr::select(BirdID,
                Brood.natal=Brood.code,
                Brood.rearing=CF.brood.no)
age2.nat.cf<-age2.nat.cf %>%
  mutate(Brood.rearing=coalesce(Brood.rearing,Brood.natal))
Tb.dat<-merge(Tb.dat,age2.nat.cf,by="BirdID")

Tb.dat$BirdID<-as.character(Tb.dat$BirdID)
Tb.dat$ind.id<-Tb.dat$BirdID

#remove missing temp values
Tb.dat.na<-Tb.dat[!is.na(Tb.dat$max.temp), ]

# need to make 2 level factor for age 2 and ages 5-12
# this is so we can fit rearing brood to only ages 5-12 and not age 2
Tb.dat.na$Age<-as.numeric(as.character(Tb.dat.na$Age))
Tb.dat.na$ages5to12<-ifelse(Tb.dat.na$Age>=5,"yes","no")

Tb.dat.na$Age<-as.factor(Tb.dat.na$Age)
Tb.dat.na$BirdID<-as.factor(Tb.dat.na$BirdID)
Tb.dat.na$ind.id<-as.factor(Tb.dat.na$ind.id)
Tb.dat.na$ages5to12<-as.factor(Tb.dat.na$ages5to12)
#-----
#ped----

Tb.ped<-read.csv("pedto19_new.csv")
colnames(Tb.ped)[1] <- "BirdID"
Tb.ped<- Tb.ped %>%
  dplyr::select(BirdID,
                Dam,
                Sire)

# #check if any missing IDs for offspring, Dams + Sires
missing.birds<-unique(Tb.dat.na$BirdID[which(!Tb.dat.na$BirdID%in%Tb.ped[,1])])
# missing.fem<-unique(Tb.ped$Dam[which(!Tb.dat$femID%in%Tb.ped[,1])])
# missing.male<-unique(Tb.ped$Sire[which(!Tb.dat$maleID%in%Tb.ped[,1])])
Tb.ped<-insertPed(Tb.ped,missing.birds)

Tb.ped.new<- prunePed(orderPed(Tb.ped[,1:3]),unique(Tb.dat.na$BirdID), make.base = TRUE)
TB.Ainv<-inverseA(Tb.ped.new)$Ainv

#----


prior_tb_age_new<-list(
  #R=list(V=1, nu=0.002), 
  R=list(R1=list(V=diag(4)*1e-6, nu=5)),
  G=list(
    Brood.rearing=list(V=diag(3), nu=1, alpha.mu=rep(0,3), alpha.V=diag(3)*100),
    Brood.rearing=list(V=diag(4), nu=1, alpha.mu=rep(0,4), alpha.V=diag(4)*100)
  )
)
a=1
mod_tb_age_new<- MCMCglmm(max.temp~Age-1+
                        Air.temp+
                        factor(Year),
                      random=~us(at.level(ages5to12,"yes"):Age):Brood.rearing+
                        us(Age):Brood.natal,
                      rcov= ~us(Age):ind.id,
                      data=Tb.dat.na,
                      nitt=23000*a,
                      thin=10*a,
                      burnin=300*a,
                      family = "gaussian",
                      prior = prior_tb_age_new, verbose=FALSE)
summary(mod_tb_age_new)
