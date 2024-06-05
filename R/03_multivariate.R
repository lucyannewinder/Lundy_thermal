rm(list=ls())

# setwd("~/github/Lundy_thermal")
#temp by age - brood rearing for ages 5-12 only

library(MCMCglmm)
library(beepr)

## --------
## Read in data
## --------
Tb.ped <- read.csv("Data/Processed/ped_for_analysis.csv")
Tb.dat.na <- read.csv("Data/Processed/data_for_analysis.csv")


# need to make 2 level factor for age 2 and ages 5-12
# this is so we can fit rearing brood to only ages 5-12 and not age 2
Tb.dat.na$Age<-as.numeric(as.character(Tb.dat.na$Age))
Tb.dat.na$ages5to12<-ifelse(Tb.dat.na$Age>=5,"yes","no")

Tb.dat.na$Age<-as.factor(Tb.dat.na$Age)
# Tb.dat.na$Year<-as.factor(Tb.dat.na$Year)
Tb.dat.na$BirdID<-as.factor(Tb.dat.na$BirdID)
Tb.dat.na$ages5to12<-as.factor(Tb.dat.na$ages5to12)

# mean center covariates used in models 
Tb.dat.na$Air.tempC <- Tb.dat.na$Air.temp - mean(Tb.dat.na$Air.temp)
Tb.dat.na$Currentage_brood_sizeC <- Tb.dat.na$Currentage_brood_size - mean(Tb.dat.na$Currentage_brood_size)
#scale(Tb.dat.na$Air.temp, scale=FALSE)

head(Tb.dat.na)


Tb.ped.new<- prunePed(Tb.ped,unique(Tb.dat.na$BirdID), make.base = TRUE)
TB.Ainv<-inverseA(Tb.ped.new)$Ainv

#----


prior_Tb<-list(
  #R=list(V=1, nu=0.002), 
  R=list(R1=list(V=diag(4)*1e-6, nu=5)),
  G=list(
    Brood.rearing=list(V=diag(3), nu=1, alpha.mu=rep(0,3), alpha.V=diag(3)*100),
    Brood.natal=list(V=diag(4), nu=1, alpha.mu=rep(0,4), alpha.V=diag(4)*100)
  )
)
a=15

### ---------
# Model 1  - multivariate model of max temp across ontogeny
### ---------

mod_multi <- MCMCglmm(max.temp~Age-1+
                            Age:(Currentage_brood_sizeC+ Air.tempC+factor(Year)),
                          random=~us(at.level(ages5to12,"yes"):Age):Brood.rearing+
                            us(Age):Brood.natal,
                          rcov= ~us(Age):BirdID,
                          data=Tb.dat.na,
                          nitt=13000*a,
                          thin=10*a,
                          burnin=300*a,
                          family = "gaussian",
                          prior = prior_Tb, verbose=FALSE)

beep(1)
summary(mod_multi)
# plot(mod_multi)



### ---------
# Model 2  - accounting for body mass
### ---------

Tb.dat.na$MassC <- Tb.dat.na$Mass - mean(Tb.dat.na$Mass)

mod_multi_mass <- MCMCglmm(max.temp~Age-1+
                        Age:(MassC + Air.tempC + factor(Year)+ Currentage_brood_sizeC),
                      random=~us(at.level(ages5to12,"yes"):Age):Brood.rearing+
                        us(Age):Brood.natal,
                      rcov= ~us(Age):BirdID,
                      data=Tb.dat.na,
                      nitt=13000*a,
                      thin=10*a,
                      burnin=3000*a,
                      family = "gaussian",
                      prior = prior_Tb, verbose=FALSE)

beep(1)
summary(mod_multi_mass)
# plot(mod_multi_mass)


### ---------
### Save all model outputs
### ---------

save(
  mod_multi,mod_multi_mass,
  file = "Data/Output/multi_model.Rdata")

