rm(list=ls())

# Proportion of variance models + plot code
# Each age modeled separately 

library(dplyr)
library(MCMCglmm)
library(tidyr)
library(survival)
library(lubridate)
library(ggplot2)
library(squidSim)


setwd("./github/Lundy_thermal")
## --------
## Read in data
## --------
load("Data/Output/uni_models.Rdata")

Tb.ped <- read.csv("Data/Processed/ped_for_analysis.csv")
Tb.dat.na <- read.csv("Data/Processed/data_for_analysis.csv")

Tb.dat.na$BirdID<-as.factor(Tb.dat.na$BirdID)
Tb.dat.na$doy<-as.factor(Tb.dat.na$doy)
Tb.dat.na$YearF <- as.factor(Tb.dat.na$Year)

## ---------
## DAY 2 model
## ---------

# subset data
Age2.tb<-subset(Tb.dat.na,Age==2)
Age2.tb<-droplevels(Age2.tb)

## remove bird not in pedigree
Age2.tb <- as.data.frame(subset(Age2.tb, BirdID != "10179"))

# mean center covariates used in models 
Age2.tb$Air.tempC <- scale(Age2.tb$Air.temp, scale=FALSE)
Age2.tb$Currentage_brood_sizeC <- scale(Age2.tb$Currentage_brood_size, scale=FALSE)

colnames(Age2.tb)<-gsub("\\.","_",colnames(Age2.tb))

# summary(mod_uni_day2)
# hist(mod_uni_day2$VCV[,1], breaks=100)

post_medians2 <-apply(mod_uni_day2$VCV,2,median)

preds2 <-model.matrix(~Air_tempC + YearF + Currentage_brood_sizeC,Age2.tb)[,-1]

nd_sims2 <- simulate_population(
	data_structure=Age2.tb[,c("BirdID", "Brood_natal", "Brood_rearing","doy","YearF")],
	parameters=list(
		intercept=as.numeric(apply(mod_uni_day2$Sol,2,median)[1]),
		# BirdID = list(vcov=0),
		doy= list(vcov=post_medians2["doy"]),
		Brood_natal= list(vcov=post_medians2["Brood.natal"]),
		residual = list(vcov=post_medians2["BirdID"]+post_medians2["units"])
		),
	# pedigree=list(BirdID=Tb.ped.day2),
	known_predictors=list(
		predictors=preds2,
		beta=apply(mod_uni_day2$Sol,2,median)[2:4]
	),
	n_pop=1000
)

nd_data2 <- get_population_data(nd_sims2,list=TRUE)

## get inverse A matrix for day 2 chicks
Tb.ped.day2<- prunePed((Tb.ped),unique(Age2.tb$BirdID), make.base = TRUE)
day2.Ainv<-inverseA(Tb.ped.day2)$Ainv

prior_day2 <- list(
  R=list(V=1, nu=0.002), 
  G=list(
    BirdID=list(V=1, nu=1, alpha.mu=0, alpha.V=1000),
    Brood.natal=list(V=1, nu=1, alpha.mu=0, alpha.V=1000),
    doy=list(V=1, nu=1, alpha.mu=0, alpha.V=1000)
  )
)

nd2 <-do.call(rbind,parallel::mclapply(1:length(nd_data2),function(x){
	a=5
	mod <- MCMCglmm(y~Air_tempC + YearF + Currentage_brood_sizeC,
   random=~BirdID + Brood_natal + doy,
   data=nd_data2[[x]],
   nitt=13000*a, thin=10*a, burnin=3000*a, 
   ginverse=list(BirdID=day2.Ainv), 
   prior=prior_day2, verbose=FALSE)
	cat(x," ")
	c(post_median=median(mod$VCV[,"BirdID"]),es=summary(mod)$Gcov["BirdID","eff.samp"])

},mc.cores=8))

save(
  mod_uni_day2,nd2,
  file = "Data/Output/uni_models_nd2.Rdata"
)



### -----
## DAY 5 model
### -----

# subset data
Age5.tb<-subset(Tb.dat.na,Age==5)
Age5.tb<-droplevels(Age5.tb)


# Remove birds not in pedigree
values_to_remove<-c("10179" ,"10821", "10827", "10829" ,"10911")
Age5.tb <- subset(Age5.tb, !BirdID %in% values_to_remove)
Age5.tb<-as.data.frame(Age5.tb)

# mean center covariates used in models 
Age5.tb$Air.tempC <- scale(Age5.tb$Air.temp, scale=FALSE)
Age5.tb$Currentage_brood_sizeC <- scale(Age5.tb$Currentage_brood_size, scale=FALSE)

colnames(Age5.tb)<-gsub("\\.","_",colnames(Age5.tb))

# simulate new datasets with no Va, adding Va to residual variance
post_medians5<-apply(mod_uni_day5$VCV,2,median)

preds5<-model.matrix(~Air_tempC + YearF + Currentage_brood_sizeC,Age5.tb)[,-1]

nd_sims5 <- simulate_population(
	data_structure=Age5.tb[,c("BirdID", "Brood_natal", "Brood_rearing","doy","YearF")],
	parameters=list(
		intercept=as.numeric(apply(mod_uni_day5$Sol,2,median)[1]),
		doy= list(vcov=post_medians5["doy"]),
		Brood_natal= list(vcov=post_medians5["Brood.natal"]),
		Brood_rearing= list(vcov=post_medians5["Brood.rearing"]),
		residual = list(vcov=post_medians5["BirdID"]+post_medians5["units"])
		),
	known_predictors=list(
		predictors=preds5,
		beta=apply(mod_uni_day5$Sol,2,median)[2:4]
	),
	n_pop=1000
)

nd_data5 <- get_population_data(nd_sims5,list=TRUE)


## get inverse A matrix for day 5 chicks
Tb.ped.day5<- prunePed(Tb.ped,unique(Age5.tb$BirdID), make.base = TRUE)
day5.Ainv<-inverseA(Tb.ped.day5)$Ainv

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

nd5 <-do.call(rbind,parallel::mclapply(1:length(nd_data5),function(x){
	a=5
	mod <- MCMCglmm(y~Air_tempC + YearF + Currentage_brood_sizeC,
   random=~BirdID + Brood_rearing + Brood_natal + doy,
   data=nd_data5[[x]],
   nitt=13000*a, thin=10*a, burnin=3000*a, 
   ginverse=list(BirdID=day5.Ainv), 
   prior=prior_day5plus, verbose=FALSE)
	cat(x," ")
	c(post_median=median(mod$VCV[,"BirdID"]),es=summary(mod)$Gcov["BirdID","eff.samp"])

},mc.cores=8))

save(
  mod_uni_day5,nd5,
  file = "Data/Output/uni_models_nd5.Rdata")



### -----
## DAY 10 model
### -----

# subset data
Age10.tb<-subset(Tb.dat.na,Age==10)
Age10.tb<-droplevels(Age10.tb)


# Remove birds not in pedigree
values_to_remove<-c("9989" , "9993" , "9998"  ,"10011" ,"10155" ,"10313", "10314" ,"10315" ,"10806" ,"10857" ,"10869", "10887", "10894")
Age10.tb <- subset(Age10.tb, !BirdID %in% values_to_remove)
Age10.tb<-as.data.frame(Age10.tb)

# mean center covariates used in models 
Age10.tb$Air.tempC <- scale(Age10.tb$Air.temp, scale=FALSE)
Age10.tb$Currentage_brood_sizeC <- scale(Age10.tb$Currentage_brood_size, scale=FALSE)

colnames(Age10.tb)<-gsub("\\.","_",colnames(Age10.tb))

# simulate new datasets with no Va, adding Va to residual variance
post_medians10<-apply(mod_uni_day10$VCV,2,median)

preds10<-model.matrix(~Air_tempC + YearF + Currentage_brood_sizeC,Age10.tb)[,-1]

nd_sims10 <- simulate_population(
	data_structure=Age10.tb[,c("BirdID", "Brood_natal", "Brood_rearing","doy","YearF")],
	parameters=list(
		intercept=as.numeric(apply(mod_uni_day10$Sol,2,median)[1]),
		doy= list(vcov=post_medians10["doy"]),
		Brood_natal= list(vcov=post_medians10["Brood.natal"]),
		Brood_rearing= list(vcov=post_medians10["Brood.rearing"]),
		residual = list(vcov=post_medians10["BirdID"]+post_medians10["units"])
		),
	known_predictors=list(
		predictors=preds10,
		beta=apply(mod_uni_day10$Sol,2,median)[2:4]
	),
	n_pop=1000
)

nd_data10 <- get_population_data(nd_sims10,list=TRUE)


## get inverse A matrix for day 5 chicks
Tb.ped.day10<- prunePed(Tb.ped,unique(Age10.tb$BirdID), make.base = TRUE)
day10.Ainv<-inverseA(Tb.ped.day10)$Ainv

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

nd10 <-do.call(rbind,parallel::mclapply(1:length(nd_data10),function(x){
	a=5
	mod <- MCMCglmm(y~Air_tempC + YearF + Currentage_brood_sizeC,
   random=~BirdID + Brood_rearing + Brood_natal + doy,
   data=nd_data10[[x]],
   nitt=13000*a, thin=10*a, burnin=3000*a, 
   ginverse=list(BirdID=day10.Ainv), 
   prior=prior_day5plus, verbose=FALSE)
	cat(x," ")
	c(post_median=median(mod$VCV[,"BirdID"]),es=summary(mod)$Gcov["BirdID","eff.samp"])

},mc.cores=8))

save(
  mod_uni_day10,nd10,
  file = "Data/Output/uni_models_nd10.Rdata")





### -----
## DAY 12 model
### -----

# subset data
Age12.tb<-subset(Tb.dat.na,Age==12)
Age12.tb<-droplevels(Age12.tb)


# Remove birds not in pedigree
values_to_remove<-c("10011", "10155", "10313" ,"10314", "10315", "10806", "10857", "10869", "10887", "10894")
Age12.tb <- subset(Age12.tb, !BirdID %in% values_to_remove)
Age12.tb<-as.data.frame(Age12.tb)

# mean center covariates used in models 
Age12.tb$Air.tempC <- scale(Age12.tb$Air.temp, scale=FALSE)
Age12.tb$Currentage_brood_sizeC <- scale(Age12.tb$Currentage_brood_size, scale=FALSE)

colnames(Age12.tb)<-gsub("\\.","_",colnames(Age12.tb))

# simulate new datasets with no Va, adding Va to residual variance
post_medians12<-apply(mod_uni_day12$VCV,2,median)

preds12<-model.matrix(~Air_tempC + YearF + Currentage_brood_sizeC,Age12.tb)[,-1]

nd_sims12 <- simulate_population(
	data_structure=Age12.tb[,c("BirdID", "Brood_natal", "Brood_rearing","doy","YearF")],
	parameters=list(
		intercept=as.numeric(apply(mod_uni_day12$Sol,2,median)[1]),
		doy= list(vcov=post_medians12["doy"]),
		Brood_natal= list(vcov=post_medians12["Brood.natal"]),
		Brood_rearing= list(vcov=post_medians12["Brood.rearing"]),
		residual = list(vcov=post_medians12["BirdID"]+post_medians12["units"])
		),
	known_predictors=list(
		predictors=preds12,
		beta=apply(mod_uni_day12$Sol,2,median)[2:4]
	),
	n_pop=1000
)

nd_data12 <- get_population_data(nd_sims12,list=TRUE)


## get inverse A matrix for day 5 chicks
Tb.ped.day12<- prunePed(Tb.ped,unique(Age12.tb$BirdID), make.base = TRUE)
day12.Ainv<-inverseA(Tb.ped.day12)$Ainv

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

nd12 <-do.call(rbind,parallel::mclapply(1:length(nd_data12),function(x){
	a=5
	mod <- MCMCglmm(y~Air_tempC + YearF + Currentage_brood_sizeC,
   random=~BirdID + Brood_rearing + Brood_natal + doy,
   data=nd_data12[[x]],
   nitt=13000*a, thin=10*a, burnin=3000*a, 
   ginverse=list(BirdID=day12.Ainv), 
   prior=prior_day5plus, verbose=FALSE)
	cat(x," ")
	c(post_median=median(mod$VCV[,"BirdID"]),es=summary(mod)$Gcov["BirdID","eff.samp"])

},mc.cores=8))

save(
  mod_uni_day12,nd12,
  file = "Data/Output/uni_models_nd12.Rdata")







%%%%%%%%%%%%%%%%%%%%%%


load("Data/Output/uni_models_nd2.Rdata")
post_medians2 <-apply(mod_uni_day2$VCV,2,median)


hist(nd2[,"post_median"], breaks=100)
abline(v=post_medians2["BirdID"],col="red")
MDE2<-quantile(nd2[,"post_median"],0.95)
MDE2/sum(post_medians2)
plot(nd2)

p_value2 <- mean(post_medians2["BirdID"]<nd2[,"post_median"])
p_value2

post_medians2[1]/sum(post_medians2)




load("Data/Output/uni_models_nd5.Rdata")
post_medians5 <-apply(mod_uni_day5$VCV,2,median)


hist(nd5[,"post_median"], breaks=100)
abline(v=post_medians5["BirdID"],col="red")
MDE5<-quantile(nd5[,"post_median"],0.95)
MDE2/sum(post_medians5)
plot(nd5)

p_value5 <- mean(post_medians5["BirdID"]<nd5[,"post_median"])
p_value5

post_medians5[1]/sum(post_medians5)