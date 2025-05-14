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

setwd("./github/Lundy_thermal")
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

## remove bird not in pedigree
Age2.tb <- as.data.frame(subset(Age2.tb, BirdID != "10179"))

# mean center covariates used in models 
Age2.tb$Air.tempC <- scale(Age2.tb$Air.temp, scale=FALSE)
Age2.tb$Currentage_brood_sizeC <- scale(Age2.tb$Currentage_brood_size, scale=FALSE)

Age2.tb$YearF <- as.factor(Age2.tb$Year)

colnames(Age2.tb)<-gsub("\\.","_",colnames(Age2.tb))

## get inverse A matrix for day 2 chicks
Tb.ped.day2<- prunePed(Tb.ped,unique(Age2.tb$BirdID), make.base = TRUE)
day2.Ainv<-inverseA(Tb.ped.day2)$Ainv

load("Data/Output/uni_models.Rdata")


summary(mod_uni_day2)
hist(mod_uni_day2$VCV[,1])

post_medians<-apply(mod_uni_day2$VCV,2,median)

library(squidSim)

preds<-model.matrix(~Air_tempC + YearF + Currentage_brood_sizeC,Age2.tb)[,-1]


nd_sims <- simulate_population(
	data_structure=Age2.tb[,c("BirdID", "Brood_natal", "Brood_rearing","YearF")],
	parameters=list(
		intercept=as.numeric(apply(mod_uni_day2$Sol,2,median)[1]),
		# BirdID = list(vcov=0),
		Brood_natal= list(vcov=post_medians["Brood.natal"]),
		residual = list(vcov=post_medians["BirdID"]+post_medians["units"])
		),
	# pedigree=list(BirdID=Tb.ped.day2),
	known_predictors=list(
		predictors=preds,
		beta=apply(mod_uni_day2$Sol,2,median)[2:4]
	),
	n_pop=1000
)



nd_data <- get_population_data(nd_sims,list=TRUE)


prior_day2 <- list(
  R=list(V=1, nu=0.002), 
  G=list(
    BirdID=list(V=1, nu=1, alpha.mu=0, alpha.V=1000),
    Brood.natal=list(V=1, nu=1, alpha.mu=0, alpha.V=1000)
  )
)
nd<-do.call(rbind,parallel::mclapply(1:length(nd_data),function(x){
	a=5
	mod <- MCMCglmm(y~Air_tempC + YearF + Currentage_brood_sizeC,
   random=~BirdID + Brood_natal,
   data=nd_data[[x]],
   nitt=13000*a, thin=10*a, burnin=3000*a, 
   ginverse=list(BirdID=day2.Ainv), 
   prior=prior_day2, verbose=FALSE)
	print(x)
	c(post_median=median(mod$VCV[,"BirdID"]),es=summary(mod)$Gcov["BirdID","eff.samp"])

},mc.cores=8))

save(
  mod_uni_day2,nd,
  file = "Data/Output/uni_models_nd.Rdata")



hist(nd[,"post_median"], breaks=100)
abline(v=post_medians["BirdID"],col="red")
MDE<-quantile(nd[,"post_median"],0.95)
MDE/sum(post_medians)
plot(nd)

p_value <- mean(post_medians["BirdID"]<nd[,"post_median"])


