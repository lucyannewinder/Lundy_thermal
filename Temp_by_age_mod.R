# Temperature by age correlations
# Also correlation with body mass interation (to remove growth effect)

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






prior_tbmass_age<-list(
  R=list(R1=list(V=diag(4)*1e-6, nu=5)),
  G=list(
    Brood.natal=list(V=diag(4), nu=1, alpha.mu=rep(0,4), alpha.V=diag(4)*100),
    Brood.rearing=list(V=diag(4), nu=1, alpha.mu=rep(0,4), alpha.V=diag(4)*100)#,
    # ind.id=list(V=diag(4), nu=1, alpha.mu=rep(0,4), alpha.V=diag(4)*100)
  )
)
a=5
mod_tbmass_age<- MCMCglmm(max.temp~-1+Age*Mass+
                            Air.temp+
                            factor(Year),
                          random=~us(Age):Brood.natal+
                            us(Age):Brood.rearing,
                          rcov=~us(Age):ind.id,
                          data=Tb.dat.na,
                          nitt=23000*a,
                          thin=10*a,
                          burnin=300*a,
                          family = "gaussian",
                          prior = prior_tbmass_age, verbose=FALSE)

summary(mod_tbmass_age)

post_vcv_tbage <- mod_tbmass_age$VCV

#posteriors
post_tbage_nat<-post_vcv_tbage[,grep("Brood.natal",colnames(post_vcv_tbage))]
post_tbage_rear<-post_vcv_tbage[,grep("Brood.rearing",colnames(post_vcv_tbage))]
post_tbage_res<-post_vcv_tbage[,grep("ind.id",colnames(post_vcv_tbage))]

## proportion of variance explains by each level (very roughly)
# get variances from these objects
var_tbage_nat <- t(apply(post_tbage_nat/rowSums(post_tbage_nat),2,function(x) mean_CI(as.mcmc(x),pMCMC=TRUE)))
var_tbage_rear <- t(apply(post_tbage_rear/rowSums(post_tbage_rear),2,function(x) mean_CI(as.mcmc(x),pMCMC=TRUE)))
var_tbage_res <- t(apply(post_tbage_res/rowSums(post_tbage_res),2,function(x) mean_CI(as.mcmc(x),pMCMC=TRUE)))


post.tb.nat<-posterior.cor(post_tbage_nat)
post.tb.rear<-posterior.cor(post_tbage_rear)
post.tb.res<-posterior.cor(post_tbage_res)

# get correlations from these objects
tb.vars.nat<-t(apply(post.tb.nat, 2, function(x) mean_CI(as.mcmc(x),pMCMC=TRUE)))
tb.vars.rear<-t(apply(post.tb.rear, 2, function(x) mean_CI(as.mcmc(x),pMCMC=TRUE)))
tb.vars.res<-t(apply(post.tb.res, 2, function(x) mean_CI(as.mcmc(x),pMCMC=TRUE)))