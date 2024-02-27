# Proportion of variance models + plot code
# Each age modeled separately 

library(dplyr)
library(MCMCglmm)
library(MasterBayes)
library(tidyr)
library(survival)
library(lubridate)
library(ggplot2)


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

tb.2018<-subset(Age2.tb,Year=="2018")
mean(tb.2018$Air.temp)
tb.2019<-subset(Age2.tb,Year=="2019")
mean(tb.2019$Age12.tb)
#ped----

Tb.ped<-read.csv("pedto19_new.csv")
colnames(Tb.ped)[1] <- "BirdID"
Tb.ped<- Tb.ped %>%
  dplyr::select(BirdID,
                Dam,
                Sire)

missing.birds<-unique(Tb.dat.na$BirdID[which(!Tb.dat.na$BirdID%in%Tb.ped[,1])])
Tb.ped<-insertPed(Tb.ped,missing.birds)

Tb.ped.new<- prunePed(orderPed(Tb.ped[,1:3]),unique(Tb.dat.na$BirdID), make.base = TRUE)
TB.Ainv<-inverseA(Tb.ped.new)$Ainv

#----
Tb.dat.na$BirdID<-as.factor(Tb.dat.na$BirdID)
Tb.dat.na$ind.id<-as.factor(Tb.dat.na$ind.id)

Age2.tb<-subset(Tb.dat.na,Age==2)
Age2.tb<-droplevels(Age2.tb)
Age5.tb<-subset(Tb.dat.na,Age==5)
Age5.tb<-droplevels(Age5.tb)
Age10.tb<-subset(Tb.dat.na,Age==10)
Age10.tb<-droplevels(Age10.tb)
Age12.tb<-subset(Tb.dat.na,Age==12)
Age12.tb<-droplevels(Age12.tb)

length(Age12.tb$max.temp)

#add brood size

Age2.tb <- Age2.tb %>%
  group_by(Brood.natal) %>%
  mutate(Natal_brood_size = n())
Age5.tb <- Age5.tb %>%
  group_by(Brood.natal) %>%
  mutate(Natal_brood_size = n())
Age5.tb <- Age5.tb %>%
  group_by(Brood.rearing) %>%
  mutate(Rearing_brood_size = n())
Age10.tb <- Age10.tb %>%
  group_by(Brood.natal) %>%
  mutate(Natal_brood_size = n())
Age10.tb <- Age10.tb %>%
  group_by(Brood.rearing) %>%
  mutate(Rearing_brood_size = n())
Age12.tb <- Age12.tb %>%
  group_by(Brood.natal) %>%
  mutate(Natal_brood_size = n())
Age12.tb <- Age12.tb %>%
  group_by(Brood.rearing) %>%
  mutate(Rearing_brood_size = n())

prior_Tbher<-list(
  R=list(V=1, nu=0.002), 
  G=list(
    BirdID=list(V=1, nu=1, alpha.mu=0, alpha.V=1000),
    Brood.natal=list(V=1, nu=1, alpha.mu=0, alpha.V=1000),
    Brood.rearing=list(V=1, nu=1, alpha.mu=0, alpha.V=1000)
  )
)
a=20


unique_to_vector2 <- setdiff(Age2.tb$BirdID, Tb.ped.new$BirdID)
Age2.tb <- subset(Age2.tb, BirdID != "10179")


Age2.tb$femID<-as.factor(Age2.tb$femID)

modM1_Tbher_2 <- MCMCglmm(max.temp~Air.temp+
                            as.factor(Year)+
                            Natal_brood_size,
                     random=~BirdID+
                       Brood.natal,
                     data=Age2.tb,
                     nitt=23000*a, thin=10*a, burnin=3000*a, 
                     ginverse=list(BirdID=TB.Ainv), prior=prior_Tbher, verbose=FALSE)
summary(modM1_Tbher_2)
plot(modM1_Tbher_2)

posther.2<-modM1_Tbher_2$VCV[,"BirdID"]/(modM1_Tbher_2$VCV[,"BirdID"]+modM1_Tbher_2$VCV[,"Brood.natal"]+modM1_Tbher_2$VCV[,"units"])
HPDinterval(posther.2,0.95)
posterior.mode(posther.2)

# posterior.mode(modM1_Tbher_2$VCV[,"units"])

plot(posther.2)

Tb.ped<-read.csv("pedto19_new.csv")
colnames(Tb.ped)[1] <- "BirdID"
Tb.ped<- Tb.ped %>%
  dplyr::select(BirdID,
                Dam,
                Sire)
missing.birds<-unique(Age5.tb$BirdID[which(!Age5.tb$BirdID%in%Tb.ped[,1])])
Tb.ped<-insertPed(Tb.ped,missing.birds)
Tb.ped.new<- prunePed(orderPed(Tb.ped[,1:3]),unique(Age5.tb$BirdID), make.base = TRUE)
TB.Ainv<-inverseA(Tb.ped.new)$Ainv

unique_to_vector2 <- setdiff(Age5.tb$BirdID, Tb.ped.new$BirdID)
values_to_remove<-c("10179" ,"10821", "10827", "10829" ,"10911")
Age5.tb <- subset(Age5.tb, !BirdID %in% values_to_remove)
Age5.tb<-as.data.frame(Age5.tb)

modM1_Tbher_5 <- MCMCglmm(max.temp~Air.temp+
                            as.factor(Year)+
                            Rearing_brood_size,
                          random=~BirdID+ 
                            Brood.rearing+
                            Brood.natal,
                          data=Age5.tb,
                          nitt=23000*a, thin=10*a, burnin=3000*a, 
                          ginverse=list(BirdID=TB.Ainv), prior=prior_Tbher, verbose=FALSE)
summary(modM1_Tbher_5)
plot(modM1_Tbher_5)

posther.5<-modM1_Tbher_5$VCV[,"BirdID"]/
  (modM1_Tbher_5$VCV[,"BirdID"]+modM1_Tbher_5$VCV[,"Brood.natal"]+modM1_Tbher_5$VCV[,"Brood.rearing"]+modM1_Tbher_5$VCV[,"units"])
HPDinterval(posther.5,0.95)
posterior.mode(posther.5)

unique_to_vector2 <- setdiff(Age10.tb$BirdID, Tb.ped.new$BirdID)
values_to_remove<-c("9989" , "9993" , "9998"  ,"10011" ,"10155" ,"10313", "10314" ,"10315" ,"10806" ,"10857" ,"10869", "10887", "10894")
Age10.tb <- subset(Age10.tb, !BirdID %in% values_to_remove)
Age10.tb<-as.data.frame(Age10.tb)

modM1_Tbher_10 <- MCMCglmm(max.temp~Air.temp+
                             as.factor(Year)+
                             Rearing_brood_size,
                          random=~BirdID+ 
                            Brood.natal+
                            Brood.rearing,
                          data=Age10.tb,
                          nitt=23000*a, thin=10*a, burnin=3000*a, 
                          ginverse=list(BirdID=TB.Ainv), prior=prior_Tbher, verbose=FALSE)
summary(modM1_Tbher_10)

posther.10<-modM1_Tbher_10$VCV[,"units"]/
  (modM1_Tbher_10$VCV[,"BirdID"]+modM1_Tbher_10$VCV[,"Brood.natal"]+modM1_Tbher_10$VCV[,"Brood.rearing"]+modM1_Tbher_10$VCV[,"units"])
HPDinterval(posther.10,0.95)
posterior.mode(posther.10)



unique_to_vector2 <- setdiff(Age12.tb$BirdID, Tb.ped.new$BirdID)
values_to_remove<-c("10011", "10155", "10313" ,"10314", "10315", "10806", "10857", "10869", "10887", "10894")
Age12.tb <- subset(Age12.tb, !BirdID %in% values_to_remove)
Age12.tb<-as.data.frame(Age12.tb)

modM1_Tbher_12 <- MCMCglmm(max.temp~Air.temp+
                             as.factor(Year)+
                             Rearing_brood_size,
                          random=~BirdID+
                            Brood.natal+
                            Brood.rearing,
                          data=Age12.tb,
                          nitt=23000*a, thin=10*a, burnin=3000*a, 
                          ginverse=list(BirdID=TB.Ainv), prior=prior_Tbher, verbose=FALSE)
summary(modM1_Tbher_12)

posther.12<-modM1_Tbher_12$VCV[,"units"]/
  (modM1_Tbher_12$VCV[,"BirdID"]+modM1_Tbher_12$VCV[,"Brood.natal"]+modM1_Tbher_12$VCV[,"Brood.rearing"]+modM1_Tbher_12$VCV[,"units"])
HPDinterval(posther.12,0.95)
posterior.mode(posther.12)

#---- temp by age plot ----
library(ggplot2)
variance.tempage<-read.csv("variance.tempbyage.csv")
variance.tempage.broodsize<-read.csv("variance.tempbyage.broodsize.csv")
str(variance.tempage)
plot_tempage<-ggplot(data = variance.tempage,aes(x=as.factor(Age),y=Proportion.of.Variance,group=Trait,fill=Trait))+
  geom_line(size=0.8,aes(linetype=Trait))+
  geom_point(size=3,shape=21,
             position=position_dodge(width=0.04))+
  geom_errorbar(aes(ymin=l.ci,ymax=u.ci,colour=Trait),width=0.05,
                position=position_dodge(width=0.04))+
  theme_light()+
  theme(text=element_text(size=20),
        legend.title=element_blank())+
  xlab("Age")+
  ylab("Proportion of Variance")+
  ggtitle("Proportion of variance")
plot_tempage

#--------
survlabs<-c("Yes","No")
surv.temp.2<-ggplot(data=Age12.tb,aes(x=as.factor(surv),y=max.temp))+
  geom_boxplot()+
  xlab("Survival")+
  ylab("Maximum head temperature (ºC)")+
  scale_x_discrete(labels=survlabs)+
  theme_classic(base_size = 12)
surv.temp.2+
  ggtitle("Age 2")+theme_classic(plot.title = element_text(face="bold"))





