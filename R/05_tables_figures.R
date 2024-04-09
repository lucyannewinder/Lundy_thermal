
post_summary <- function(x) {
  out <- c(mean(x),quantile(x, c(0.5,0.025,0.975)))
  names(out) <- c("mean","median","lCI","uCI")
  out
}


## Univariate
load("Data/Output/uni_models.Rdata")



uni_day2 <- t(apply(mod_uni_day2$VCV/rowSums(mod_uni_day2$VCV),2,post_summary))

uni_day5 <- t(apply(mod_uni_day5$VCV/rowSums(mod_uni_day5$VCV),2,post_summary))

uni_day10 <- t(apply(mod_uni_day10$VCV/rowSums(mod_uni_day10$VCV),2,post_summary))

uni_day12 <- t(apply(mod_uni_day12$VCV/rowSums(mod_uni_day12$VCV),2,post_summary))



var_tbage_res <- t(apply(post_tbage_res/rowSums(post_tbage_res),2,function(x) mean_CI(as.mcmc(x),pMCMC=TRUE)))









cor(modM1_Tbher_2$VCV)
cor(modM1_Tbher_5$VCV)
cor(modM1_Tbher_10$VCV)
cor(modM1_Tbher_12$VCV)




## ------
#---- temp by age plot ----
## ------

library(ggplot2)
variance.tempage<-read.csv("Data/Output/variance.tempbyage.csv")
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


load("Data/Output/multi_model.Rdata")


post_vcv_tbage <- mod_tb_age_new$VCV

#posteriors
post_tbage_nat<-post_vcv_tbage[,grep("Brood.natal",colnames(post_vcv_tbage))]
post_tbage_rear<-post_vcv_tbage[,grep("Brood.rearing",colnames(post_vcv_tbage))]
post_tbage_res<-post_vcv_tbage[,grep("BirdID",colnames(post_vcv_tbage))]

## proportion of variance explains by each level (very roughly)
# get variances from these objects
var_tbage_nat <- t(apply(post_tbage_nat,2,function(x) mean_CI(as.mcmc(x),pMCMC=TRUE)))
var_tbage_rear <- t(apply(post_tbage_rear,2,function(x) mean_CI(as.mcmc(x),pMCMC=TRUE)))
var_tbage_res <- t(apply(post_tbage_res,2,function(x) mean_CI(as.mcmc(x),pMCMC=TRUE)))

post.tb.nat<-posterior.cor(post_tbage_nat)
post.tb.rear<-posterior.cor(post_tbage_rear)
post.tb.res<-posterior.cor(post_tbage_res)

# get correlations from these objects
tb.vars.nat<-t(apply(post.tb.nat, 2, function(x) mean_CI(as.mcmc(x),pMCMC=TRUE)))
tb.vars.rear<-t(apply(post.tb.rear, 2, function(x) mean_CI(as.mcmc(x),pMCMC=TRUE)))
tb.vars.res<-t(apply(post.tb.res, 2, function(x) mean_CI(as.mcmc(x),pMCMC=TRUE)))












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

