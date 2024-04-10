rm(list=ls())

library(stringr)
library(MCMCglmm)

# setwd("~/github/Lundy_thermal")

##function to get 'pvalues' from MCMC chain
pMCMC <- function(x) 2*pmax(0.5/length(x), pmin(sum(x > 0)/length(x), 1 - sum(x > 0)/length(x)))

## function to summarise a posterior distribution
post_summary <- function(x,pMCMC=FALSE) {
  y <- c(mean(x),quantile(x, c(0.5,0.025,0.975)))
  names(y) <- c("mean","median","lCI","uCI")
  if(pMCMC) y <- c(y,pMCMC=pMCMC(x))
  y
}

## ------
## Univariate Models 
## ------

load("Data/Output/uni_models.Rdata")

uni_day2 <- cbind(t(apply(mod_uni_day2$VCV/rowSums(mod_uni_day2$VCV),2,post_summary)),Age=2)

uni_day5 <- cbind(t(apply(mod_uni_day5$VCV/rowSums(mod_uni_day5$VCV),2,post_summary)),Age=5)

uni_day10 <- cbind(t(apply(mod_uni_day10$VCV/rowSums(mod_uni_day10$VCV),2,post_summary)),Age=10)

uni_day12 <- cbind(t(apply(mod_uni_day12$VCV/rowSums(mod_uni_day12$VCV),2,post_summary)),Age=12)

prop_uni_all <- rbind(uni_day2,uni_day5,uni_day10,uni_day12)
 
trait_names <- ifelse(rownames(prop_uni_all)=="BirdID","genetic",
    ifelse(rownames(prop_uni_all)=="Brood.rearing","Rearing brood",
      ifelse(rownames(prop_uni_all)=="Brood.natal","Natal brood",
        ifelse(rownames(prop_uni_all)=="units","Residual",
    NA))))
variance.tempage <- data.frame(prop_uni_all,Trait=trait_names)

## ------
#---- temp by age plot ----
## ------

library(ggplot2)
# variance.tempage<-read.csv("Data/Output/variance.tempbyage.csv")
# str(variance.tempage)
plot_tempage<-ggplot(data = variance.tempage,aes(x=as.factor(Age),y=median,group=Trait,fill=Trait))+
  geom_line(size=0.8,aes(linetype=Trait))+
  geom_point(size=3,shape=21,
             position=position_dodge(width=0.04))+
  geom_errorbar(aes(ymin=lCI,ymax=uCI,colour=Trait),width=0.05,
                position=position_dodge(width=0.04))+
  theme_light()+
  theme(text=element_text(size=20),
        legend.title=element_blank())+
  xlab("Age")+
  ylab("Proportion of Variance")+
  ggtitle("Proportion of variance")
plot_tempage



## ------
## Multivariate Models 
## ------

load("Data/Output/multi_model.Rdata")


post_multi <- mod_multi$VCV
colnames(post_multi) <- gsub('at.level\\(ages5to12, \\"yes\\"\\):',"",colnames(post_multi))

## posteriors grouped by age
post_multi_age2 <- post_multi[,which(str_count(colnames(post_multi),"Age2")==2)]
post_multi_age5 <- post_multi[,which(str_count(colnames(post_multi),"Age5")==2)]
post_multi_age10 <- post_multi[,which(str_count(colnames(post_multi),"Age10")==2)]
post_multi_age12 <- post_multi[,which(str_count(colnames(post_multi),"Age12")==2)]

#posteriors grouped by component
post_multi_nat<-post_multi[,grep("Brood.natal",colnames(post_multi))]
post_multi_rear<-post_multi[,grep("Brood.rearing",colnames(post_multi))]
post_multi_res<-post_multi[,grep("BirdID",colnames(post_multi))]

### function to make covariance-correlation matrices
cov_cor_matrix <- function(post_cov, round=3){
  
  matrix_size <-sqrt(ncol(post_cov))
  
  # get sunmmaries of (co)variances from posterior
  cov_sum <- t(apply(post_cov,2,function(x) post_summary(as.mcmc(x),pMCMC=TRUE)))
  covs_out <- matrix(paste0(round(cov_sum[,"median"],round)," (",round(cov_sum[,"lCI"],round),", ",round(cov_sum[,"uCI"],round),")"),ncol=matrix_size)
  
  # generate correlations and get summaries
  post_cor <- posterior.cor(post_cov)
  cor_sum <- t(apply(post_cor, 2, function(x) post_summary(as.mcmc(x),pMCMC=TRUE)))
  cors_out <- matrix(paste0(round(cor_sum[,"median"],round)," (",round(cor_sum[,"lCI"],round),", ",round(cor_sum[,"uCI"],round),")"),ncol=matrix_size)

  # make a covariance (lower Tri)/ correlation (upper tri) matrix, with variances on the diagonal
  out <- matrix(NA,matrix_size,matrix_size)
  out[upper.tri(out)] <- cors_out[upper.tri(cors_out)]
  out[lower.tri(out,diag=TRUE)] <- covs_out[lower.tri(covs_out,diag=TRUE)]
  out

}
cov_cor_matrix(post_multi_nat)
cov_cor_matrix(post_multi_rear)
cov_cor_matrix(post_multi_res)










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

