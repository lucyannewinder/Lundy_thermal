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
 
trait_names <- ifelse(rownames(prop_uni_all)=="BirdID","Genetic",
    ifelse(rownames(prop_uni_all)=="Brood.rearing","Rearing brood",
      ifelse(rownames(prop_uni_all)=="Brood.natal","Natal brood",
        ifelse(rownames(prop_uni_all)=="units","Residual",
               ifelse(rownames(prop_uni_all)=="doy","Date",
    NA)))))
variance.tempage <- data.frame(prop_uni_all,Trait=trait_names)

## ------
#---- temp by age plot ----
## ------

library(ggplot2)
# variance.tempage<-read.csv("Data/Output/variance.tempbyage.csv")
# str(variance.tempage)

variance.tempage$Trait<-ifelse(variance.tempage$Age==2 & variance.tempage$Trait == "Natal brood", "Brood", variance.tempage$Trait)

# version 1 - combined plot
plot_tempage<-ggplot(data = variance.tempage,aes(x=as.factor(Age),y=median,group=Trait,fill=Trait))+
  geom_line(aes(linetype=Trait))+
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

#version 2 - faceted
plot_tempage_facet<-ggplot(data = variance.tempage,aes(x=as.factor(Age),y=median,group=Trait,fill=Trait))+
  geom_point(size=3,shape=21,
             position=position_dodge(width=4))+
  geom_errorbar(aes(ymin=lCI,ymax=uCI,colour=Trait),width=0.05,
                position=position_dodge(width=4))+
  scale_fill_manual(values=c("goldenrod","lightpink2","palegreen3","steelblue1","royalblue3","indianred"))+
  scale_colour_manual(values=c("goldenrod","lightpink2","palegreen3","steelblue1","royalblue3","indianred"))+
  theme_linedraw()+
  theme(text=element_text(size=20),
        legend.title=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_blank())+
  ylab("Proportion of Variance")+
  xlab("")+
  facet_wrap("Age",nrow=1)+
  geom_hline(yintercept = c(0,0.25,0.5,0.75), linetype=2,alpha=0.2)
plot_tempage_facet

## ------
## Multivariate Models 
## ------

load("Data/Output/multi_model.Rdata")

#model options mod_multi or mod_multi_mass - former without body mass effect and latter with
# post_multi <- mod_multi$VCV
post_multi <- mod_multi_mass$VCV

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





## ------
## Survival Models 
## ------

load("Data/Output/surv_models.Rdata")

# Function to create survival plots
make_survival_plot <- function(data, model, age_label = "", plot_label = "") {
  years <- levels(data$Year) #colour + predictions for each year separately in plot
  
  # Generate model predictions for predicted survival curve
  newdata_list <- lapply(years, function(y) {
    df <- data.frame(
      max.tempC = seq(min(data$max.tempC, na.rm = TRUE),
                      max(data$max.tempC, na.rm = TRUE),
                      length.out = 100),
      MassC = mean(data$MassC, na.rm = TRUE),
      Year = factor(y, levels = levels(data$Year)),
      Brood.natal = data$Brood.natal[1],
      doy = data$doy[1]
    )
    
    # Create matrix of predictions from model
    X <- model.matrix(~ max.tempC + MassC + as.factor(Year), data = df)
    betas <- fixef(model)
    vcov_mat <- vcov(model)
    df$logit <- X %*% betas
    df$se <- sqrt(diag(X %*% vcov_mat %*% t(X)))
    df$logit.lower <- df$logit - 1.96 * df$se
    df$logit.upper <- df$logit + 1.96 * df$se
    
    # Convert from logit to predicted probabilities
    invlogit <- function(x) exp(x) / (1 + exp(x))
    df$fit <- invlogit(df$logit)
    df$lower <- invlogit(df$logit.lower)
    df$upper <- invlogit(df$logit.upper)
    return(df)
  })
  
  newdata_all <- do.call(rbind, newdata_list)
  
  # plot
  
  ggplot() +
    geom_point(data = data,
               aes(x = max.tempC, y = survnew, colour = Year),
               position = position_jitter(width = 0, height = 0.08),
               alpha = 0.4, size = 2) +
    geom_line(data = newdata_all,
              aes(x = max.tempC, y = fit, colour = Year),
              size = 0.8) +
    geom_ribbon(data = newdata_all,
                aes(x = max.tempC, ymin = lower, ymax = upper, fill = Year),
                alpha = 0.2) +
    ylab("Survival") +
    xlab("Mean centred body temperature (ºC)") +
    ggtitle(paste0(plot_label, " Age ", age_label)) +
    theme_light(base_size = 18) +
    scale_y_continuous(breaks = c(0, 1), labels = c("0", "1")) +
    scale_color_brewer(palette = "Set1") +
    scale_fill_brewer(palette = "Set1")
}

survplot_age2 <- make_survival_plot(Age2.tbsurv, mod_age2_tbsurv, age_label = "2", plot_label = "A) ")
survplot_age5 <- make_survival_plot(Age5.tbsurv, mod_age5_tbsurv, age_label = "5", plot_label = "B) ")
survplot_age10 <- make_survival_plot(Age10.tbsurv, mod_age10_tbsurv, age_label = "10", plot_label = "C) ")
survplot_age12 <- make_survival_plot(Age12.tbsurv, mod_age12_tbsurv, age_label = "12", plot_label = "D) ")


survplot_age2_noleg<-survplot_age2  +theme(legend.position = "none")
survplot_age5_noleg<-survplot_age5  +theme(legend.position = "none")
survplot_age10_noleg<-survplot_age10 +theme(legend.position = "none")
survplot_age12_noleg<-survplot_age12 +theme(legend.position = "none")

legend_survplot <- get_legend(survplot_age2 + theme(legend.position = "right"))

survplot_combined <- plot_grid(
  survplot_age2_noleg, survplot_age5_noleg,
  survplot_age10_noleg, survplot_age12_noleg,
  labels = NULL,
  ncol = 2
)

survplot_combinedleg <- plot_grid(
  survplot_combined,
  legend_survplot,
  ncol = 2,
  rel_heights = c(1, 0.1),
  rel_widths = c(1,0.1)
)
survplot_combinedleg


