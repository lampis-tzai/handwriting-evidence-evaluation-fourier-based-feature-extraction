setwd("C:/Users/ltzai/Desktop/PhD/Handwritten_Loop_characters/Handwriting_Multivariate_Approach")
library(readxl)
library(dplyr)
library(MASS)
library(matrixcalc)
library(Matrix)
library(abind)
library(matlib)
library(parallel)
library(writexl)
library(rstan)
library(bridgesampling)
source('Stan_code/Stan_BF_calculation.R')



set.seed(2)
adoq_data <- read_excel("Data/adoq colonnes.xls")
adoq_data = as.data.frame(adoq_data)

#coefficients
deg2rad <- function(deg) {deg * pi/180}
coef_data = data.frame(Surface = adoq_data$Surface)

for(h in 1:4){
  ampl = paste0('Ampl',h)
  phase = paste0('Phase',h)
  a_h = adoq_data[,ampl]*cos(deg2rad(adoq_data[,phase]))
  b_h = adoq_data[,ampl]*sin(deg2rad(adoq_data[,phase]))
  har_coef = data.frame(a_h,b_h)
  colnames(har_coef) = c(paste0('a_',h),paste0('b_',h))
  coef_data = cbind(coef_data,har_coef)
}

adoq_data = cbind(adoq_data[,1:4],scale(coef_data))

background_statistics_niw <- function(background_data){
  p=9
  nw.min =  p + 2
  nw_hat = nw.min
  
  mu_hat=matrix(colMeans(background_data[,5:ncol(background_data)]),
                nrow = 1)
  
  B_hat = cov(background_data[,5:ncol(background_data)])
  if (!is.positive.definite(B_hat)){B_hat = as.matrix(nearPD(B_hat)$mat)}
  
  Sw = 0
  for (w in unique(background_data$N)){
    df_writer = background_data[(background_data$N==w),]
    var_data = unname(as.matrix(df_writer[,5:ncol(df_writer)]))
    Cov.this = cov(var_data)*(nrow(var_data)-1) 
    Sw <- Sw + Cov.this
  } 
  
  W_hat <- Sw/(nrow(background_data) - length(unique(background_data$N)))
  U_hat <- W_hat * (nw_hat - p - 1) 
  
  loc <- mean(log(diag(W_hat)))
  sc <- sd(log(diag(W_hat)))
  eta=1
  
  return(list(mu_hat,B_hat,nw_hat,U_hat,loc,sc,eta))
}


background_statistics_br <- function(background_data){
  p=9
  l = 4
  nw.min =  p + 2
  nw_hat = nw.min
  
  a_data = background_data[(background_data$Lettre==1),5:ncol(background_data)]
  mu_hat=matrix(colMeans(a_data),nrow = 1)
  
  B_hat = cov(a_data)
  if (!is.positive.definite(B_hat)){B_hat = as.matrix(nearPD(B_hat)$mat)}
  
  beta_mu = array(0, dim=c(l,p))
  beta_cov = array(0, dim=c(p,p,l))
  for (l_id in 1:l){
    letter_data = as.matrix(unname(background_data[(
      background_data$Lettre==l_id),5:ncol(background_data)]))
    
    letter_diff = letter_data - matrix(mu_hat[col(letter_data)],ncol = p)
    beta_l = colMeans(letter_diff)
    beta_mu[l_id,] = beta_l
    
    B_hat_l = cov(letter_diff)
    if (!is.positive.definite(B_hat_l)){B_hat_l = as.matrix(nearPD(B_hat_l)$mat)}
    beta_cov[,,l_id] = B_hat_l
  }
  
  
  Sw = 0
  for (w in unique(background_data$N)){
    df_writer = background_data[(background_data$N==w),]
    var_data = unname(as.matrix(df_writer[,5:ncol(df_writer)]))
    Cov.this = cov(var_data)*(nrow(df_writer)-1)
    Sw <- Sw + Cov.this
    
  }
  
  W_hat <- Sw/(nrow(background_data) - length(unique(background_data$N)))
  U_hat <- W_hat * (nw_hat - p -1)
  
  loc <- mean(log(diag(W_hat)))
  sc <- sd(log(diag(W_hat)))
  eta=1
  
  return(list(mu_hat,B_hat,beta_mu,beta_cov,nw_hat,U_hat,loc,sc,eta))
}

stan_model_nlkj <- stan_model(file = "Stan_code/normal_lkj_model.stan", model_name = "normal_lkj_model")
stan_model_manova_lkj <- stan_model(file = "Stan_code/MANOVA_lkj_model.stan", model_name = "MANOVA_lkj")

write_xlsx(data.frame(),"Stan_code/different_eta_per_comparison_iter.xlsx")

different_source_eta_def <- function(character_data,composition,w){
  df_all=data.frame()
  
  writer_data_1 = character_data[(character_data$N == composition[w,1]),]
  
  writer_data_2 = character_data[(character_data$N == composition[w,2]),]
  
  background_data_all = character_data[!(character_data$N %in% c(composition[w,1],
                                                                 composition[w,2])),]
  
  questioned_data = data.frame()
  suspect_data = data.frame()
  
  random_percentage_list = c()
  for (c in 1:4){
    writer_data_1_c = writer_data_1[(writer_data_1$Lettre==c),]
    random_percentage = runif(1,0.35,0.65)
    random_percentage_list = c(random_percentage_list,random_percentage)
    smp_size <- round(random_percentage  * nrow(writer_data_1_c))
    ind <- sample(seq_len(nrow(writer_data_1_c)), size = smp_size)
    questioned_data = rbind(questioned_data,writer_data_1_c[ind, ])
    
    writer_data_2_c = writer_data_2[(writer_data_2$Lettre==c),]
    smp_size <- round((1-random_percentage) * nrow(writer_data_2_c))
    ind <- sample(seq_len(nrow(writer_data_2_c)), size = smp_size)
    suspect_data = rbind(suspect_data,writer_data_2_c[ind, ])
  }
  
  background_stats_niw = background_statistics_niw(background_data_all)
  background_stats_br = background_statistics_br(background_data_all)
  
  p = nrow(background_stats_niw[[2]])
  
  for (eta in c(1,2,5,10,20)){
    
    background_stats_niw[[7]]<- eta
    background_stats_br[[9]]<- eta
    
    BF_normal_lkj = normal_lkj(questioned_data,
                               suspect_data,
                               background_stats_niw)
    
    
    BF_manova_lkj = MANOVA_LKJ(questioned_data,
                            suspect_data,
                            background_stats_br)
    
    
    df_new = cbind(data.frame(writer_1=composition[w,1],
                              writer_2=composition[w,2],
                              lkj_eta = eta,
                              BF_normal_lkj,BF_manova_lkj))
    
    print(paste0("Writer1: ",composition[w,1]," vs Writer2: ",
                 composition[w,2],", eta:",eta))
    
    df_all = rbind(df_all,df_new)
    
  }
  
  ssr_i <- read_excel("Stan_code/different_eta_per_comparison_iter.xlsx")
  ssr_i = rbind(ssr_i,df_all)
  write_xlsx(ssr_i,"Stan_code/different_eta_per_comparison_iter.xlsx")
  return(df_all)
}

detectCores()
cl <- makeCluster(6,
                  outfile="C:/Users/ltzai/Desktop/PhD/Handwritten_Loop_characters/Handwriting_Multivariate_Approach/Stan_code/log.txt")
clusterExport(cl,
              list("background_statistics_niw","background_statistics_br","abind",
                   "is.positive.definite","nearPD","fitdistr","lmvgamma", "inv",
                   "stan_model_nlkj","stan_model_manova_lkj",
                   "normal_lkj","MANOVA_LKJ",
                   "assess_BF","sampling","extract","bridge_sampler",
                   "marginal_likelihood_manova_conjugate","MANOVA_conjugate",
                   "MANOVA_iw","MANOVA_LKJ",
                   "read_excel","write_xlsx"),
              envir=globalenv())



comp_writers = t(combn(unique(adoq_data$N), 2))
#comp_writers_cons <- comp_writers[c(51,52,58,sample(c(1:50,53:57,59:78),10)),]
w.list <- sapply(1:nrow(comp_writers), list)

system.time({saves = parLapply(cl, w.list,
                               different_source_eta_def,
                               character_data = adoq_data,
                               composition = comp_writers)})

stopCluster(cl)

df_all <- do.call("rbind", saves)

write_xlsx(df_all,"Stan_code/different_eta_per_comparison.xlsx")


dsr <- read_excel("Stan_code/Experimental_data/different_eta_per_comparison.xlsx")



dsr = as.data.frame(dsr)



#dsr = dsr[(dsr$writer_1==7) & (dsr$writer_2==8),]


library(reshape2)
melt_dsr_df <- melt(dsr, id = c("writer_1","writer_2","lkj_eta"), 
                    variable.name = 'Model') 



library(stringr)
split_data = str_split_fixed(melt_dsr_df$Model, "_", 3)

melt_dsr_df$Model = split_data[,2]
melt_dsr_df['BF_approach'] = split_data[,3]

library(ggplot2)
library(plyr)
data_summary <- function(data, varname, groupnames){
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}

melt_dsr_df2 <- data_summary(melt_dsr_df, varname="value", 
                             groupnames=c("lkj_eta", "BF_approach","Model"))



melt_dsr_df2$lkj_eta <- factor(melt_dsr_df2$lkj_eta, 
                           levels = unique(melt_dsr_df2$lkj_eta))

melt_dsr_df2$Model <- factor(melt_dsr_df2$Model, 
                             levels = unique(melt_dsr_df2$Model))
levels(melt_dsr_df2$Model)<- c('MANOVA','Normal')

melt_dsr_df2$BF_approach <- factor(melt_dsr_df2$BF_approach, 
                             levels = unique(melt_dsr_df2$BF_approach))

levels(melt_dsr_df2$BF_approach)<- c('Normal-LogNormal-LKJ')


melt_dsr_df2$value = as.numeric(melt_dsr_df2$value)


library(latex2exp)
p<- ggplot(melt_dsr_df2, aes(x=lkj_eta, y=value, group=Model, color=Model)) + 
  #geom_bar(stat="identity",color="black",
  #         position=position_dodge()) +
  geom_line() +
  geom_point(aes(shape=Model), size=2)+
  #geom_errorbar(aes(ymin=value-sd, ymax=value+sd), width=.2,
  #              position=position_dodge(0.05))+
  facet_grid(~BF_approach) +
  scale_y_continuous(name = TeX(r"(\textbf{LogBF})")) +
  scale_x_discrete(labels = abbreviate, name = TeX(r"($\eta$)"))+
  labs(title="Different Source Comparisons for Different values of LKJ parameter", fill = "Model") + # for \n Writers 7 and 8
  theme(plot.title = element_text(hjust = 0.5,face="bold"),
        axis.title=element_text(size=11,face="bold"),
        legend.title = element_text(size=15,face="bold"))+
  scale_color_brewer(palette="Accent")

p

jpeg("Stan_code/plots/ds_diffent_eta.jpg",width=2920, height=1080, res=300)
p
dev.off()


#####################################
## Different degrees of freedom old
####################################

dsr <- read_excel("Rjags_code/Experimental_data/different_source_results_DoF.xlsx")


dsr = as.data.frame(dsr)

#exclude naive
dsr = dsr[,c(1:7,9:12,14:17,19:22,24:27,29:32,34:37)]

colnames(dsr) = c("writer_1","writer_2", "a_questioned_per", "d_questioned_per",
                  "o_questioned_per", "q_questioned_per",
                  "Normal-inverse-Wishart_all_conjugate_approach",
                  "Normal-inverse-Wishart_all_generalized_harmonic_mean",
                  "Normal-inverse-Wishart_all_laplace_metropolis",
                  "Normal-inverse-Wishart_all_bridge_sampling",
                  "Normal-inverse-Wishart_a_conjugate_approach",
                  "Normal-inverse-Wishart_a_generalized_harmonic_mean",
                  "Normal-inverse-Wishart_a_laplace_metropolis",
                  "Normal-inverse-Wishart_a_bridge_sampling",
                  "Normal-inverse-Wishart_d_conjugate_approach",
                  "Normal-inverse-Wishart_d_generalized_harmonic_mean",
                  "Normal-inverse-Wishart_d_laplace_metropolis",
                  "Normal-inverse-Wishart_d_bridge_sampling",
                  "Normal-inverse-Wishart_o_conjugate_approach",
                  "Normal-inverse-Wishart_o_generalized_harmonic_mean",
                  "Normal-inverse-Wishart_o_laplace_metropolis",
                  "Normal-inverse-Wishart_o_bridge_sampling",
                  "Normal-inverse-Wishart_q_conjugate_approach",
                  "Normal-inverse-Wishart_q_generalized_harmonic_mean",
                  "Normal-inverse-Wishart_q_laplace_metropolis",
                  "Normal-inverse-Wishart_q_bridge_sampling",
                  "Bayesian_MANOVA_conjugate_approach",
                  "Bayesian_MANOVA_generalized_harmonic_mean",
                  "Bayesian_MANOVA_laplace_metropolis",
                  "Bayesian_MANOVA_bridge_sampling","DoF")

dsr[is.na(dsr)] = 0



dsr = dsr[, c( "writer_1","writer_2", "a_questioned_per", "d_questioned_per",
               "o_questioned_per", "q_questioned_per",
               "Normal-inverse-Wishart_all_conjugate_approach",
               "Normal-inverse-Wishart_all_bridge_sampling",
               "Bayesian_MANOVA_conjugate_approach",
               "Bayesian_MANOVA_bridge_sampling","DoF")]

#dsr = dsr[(dsr$writer_1==7) & (dsr$writer_2==8),]


library(reshape2)
melt_dsr_df <- melt(dsr, id = c("writer_1","writer_2","a_questioned_per",
                                "d_questioned_per","o_questioned_per",
                                "q_questioned_per","DoF"), 
                    variable.name = 'Model') 



library(stringr)
split_data = str_split_fixed(melt_dsr_df$Model, "_", 3)

melt_dsr_df$Model = str_c(split_data[,1],'_',split_data[,2])
melt_dsr_df['BF_approach'] = split_data[,3]

library(ggplot2)
library(plyr)
data_summary <- function(data, varname, groupnames){
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}

melt_dsr_df2 <- data_summary(melt_dsr_df, varname="value", 
                             groupnames=c("DoF", "BF_approach","Model"))


niw_lm = melt_dsr_df2[melt_dsr_df2$Model=="Normal-inverse-Wishart_all",]

summary(lm(value~ DoF, data = niw_lm[niw_lm$BF_approach=='bridge_sampling',]))


manova_lm = melt_dsr_df2[melt_dsr_df2$Model=="Bayesian_MANOVA",]
summary(lm(value~ DoF, data = manova_lm[manova_lm$BF_approach=='bridge_sampling',]))



melt_dsr_df2$DoF <- factor(melt_dsr_df2$DoF, 
                           levels = unique(melt_dsr_df2$DoF))

melt_dsr_df2$Model <- factor(melt_dsr_df2$Model, 
                             levels = unique(melt_dsr_df2$Model))

levels(melt_dsr_df2$Model) <- c("MANOVA","Normal")

melt_dsr_df2$BF_approach <- factor(melt_dsr_df2$BF_approach, 
                                   levels = c("conjugate_approach",
                                              "bridge_sampling"))

levels(melt_dsr_df2$BF_approach) <- c("NIW Conjugate",
                                      "NIW Hierarchical")
melt_dsr_df2$value = as.numeric(melt_dsr_df2$value)


library(latex2exp)
p<- ggplot(melt_dsr_df2, aes(x=DoF, y=value, group=Model, color=Model)) + 
  #geom_bar(stat="identity",color="black",
  #         position=position_dodge()) +
  geom_line() +
  geom_point(aes(shape=Model), size=2)+
  #geom_errorbar(aes(ymin=value-sd, ymax=value+sd), width=.2,
  #              position=position_dodge(0.05))+
  facet_grid(~BF_approach) +
  scale_y_continuous(name = TeX(r"(\textbf{LogBF})")) +
  scale_x_discrete(labels = abbreviate, name = "Degrees of Freedom")+
  labs(title="Different Source Comparisons for Different Degrees of Freedom", fill = "Model") + # for \n Writers 7 and 8
  theme(plot.title = element_text(hjust = 0.5,face="bold"),
        axis.title=element_text(size=11,face="bold"),
        legend.title = element_text(size=15,face="bold"))+
  scale_color_brewer(palette="Accent")

p

jpeg("Stan_code/plots/niw_ds_dof.jpg",width=2920, height=1080, res=300)
p
dev.off()

