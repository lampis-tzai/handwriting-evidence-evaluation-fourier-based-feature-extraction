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

comp_writers = t(combn(unique(adoq_data$N), 2))


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

stan_model_niw <- stan_model(file = "Stan_code/niw.stan", model_name = "niw")
stan_model_nlkj <- stan_model(file = "Stan_code/normal_lkj_model.stan", model_name = "normal_lkj_model")
stan_model_manova_iw <- stan_model(file = "Stan_code/MANOVA_iw_model.stan", model_name = "MANOVA_iw")
stan_model_manova_lkj <- stan_model(file = "Stan_code/MANOVA_lkj_model.stan", model_name = "MANOVA_lkj")


write_xlsx(data.frame(),"Stan_code/different_source_results_iter.xlsx")

different_source_def <- function(character_data,composition,w){
  
  df_all=data.frame()
  
  writer_data_1 = character_data[(character_data$N == composition[w,1]),]
  
  writer_data_2 = character_data[(character_data$N == composition[w,2]),]
  
  background_data = character_data[!(character_data$N %in% c(composition[w,1],
                                                             composition[w,2])),]
  
  background_stats_niw_all = background_statistics_niw(background_data)
  background_stats_niw_a = background_statistics_niw(background_data[(background_data$Lettre==1),])
  background_stats_niw_d = background_statistics_niw(background_data[(background_data$Lettre==2),])
  background_stats_niw_o = background_statistics_niw(background_data[(background_data$Lettre==3),])
  background_stats_niw_q = background_statistics_niw(background_data[(background_data$Lettre==4),])
  background_stats_br = background_statistics_br(background_data)
  
  for (iter_for_eval in (1)){   
    
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
    
    niw_conjugate_a = niw_conjugate(questioned_data[(questioned_data$Lettre==1),],
                                    suspect_data[(suspect_data$Lettre==1),],
                                    background_stats_niw_a)
    
    niw_a = normal_iW(questioned_data[(questioned_data$Lettre==1),],
                      suspect_data[(suspect_data$Lettre==1),],
                      background_stats_niw_a)
    
    nlkj_a = normal_lkj(questioned_data[(questioned_data$Lettre==1),],
                        suspect_data[(suspect_data$Lettre==1),],
                        background_stats_niw_a)
    
    df_niw_a <- data.frame(niw_a_conjugate  = niw_conjugate_a,
                           niw_a  = niw_a,
                           nlkj_a = nlkj_a)
    
    
    niw_conjugate_d = niw_conjugate(questioned_data[(questioned_data$Lettre==2),],
                                    suspect_data[(suspect_data$Lettre==2),],
                                    background_stats_niw_d)
    
    niw_d = normal_iW(questioned_data[(questioned_data$Lettre==2),],
                      suspect_data[(suspect_data$Lettre==2),],
                      background_stats_niw_d)
    
    nlkj_d = normal_lkj(questioned_data[(questioned_data$Lettre==2),],
                        suspect_data[(suspect_data$Lettre==2),],
                        background_stats_niw_d)
    
    df_niw_d <- data.frame(niw_d_conjugate  = niw_conjugate_d,
                           niw_d  = niw_d,
                           nlkj_d = nlkj_d)
    
    niw_conjugate_o = niw_conjugate(questioned_data[(questioned_data$Lettre==3),],
                                    suspect_data[(suspect_data$Lettre==3),],
                                    background_stats_niw_o)
    
    niw_o = normal_iW(questioned_data[(questioned_data$Lettre==3),],
                      suspect_data[(suspect_data$Lettre==3),],
                      background_stats_niw_o)
    
    nlkj_o = normal_lkj(questioned_data[(questioned_data$Lettre==3),],
                        suspect_data[(suspect_data$Lettre==3),],
                        background_stats_niw_o)
    
    df_niw_o <- data.frame(niw_o_conjugate  = niw_conjugate_o,
                           niw_o  = niw_o,
                           nlkj_o = nlkj_o)
    
    niw_conjugate_q = niw_conjugate(questioned_data[(questioned_data$Lettre==4),],
                                    suspect_data[(suspect_data$Lettre==4),],
                                    background_stats_niw_q)
    
    niw_q = normal_iW(questioned_data[(questioned_data$Lettre==4),],
                      suspect_data[(suspect_data$Lettre==4),],
                      background_stats_niw_q)
    
    nlkj_q = normal_lkj(questioned_data[(questioned_data$Lettre==4),],
                        suspect_data[(suspect_data$Lettre==4),],
                        background_stats_niw_q)
    
    df_niw_q <- data.frame(niw_q_conjugate  = niw_conjugate_q,
                           niw_q  = niw_q,
                           nlkj_q = nlkj_q)
    
    niw_conjugate_all = niw_conjugate(questioned_data,
                                      suspect_data,
                                      background_stats_niw_all)
    
    niw_all = normal_iW(questioned_data,
                        suspect_data,
                        background_stats_niw_all)
    
    nlkj = normal_lkj(questioned_data,
                      suspect_data,
                      background_stats_niw_all)
    
    df_niw_all <- data.frame(niw_all_conjugate  = niw_conjugate_all,
                             niw_all  = niw_all,
                             nlkj_all = nlkj)
    
    
    manova_conjugate = MANOVA_conjugate(questioned_data,
                                        suspect_data,
                                        background_stats_br)
    
    manova_iw = MANOVA_iw(questioned_data,
                          suspect_data,
                          background_stats_br)
    
    manova_lkj = MANOVA_LKJ(questioned_data,
                            suspect_data,
                            background_stats_br)
    
    df_br <- data.frame(manova_conjugate  = manova_conjugate,
                        manova_iw  = manova_iw,
                        manova_lkj = manova_lkj)
    
    df_new = cbind(data.frame(writer_1=composition[w,1],
                              writer_2=composition[w,2],
                              a_questioned_per = random_percentage_list[1],
                              d_questioned_per = random_percentage_list[2],
                              o_questioned_per = random_percentage_list[3],
                              q_questioned_per = random_percentage_list[4]),
                   df_niw_a,df_niw_d,df_niw_o,df_niw_q,
                   df_niw_all,df_br)
    
    df_all = rbind(df_all,df_new)
    print(paste0("Writer1: ",composition[w,1]," vs Writer2: ",
                 composition[w,2],", iteration:",iter_for_eval))
  }
  dsr_i <- read_excel("Stan_code/different_source_results_iter.xlsx")
  dsr_i = rbind(dsr_i,df_all)
  write_xlsx(dsr_i,"Stan_code/different_source_results_iter.xlsx")

  return(df_all)
}

cl <- makeCluster(6,
                  outfile="C:/Users/ltzai/Desktop/PhD/Handwritten_Loop_characters/Handwriting_Multivariate_Approach/Stan_code/log.txt")
clusterExport(cl,
              list("background_statistics_br","background_statistics_niw","abind",
                   "is.positive.definite","nearPD","fitdistr","lmvgamma","inv",
                   "stan_model_niw", "stan_model_nlkj",
                   "stan_model_manova_iw","stan_model_manova_lkj",
                   "assess_BF","sampling","extract","bridge_sampler",
                   "marginal_likelihood_niw_conjugate", "niw_conjugate",
                   "marginal_likelihood_manova_conjugate","MANOVA_conjugate",
                   "normal_iW","normal_lkj","MANOVA_iw","MANOVA_LKJ",
                   "read_excel","write_xlsx"),
              envir=globalenv())

w.list <- sapply(1:nrow(comp_writers), list)

system.time({saves = parLapply(cl, w.list,
                               different_source_def,
                               character_data = adoq_data,
                               composition = comp_writers)})

stopCluster(cl)
df_all <- do.call("rbind", saves)


write_xlsx(df_all,"Stan_code/different_source_results.xlsx")


dsr <- read_excel("Stan_code/Experimental_data/different_source_results.xlsx")


dsr = as.data.frame(dsr)


indx <- apply(dsr, 2, function(x) any(is.na(x) | is.infinite(x)))
colnames(dsr)[indx]
dsr[sapply(dsr, is.infinite)] <- NA
dsr[is.na(dsr)] = 0

dsr_binary = dsr
dsr_binary[,7:ncol(dsr_binary)] = dsr_binary[,7:ncol(dsr_binary)]>0

View(as.data.frame(colSums(dsr_binary[,7:ncol(dsr_binary)])))
View(round(as.data.frame(colMeans(dsr_binary[,7:ncol(dsr_binary)])),3)*100)

dsr_grouped = dsr_binary %>%
  group_by(writer_1,writer_2) %>%
  summarise_all("sum")

dsr_grouped = as.data.frame(dsr_grouped)
View(dsr_grouped)


colnames(dsr) = c("writer_1","writer_2", "a_questioned_per", "d_questioned_per",
                  "o_questioned_per", "q_questioned_per",
                  "Normal-conjugate_a",
                  "Normal-inverse-Wishart_a",
                  "Normal-LN-LKJ_a",
                  "Normal-conjugate_d",
                  "Normal-inverse-Wishart_d",
                  "Normal-LN-LKJ_d",
                  "Normal-conjugate_o",
                  "Normal-inverse-Wishart_o",
                  "Normal-LN-LKJ_o",
                  "Normal-conjugate_q",
                  "Normal-inverse-Wishart_q",
                  "Normal-LN-LKJ_q",
                  "Normal-conjugate_all",
                  "Normal-inverse-Wishart_all",
                  "Normal-LN-LKJ_all",
                  "MANOVA-conjugate",
                  "MANOVA-inverse-Wishart",
                  "MANOVA-LN-LKJ")


library(reshape2)
melt_dsr_df <- melt(dsr, id = c("writer_1","writer_2","a_questioned_per",
                                "d_questioned_per","o_questioned_per",
                                "q_questioned_per"), 
                    variable.name = 'model') 



library(stringr)
split_data = str_split_fixed(melt_dsr_df$model, "-", 2)
split_data2 = str_split_fixed(split_data[,2],"_",2)
melt_dsr_df$model = paste0(split_data[,1],'_',split_data2[,2])
melt_dsr_df['Prior_approach'] = split_data2[,1]

library(ggplot2)

melt_dsr_df$Prior_approach <- factor(melt_dsr_df$Prior_approach, 
                                  levels = unique(melt_dsr_df$Prior_approach))




melt_dsr_df$model <- factor(melt_dsr_df$model, 
                            levels = c("Normal_a",
                                       "Normal_d",
                                       "Normal_o",
                                       "Normal_q",
                                       "Normal_all",
                                       "MANOVA_"))

levels(melt_dsr_df$model) <- c("Normal a","Normal d","Normal o","Normal q",
                              "Normal all", "MANOVA")


levels(melt_dsr_df$Prior_approach) <- c("(1) NIW Conjugate",
                                        "(2) NIW Hierarchical",
                                        "(3) Normal-LogNormal-LKJ")

melt_dsr_df$value = as.numeric(melt_dsr_df$value)

library(latex2exp)
plot = ggplot(melt_dsr_df,
              aes(x = Prior_approach, y = value, fill = Prior_approach)) +
  geom_boxplot() +
  facet_wrap(~model,ncol = 6) +
  scale_y_continuous(name = TeX(r"(\textbf{LogBF})"), limits = c(-700, 100)) +
  scale_x_discrete(labels = c("(1)","(2)","(3)"), name = "Models")+
  #labs(title="Logarithmic Bayes Factors for \n Different Source Comparisons") + 
  #theme(plot.title = element_text(hjust = 0.5))+
  scale_fill_brewer(palette="Set2")+
  geom_hline(yintercept = 0, color = 'brown',lty='dashed')+ 
  labs(fill = "Prior approach") +
  theme(#plot.title = element_text(hjust = 0.5),
    #legend.spacing.y = unit(0.5, 'cm'),
    strip.text = element_text(size = 12,face="bold"),
    axis.title=element_text(size=11,face="bold"),
    legend.text = element_text(size=10),
    #legend.title=element_blank())+
  legend.title = element_text(size=15,face="bold"))+
  guides(fill = guide_legend(byrow = TRUE),
         colour = guide_legend(override.aes = list(size=5)))

jpeg("Stan_code/plots/ds_boxplot.jpg",width=3920, height=2000, res=300)
plot
dev.off()


