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

stan_model_manova_iw <- stan_model(file = "Stan_code/MANOVA_iw_model.stan", model_name = "MANOVA_iw")
stan_model_manova_lkj <- stan_model(file = "Stan_code/MANOVA_lkj_model.stan", model_name = "MANOVA_lkj")


write_xlsx(data.frame(),"Stan_code/different_source_results_background_bootsrap_iter.xlsx")

different_source_bootstrap_def <- function(character_data,composition,w){
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
  
  background_stats_br = background_statistics_br(background_data_all)
  
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
                 bootstrap_id = "all",
                 df_br)
  
  df_all = rbind(df_all,df_new)
  
  print(paste0("Writer1: ",composition[w,1]," vs Writer2: ",
               composition[w,2],", iteration:","all"))
  
  for (iter_for_eval in (1:3)){
    
    n_subsample <- length(unique(background_data_all$Writer))#ceiling(0.5 * length(unique(background_data_all$Writer)))
    
    # Step 1: Subsample writers
    sampled_writers <- sample(unique(background_data_all$Writer), size = n_subsample, replace = FALSE)
    
    background_data <- do.call(rbind, lapply(sampled_writers, function(w) {
      writer_data_db <- background_data_all[background_data_all$Writer == w, ]
      n_subsample <- ceiling(0.5 * nrow(writer_data_db))
      writer_data_db[sample(nrow(writer_data_db), size =n_subsample, replace = TRUE), ]
    }))
    
    background_stats_br = background_statistics_br(background_data)
    
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
                   bootstrap_id = iter_for_eval,
                   df_br)
    
    print(paste0("Writer1: ",composition[w,1]," vs Writer2: ",
                 composition[w,2],", iteration:",iter_for_eval))
    
    df_all = rbind(df_all,df_new)
  }
  dsr_i <- read_excel("Stan_code/different_source_results_background_bootsrap_iter.xlsx")
  dsr_i = rbind(dsr_i,df_all)
  write_xlsx(dsr_i,"Stan_code/different_source_results_background_bootsrap_iter.xlsx")
  return(df_all)
}

detectCores()
cl <- makeCluster(6,
                  outfile="C:/Users/ltzai/Desktop/PhD/Handwritten_Loop_characters/Handwriting_Multivariate_Approach/Stan_code/log.txt")
clusterExport(cl,
              list("background_statistics_br","abind",
                   "is.positive.definite","nearPD","fitdistr","lmvgamma", "inv",
                   "stan_model_manova_iw","stan_model_manova_lkj",
                   "assess_BF","sampling","extract","bridge_sampler",
                   "marginal_likelihood_manova_conjugate","MANOVA_conjugate",
                   "MANOVA_iw","MANOVA_LKJ",
                   "read_excel","write_xlsx"),
              envir=globalenv())



comp_writers = t(combn(unique(adoq_data$N), 2))
comp_writers_cons <- comp_writers[c(51,52,58,sample(c(1:50,53:57,59:78),10)),]
w.list <- sapply(1:nrow(comp_writers_cons), list)

system.time({saves = parLapply(cl, w.list,
                               different_source_bootstrap_def,
                               character_data = adoq_data,
                               composition = comp_writers_cons)})

stopCluster(cl)

df_all <- do.call("rbind", saves)

write_xlsx(df_all,"Stan_code/different_source_results_background_bootsrap.xlsx")


dsr <- read_excel("Stan_code/Experimental_data/different_source_results_background_bootsrap.xlsx")



dsr = as.data.frame(dsr)


#indx <- apply(dsr, 2, function(x) any(is.na(x) | is.infinite(x)))
#colnames(dsr)[indx]
#dsr[sapply(dsr, is.infinite)] <- NA
#dsr[is.na(dsr)] = -10000

dsr_binary = dsr
dsr_binary[,8:ncol(dsr_binary)] = dsr_binary[,8:ncol(dsr_binary)]>0

View(as.data.frame(colSums(dsr_binary[,8:ncol(dsr_binary)])))
View(round(as.data.frame(colMeans(dsr_binary[,8:ncol(dsr_binary)])),3))

dsr_grouped = dsr_binary %>%
  group_by(writer_1,writer_2) %>%
  summarise(across(where(~ is.numeric(.x) | is.logical(.x)), sum), .groups = "drop")


dsr_grouped = as.data.frame(dsr_grouped)
View(dsr_grouped)



colnames(dsr) = c("writer_1","writer_2", "a_questioned_per", "d_questioned_per",
                  "o_questioned_per", "q_questioned_per",
                  "bootstrap_id",
                  "MANOVA-conjugate",
                  "MANOVA-inverse-Wishart",
                  "MANOVA-LN-LKJ")


library(reshape2)
melt_dsr_df <- melt(dsr, id = c("writer_1","writer_2","a_questioned_per",
                                "d_questioned_per","o_questioned_per",
                                "q_questioned_per","bootstrap_id"), 
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
                            levels = c("MANOVA_"))

levels(melt_dsr_df$model) <- c("MANOVA")


levels(melt_dsr_df$Prior_approach) <- c("(1) NIW Conjugate",
                                        "(2) NIW Hierarchical",
                                        "(3) Normal-LogNormal-LKJ")

melt_dsr_df$value = as.numeric(melt_dsr_df$value)

melt_dsr_df['is_reference'] <- melt_dsr_df$bootstrap_id=='all'

melt_dsr_df_less = melt_dsr_df[!((melt_dsr_df$writer_1 %in% 6:8) | (melt_dsr_df$writer_2 %in% 6:8)),]

melt_dsr_df_less$writer_1 <- paste0('W', melt_dsr_df_less$writer_1)
melt_dsr_df_less$writer_2 <- paste0('W', melt_dsr_df_less$writer_2)

melt_dsr_df_less['writer']<- paste0(melt_dsr_df_less$writer_1,'-',melt_dsr_df_less$writer_2)

library(latex2exp)
plot = ggplot(melt_dsr_df_less,
              aes(x = Prior_approach, y = value, fill = Prior_approach)) +
  geom_boxplot(alpha = 0.3) +
  geom_point(
    data = melt_dsr_df_less %>% filter(is_reference),
    aes(x = Prior_approach, y = value, color = Prior_approach),
    shape = 8, size = 2, stroke = 1, show.legend = FALSE
  ) +
  facet_wrap(~writer,ncol = 13, scales = "free_y") +
  scale_y_continuous(
    name = TeX(r"(\textbf{LogBF})"),
    limits = function(lims) {
      ymin <- floor(lims[1])                    # round down the min
      ymax <- ymin + 50                         # set max = min + 20
      c(ymin, ymax)
    },
    breaks = function(lims) {
      ymin <- floor(lims[1])
      ymax <- ymin + 50
      seq(ymin, ymax, by = 15)                  # tick marks every 5 units
    }
  ) +
  scale_x_discrete(labels = c("(1)","(2)","(3)"), name = "Models")+
  #labs(title="Logarithmic Bayes Factors for \n Different Source Comparisons") + 
  #theme(plot.title = element_text(hjust = 0.5))+
  scale_fill_brewer(palette = "Set2") +
  scale_color_brewer(palette = "Set2") +
  #geom_hline(yintercept = 0, color = 'brown',lty='dashed')+ 
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

jpeg("Stan_code/plots/ds_boxplot_bootstrap2.jpg",width=3920, height=2000, res=300)
plot
dev.off()
