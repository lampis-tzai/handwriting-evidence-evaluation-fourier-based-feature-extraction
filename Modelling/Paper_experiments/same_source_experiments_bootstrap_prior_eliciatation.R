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


write_xlsx(data.frame(),"Stan_code/same_source_results_background_bootsrap_iter.xlsx")

same_source_bootstrap_def <- function(character_data,w){
  df_all=data.frame()
  
  writer_data = character_data[(character_data$N==w),]
  background_data_all = character_data[(character_data$N!=w),]
  
  questioned_data = data.frame()
  suspect_data = data.frame()
  random_percentage_list = c()
  for (c in 1:4){
    writer_data_c = writer_data[(writer_data$Lettre==c),]
    random_percentage = runif(1,0.35,0.65)
    random_percentage_list = c(random_percentage_list,random_percentage)
    smp_size <- round(random_percentage * nrow(writer_data_c))
    suspect_ind <- sample(seq_len(nrow(writer_data_c)), size = smp_size)
    
    questioned_data = rbind(questioned_data,writer_data_c[suspect_ind, ])
    suspect_data = rbind(suspect_data, writer_data_c[-suspect_ind, ])
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
  
  df_new = cbind(data.frame(writer=w,
                            a_questioned_per = random_percentage_list[1],
                            d_questioned_per = random_percentage_list[2],
                            o_questioned_per = random_percentage_list[3],
                            q_questioned_per = random_percentage_list[4]),
                 bootstrap_id = "all",
                 df_br)
  
  df_all = rbind(df_all,df_new)
  
  for (iter_for_eval in (1:30)){
    
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
    
    df_new = cbind(data.frame(writer=w,
                              a_questioned_per = random_percentage_list[1],
                              d_questioned_per = random_percentage_list[2],
                              o_questioned_per = random_percentage_list[3],
                              q_questioned_per = random_percentage_list[4]),
                   bootstrap_id = iter_for_eval,
                   df_br)
    
    print(paste0("Writer:",w,", iteration:",iter_for_eval))
    
    df_all = rbind(df_all,df_new)
  }
  ssr_i <- read_excel("Stan_code/same_source_results_background_bootsrap_iter.xlsx")
  ssr_i = rbind(ssr_i,df_all)
  write_xlsx(ssr_i,"Stan_code/same_source_results_background_bootsrap_iter.xlsx")
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

w.list <- sapply(1:length(unique(adoq_data$N)), list)

system.time({saves = parLapply(cl, w.list,
                               same_source_bootstrap_def,
                               character_data = adoq_data)})

stopCluster(cl)

df_all <- do.call("rbind", saves)

write_xlsx(df_all,"Stan_code/same_source_results_background_bootsrap.xlsx")


ssr <- read_excel("Stan_code/Experimental_data/same_source_results_background_bootsrap.xlsx")



ssr = as.data.frame(ssr)


# indx <- apply(ssr, 2, function(x) any(is.na(x) | is.infinite(x)))
# colnames(ssr)[indx]
# ssr[sapply(ssr, is.infinite)] <- NA
# ssr[is.na(ssr)] = -10000

ssr_binary = ssr
ssr_binary[,6:ncol(ssr_binary)] = ssr_binary[,6:ncol(ssr_binary)]<0

View(as.data.frame(colSums(ssr_binary[,6:ncol(ssr_binary)])))

View(round(as.data.frame(colMeans(ssr_binary[,6:ncol(ssr_binary)])),3))

ssr_grouped = ssr_binary %>%
  group_by(writer) %>%
  summarise_all("mean")

ssr_grouped = as.data.frame(ssr_grouped)
View(ssr_grouped)



colnames(ssr) = c("writer", "a_questioned_per", "d_questioned_per",
                  "o_questioned_per", "q_questioned_per",
                  "bootstrap_id",
                  "MANOVA-conjugate",
                  "MANOVA-inverse-Wishart",
                  "MANOVA-LN-LKJ")


library(reshape2)
melt_ssr_df <- melt(ssr, id = c("writer","a_questioned_per",
                                "d_questioned_per","o_questioned_per",
                                "q_questioned_per","bootstrap_id"), 
                    variable.name = 'model') 


library(stringr)
split_data = str_split_fixed(melt_ssr_df$model, "-", 2)
split_data2 = str_split_fixed(split_data[,2],"_",2)
melt_ssr_df$model = paste0(split_data[,1],'_',split_data2[,2])
melt_ssr_df['Prior_approach'] = split_data2[,1]

library(ggplot2)

melt_ssr_df$Prior_approach <- factor(melt_ssr_df$Prior_approach, 
                                     levels = unique(melt_ssr_df$Prior_approach))




melt_ssr_df$model <- factor(melt_ssr_df$model, 
                            levels = c("MANOVA_"))

levels(melt_ssr_df$model) <- c("MANOVA")


levels(melt_ssr_df$Prior_approach) <-c("(1) NIW Conjugate",
                                       "(2) NIW Hierarchical",
                                       "(3) Normal-LogNormal-LKJ")

melt_ssr_df$value = as.numeric(melt_ssr_df$value)

melt_ssr_df['is_reference'] <- melt_ssr_df$bootstrap_id=='all'

melt_ssr_df_less <- melt_ssr_df[melt_ssr_df$writer %in% 7:13,]

melt_ssr_df_less$writer <- factor(melt_ssr_df_less$writer, 
                             levels = sort(unique(melt_ssr_df_less$writer)))
levels(melt_ssr_df_less$writer) <- paste0('W', levels(melt_ssr_df_less$writer))

library(latex2exp)
library(scales)

plot = ggplot(melt_ssr_df_less,
              aes(x = Prior_approach, y = value, fill = Prior_approach)) +
  geom_boxplot(alpha = 0.2) +
  geom_point(
    data = melt_ssr_df_less %>% filter(is_reference),
    aes(x = Prior_approach, y = value, color = Prior_approach),
    shape = 8, size = 2, stroke = 1, show.legend = FALSE
  ) +
  facet_wrap(~writer,ncol = 13, scale="free_y") +
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
      seq(ymin, ymax, by = 15)                  # tick marks every 15 units
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

jpeg("Stan_code/plots/ss_boxplot_bootstrap2.jpg",width=3920, height=2000, res=300)
plot
dev.off()
