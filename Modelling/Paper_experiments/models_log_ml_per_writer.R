setwd("C:/Users/lampis/Desktop/PhD/Handwritten_Loop_characters/handwriting-evidence-evaluation-fourier-based-feature-extraction/Modelling")
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
source('Paper_experiments/Stan_BF_calculation.R')



set.seed(2)
db <- read_excel("IAM_fourier_features_dataset/DB_loop_handwriting_ls.xlsx")
db = as.data.frame(db)

db[,1] = log(db[,1])

db = cbind(scale(db[,1:9]),db[,10:ncol(db)])

background_statistics_niw <- function(background_data){
  
  
  p=9
  nw.min = p + 2
  nw_hat = nw.min
  
  mu_hat=matrix(colMeans(do.call(rbind, lapply(unique(background_data$writer_id), function(w)
    colMeans(background_data[background_data$writer_id == w, 1:p])))), nrow = 1)
  
  S = 0
  Sw = 0
  for (w in unique(background_data$writer_id)){
    df_writer = background_data[(background_data$writer_id==w),]
    if (nrow(df_writer)>2){
      var_data = unname(as.matrix(df_writer[,1:p]))
      theta_w = matrix(colMeans(var_data), nrow = 1)
      S.this <- (t(theta_w - mu_hat) %*% (theta_w - mu_hat))
      S <- S + S.this
      Cov.this = cov(var_data)*(nrow(df_writer)-1) 
      Sw <- Sw + Cov.this
    }
  } 
  
  B_hat = S/(length(unique(background_data$writer_id)) - 1)
  #B_hat = cov(background_data[,1:p])
  if (any(is.na(B_hat)) || any(is.nan(B_hat))) {
    B_hat[is.na(B_hat) | is.nan(B_hat)] <- 0
    diag(B_hat) <- pmax(diag(B_hat), 1e-6)
  }
  if (!is.positive.definite(B_hat)) { B_hat <- as.matrix(nearPD(B_hat)$mat) }
  
  W_hat <- Sw/(nrow(background_data) - length(unique(background_data$writer_id)))
  U_hat <- W_hat*(nw_hat-p-1)
  
  
  
  log_sd_mat <- sapply(unique(background_data$writer_id), function(w) {
    df_w <- background_data[background_data$writer_id == w, 1:p]
    if (nrow(df_w) <= 3) return(rep(NA_real_, p))  # skip degenerate writers
    log(sqrt(diag(cov(as.matrix(df_w)))))
  })
  
  # Remove degenerate writers (columns with NA)
  log_sd_mat <- log_sd_mat[, colSums(is.na(log_sd_mat)) == 0]
  
  loc <- rowMeans(log_sd_mat)
  sc  <- apply(log_sd_mat, 1, sd)
  
  eta <- 1
  
  return(list(mu_hat,B_hat,nw_hat,U_hat,loc,sc,eta))
}


background_statistics_br <- function(background_data){
  
  p=9
  l = length(unique(background_data$character))
  nw.min = p + 2
  nw_hat = nw.min
  
  a_data = background_data[(background_data$character==1),]
  mu_hat=matrix(colMeans(do.call(rbind, lapply(unique(a_data$writer_id), function(w)
    colMeans(a_data[a_data$writer_id == w, 1:p])))), nrow = 1)
  
  S <- 0
  n_writers <- 0
  for (w in unique(background_data$writer_id)){
    df_writer = background_data[(
      background_data$character==1)& (background_data$writer_id==w),]
    if (nrow(df_writer)>2){
      n_writers <- n_writers+1
      var_data <- unname(as.matrix(df_writer[,1:p]))
      theta_w <- matrix(colMeans(var_data), nrow = 1)
      S.this <- (t(theta_w - mu_hat) %*% (theta_w - mu_hat))
      S <- S + S.this
    }
  } 
  
  B_hat = S/(n_writers - 1)
  #B_hat = cov(background_data[,1:p])
  if (!is.positive.definite(B_hat)){B_hat = as.matrix(nearPD(B_hat)$mat)}
  
  
  
  beta_mu = array(0, dim=c(l,p))
  beta_cov = array(0, dim=c(p,p,l))
  for (l_id in 1:l){
    letter_data = as.matrix(unname(background_data[(background_data$character==l_id),1:p]))
    
    letter_diff = letter_data - matrix(mu_hat[col(letter_data)], ncol = p)
    beta_l = colMeans(letter_diff)
    beta_mu[l_id,] = beta_l
    
    S = matrix(0, nrow = p, ncol = p)
    n_writers_letter <- 0
    
    for (w in unique(background_data$writer_id)){
      letter_writer = background_data[
        (background_data$character==l_id) & 
          (background_data$writer_id==w), 1:p, drop = FALSE
      ]
      
      if (nrow(letter_writer)>2){
        a_data_writer = background_data[
          (background_data$character==1) & 
            (background_data$writer_id==w), 1:p, drop = FALSE
        ]
        
        if (nrow(a_data_writer)>2){
          n_writers_letter <- n_writers_letter + 1
          
          mu_hat_writer = matrix(colMeans(a_data_writer), nrow = 1)
          letter_diff_writer = as.matrix(letter_writer) - 
            matrix(mu_hat_writer[col(as.matrix(letter_writer))], ncol = p)
          
          beta_w = matrix(colMeans(letter_diff_writer), nrow = 1)
          S.this <- t(beta_w - beta_l) %*% (beta_w - beta_l)
          S <- S + S.this
        }
      }
    }
    
    if (n_writers_letter > 1){
      B_hat_l = S/(n_writers_letter - 1)
    } else {
      B_hat_l = diag(1e-6, p)
    }
    
    if (!is.positive.definite(B_hat_l)){B_hat_l = as.matrix(nearPD(B_hat_l)$mat)}
    beta_cov[,,l_id] = B_hat_l
  }
  
  
  Sw = 0
  n_rows <- 0
  n_writers <- 0
  for (w in unique(background_data$writer_id)){
    df_writer = background_data[(background_data$writer_id==w),]
    if (nrow(df_writer)>2){
      n_rows <- n_rows + nrow(df_writer)
      n_writers <- n_writers+1
      Cov.this = cov(df_writer[,1:p])*(nrow(df_writer)-1)
      Sw <- Sw + Cov.this
    }
  }
  W_hat <- Sw/(n_rows - n_writers)
  U_hat <- W_hat * (nw_hat - p  -1)
  
  
  log_sd_mat <- sapply(unique(background_data$writer_id), function(w) {
    df_w <- background_data[background_data$writer_id == w, 1:p]
    if (nrow(df_w) <= 3) return(rep(NA_real_, p))  # skip degenerate writers
    log(sqrt(diag(cov(as.matrix(df_w)))))
  })
  
  # Remove degenerate writers (columns with NA)
  log_sd_mat <- log_sd_mat[, colSums(is.na(log_sd_mat)) == 0]
  
  loc <- rowMeans(log_sd_mat)
  sc  <- apply(log_sd_mat, 1, sd)
  
  eta <- 1
  
  return(list(mu_hat,B_hat,beta_mu,beta_cov,nw_hat,U_hat,loc,sc,eta))
}





stan_model_niw <- stan_model(file = "Stan_models/niw.stan", model_name = "niw")
stan_model_nlkj <- stan_model(file = "Stan_models/normal_lkj_model.stan", model_name = "normal_lkj_model")
stan_model_manova_iw <- stan_model(file = "Stan_models/MANOVA_iw_model.stan", model_name = "MANOVA_iw")
stan_model_manova_lkj <- stan_model(file = "Stan_models/MANOVA_lkj_model.stan", model_name = "MANOVA_lkj")


write_xlsx(data.frame(),"Paper_experiments/models_ml_per_writer_iter.xlsx")

ml_per_writer_def <- function(character_data,w){
  df_all=data.frame()
  
  writer_data_all = character_data[(character_data$writer_id==w),]
  
  sample_size <- min(100, nrow(writer_data_all))
  
  writer_data <- writer_data_all %>%
    add_count(character, name = "char_freq") %>%  # add frequency column
    slice_sample(
      n = sample_size,       # now it's a constant
      weight_by = char_freq, # weighted sampling
      replace = FALSE
    )
  
  int_characters = unique(writer_data$character)
  l = length(int_characters)
  
  alphabet_map <- setNames(seq_along(int_characters), int_characters)
  
  writer_data$character <- as.numeric(alphabet_map[writer_data$character])
  
  background_data = character_data[((character_data$writer_id!=w) & (character_data$character %in% int_characters)),]
  
  background_data$character <- as.numeric(alphabet_map[background_data$character])
  
  
  background_stats_niw = background_statistics_niw(background_data)
  background_stats_br = background_statistics_br(background_data)
  
  p = nrow(background_stats_niw[[2]])
  
  stan_data_normal <- list(N = nrow(writer_data), 
                       P = p, 
                       y = unname(as.matrix(writer_data[,1:p])), 
                       mu = as.vector(background_stats_niw[[1]]), 
                       B = background_stats_niw[[2]],
                       U = background_stats_niw[[4]],
                       nu = background_stats_niw[[3]],
                       loc = background_stats_niw[[5]],
                       sc = background_stats_niw[[6]],
                       eta = background_stats_niw[[7]]
  )

  
  
  beta_cov <-background_stats_br[[4]]
  beta_cov_list <- vector("list", dim(beta_cov)[3])
  for (i in 1:dim(beta_cov)[3]) {
    beta_cov_list[[i]] <- beta_cov[, , i]
  }
  
  stan_data_manova <- list(N = nrow(writer_data), 
                       P = p, 
                       L = l,
                       letters = writer_data$character,
                       y = unname(as.matrix(writer_data[,1:p])),
                       mu = as.vector(background_stats_br[[1]]), 
                       B = as.matrix(background_stats_br[[2]]), 
                       beta_mu=background_stats_br[[3]], 
                       beta_cov=beta_cov_list, 
                       U = background_stats_br[[6]],
                       nu = background_stats_br[[5]],
                       loc = background_stats_br[[7]], 
                       sc = background_stats_br[[8]], 
                       eta = background_stats_br[[9]])
  
  
    
  loglik_niw_conjugate = marginal_likelihood_niw_conjugate(stan_data_normal)
    
  niw = sampling(stan_model_niw, data = stan_data_normal, iter = 2000, chains = 1, cores=1, refresh = 0)
  loglik_niw <- bridge_sampler(niw, method = 'warp3', silent = TRUE)$logml
    
  nlkj = sampling(stan_model_nlkj, data = stan_data_normal, iter = 2000, chains = 1, cores=1, refresh = 0)
  loglik_nlkj <- bridge_sampler(nlkj, method = 'warp3', silent = TRUE)$logml
  
    
    
  loglik_manova_conjugate = marginal_likelihood_manova_conjugate(writer_data, stan_data_manova)
    
  manova_iw = sampling(stan_model_manova_iw, data = stan_data_manova, iter = 2000, chains = 1, cores=1, refresh = 0)
  loglik_manova_iw <- bridge_sampler(manova_iw, method = 'warp3', silent = TRUE)$logml
    
  manova_lkj = sampling(stan_model_manova_lkj, data = stan_data_manova, iter = 2000, chains = 1, cores=1, refresh = 0)
  loglik_manova_lkj <- bridge_sampler(manova_lkj, method = 'warp3', silent = TRUE)$logml
    
    
  df_new = cbind(data.frame(writer=w,
                            loglik_niw_conjugate,loglik_niw,loglik_nlkj,
                            loglik_manova_conjugate,loglik_manova_iw,loglik_manova_lkj
                            ))
    
  print(paste0("Writer:",w))
    
  df_all = rbind(df_all,df_new)

  ssr_i <- read_excel("Paper_experiments/models_ml_per_writer_iter.xlsx")
  ssr_i = rbind(ssr_i,df_all)
  write_xlsx(ssr_i,"Paper_experiments/models_ml_per_writer_iter.xlsx")
  return(df_all)
}


#ml_per_writer_def(db,'123')

unique(db$writer_id)

detectCores()
cl <- makeCluster(12,
                  outfile="C:/Users/lampis/Desktop/PhD/Handwritten_Loop_characters/handwriting-evidence-evaluation-fourier-based-feature-extraction/Modelling/Paper_experiments/log.txt")

clusterEvalQ(cl, {
  library(dplyr)
  library(Matrix)
})

clusterExport(cl,
              list("background_statistics_br","background_statistics_niw","abind",
                   "is.positive.definite","nearPD","fitdistr","lmvgamma", "inv",
                   "stan_model_niw", "stan_model_nlkj",
                   "stan_model_manova_iw","stan_model_manova_lkj",
                   "assess_BF","sampling","extract","bridge_sampler",
                   "marginal_likelihood_niw_conjugate", "niw_conjugate",
                   "marginal_likelihood_manova_conjugate","MANOVA_conjugate",
                   "normal_iW","normal_lkj","MANOVA_iw","MANOVA_LKJ",
                   "read_excel","write_xlsx"),
              envir=globalenv())

w.list <- sapply(unique(db$writer_id), list)

system.time({saves = parLapply(cl, w.list,
                               ml_per_writer_def,
                               character_data = db)})

stopCluster(cl)

df_all <- do.call("rbind", saves)

write_xlsx(df_all,"Paper_experiments/models_ml_per_writer.xlsx")


ssr <- read_excel("Paper_experiments/models_ml_per_writer.xlsx")


ssr = as.data.frame(ssr)

log_BF_manova_vs_normal <- data.frame(Writer = ssr$writer) 

log_BF_manova_vs_normal['BF_4_1'] = ssr$loglik_manova_conjugate-ssr$loglik_niw_conjugate

log_BF_manova_vs_normal['BF_5_2'] = ssr$loglik_manova_iw-ssr$loglik_niw

log_BF_manova_vs_normal['BF_6_3'] = ssr$loglik_manova_lkj-ssr$loglik_nlkj

log_BF_manova_vs_normal

colMeans(log_BF_manova_vs_normal>0, na.rm=T)


log_BF_conjugate_vs_non_iW <- data.frame(Writer = ssr$writer) 

log_BF_conjugate_vs_non_iW['BF_1_2'] = ssr$loglik_niw_conjugate-ssr$loglik_niw

log_BF_conjugate_vs_non_iW['BF_4_5'] = ssr$loglik_manova_conjugate-ssr$loglik_manova_iw


log_BF_conjugate_vs_non_iW

colMeans(log_BF_conjugate_vs_non_iW>0, na.rm=T)



log_BF_conjugate_vs_lkj <- data.frame(Writer = ssr$writer) 

log_BF_conjugate_vs_lkj['BF_1_3'] = ssr$loglik_niw_conjugate-ssr$loglik_nlkj

log_BF_conjugate_vs_lkj['BF_4_6'] = ssr$loglik_manova_conjugate-ssr$loglik_manova_lkj


log_BF_conjugate_vs_lkj
colMeans(log_BF_conjugate_vs_lkj>0, na.rm=T)



log_BF_iw_vs_lkj <- data.frame(Writer = ssr$writer) 

log_BF_iw_vs_lkj['BF_2_3'] = ssr$loglik_niw-ssr$loglik_nlkj

log_BF_iw_vs_lkj['BF_5_6'] = ssr$loglik_manova_iw-ssr$loglik_manova_lkj

log_BF_iw_vs_lkj
colMeans(log_BF_iw_vs_lkj>0, na.rm=T)

# 
# 
# log_BF_normal_all <- data.frame(Writer = ssr$writer) 
# 
# log_BF_normal_all['BF_1_2'] = ssr$loglik_niw_conjugate-ssr$loglik_niw
# 
# log_BF_normal_all['BF_1_3'] = ssr$loglik_niw_conjugate-ssr$loglik_nlkj
# 
# log_BF_normal_all['BF_2_3'] = ssr$loglik_niw-ssr$loglik_nlkj
# 
# log_BF_normal_all
# 
# 
# 
# log_BF_manova_all <- data.frame(Writer = ssr$writer) 
# 
# log_BF_manova_all['BF_4_5'] = ssr$loglik_manova_conjugate-ssr$loglik_manova_iw
# 
# log_BF_manova_all['BF_4_6'] = ssr$loglik_manova_conjugate-ssr$loglik_manova_lkj
# 
# log_BF_manova_all['BF_5_6'] = ssr$loglik_manova_iw-ssr$loglik_manova_lkj
# 
# log_BF_manova_all
