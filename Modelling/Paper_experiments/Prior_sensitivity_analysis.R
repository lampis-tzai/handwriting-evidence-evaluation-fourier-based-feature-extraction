setwd("C:/Users/ltzai/Desktop/PhD/Handwritten_Loop_characters/Handwriting_Multivariate_Approach")
library(readxl)
library(dplyr)
library(MASS)
library(matrixcalc)
library(Matrix)
library(matlib)

adoq_data <- read_excel("Data/adoq colonnes.xls")
adoq_data = as.data.frame(adoq_data)

names(adoq_data)

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

writer_data = adoq_data[(adoq_data$N==1),]
background_data = adoq_data[(adoq_data$N!=1),]

questioned_data = data.frame()
suspect_data = data.frame()
for (c in 1:4){
  writer_data_c = writer_data[(writer_data$Lettre==c),]
  random_percentage = runif(1,0.35,0.65)
  smp_size <- floor(random_percentage * nrow(writer_data_c))
  suspect_ind <- sample(seq_len(nrow(writer_data_c)), size = smp_size)
  
  questioned_data = rbind(questioned_data,writer_data_c[suspect_ind, ])
  suspect_data = rbind(suspect_data, writer_data_c[-suspect_ind, ])
}

writer_data = rbind(suspect_data,questioned_data)

writer_data_1 = adoq_data[(adoq_data$N == 7),]

writer_data_2 = adoq_data[(adoq_data$N == 8),]

background_data_all = adoq_data[!(adoq_data$N %in% c(7,8)),]

#for (i in 1:100){
questioned_data = data.frame()
suspect_data = data.frame()
for (c in 1:4){
  writer_data_1_c = writer_data_1[(writer_data_1$Lettre==c),]
  random_percentage1 = runif(1,0.35,0.65)
  smp_size <- round(random_percentage1 * nrow(writer_data_1_c))
  ind <- sample(seq_len(nrow(writer_data_1_c)), size = smp_size)
  questioned_data = rbind(questioned_data,writer_data_1_c[ind, ])
  
  writer_data_2_c = writer_data_2[(writer_data_2$Lettre==c),]
  smp_size <- round((1-random_percentage1) * nrow(writer_data_2_c))
  ind <- sample(seq_len(nrow(writer_data_2_c)), size = smp_size)
  print(dim(writer_data_2_c[ind, ]))
  suspect_data = rbind(suspect_data,writer_data_2_c[ind, ])
}

writer_data = rbind(suspect_data,questioned_data)


n_subsample <- length(unique(background_data_all$Writer))#ceiling(0.5 * length(unique(background_data_all$Writer)))

# Step 1: Subsample writers
sampled_writers <- sample(unique(background_data_all$Writer), size = n_subsample, replace = FALSE)
  
# Step 2: Bootstrap within each writer
background_data <- do.call(rbind, lapply(sampled_writers, function(w) {
  writer_data_db <- background_data_all[background_data_all$Writer == w, ]
  n_subsample <- ceiling(0.5 * nrow(writer_data_db))
  writer_data_db[sample(nrow(writer_data_db), size =n_subsample, replace = TRUE), ]
}))

#background_data<-background_data_all

p=9
l = length(unique(writer_data$Lettre))
nw.min = p + 2
nw_hat = nw.min

a_data = background_data[(background_data$Lettre==1),5:ncol(background_data)]
mu_hat=matrix(colMeans(a_data),nrow = 1)

# S = 0
# for (w in unique(background_data$N)){
#   df_writer = background_data[(
#     background_data$Lettre==1)& (background_data$N==w),]
#   
#   theta_w = matrix(colMeans(df_writer[,5:ncol(df_writer)]), nrow = 1)
#   S.this <- (t(theta_w - mu_hat) %*% (theta_w - mu_hat))
#   S <- S + S.this
# } 
# 
# B_hat = S/(length(unique(background_data$N)) - 1)
B_hat = cov(a_data)

if (!is.positive.definite(B_hat)){B_hat = as.matrix(nearPD(B_hat)$mat)}

library(Matrix)

beta_mu = array(0, dim=c(l,p))
beta_cov = array(0, dim=c(p,p,l))
for (l_id in 1:l){
  letter_data = as.matrix(unname(background_data[(
    background_data$Lettre==l_id),5:ncol(background_data)]))
  
  letter_diff = letter_data - matrix(mu_hat[col(letter_data)],ncol = p)
  beta_l = colMeans(letter_diff)
  beta_mu[l_id,] = beta_l
  
  # S = 0
  # for (w in unique(background_data$N)){
  #   
  #   letter_writer = background_data[(
  #     background_data$Lettre==l_id)& (background_data$N==w),
  #     5:ncol(background_data)]
  #   
  #   a_data_writer = background_data[(
  #     background_data$Lettre==1)& (background_data$N==w),
  #     5:ncol(background_data)]
  #   
  #   mu_hat_writer=matrix(colMeans(a_data_writer),nrow = 1)
  #   
  #   letter_diff_writer = letter_writer - matrix(mu_hat_writer[col(letter_writer)],ncol = p)
  #   
  #   beta_w = matrix(colMeans(letter_diff_writer), nrow = 1)
  #   S.this <- (t(beta_w - beta_l) %*% (beta_w - beta_l))
  #   S <- S + S.this
  # }
  # B_hat_l = S/(length(unique(background_data$N)) - 1)
  B_hat_l = cov(letter_data)
  if (!is.positive.definite(B_hat_l)){B_hat_l = as.matrix(nearPD(B_hat_l)$mat)}
  beta_cov[,,l_id] = B_hat_l
}


Sw = 0
for (w in unique(background_data$N)){
  df_writer = background_data[(background_data$N==w),]
  Cov.this = cov(df_writer[,5:ncol(df_writer)])*(nrow(df_writer)-1)
  Sw <- Sw + Cov.this
}
W_hat <- Sw/(nrow(background_data) - length(unique(background_data$N)))
U_hat <- W_hat * (nw_hat - p  -1)


beta_cov_list <- vector("list", dim(beta_cov)[3])
for (i in 1:dim(beta_cov)[3]) {
  beta_cov_list[[i]] <- beta_cov[, , i]
}



mle_eta=1

#fit <- fitdistr(all_diagonals, "cauchy")
#print(fit)
#hist(rcauchy(1000, location = fit$estimate[1], scale = fit$estimate[2]))

loc <- mean(log(diag(W_hat)))
sc <- sd(log(diag(W_hat)))

# Stan data list
stan_data_H0 <- list(N = nrow(writer_data), 
                     P = p, 
                     L = l,
                     letters = as.numeric(writer_data$Lettre),
                     y = unname(as.matrix(writer_data[,5:ncol(writer_data)])),
                     mu = as.vector(mu_hat),
                     B = as.matrix(B_hat),
                     beta_mu=beta_mu,
                     beta_cov=beta_cov_list,
                     loc = loc,
                     sc = sc,
                     eta = mle_eta,
                     U = U_hat,
                     nu = nw_hat)

# mu = rep(0,p),
# B = diag(1000,p),
# beta_mu=array(0,c(l,p)),
# beta_cov=rep(list(diag(1000,p)),4),
# loc = 0,
# sc = 5,
# eta = 1,
# U = diag(1,p),
# nu = nw_hat

stan_data_H1_1 <- list(N = nrow(questioned_data), 
                       P = p, 
                       L = l,
                       letters = as.numeric(questioned_data$Lettre),
                       y = unname(as.matrix(questioned_data[,5:ncol(questioned_data)])),
                       mu = as.vector(mu_hat),
                       B = as.matrix(B_hat),
                       beta_mu=beta_mu,
                       beta_cov=beta_cov_list,
                       loc = loc,
                       sc = sc,
                       eta = mle_eta,
                       U = U_hat,
                       nu = nw_hat)

stan_data_H1_2 <- list(N = nrow(suspect_data), 
                       P = p, 
                       L = l,
                       letters =  as.numeric(suspect_data$Lettre),
                       y = unname(as.matrix(suspect_data[,5:ncol(suspect_data)])),
                       mu = as.vector(mu_hat),
                       B = as.matrix(B_hat),
                       beta_mu=beta_mu,
                       beta_cov=beta_cov_list,
                       loc = loc,
                       sc = sc,
                       eta = mle_eta,
                       U = U_hat,
                       nu = nw_hat)

library(rstan)
library(bridgesampling)
stan_model_manova_iw <- stan_model(file = "Stan_code/MANOVA_iw_model.stan", model_name = "MANOVA_iw")
stan_model__manova_lkj <- stan_model(file = "Stan_code/MANOVA_lkj_model.stan", model_name = "MANOVA_lkj")

assess_BF <- function(stan_model,stan_data_H0,stan_data_H1_1,stan_data_H1_2){
  fit_H0 <- sampling(stan_model, data = stan_data_H0, iter = 2000, chains = 1, cores=1)
  fit_H1_1 <- sampling(stan_model, data = stan_data_H1_1, iter = 2000, chains = 1, cores=1)
  fit_H1_2 <- sampling(stan_model, data = stan_data_H1_2, iter = 2000, chains = 1, cores=1)
  
  samples_H0 <- extract(fit_H0)
  samples_H1_1 <- extract(fit_H1_1)
  samples_H1_2 <- extract(fit_H1_2)
  
  # Compute the log marginal likelihood
  bs_lik_H0 <- bridge_sampler(fit_H0, method = 'warp3')$logml
  bs_lik_H1_1 <- bridge_sampler(fit_H1_1, method = 'warp3')$logml
  bs_lik_H1_2 <- bridge_sampler(fit_H1_2, method = 'warp3')$logml
  
  return(bs_lik_H0 - bs_lik_H1_1 - bs_lik_H1_2)
  
}


BF_manova_iw<- assess_BF(stan_model_manova_iw,stan_data_H0,stan_data_H1_1,stan_data_H1_2)
BF_manova_lkj<- assess_BF(stan_model__manova_lkj,stan_data_H0,stan_data_H1_1,stan_data_H1_2)



marginal_likelihood_manova_conjugate<- function(w_data, stan_data){
  
  
  l <- stan_data$L
  m0 = stan_data$beta_mu
  m0[1,] = stan_data$mu
  v0 = stan_data$nu
  k0 = diag(0.5,l,l)
  U0 = stan_data$U
  
  
  
  writer_x = as.matrix(data.frame(a = rep(1,nrow(w_data)),
                                  d = ifelse(w_data$Lettre == 2, 1, 0),
                                  o = ifelse(w_data$Lettre == 3, 1, 0),
                                  q = ifelse(w_data$Lettre == 4, 1, 0)))
  
  writer_y = as.matrix(w_data[,5:ncol(w_data)])
  
  n = nrow(w_data)
  d = ncol(writer_y)
  vn = v0 + n 
  
  kn = t(writer_x) %*% writer_x + k0
  mn = inv(kn)%*%(t(writer_x)%*%writer_y + k0%*%m0)
  Un = U0 + t(writer_y) %*% writer_y + t(m0)%*%k0%*%m0 - t(mn)%*%kn%*%mn
  
  
  logml = - ((d*n)/2)*log(2*pi) + (d/2)*determinant(k0,logarithm = TRUE)$modulus[1] -
    (d/2)*determinant(kn,logarithm = TRUE)$modulus[1] + 
    (v0/2)*determinant(U0/2,logarithm = TRUE)$modulus[1] - 
    (vn/2)*determinant(Un/2,logarithm = TRUE)$modulus[1] +
    lmvgamma(vn/2, d) - lmvgamma(v0/2, d)
  
  return(logml)
  
}

library(CholWishart)

log_lik_H0 <- marginal_likelihood_manova_conjugate(writer_data, stan_data_H0)
log_lik_H1_1 <- marginal_likelihood_manova_conjugate(questioned_data, stan_data_H1_1)
log_lik_H1_2 <- marginal_likelihood_manova_conjugate(suspect_data, stan_data_H1_2)


BF_manova_conjugate <- log_lik_H0-log_lik_H1_1-log_lik_H1_2


print(BF_manova_conjugate)
print(BF_manova_iw)
print(BF_manova_lkj)
 
