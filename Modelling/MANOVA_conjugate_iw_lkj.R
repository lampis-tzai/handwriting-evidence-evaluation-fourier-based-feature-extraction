setwd("C:/Users/ltzai/Desktop/PhD/Handwritten_Loop_characters/handwriting-evidence-evaluation-fourier-based-feature-extraction/Modelling")
library(readxl)
library(dplyr)
library(MASS)
library(matrixcalc)
library(Matrix)
library(matlib)
library(CholWishart)

set.seed(2)

#writers ids: "93"  "62"  "87"  "118" "332" "151" "123" "26"  "112" "334" "60"  "90"  "92"  "61"  "37"  "111" "109" "113" "114" "25" 

IAM_data <- read_excel("IAM_fourier_features_dataset/DB_loop_handwriting.xlsx")
IAM_data = as.data.frame(IAM_data)

IAM_data[,2:9] = IAM_data[,2:9]/sqrt(IAM_data$area)

IAM_data = cbind(scale(IAM_data[,1:9]),IAM_data[,10:ncol(IAM_data)])

writers_ids <- unique(IAM_data$writer_id)


#same person

writer_data_all = IAM_data[(IAM_data$writer_id==writers_ids[1]),]
background_data = IAM_data[(IAM_data$writer_id!=writers_ids[1]),]

sample_size <- min(200, nrow(writer_data_all))

writer_data <- writer_data_all %>%
  add_count(character, name = "char_freq") %>%  # add frequency column
  slice_sample(
    n = sample_size,       # now it's a constant
    weight_by = char_freq, # weighted sampling
    replace = FALSE
  )

questioned_data <- writer_data[1:100,]
suspect_data <- writer_data[101:200,]

# intersect characters
#int_characters <- sort(intersect(questioned_data$character,suspect_data$character))

#questioned_data$character <- questioned_data[questioned_data$character %in% int_characters, ]
#suspect_data$character <- suspect_data[suspect_data$character %in% int_characters, ]

#alphabet_map <- setNames(seq_along(int_characters), int_characters) 

#background_data <- background_data_all[background_data_all$character %in% int_characters,]

int_characters <- sort(unique(IAM_data$character))
l <- length(int_characters)

questioned_data$character <- factor(questioned_data$character, levels = int_characters)
suspect_data$character <- factor(suspect_data$character, levels = int_characters)

table(questioned_data$character)
table(suspect_data$character)

alphabet_map <- setNames(seq_along(int_characters), int_characters) 


# different writers
writer_data_1 <- IAM_data[IAM_data$writer_id == "88", ]

writer_data_2 <- IAM_data[IAM_data$writer_id == "152", ]

background_data <- IAM_data[!(IAM_data$writer_id %in% c("88", "152")), ]

sample_size <- min(100, nrow(writer_data_1))

questioned_data <- writer_data_1 %>%
  add_count(character, name = "char_freq") %>%  # add frequency column
  slice_sample(
    n = sample_size,       # now it's a constant
    weight_by = char_freq, # weighted sampling
    replace = FALSE
  )

sample_size <- min(100, nrow(writer_data_2))
suspect_data <- writer_data_2 %>%
  add_count(character, name = "char_freq") %>%  # add frequency column
  slice_sample(
    n = sample_size,       # now it's a constant
    weight_by = char_freq, # weighted sampling
    replace = FALSE
  )

#questioned_data_old = questioned_data
#suspect_data_old = suspect_data

# questioned_data = questioned_data_old
# suspect_data = suspect_data_old

int_characters <- sort(unique(IAM_data$character))
l <- length(int_characters)

questioned_data$character <- factor(questioned_data$character, levels = int_characters)
suspect_data$character <- factor(suspect_data$character, levels = int_characters)

alphabet_map <- setNames(seq_along(int_characters), int_characters)

# int_characters <- sort(intersect(questioned_data$character,suspect_data$character))
# 
# l <- length(int_characters)
# 
# alphabet_map <- setNames(seq_along(int_characters), int_characters)
# 
# 
# questioned_data$character <- ifelse(questioned_data$character %in% int_characters,
#                                    questioned_data$character,
#                                    'o')
# 
# suspect_data$character <- ifelse(suspect_data$character %in% int_characters,
#                                    suspect_data$character,
#                                    'o')
# 
# background_data$character <- ifelse(background_data$character %in% int_characters,
#                                    background_data$character,
#                                 'o')

table(questioned_data$character)
table(suspect_data$character)


writer_data = rbind(suspect_data,questioned_data)

# Hyperparameter Elicitation

p=9
nw.min = p + 2
nw_hat = nw.min

a_data = background_data[(background_data$character=='a'),1:p]
mu_hat=matrix(colMeans(a_data),nrow = 1)

S = 0
for (w in unique(background_data$writer_id)){
  df_writer = background_data[(
    background_data$character=='a')& (background_data$writer_id==w),]

  theta_w = matrix(colMeans(df_writer[,1:p]), nrow = 1)
  S.this <- (t(theta_w - mu_hat) %*% (theta_w - mu_hat))
  S <- S + S.this
}

B_hat = S/(length(unique(background_data$writer_id)) - 1)
#B_hat = cov(a_data)

if (!is.positive.definite(B_hat)){B_hat = as.matrix(nearPD(B_hat)$mat)}


beta_mu = array(0, dim=c(l,p))
beta_cov = array(0, dim=c(p,p,l))
for (l_id in 1:l){
  letter_data = as.matrix(unname(background_data[(background_data$character==int_characters[l_id]),1:p]))
  
  letter_diff = letter_data - matrix(mu_hat[col(letter_data)],ncol = p)
  beta_l = colMeans(letter_diff)
  beta_mu[l_id,] = beta_l
  S = 0
  for (w in unique(background_data$writer_id)){
    letter_writer = background_data[(
      background_data$character==int_characters[l_id])& (background_data$writer_id==w),1:p]
    if (nrow(letter_writer)>3){
      a_data_writer = background_data[(
        background_data$character=='a')& (background_data$writer_id==w),1:p]

      mu_hat_writer=matrix(colMeans(a_data_writer),nrow = 1)

      letter_diff_writer = letter_writer - matrix(mu_hat_writer[col(letter_writer)],ncol = p)

      beta_w = matrix(colMeans(letter_diff_writer), nrow = 1)
      S.this <- (t(beta_w - beta_l) %*% (beta_w - beta_l))
      S <- S + S.this
    }
  }
  B_hat_l = S/(length(unique(background_data$writer_id)) - 1)
  #B_hat_l = cov(letter_data)
  if (!is.positive.definite(B_hat_l)){B_hat_l = as.matrix(nearPD(B_hat_l)$mat)}
  beta_cov[,,l_id] = B_hat_l
}


Sw = 0
for (w in unique(background_data$writer_id)){
  df_writer = background_data[(background_data$writer_id==w),]
  Cov.this = cov(df_writer[,1:p])*(nrow(df_writer)-1)
  Sw <- Sw + Cov.this
}
W_hat <- Sw/(nrow(background_data) - length(unique(background_data$writer_id)))
U_hat <- W_hat * (nw_hat - p  -1)


beta_cov_list <- vector("list", dim(beta_cov)[3])
for (i in 1:dim(beta_cov)[3]) {
  beta_cov_list[[i]] <- beta_cov[, , i]
}

eta=1

#fit <- fitdistr(all_diagonals, "cauchy")
#print(fit)
#hist(rcauchy(1000, location = fit$estimate[1], scale = fit$estimate[2]))

loc <- mean(log(diag(W_hat)))
sc <- sd(log(diag(W_hat)))




# Create mapping: letter -> number based on alphabetical order



# Stan data list
stan_data_H0 <- list(N = nrow(writer_data), 
                     P = p, 
                     L = l,
                     mu = as.vector(mu_hat), 
                     B = as.matrix(B_hat),
                     beta_mu=beta_mu,
                     beta_cov=beta_cov_list,
                     letters = as.numeric(alphabet_map[writer_data$character]),
                     y = unname(as.matrix(writer_data[,1:p])),
                     loc =loc,
                     sc = sc,
                     eta = eta,
                     U = U_hat,
                     nu = nw_hat)

stan_data_H1_1 <- list(N = nrow(questioned_data), 
                       P = p, 
                       L = l,
                       mu = as.vector(mu_hat), 
                       B = as.matrix(B_hat),
                       beta_mu=beta_mu,
                       beta_cov=beta_cov_list,
                       letters = as.numeric(alphabet_map[questioned_data$character]),
                       y = unname(as.matrix(questioned_data[,1:p])),
                       loc =loc,
                       sc = sc,
                       eta = eta,
                       U = U_hat,
                       nu = nw_hat)

stan_data_H1_2 <- list(N = nrow(suspect_data), 
                       P = p, 
                       L = l,
                       mu = as.vector(mu_hat), 
                       B = as.matrix(B_hat),
                       beta_mu=beta_mu,
                       beta_cov=beta_cov_list,
                       letters =  as.numeric(alphabet_map[suspect_data$character]),
                       y = unname(as.matrix(suspect_data[,1:p])),
                       loc =loc,
                       sc = sc,
                       eta = eta,
                       U = U_hat,
                       nu = nw_hat)

library(rstan)
library(bridgesampling)
stan_model_manova_iw <- stan_model(file = "Stan_Models/MANOVA_iw_model.stan", model_name = "MANOVA_iw")
stan_model__manova_lkj <- stan_model(file = "Stan_Models/MANOVA_lkj_model.stan", model_name = "MANOVA_lkj")

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



marginal_likelihood_manova_conjugate<- function(stan_data){
  

  l <- stan_data$L
  m0 = stan_data$beta_mu
  m0[1,] = stan_data$mu
  v0 = stan_data$nu
  k0 = diag(0.5,l,l)
  U0 = stan_data$U
  
  writer_x = unname(model.matrix(~factor(stan_data$letters,levels = 1:l)))
  
  writer_y = as.matrix(stan_data$y)
  
  n = nrow(writer_y)
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

log_lik_H0 <- marginal_likelihood_manova_conjugate(stan_data_H0)
log_lik_H1_1 <- marginal_likelihood_manova_conjugate(stan_data_H1_1)
log_lik_H1_2 <- marginal_likelihood_manova_conjugate(stan_data_H1_2)


BF_manova_conjugate <- log_lik_H0-log_lik_H1_1-log_lik_H1_2


print(BF_manova_conjugate)
print(BF_manova_iw)
print(BF_manova_lkj)

