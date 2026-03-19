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

p <- 9

IAM_data = cbind(scale(IAM_data[,1:p]),IAM_data[,(p+1):ncol(IAM_data)])

writers_ids <- unique(IAM_data$writer_id)

#same person

writer_data = IAM_data[(IAM_data$writer_id==writers_ids[1]),]
background_data = IAM_data[(IAM_data$writer_id!=writers_ids[1]),]

int_characters = sort(intersect(writer_data$character,background_data$character))
l = length(int_characters)

questioned_data = data.frame()
suspect_data = data.frame()
for (c in int_characters){
  writer_data_c = writer_data[(writer_data$character==c),]
  random_percentage = runif(1,0.35,0.65)
  smp_size <- floor(random_percentage * nrow(writer_data_c))
  suspect_ind <- sample(seq_len(nrow(writer_data_c)), size = smp_size)
  
  questioned_data = rbind(questioned_data,writer_data_c[suspect_ind, ])
  suspect_data = rbind(suspect_data, writer_data_c[-suspect_ind, ])
}

writer_data = rbind(suspect_data,questioned_data)


# different writers
writer_data_1 <- IAM_data[IAM_data$writer_id == "60", ]

writer_data_2 <- IAM_data[IAM_data$writer_id == "90", ]

background_data <- IAM_data[!(IAM_data$writer_id %in% c("60", "90")), ]

# Function to find low-frequency characters
get_small_groups <- function(data, threshold) {
  tab <- table(data$character)
  names(tab[tab < threshold])
}

# Find rare characters in each dataset
small_1 <- get_small_groups(writer_data_1, p)
small_2 <- get_small_groups(writer_data_2, p)
small_bg <- get_small_groups(background_data, p)

# Determine reference based on who has the most rare characters
rare_counts <- c(length(small_1), length(small_2), length(small_bg))
ref_index <- which.max(rare_counts)
reference_rare_chars <- list(small_1, small_2, small_bg)[[ref_index]]

# Replace rare characters with "OTHERS"
replace_rare <- function(data, rare_chars) {
  data$character <- ifelse(data$character %in% rare_chars, "z-others", data$character)
  return(data)
}

# Apply to all datasets
# writer_data_1 <- replace_rare(writer_data_1, reference_rare_chars)
# writer_data_2 <- replace_rare(writer_data_2, reference_rare_chars)
# background_data <- replace_rare(background_data, reference_rare_chars)

# writer_data_1 <- writer_data_1[!(writer_data_1$character %in% reference_rare_chars),]
# writer_data_2 <- writer_data_2[!(writer_data_2$character %in% reference_rare_chars),]
# background_data <- background_data[!(background_data$character %in% reference_rare_chars),]

# Optional: restrict all datasets to shared characters only
int_characters <- sort(intersect(background_data$character,
                                 intersect(writer_data_1$character, writer_data_2$character)))

l = length(int_characters)

questioned_data = data.frame()
suspect_data = data.frame()
for (c in int_characters){
  writer_data_1_c = writer_data_1[(writer_data_1$character == c),]
  if (nrow(writer_data_1_c)<5){
    questioned_data = rbind(questioned_data,writer_data_1_c)
  } else{
    random_percentage1 = runif(1,0.35,0.65)
    smp_size <- round(random_percentage1 * nrow(writer_data_1_c))
    ind <- sample(seq_len(nrow(writer_data_1_c)), size = smp_size)
    questioned_data = rbind(questioned_data,writer_data_1_c[ind, ])
  }
  
  writer_data_2_c = writer_data_2[(writer_data_2$character==c),]
  if (nrow(writer_data_2_c)<5){
    suspect_data = rbind(suspect_data,writer_data_2_c)
  } else{
    smp_size <- round((1-random_percentage1) * nrow(writer_data_2_c))
    ind <- sample(seq_len(nrow(writer_data_2_c)), size = smp_size)
    suspect_data = rbind(suspect_data,writer_data_2_c[ind, ])
  }
  print(dim(writer_data_2_c[ind, ]))
}

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

eta=4

#fit <- fitdistr(all_diagonals, "cauchy")
#print(fit)
#hist(rcauchy(1000, location = fit$estimate[1], scale = fit$estimate[2]))

loc <- mean(log(diag(W_hat)))
sc <- sd(log(diag(W_hat)))




# Create mapping: letter -> number based on alphabetical order

alphabet_map <- setNames(seq_along(int_characters), int_characters)

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
  
  writer_x = unname(model.matrix(~factor(stan_data$letters)))
  
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

