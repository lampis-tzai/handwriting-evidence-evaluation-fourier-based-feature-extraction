library(readxl)
library(dplyr)
library(MASS)
library(Matrix)
library(matrixcalc)
library(CholWishart)


set.seed(2)

IAM_data <- read_excel("IAM_fourier_features_dataset/DB_loop_handwriting_ls.xlsx")
IAM_data = as.data.frame(IAM_data)

# IAM_data[,2:9] = IAM_data[,2:9]/sqrt(IAM_data$area)
# IAM_data[,1] = log(IAM_data[,1])

IAM_data = cbind(scale(IAM_data[,1:9]),IAM_data[,10:ncol(IAM_data)])

writers_ids <- unique(IAM_data$writer_id)

#same person

writer_data_all = IAM_data[(IAM_data$writer_id==writers_ids[1]),]
background_data = IAM_data[(IAM_data$writer_id!=writers_ids[1]),]

sample_size <- min(100, nrow(writer_data_all))

writer_data <- writer_data_all %>%
  add_count(character, name = "char_freq") %>%  # add frequency column
  slice_sample(
    n = sample_size,       # now it's a constant
    weight_by = char_freq, # weighted sampling
    replace = FALSE
  )

questioned_data <- writer_data[1:50,]
suspect_data <- writer_data[51:100,]

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

writer_data_2 <- IAM_data[IAM_data$writer_id == "92", ]

background_data <- IAM_data[!(IAM_data$writer_id %in% c("88", "92")), ]

sample_size <-  floor(nrow(writer_data_1)/2)#min(50, nrow(writer_data_1))

questioned_data <- writer_data_1 %>%
  add_count(character, name = "char_freq") %>%  # add frequency column
  slice_sample(
    n = sample_size,       # now it's a constant
    weight_by = char_freq, # weighted sampling
    replace = FALSE
  )

sample_size <-  floor(nrow(writer_data_2)/2)#min(50, nrow(writer_data_2))
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

library(Hotelling)
htest<- hotelling.test(x = questioned_data[,1:9], y = suspect_data[,1:9])
htest

# Hyperparameter Elicitation

p=9
nw.min = p + 2


# 
# function_NR <- function(data,df){
#   sum_data <- Reduce('+', data)
#   M = length(data)
#   V <- sum_data/(df*M)
#   model_data <- abind(data, along = 3)
#   lh = sum(dWishart(model_data,df,V,log = T))
#   return(-lh)
# }
# 
# model_data_list <- list()
# i=1
# for (w in unique(background_data$writer_id)){
#   df_writer = background_data[(background_data$writer_id==w),]
#   var_data = unname(as.matrix(df_writer[,1:p]))
#   model_data_list[[i]] = cov(var_data)
#   i=i+1
# }
# 
# nlm(function_NR,10,data=model_data_list)

nw_hat = nw.min

mu_hat=matrix(colMeans(do.call(rbind, lapply(unique(background_data$writer_id), function(w)
       colMeans(background_data[background_data$writer_id == w, 1:p])))), nrow = 1)

S = 0
Sw = 0
for (w in unique(background_data$writer_id)){
  df_writer = background_data[(background_data$writer_id==w),]
  var_data = unname(as.matrix(df_writer[,1:p]))
  theta_w = matrix(colMeans(var_data), nrow = 1)
  S.this <- (t(theta_w - mu_hat) %*% (theta_w - mu_hat))
  S <- S + S.this
  Cov.this = cov(var_data)*(nrow(df_writer)-1) 
  Sw <- Sw + Cov.this
} 

B_hat = S/(length(unique(background_data$writer_id)) - 1)
#B_hat = cov(background_data[,1:p])
if (!is.positive.definite(B_hat)){B_hat = as.matrix(nearPD(B_hat)$mat)}


W_hat <- Sw/(nrow(background_data) - length(unique(background_data$writer_id)))
U_hat <- W_hat*(nw_hat-p-1)


# sd_w <- sqrt(diag(W_hat))
# loc <- mean(log(sd_w))
# sc <- sd(log(sd_w))
log_sd_mat <- sapply(unique(background_data$writer_id), function(w) {
  df_w <- background_data[background_data$writer_id == w, 1:p]
  if (nrow(df_w) <= p) return(rep(NA_real_, p))  # skip degenerate writers
  log(sqrt(diag(cov(as.matrix(df_w)))))
})

# Remove degenerate writers (columns with NA)
log_sd_mat <- log_sd_mat[, colSums(is.na(log_sd_mat)) == 0]

loc <- rowMeans(log_sd_mat)
sc  <- apply(log_sd_mat, 1, sd)

# loc <- mean(loc)
# sc  <- mean(sc)


# log_sds <- do.call(c, lapply(unique(background_data$writer_id), function(w) {
#   df_w <- background_data[background_data$writer_id == w, 1:p]
#   log(apply(df_w, 2, sd))  # log-SD per feature per writer
# }))

# log_sds<-log(apply(background_data[,1:p],2,sd))
# 
# loc <- mean(log_sds)
# sc  <- sd(log_sds)


# per writer covariance and correlation matrices
# writers <- unique(background_data$writer_id)
# Rs <- lapply(writers, function(w) {
#   df_w <- as.matrix(background_data[background_data$writer_id == w, 1:p])
#   S_w <- cov(df_w)
#   D_half <- diag(1 / sqrt(diag(S_w)))
#   R_w <- D_half %*% S_w %*% D_half
#   R_w
# })
# 
# logdet_R <- vapply(Rs, function(R) determinant(R, logarithm = TRUE)$modulus, numeric(1))
# 
# # log normalizing constant for LKJ_K(eta)
# logZ_lkj <- function(eta, K) {
#   # Using Lewandowski–Kurowicka–Joe form: product of Beta functions
#   # See e.g. Stan / TFP implementation for exact formula
#   # Skeleton (you need to fill-in carefully from a reference):
#   s <- 0
#   for (i in 1:(K - 1)) {
#     a <- 0.5 * (K - i + 1)
#     b <- eta + 0.5 * (K - i - 1)
#     s <- s + lbeta(a, b)
#   }
#   -s
# }
# 
# loglik_eta <- function(eta, logdet_R, K) {
#   if (eta <= 0) return(-Inf)
#   M <- length(logdet_R)
#   (eta - 1) * sum(logdet_R) + M * logZ_lkj(eta, K)
# }
# 
# # Maximize numerically (e.g. nlm or optim with 1D search)
# eta_hat <- optimize(
#   f = function(e) -loglik_eta(e, logdet_R, K = p),
#   interval = c(1e-3, 100),
#   maximum = FALSE
# )$minimum
# eta_hat

eta <- 1



# Stan data list
stan_data_H0 <- list(N = nrow(writer_data), 
                     P = p, 
                     y = unname(as.matrix(writer_data[,1:p])), 
                     mu = as.vector(mu_hat), 
                     B = as.matrix(B_hat),
                     loc = loc,
                     sc = sc,
                     eta = eta,
                     U = U_hat,
                     nu = nw_hat
)

stan_data_H1_1 <- list(N = nrow(questioned_data), 
                       P = p, 
                       y = unname(as.matrix(questioned_data[,1:p])), 
                       mu = as.vector(mu_hat), 
                       B = as.matrix(B_hat),
                       loc = loc,
                       sc = sc,
                       eta = eta,
                       U = U_hat,
                       nu = nw_hat
)

stan_data_H1_2 <- list(N = nrow(suspect_data), 
                       P = p, 
                       y = unname(as.matrix(suspect_data[,1:p])), 
                       mu = as.vector(mu_hat), 
                       B = as.matrix(B_hat),
                       loc = loc,
                       sc = sc,
                       eta = eta,
                       U = U_hat,
                       nu = nw_hat
)



# Fit the model
library(rstan)
library(bridgesampling)
stan_model_niw <- stan_model(file = "Stan_models/niw.stan", model_name = "niw")
stan_model_nlkj <- stan_model(file = "Stan_models/normal_lkj_model.stan", model_name = "normal_lkj_model")

assess_BF <- function(stan_model,stan_data_H0,stan_data_H1_1,stan_data_H1_2){
  fit_H0 <- sampling(stan_model, data = stan_data_H0, iter = 2000, chains =1, cores=1)
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


BF_niw<- assess_BF(stan_model_niw,stan_data_H0,stan_data_H1_1,stan_data_H1_2)
BF_nlkj<- assess_BF(stan_model_nlkj,stan_data_H0,stan_data_H1_1,stan_data_H1_2)


marginal_likelihood_niw_conjugate<- function(stan_data){
  w_data = stan_data$y
  m0 = matrix(stan_data$mu,nrow = 1, byrow = TRUE)
  v0 = stan_data$nu
  k0 = 0.5
  U0 = stan_data$U
  n = nrow(w_data)
  d = ncol(w_data)
  κn = k0 + n
  vn = v0 + n
  S = cov(w_data)*(nrow(w_data)-1)
  Wn = U0 + S + (k0*n)/(k0+n)*(t((colMeans(w_data)-m0))%*%(colMeans(w_data)-m0))

  logml = - ((d*n)/2)*log(pi) + lmvgamma(vn/2, d) - lmvgamma(v0/2, d) +
    (v0/2)*determinant(U0,logarithm = TRUE)$modulus[1] -
    (vn/2)*determinant(Wn,logarithm = TRUE)$modulus[1]+
    (d/2)*(log(k0)-log(κn))

  return(logml)

}

log_lik_H0 <- marginal_likelihood_niw_conjugate(stan_data_H0)
log_lik_H1_1 <- marginal_likelihood_niw_conjugate(stan_data_H1_1)
log_lik_H1_2 <- marginal_likelihood_niw_conjugate(stan_data_H1_2)


BF_niw_conjugate <- log_lik_H0-log_lik_H1_1-log_lik_H1_2


print(BF_niw_conjugate)
print(BF_niw)
print(BF_nlkj)

