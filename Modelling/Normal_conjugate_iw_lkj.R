setwd("C:/Users/ltzai/Desktop/PhD/Handwritten_Loop_characters/handwriting-evidence-evaluation-fourier-based-feature-extraction/Modelling")
library(readxl)
library(dplyr)
library(MASS)
library(Matrix)
library(matrixcalc)
library(CholWishart)


set.seed(2)

IAM_data <- read_excel("IAM_fourier_features_dataset/DB_loop_handwriting.xlsx")
IAM_data = as.data.frame(IAM_data)

IAM_data = cbind(scale(IAM_data[,1:9]),IAM_data[,10:ncol(IAM_data)])

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
writer_data_1 = IAM_data[(IAM_data$writer_id==writers_ids[1]),]

writer_data_2 = IAM_data[(IAM_data$writer_id==writers_ids[2]),]

background_data = IAM_data[!(IAM_data$writer_id %in% c(writers_ids[1],writers_ids[2])),]


int_characters = sort(intersect(background_data$character,
                                intersect(writer_data_1$character,writer_data_2$character)))
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

mu_hat=matrix(colMeans(background_data[,1:p]),nrow = 1)

S = 0
Sw = 0
for (w in unique(background_data$writer_id)){
  df_writer = background_data[(background_data$writer_id==w),]
  var_data = unname(as.matrix(df_writer[,1:p]))
  #theta_w = matrix(colMeans(var_data), nrow = 1)
  #S.this <- (t(theta_w - mu_hat) %*% (theta_w - mu_hat))
  #S <- S + S.this
  Cov.this = cov(var_data)*(nrow(df_writer)-1) 
  Sw <- Sw + Cov.this
} 

#B_hat = S/(length(unique(background_data$writer_id)) - 1)
B_hat = cov(background_data[,1:p])
if (!is.positive.definite(B_hat)){B_hat = as.matrix(nearPD(B_hat)$mat)}


W_hat <- Sw/(nrow(background_data) - length(unique(background_data$writer_id)))
U_hat <- W_hat*(nw_hat-p-1)

eta=1

#fit <- fitdistr(diag(W_hat), "cauchy")
#print(fit)
#hist(rcauchy(1000, location = fit$estimate[1], scale = fit$estimate[2]))

loc <- mean(log(diag(W_hat)))
sc <- sd(log(diag(W_hat))) # biased estimation: (sum((log(diag(W_hat))-loc)^2)/p)^(1/2)

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
stan_model_niw <- stan_model(file = "niw.stan", model_name = "niw")
stan_model_nlkj <- stan_model(file = "normal_lkj_model.stan", model_name = "normal_lkj_model")

assess_BF <- function(stan_model,stan_data_H0,stan_data_H1_1,stan_data_H1_2){
  fit_H0 <- sampling(stan_model, data = stan_data_H0, iter = 2000, chains = 2, cores=2)
  fit_H1_1 <- sampling(stan_model, data = stan_data_H1_1, iter = 2000, chains = 2, cores=2)
  fit_H1_2 <- sampling(stan_model, data = stan_data_H1_2, iter = 2000, chains = 2, cores=2)
  
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
