library(CholWishart)





assess_BF <- function(stan_model,stan_data_H0,stan_data_H1_1,stan_data_H1_2){
  
  fit_H0 <- sampling(stan_model, data = stan_data_H0, iter = 2000, chains = 1, cores=1, refresh = 0)
  fit_H1_1 <- sampling(stan_model, data = stan_data_H1_1, iter = 2000, chains = 1, cores=1, refresh = 0)
  fit_H1_2 <- sampling(stan_model, data = stan_data_H1_2, iter = 2000, chains = 1, cores=1, refresh = 0)
  
  # samples_H0 <- extract(fit_H0)
  # samples_H1_1 <- extract(fit_H1_1)
  # samples_H1_2 <- extract(fit_H1_2)
  # 
  # Compute the log marginal likelihood
  bs_lik_H0 <- bridge_sampler(fit_H0, method = 'warp3', silent = TRUE)$logml
  bs_lik_H1_1 <- bridge_sampler(fit_H1_1, method = 'warp3', silent = TRUE)$logml
  bs_lik_H1_2 <- bridge_sampler(fit_H1_2, method = 'warp3', silent = TRUE)$logml
  
  return(bs_lik_H0 - bs_lik_H1_1 - bs_lik_H1_2)
  
}



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



niw_conjugate <- function(questioned_data, suspect_data, background_stats){
  
  p = nrow(background_stats[[2]])
  
  writer_data <- rbind(questioned_data,suspect_data)
  
  stan_data_H0 <- list(N = nrow(writer_data), 
                       P = p, 
                       y = unname(as.matrix(writer_data[,1:p])), 
                       mu = background_stats[[1]], 
                       B = background_stats[[2]],
                       U = background_stats[[4]],
                       nu = background_stats[[3]]
  )
  
  stan_data_H1_1 <- list(N = nrow(questioned_data), 
                         P = p, 
                         y = unname(as.matrix(questioned_data[,1:p])), 
                         mu = background_stats[[1]], 
                         B = background_stats[[2]],
                         U = background_stats[[4]],
                         nu = background_stats[[3]]
  )
  
  stan_data_H1_2 <- list(N = nrow(suspect_data), 
                         P = p, 
                         y = unname(as.matrix(suspect_data[,1:p])), 
                         mu = background_stats[[1]], 
                         B = background_stats[[2]],
                         U = background_stats[[4]],
                         nu = background_stats[[3]]
  )
  
  log_ml_H0 <- marginal_likelihood_niw_conjugate(stan_data_H0)
  log_ml_H1_1 <- marginal_likelihood_niw_conjugate(stan_data_H1_1)
  log_ml_H1_2 <- marginal_likelihood_niw_conjugate(stan_data_H1_2)
  
  logBF <- log_ml_H0 - log_ml_H1_1 - log_ml_H1_2
  return(logBF)
}



normal_iW <- function(questioned_data, suspect_data, background_stats){
  
  p = nrow(background_stats[[2]])
  
  writer_data <- rbind(questioned_data,suspect_data)
  
  stan_data_H0 <- list(N = nrow(writer_data), 
                       P = p, 
                       y = unname(as.matrix(writer_data[,1:p])), 
                       mu = as.vector(background_stats[[1]]), 
                       B = background_stats[[2]],
                       U = background_stats[[4]],
                       nu = background_stats[[3]]
  )
  
  stan_data_H1_1 <- list(N = nrow(questioned_data), 
                         P = p, 
                         y = unname(as.matrix(questioned_data[,1:p])), 
                         mu =as.vector(background_stats[[1]]), 
                         B = background_stats[[2]],
                         U = background_stats[[4]],
                         nu = background_stats[[3]]
  )
  
  stan_data_H1_2 <- list(N = nrow(suspect_data), 
                         P = p, 
                         y = unname(as.matrix(suspect_data[,1:p])), 
                         mu = as.vector(background_stats[[1]]), 
                         B = background_stats[[2]],
                         U = background_stats[[4]],
                         nu = background_stats[[3]]
  )
  
  
  logBF <- assess_BF(stan_model_niw,stan_data_H0,stan_data_H1_1,stan_data_H1_2)
  return(logBF)
}

normal_lkj <- function(questioned_data, suspect_data, background_stats){
  
  p = nrow(background_stats[[2]])
  
  writer_data <- rbind(questioned_data,suspect_data)
  
  stan_data_H0 <- list(N = nrow(writer_data), 
                       P = p, 
                       y = unname(as.matrix(writer_data[,1:p])), 
                       mu = as.vector(background_stats[[1]]), 
                       B = background_stats[[2]],
                       loc = background_stats[[5]],
                       sc = background_stats[[6]],
                       eta = background_stats[[7]]
  )
  
  stan_data_H1_1 <- list(N = nrow(questioned_data), 
                         P = p, 
                         y = unname(as.matrix(questioned_data[,1:p])), 
                         mu = as.vector(background_stats[[1]]), 
                         B = background_stats[[2]],
                         loc = background_stats[[5]],
                         sc = background_stats[[6]],
                         eta = background_stats[[7]]
  )
  
  stan_data_H1_2 <- list(N = nrow(suspect_data), 
                         P = p, 
                         y = unname(as.matrix(suspect_data[,1:p])), 
                         mu = as.vector(background_stats[[1]]), 
                         B = background_stats[[2]],
                         loc = background_stats[[5]],
                         sc = background_stats[[6]],
                         eta = background_stats[[7]]
  )
  
  logBF <- assess_BF(stan_model_nlkj,stan_data_H0,stan_data_H1_1,stan_data_H1_2)
  return(logBF)
}


marginal_likelihood_manova_conjugate<- function(w_data,stan_data){
  
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

MANOVA_conjugate <- function(questioned_data, suspect_data, background_stats){
  
  p = nrow(background_stats[[2]])
  
  writer_data = rbind(questioned_data,suspect_data)
  
  l = dim(background_stats[[4]])[3]
  
  
  
  stan_data_H0 <- list(N = nrow(writer_data), 
                       P = p, 
                       L = l,
                       letters = writer_data$character,
                       y = unname(as.matrix(writer_data[,1:p])),
                       mu = background_stats[[1]], 
                       B = background_stats[[2]], 
                       beta_mu=background_stats[[3]], 
                       beta_cov=background_stats[[4]], 
                       U = background_stats[[6]],
                       nu = background_stats[[5]])
  
  stan_data_H1_1 <- list(N = nrow(questioned_data), 
                         P = p, 
                         L = l,
                         letters = questioned_data$character,
                         y = unname(as.matrix(questioned_data[,1:p])),
                         mu = background_stats[[1]], 
                         B = background_stats[[2]], 
                         beta_mu=background_stats[[3]], 
                         beta_cov=background_stats[[4]], 
                         U = background_stats[[6]],
                         nu = background_stats[[5]])
  
  stan_data_H1_2 <- list(N = nrow(suspect_data), 
                         P = p, 
                         L = l,
                         letters =  suspect_data$character,
                         y = unname(as.matrix(suspect_data[,1:p])),
                         mu = background_stats[[1]], 
                         B = background_stats[[2]], 
                         beta_mu=background_stats[[3]], 
                         beta_cov=background_stats[[4]], 
                         U = background_stats[[6]],
                         nu = background_stats[[5]])
  
  
  log_ml_H0 <- marginal_likelihood_manova_conjugate(writer_data,stan_data_H0)
  log_ml_H1_1 <- marginal_likelihood_manova_conjugate(questioned_data,stan_data_H1_1)
  log_ml_H1_2 <- marginal_likelihood_manova_conjugate(suspect_data,stan_data_H1_2)
  
  logBF <- log_ml_H0 - log_ml_H1_1 - log_ml_H1_2
  return(logBF)
  
}


MANOVA_iw <- function(questioned_data, suspect_data, background_stats){
  
  p = nrow(background_stats[[2]])
  
  writer_data <- rbind(questioned_data,suspect_data)
  
  l = dim(background_stats[[4]])[3]
  
  beta_cov <-background_stats[[4]]
  beta_cov_list <- vector("list", dim(beta_cov)[3])
  for (i in 1:dim(beta_cov)[3]) {
    beta_cov_list[[i]] <- beta_cov[, , i]
  }
  
  stan_data_H0 <- list(N = nrow(writer_data), 
                       P = p, 
                       L = l,
                       letters = writer_data$character,
                       y = unname(as.matrix(writer_data[,1:p])),
                       mu = as.vector(background_stats[[1]]), 
                       B = as.matrix(background_stats[[2]]), 
                       beta_mu=background_stats[[3]], 
                       beta_cov=beta_cov_list, 
                       U = background_stats[[6]],
                       nu = background_stats[[5]])
  
  stan_data_H1_1 <- list(N = nrow(questioned_data), 
                         P = p, 
                         L = l,
                         letters = questioned_data$character,
                         y = unname(as.matrix(questioned_data[,1:p])),
                         mu = as.vector(background_stats[[1]]), 
                         B = as.matrix(background_stats[[2]]), 
                         beta_mu=background_stats[[3]], 
                         beta_cov=beta_cov_list, 
                         U = background_stats[[6]],
                         nu = background_stats[[5]])
  
  stan_data_H1_2 <- list(N = nrow(suspect_data), 
                         P = p, 
                         L = l,
                         letters =  suspect_data$character,
                         y = unname(as.matrix(suspect_data[,1:p])),
                         mu = as.vector(background_stats[[1]]), 
                         B = as.matrix(background_stats[[2]]), 
                         beta_mu=background_stats[[3]], 
                         beta_cov=beta_cov_list, 
                         U = background_stats[[6]],
                         nu = background_stats[[5]])
  
  
  logBF <- assess_BF(stan_model_manova_iw,stan_data_H0,stan_data_H1_1,stan_data_H1_2)
  return(logBF)
  
}

MANOVA_LKJ <- function(questioned_data, suspect_data, background_stats){
  
  p = nrow(background_stats[[2]])
  
  writer_data <- rbind(questioned_data,suspect_data)
  
  l = dim(background_stats[[4]])[3]
  
  beta_cov <-background_stats[[4]]
  beta_cov_list <- vector("list", dim(beta_cov)[3])
  for (i in 1:dim(beta_cov)[3]) {
    beta_cov_list[[i]] <- beta_cov[, , i]
  }
  
  stan_data_H0 <- list(N = nrow(writer_data), 
                       P = p, 
                       L = l,
                       letters = writer_data$character,
                       y = unname(as.matrix(writer_data[,1:p])),
                       mu = as.vector(background_stats[[1]]), 
                       B = as.matrix(background_stats[[2]]), 
                       beta_mu=background_stats[[3]], 
                       beta_cov=beta_cov_list, 
                       loc = background_stats[[7]], 
                       sc = background_stats[[8]], 
                       eta = background_stats[[9]])
  
  stan_data_H1_1 <- list(N = nrow(questioned_data), 
                         P = p, 
                         L = l,
                         letters = questioned_data$character,
                         y = unname(as.matrix(questioned_data[,1:p])),
                         mu = as.vector(background_stats[[1]]), 
                         B = as.matrix(background_stats[[2]]), 
                         beta_mu=background_stats[[3]], 
                         beta_cov=beta_cov_list, 
                         loc = background_stats[[7]], 
                         sc = background_stats[[8]], 
                         eta = background_stats[[9]])
  
  stan_data_H1_2 <- list(N = nrow(suspect_data), 
                         P = p, 
                         L = l,
                         letters =  suspect_data$character,
                         y = unname(as.matrix(suspect_data[,1:p])),
                         mu = as.vector(background_stats[[1]]), 
                         B = as.matrix(background_stats[[2]]), 
                         beta_mu=background_stats[[3]], 
                         beta_cov=beta_cov_list, 
                         loc = background_stats[[7]], 
                         sc = background_stats[[8]], 
                         eta = background_stats[[9]])
  
  logBF <- assess_BF(stan_model_manova_lkj,stan_data_H0,stan_data_H1_1,stan_data_H1_2)
  return(logBF)
  
}