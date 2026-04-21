setwd("C:/Users/ltzai/Desktop/PhD/Handwritten_Loop_characters/handwriting-evidence-evaluation-fourier-based-feature-extraction/Modelling")
library(readxl)
library(dplyr)
library(MASS)
library(matrixcalc)
library(Matrix)
library(matlib)
library(CholWishart)

set.seed(2)

IAM_data <- read_excel("IAM_fourier_features_dataset/DB_loop_handwriting_ls.xlsx")
IAM_data = as.data.frame(IAM_data)

#IAM_data[,2:9] = IAM_data[,2:9]/sqrt(IAM_data$area)
#IAM_data[,1] = log(IAM_data[,1])

IAM_data = cbind(scale(IAM_data[,1:9]),IAM_data[,10:ncol(IAM_data)])

writers_ids <- unique(IAM_data$writer_id)


#library(PerformanceAnalytics)
#chart.Correlation(IAM_data[,2:9], histogram = TRUE, method = "pearson")


library(nortest)
p_values <- array(0, dim = c(length(writers_ids),9))
i=1
for (w in writers_ids){

  writer_data_all = IAM_data[(IAM_data$writer_id==w),1:9]

  for (f in 1:9){
    p_values[i,f] <- lillie.test(writer_data_all[,f])$p.value
  }
  i=i+1
}
colMeans(p_values>0.01)



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
writer_data_1 <- IAM_data[IAM_data$writer_id == "62", ]

writer_data_2 <- IAM_data[IAM_data$writer_id == "112", ]

background_data <- IAM_data[!(IAM_data$writer_id %in% c("62", "112")), ]

sample_size <- floor(nrow(writer_data_1)/2)#min(50, nrow(writer_data_1))

questioned_data <- writer_data_1 %>%
  add_count(character, name = "char_freq") %>%  # add frequency column
  slice_sample(
    n = sample_size,       # now it's a constant
    weight_by = char_freq, # weighted sampling
    replace = FALSE
  )

sample_size <- floor(nrow(writer_data_2)/2)#min(50, nrow(writer_data_2))
suspect_data <- writer_data_2 %>%
  add_count(character, name = "char_freq") %>%  # add frequency column
  slice_sample(
    n = sample_size,       # now it's a constant
    weight_by = char_freq, # weighted sampling
    replace = FALSE
  )



int_characters <- sort(unique(IAM_data$character))
l <- length(int_characters)

questioned_data$character <- factor(questioned_data$character, levels = int_characters)
suspect_data$character <- factor(suspect_data$character, levels = int_characters)

alphabet_map <- setNames(seq_along(int_characters), int_characters)
# 
# int_characters <- sort(intersect(questioned_data$character,suspect_data$character))
# 
# l <- length(int_characters)
# 
# alphabet_map <- setNames(seq_along(int_characters), int_characters)
# 
# 
# questioned_data$character <- ifelse(questioned_data$character %in% int_characters,
#                                    questioned_data$character,
#                                    'a')
# 
# suspect_data$character <- ifelse(suspect_data$character %in% int_characters,
#                                    suspect_data$character,
#                                    'a')
# 
# background_data$character <- ifelse(background_data$character %in% int_characters,
#                                    background_data$character,
#                                 'a')

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
#~27.5

nw_hat = nw.min

a_data = background_data[(background_data$character=='a'),]
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
  letter_data = as.matrix(unname(background_data[(background_data$character==int_characters[l_id]),1:p]))
  
  letter_diff = letter_data - matrix(mu_hat[col(letter_data)],ncol = p)
  beta_l = colMeans(letter_diff)
  beta_mu[l_id,] = beta_l
  S = 0
  n_rows <- 0
  for (w in unique(background_data$writer_id)){
    letter_writer = background_data[(
      background_data$character==int_characters[l_id])& (background_data$writer_id==w),1:p]
    if (nrow(letter_writer)>3){
      n_rows <- n_rows + nrow(letter_writer)
      a_data_writer = background_data[(
        background_data$character=='a')& (background_data$writer_id==w),1:p]

      mu_hat_writer=matrix(colMeans(a_data_writer),nrow = 1)

      letter_diff_writer = letter_writer - matrix(mu_hat_writer[col(letter_writer)],ncol = p)

      beta_w = matrix(colMeans(letter_diff_writer), nrow = 1)
      S.this <- (t(beta_w - beta_l) %*% (beta_w - beta_l))#*nrow(letter_writer)
      S <- S + S.this
    }
  }
  B_hat_l = S/(n_rows - 1)
  #B_hat_l = cov(letter_data)
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


#loc <- mean(log(diag(W_hat))/2)
#sc <-  sd(log(diag(W_hat))/2)



log_sd_mat <- sapply(unique(background_data$writer_id), function(w) {
  df_w <- background_data[background_data$writer_id == w, 1:p]
  if (nrow(df_w) <= 3) return(rep(NA_real_, p))  # skip degenerate writers
  log(sqrt(diag(cov(as.matrix(df_w)))))
})

# Remove degenerate writers (columns with NA)
log_sd_mat <- log_sd_mat[, colSums(is.na(log_sd_mat)) == 0]

loc <- rowMeans(log_sd_mat)
sc  <- apply(log_sd_mat, 1, sd)

# loc <- mean(loc)
# sc  <- mean(sc)

# sds <- do.call(c, lapply(unique(background_data$writer_id), function(w) {
#   df_w <- background_data[background_data$writer_id == w, 1:p]
#   apply(df_w, 2, sd)  # log-SD per feature per writer
# }))
# 
# # log_sds<-log(apply(background_data[,1:p],2,sd))
# #
#  loc <- mean(log(sds))
#  sc  <- sd(log(sds))




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
# ~9.9
eta <- 1



beta_cov_list <- vector("list", dim(beta_cov)[3])
for (i in 1:dim(beta_cov)[3]) {
  beta_cov_list[[i]] <- beta_cov[, , i]
}





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





#################################
#posterior plots
#################################

fit_H1_1 <- sampling(stan_model__manova_lkj, data = stan_data_H1_1, iter = 2000, chains =1, cores=1)



post_samples = extract(fit_H1_1)$theta

post_samples = abind::abind(post_samples, extract(fit_H1_1)$beta_raw, along = 3) 

i = 1
questioned_posterior_list = list()
for (character in int_characters){
  questioned_posterior_character = as.data.frame(post_samples[,,alphabet_map[[character]]])
  questioned_posterior_character['character'] =character
  names(questioned_posterior_character) = c( "area","a1", "b1", "a2", "b2", "a3","b3", "a4","b4","character")
  questioned_posterior_list[[i]] = questioned_posterior_character
  i = i + 1
}

questioned_posterior <- do.call("rbind", questioned_posterior_list)
questioned_posterior['Source'] = "Writer 152" 


fit_H1_2 <- sampling(stan_model__manova_lkj, data = stan_data_H1_2, iter = 2000, chains = 1, cores=1)

post_samples = extract(fit_H1_2)$theta

post_samples = abind::abind(post_samples, extract(fit_H1_2)$beta_raw, along = 3) 

i = 1
suspect_posterior_list = list()
for (character in int_characters){
  suspect_posterior_character = as.data.frame(post_samples[,,alphabet_map[[character]]])
  suspect_posterior_character['character'] =character
  names(suspect_posterior_character) = c("area", "a1", "b1", "a2", "b2", "a3","b3", "a4","b4","character")
  suspect_posterior_list[[i]] = suspect_posterior_character
  i = i + 1
}

suspect_posterior <- do.call("rbind", suspect_posterior_list)

suspect_posterior['Source'] = 'Writer 88'

posterior_data = rbind(questioned_posterior,suspect_posterior)

library(reshape2)
posterior_data_melt = melt(posterior_data, id = c("character","Source"), 
                           variable.name = "Coefficient")

posterior_data_melt$Source <- factor(posterior_data_melt$Source, levels = c("Writer 152",'Writer 88'))

#posterior_data_melt_abde = posterior_data_melt[posterior_data_melt$character %in% c('a','b','d','e'),]

library(ggplot2)
p1 = ggplot(posterior_data_melt, aes(x = value, fill = Source)) + 
  geom_histogram(aes(y = ..density..),alpha = 0.5, position = "identity")+
  facet_grid(character~Coefficient,scales="free_y")+ 
  #labs(title="Posterior distribution of Fourier coefficients") + 
  theme(strip.text = element_text(size = 20),
        axis.title=element_blank(),#element_text(size=20,face="bold"),
        plot.title = element_text(hjust = 0.5,size =20,face="bold"))

p1
jpeg("Paper_experiments/questioned_suspect_posterior.jpg",width=3920, height=2000, res=300)
p1
dev.off()



post_samples =  extract(fit_H1_1)$theta


questioned_posterior_list = list()
questioned_posterior_character = as.data.frame(post_samples)
values_a = questioned_posterior_character
questioned_posterior_character['character'] =int_characters[1]
names(questioned_posterior_character) = c( "a1", "b1", "a2", "b2", "a3","b3", "a4","b4","character")
questioned_posterior_list[[1]] = questioned_posterior_character


post_samples = abind::abind(post_samples, extract(fit_H1_2)$beta_raw, along = 3) 

i = 2
for (character in int_characters[-1]){
  questioned_posterior_character = as.data.frame(post_samples[,,alphabet_map[[character]]])+values_a
  questioned_posterior_character['character'] =character
  names(questioned_posterior_character) = c("a1", "b1", "a2", "b2", "a3","b3", "a4","b4","character")
  questioned_posterior_list[[i]] = questioned_posterior_character
  i = i + 1
}

questioned_posterior <- do.call("rbind", questioned_posterior_list)
questioned_posterior['Source'] = 'Questioned'



post_samples = extract(fit_H1_2)$theta

suspect_posterior_list = list()
suspect_posterior_character = as.data.frame(post_samples)
values_a = suspect_posterior_character
suspect_posterior_character['character'] =int_characters[1]
names(suspect_posterior_character) = c( "a1", "b1", "a2","b2", "a3","b3", "a4","b4","character")
suspect_posterior_list[[1]] = suspect_posterior_character

post_samples = abind::abind(post_samples, extract(fit_H1_2)$beta_raw, along = 3) 

i = 2
for (character in int_characters[-1]){
  suspect_posterior_character = as.data.frame(post_samples[,,alphabet_map[[character]]])+values_a
  suspect_posterior_character['character'] =character
  names(suspect_posterior_character) = c( "a1", "b1", "a2", "b2", "a3","b3", "a4","b4","character")
  suspect_posterior_list[[i]] = suspect_posterior_character
  i = i + 1
}


suspect_posterior <- do.call("rbind", suspect_posterior_list)

suspect_posterior['Source'] = 'Writer 36'#'P.o.i.'#

posterior_data = rbind(questioned_posterior,suspect_posterior)

library(reshape2)
posterior_data_melt = melt(posterior_data, id = c("character","Source"),
                           variable.name = "Coefficient")

posterior_data_melt$Source <- factor(posterior_data_melt$Source, levels = c("Questioned","Writer 36"))#'P.o.i.'))#


p2 = ggplot(posterior_data_melt, 
            aes(x = character , y = value, fill = character)) + 
  geom_boxplot(alpha = 0.5) + 
  facet_grid(Source~Coefficient, axes = "all", axis.labels = "all") + 
  #labs(title = "Comparison of Fourier Coefficients by Character") + 
  theme(strip.text = element_text(size = 20),
        axis.title = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
        legend.position = "none")

p2
# jpeg("Stan_code/Plots/character_difference.jpg",width=3920, height=3000, res=300)
# p2
# dev.off()

library(gridExtra)

plots <- list()

# Loop through unique combinations of Source and Coefficient
for(src in unique(posterior_data_melt$Source)) {
  for(coef in unique(posterior_data_melt$Coefficient)) {
    # Subset the data
    subset_data <- subset(posterior_data_melt, Source == src & Coefficient == coef)
    
    # Create individual plot
    p <- ggplot(subset_data, aes(x = character, y = value, fill = character)) + 
      geom_boxplot(alpha = 0.5) + 
      ggtitle(paste(src, "\n", coef)) +  # Simplified title
      theme(strip.text = element_text(size = 20),
            axis.title = element_blank(),
            plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
            legend.position = "none")
    
    # Add plot to list
    plots[[length(plots) + 1]] <- p
  }
}

jpeg("Stan_code/Plots/character_difference_writer_36.jpg",width=5000, height=3000, res=300)
grid.arrange(grobs = plots, ncol = length(unique(posterior_data_melt$Coefficient)))
dev.off()

# posterior_data_melt_qop = posterior_data_melt[posterior_data_melt$character %in% c('g','o','p'),]
# 
# 
# p2 = ggplot(posterior_data_melt_qop, aes(x = value, fill = indicator)) + 
#   geom_histogram(aes(y = ..density..), alpha = 0.5, position = "identity")+
#   facet_grid(character~Coefficient,scales="free")+ 
#   labs(title="Posterior distribution of Fourier coefficients parameters \n per character") + 
#   theme(strip.text = element_text(size = 20),
#         axis.title=element_text(size=20,face="bold"),
#         plot.title = element_text(hjust = 0.5,size =20,face="bold"))
# 
# p2
# jpeg("questioned_suspect_posterior_gop.jpg",width=3920, height=2000, res=300)
# p2
# dev.off()




p3 = ggplot(posterior_data_melt, aes(x = value, fill = Source)) + 
  geom_histogram(aes(y = ..density..), alpha = 0.5, position = "identity")+
  facet_grid(character~Coefficient,scales="free")+ 
  labs(title="Posterior distribution of Fourier coefficients parameters \n per character") + 
  theme(strip.text = element_text(size = 20),
        axis.title=element_text(size=20,face="bold"),
        plot.title = element_text(hjust = 0.5,size =20,face="bold"))

jpeg("Stan_code/Plots/questioned_suspect_posterior.jpg",width=3920, height=2000, res=300)
p3
dev.off()



