setwd("C:/Users/ltzai/Desktop/PhD/Handwritten_Loop_characters/handwriting-evidence-evaluation-fourier-based-feature-extraction/Modelling")
library(readxl)
library(dplyr)
library(MASS)
library(matrixcalc)
library(Matrix)
library(abind)
library(matlib)
library(parallel)
library(writexl)

set.seed(2)

IAM_data <- read_excel("IAM_fourier_features_dataset/DB_loop_handwriting.xlsx")
IAM_data = as.data.frame(IAM_data)

#ALREADY /sqrt(area)
#IAM_data[,2:9] = IAM_data[,2:9]/sqrt(IAM_data$area)
#IAM_data[,1] = log(IAM_data[,1])
#IAM_data = cbind(scale(IAM_data[,1:9]),IAM_data[,10:ncol(IAM_data)])

writers_ids <- unique(IAM_data$writer_id)



count_ch = as.data.frame(IAM_data %>% group_by(writer_id,character) %>% 
                           count()%>% 
                           group_by(character) %>% 
                           summarize(count = sum(n>18)))

count_ch


table(IAM_data$character)

log_marginal_referee <- function(z, alpha = 1, beta = 1) {
  # z: vector of squared amplitudes |A_h|^2 for one harmonic
  n <- length(z)
  lgamma(alpha + n) - lgamma(alpha) +
    alpha * log(beta) -
    (alpha + n) * log(beta + sum(z))
}

# BF under referee's model (across all harmonics, treated independently)
BF_referee <- function(y_coef, x_coef, alpha = 1, beta) {
  # y_coef, x_coef: matrices of shape [n_loops x 4] of squared amplitudes
  # columns = harmonics h=1,...,4
  
  H <- ncol(y_coef)
  log_BF <- 0
  
  for (h in 1:H) {
    z_y   <- y_coef[, h]
    z_x   <- x_coef[, h]
    z_all <- c(z_y, z_x)
    
    log_BF <- log_BF +
      log_marginal_referee(z_all, alpha, beta[h]) -
      log_marginal_referee(z_y,   alpha, beta[h]) -
      log_marginal_referee(z_x,   alpha, beta[h])
  }
  
  return(list(log_BF = log_BF, BF = exp(log_BF)))
}


# # Log marginal under posterior updated on training data z_train
# log_pred_LR <- function(z_test, z_train, alpha = 1, beta_h) {
#   # Numerator: predict z_test given z_train (same writer)
#   alpha_n <- alpha + length(z_train)
#   beta_n  <- beta_h + sum(z_train)
#   log_num <- log_marginal_referee(z_test, alpha = alpha_n, beta = beta_n)
#   
#   # Denominator: predict z_test from background only (different writer)
#   log_den <- log_marginal_referee(z_test, alpha = alpha, beta = beta_h)
#   
#   return(log_num - log_den)
# }
# 
# LR_predictive <- function(y_coef, x_coef, alpha = 1, beta) {
#   H <- ncol(y_coef)
#   log_LR <- 0
#   for (h in 1:H) {
#     log_LR <- log_LR +
#       log_pred_LR(x_coef[,h], y_coef[,h], alpha, beta[h])
#   }
#   return(list(log_BF = log_LR, LR = exp(log_LR)))
# }




same_source_def <- function(character_data,w){
  
  all_chars <- sort(unique(character_data$character))
  l <- length(all_chars)
  
  df_all<-data.frame()
  
  writer_data_all <- character_data[(character_data$writer_id==w),]
  background_data <- character_data[(character_data$writer_id!=w),]
  
  for (iter_for_eval in (1:10)){
    
    sample_size <- min(100, nrow(writer_data_all))
    
    writer_data <- writer_data_all %>%
      add_count(character, name = "char_freq") %>%  # add frequency column
      slice_sample(
        n = sample_size,       # now it's a constant
        weight_by = char_freq, # weighted sampling
        replace = FALSE
      )
    
    questioned_data <- writer_data[1:floor(sample_size/2),]
    suspect_data <- writer_data[(floor(sample_size/2)+1):sample_size,]
    
    questioned_data$character <- factor(questioned_data$character, levels = all_chars)
    suspect_data$character <- factor(suspect_data$character, levels = all_chars)
    
    alphabet_map <- setNames(seq_along(all_chars), all_chars) 
    
    
    Ampl1 <- background_data$a1^2 + background_data$b1^2
    Ampl2 <- background_data$a2^2 + background_data$b2^2
    Ampl3 <- background_data$a3^2 + background_data$b3^2
    Ampl4 <- background_data$a4^2 + background_data$b4^2
    
    
    back_fourier <- cbind(Ampl1, Ampl2, Ampl3, Ampl4)
    
    beta_h <- colMeans(back_fourier)
    
    Ampl1 <- questioned_data$a1^2 + questioned_data$b1^2
    Ampl2 <- questioned_data$a2^2 + questioned_data$b2^2
    Ampl3 <- questioned_data$a3^2 + questioned_data$b3^2
    Ampl4 <- questioned_data$a4^2 + questioned_data$b4^2
    
    
    y_fourier <- cbind(Ampl1, Ampl2, Ampl3, Ampl4)
    
    Ampl1 <- suspect_data$a1^2 + suspect_data$b1^2
    Ampl2 <- suspect_data$a2^2 + suspect_data$b2^2
    Ampl3 <- suspect_data$a3^2 + suspect_data$b3^2
    Ampl4 <- suspect_data$a4^2 + suspect_data$b4^2
    
    x_fourier <- cbind(Ampl1, Ampl2, Ampl3, Ampl4)
    
    BF <- BF_referee(y_fourier, x_fourier, alpha = 1, beta = beta_h)
      
      df_new <- data.frame(
        writer = w,
        character = 'all',
        model     = c("x2"),
        BF        = c(BF$log_BF)
      )
    
    df_all = rbind(df_all,df_new)
  }
  return(df_all)
}


w.list <- sapply(unique(IAM_data$writer_id), list)

df_conjugate_all_list = list()
for (w in w.list){
  df_conjugate_all_list[[w]] <- same_source_def(IAM_data,w)
  print(w)
}

df_conjugate_all <- do.call(rbind, df_conjugate_all_list)

mean(df_conjugate_all$BF<0)

df_conjugate_all[df_conjugate_all$BF<0,]


different_source_def <- function(character_data,composition,w){
  
  all_chars <- sort(unique(character_data$character))
  l <- length(all_chars)
  
  df_all=data.frame()
  
  writer_data_1 = character_data[(character_data$writer_id == composition[w,1]),]
  
  writer_data_2 = character_data[(character_data$writer_id == composition[w,2]),]
  
  writer_data_all <- rbind(writer_data_1,writer_data_2)
  
  background_data = character_data[!(character_data$writer_id %in% c(composition[w,1],
                                                                     composition[w,2])),]
  
  
  for (iter_for_eval in (1:5)){   
    
    
    sample_size <- min(50, nrow(writer_data_1))
    
    questioned_data <- writer_data_1 %>%
      add_count(character, name = "char_freq") %>%  # add frequency column
      slice_sample(
        n = sample_size,       # now it's a constant
        weight_by = char_freq, # weighted sampling
        replace = FALSE
      )
    
    sample_size <- min(50, nrow(writer_data_2))
    suspect_data <- writer_data_2 %>%
      add_count(character, name = "char_freq") %>%  # add frequency column
      slice_sample(
        n = sample_size,       # now it's a constant
        weight_by = char_freq, # weighted sampling
        replace = FALSE
      )
    
    # intersect characters
    #int_characters <- sort(intersect(questioned_data$character,suspect_data$character))
    
    #questioned_data$character <- questioned_data[questioned_data$character %in% int_characters, ]
    #suspect_data$character <- suspect_data[suspect_data$character %in% int_characters, ]
    
    #alphabet_map <- setNames(seq_along(int_characters), int_characters) 
    
    #background_data <- background_data_all[background_data_all$character %in% int_characters,]
    
    questioned_data$character <- factor(questioned_data$character, levels = all_chars)
    suspect_data$character <- factor(suspect_data$character, levels = all_chars)
    
    Ampl1 <- background_data$a1^2 + background_data$b1^2
    Ampl2 <- background_data$a2^2 + background_data$b2^2
    Ampl3 <- background_data$a3^2 + background_data$b3^2
    Ampl4 <- background_data$a4^2 + background_data$b4^2
    
    
    back_fourier <- cbind(Ampl1, Ampl2, Ampl3, Ampl4)
    
    beta_h <- colMeans(back_fourier)
    

    
    Ampl1 <- questioned_data$a1^2 + questioned_data$b1^2
    Ampl2 <- questioned_data$a2^2 + questioned_data$b2^2
    Ampl3 <- questioned_data$a3^2 + questioned_data$b3^2
    Ampl4 <- questioned_data$a4^2 + questioned_data$b4^2
    
    
    y_fourier <- cbind(Ampl1, Ampl2, Ampl3, Ampl4)
    
    Ampl1 <- suspect_data$a1^2 + suspect_data$b1^2
    Ampl2 <- suspect_data$a2^2 + suspect_data$b2^2
    Ampl3 <- suspect_data$a3^2 + suspect_data$b3^2
    Ampl4 <- suspect_data$a4^2 + suspect_data$b4^2
    
    x_fourier <- cbind(Ampl1, Ampl2, Ampl3, Ampl4)
    
    BF <- BF_referee(y_fourier, x_fourier, alpha = 1, beta = beta_h)
    
    df_new <- data.frame(
      writer1 = composition[w,1],
      writer2 = composition[w,2],
      character = 'all',
      model     = c("x2"),
      BF        = c(BF$log_BF)
    )
    
    
    df_all = rbind(df_all,df_new)
  }
  return(df_all)
}

comp_writers = t(combn(unique(IAM_data$writer_id), 2))

w.list <- sapply(1:nrow(comp_writers), list)

df_conjugate_all_list = list()
for (w in 1:nrow(comp_writers)){
  df_conjugate_all_list[[w]] <- different_source_def(IAM_data,comp_writers,w)
  print(w)
}

df_conjugate_all_ds <- do.call(rbind, df_conjugate_all_list)

mean(df_conjugate_all_ds$BF>0)

df_conjugate_all_ds[df_conjugate_all_ds$BF>0,]

