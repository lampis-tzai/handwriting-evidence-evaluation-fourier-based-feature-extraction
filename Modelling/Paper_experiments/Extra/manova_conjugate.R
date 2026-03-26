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
library(rstan)
library(bridgesampling)
source('Paper_experiments/Stan_BF_calculation.R')

set.seed(2)

IAM_data <- read_excel("IAM_fourier_features_dataset/DB_loop_handwriting.xlsx")
IAM_data = as.data.frame(IAM_data)


IAM_data[,2:9] = IAM_data[,2:9]/sqrt(IAM_data$area)
IAM_data[,1] = log(IAM_data[,1])
IAM_data = cbind(scale(IAM_data[,1:9]),IAM_data[,10:ncol(IAM_data)])

writers_ids <- unique(IAM_data$writer_id)



count_ch = as.data.frame(IAM_data %>% group_by(writer_id,character) %>% 
                           count()%>% 
                           group_by(character) %>% 
                           summarise(count = sum(n>18)))

count_ch


table(IAM_data$character)


background_statistics_br <- function(background_data){
  
  p=9
  l = length(unique(background_data$character))
  nw.min = p + 2
  nw_hat = nw.min
  
  a_data = background_data[(background_data$character==1),1:p]
  mu_hat=matrix(colMeans(a_data),nrow = 1)
  
  S = 0
  for (w in unique(background_data$writer_id)){
    df_writer = background_data[(
      background_data$character==1)& (background_data$writer_id==w),]
    if (nrow(df_writer)>2){
      theta_w = matrix(colMeans(df_writer[,1:p]), nrow = 1)
      S.this <- (t(theta_w - mu_hat) %*% (theta_w - mu_hat))
      S <- S + S.this
    }
  }
  
  B_hat = S/(length(unique(background_data$writer_id)) - 1)
  #B_hat = cov(a_data)
  
  if (!is.positive.definite(B_hat)){B_hat = as.matrix(nearPD(B_hat)$mat)}
  
  beta_mu = array(0, dim=c(l,p))
  beta_cov = array(0, dim=c(p,p,l))
  for (l_id in 1:l){
    letter_data = as.matrix(unname(background_data[(
      background_data$character==l_id),1:p]))
    if (nrow(letter_data)>2){
      letter_diff = letter_data - matrix(mu_hat[col(letter_data)],ncol = p)
      beta_l = colMeans(letter_diff)
      beta_mu[l_id,] = beta_l
      S = 0
      for (w in unique(background_data$writer_id)){
        letter_writer = background_data[(
          background_data$character==l_id)& (background_data$writer_id==w),1:p]
        if (nrow(letter_writer)>2){
          a_data_writer = background_data[(
            background_data$character==1)& (background_data$writer_id==w),1:p]
          
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
  }
  
  
  Sw = 0
  for (w in unique(background_data$writer_id)){
    df_writer = background_data[(background_data$writer_id==w),]
    if (nrow(df_writer)>2){
      Cov.this = cov(df_writer[,1:p])*(nrow(df_writer)-1)
      Sw <- Sw + Cov.this
    }
  }
  W_hat <- Sw/(nrow(background_data) - length(unique(background_data$writer_id)))
  U_hat <- W_hat * (nw_hat - p  -1)
  
  eta <- 4
  
  
  loc <- mean(log(0.5 * diag(W_hat)))
  sc  <- sd(log(0.5 * diag(W_hat)))
  
  return(list(mu_hat,B_hat,beta_mu,beta_cov,nw_hat,U_hat,loc,sc,eta))
}



same_source_def <- function(character_data,w){
  
  all_chars <- sort(unique(character_data$character))
  l <- length(all_chars)
  
  df_all<-data.frame()
  
  writer_data_all <- character_data[(character_data$writer_id==w),]
  background_data <- character_data[(character_data$writer_id!=w),]
  
  for (iter_for_eval in (1:1)){
    
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
    
    
    questioned_data$character <- as.numeric(alphabet_map[questioned_data$character]) 
    suspect_data$character <- as.numeric(alphabet_map[suspect_data$character])
    
    
    background_data$character <- as.numeric(alphabet_map[background_data$character])
    
    background_stats_br <- background_statistics_br(background_data)
    
    
    
    manova_conjugate <- MANOVA_conjugate(questioned_data,
                                         suspect_data,
                                         background_stats_br)
    
   
    
    df_new <- data.frame(
      writer = w,
      character = 'all',
      model     = c("manova_conjugate"),
      BF        = c(manova_conjugate)
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
  
  
  for (iter_for_eval in (1)){   
    
    
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
    
    #table(questioned_data$character)
    #table(suspect_data$character)
    
    alphabet_map <- setNames(seq_along(all_chars), all_chars) 
    
    
    questioned_data$character <- as.numeric(alphabet_map[questioned_data$character]) 
    suspect_data$character <- as.numeric(alphabet_map[suspect_data$character])
    
    
    background_data$character <- as.numeric(alphabet_map[background_data$character])

    
    background_stats_br <- background_statistics_br(background_data)
    
    
    
    manova_conjugate <- MANOVA_conjugate(questioned_data,
                                         suspect_data,
                                         background_stats_br)
    
    df_new <- data.frame(
      writer1 = composition[w,1],
      writer2 = composition[w,2],
      character = 'all',
      model     = c("manova_conjugate"),
      BF        = c(manova_conjugate)
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
