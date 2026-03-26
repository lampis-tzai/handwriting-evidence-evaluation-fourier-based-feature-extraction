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
                           summarize(count = sum(n>18)))

count_ch


table(IAM_data$character)


background_statistics_niw <- function(background_data){
  
  
  p=9
  nw.min = p + 2
  nw_hat = 25
  
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
  if (!is.positive.definite(B_hat)){B_hat = as.matrix(nearPD(B_hat)$mat)}
  
  W_hat <- Sw/(nrow(background_data) - length(unique(background_data$writer_id)))
  U_hat <- W_hat*(nw_hat-p-1)
  
  log_sds <- do.call(c, lapply(unique(background_data$writer_id), function(w) {
    df_w <- background_data[background_data$writer_id == w, 1:p]
    log(apply(df_w, 2, sd))  # log-SD per feature per writer
  }))
  
  loc <- mean(log_sds)
  sc  <- sd(log_sds)
  eta=4
  
  return(list(mu_hat,B_hat,nw_hat,U_hat,loc,sc,eta))
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
    
    background_stats_niw_all <- background_statistics_niw(background_data)
    
    niw_conjugate_all <- niw_conjugate(questioned_data,
                                       suspect_data,
                                       background_stats_niw_all)
    
    
    
    df_new <- data.frame(
      writer = w,
      character = 'all',
      model     = c("niw_conjugate"),
      BF        = c(niw_conjugate_all)
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
    
    
    background_stats_niw_all <- background_statistics_niw(background_data)
    
    niw_conjugate_all <- niw_conjugate(questioned_data,
                                       suspect_data,
                                       background_stats_niw_all)
    
    
    df_new <- data.frame(
      writer1 = composition[w,1],
      writer2 = composition[w,2],
      character = 'all',
      model     = c("niw_conjugate"),
      BF        = c(niw_conjugate_all)
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

