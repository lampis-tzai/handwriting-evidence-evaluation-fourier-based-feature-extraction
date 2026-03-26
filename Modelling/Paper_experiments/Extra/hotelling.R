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
library(Hotelling)

set.seed(2)

IAM_data <- read_excel("IAM_fourier_features_dataset/DB_loop_handwriting.xlsx")
IAM_data = as.data.frame(IAM_data)


IAM_data[,2:9] = IAM_data[,2:9]#/sqrt(IAM_data$area)
IAM_data[,1] = log(IAM_data[,1])
IAM_data = cbind(scale(IAM_data[,1:9]),IAM_data[,10:ncol(IAM_data)])

writers_ids <- unique(IAM_data$writer_id)


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
    
    htest<- hotelling.test(x = questioned_data[,1:9], y = suspect_data[,1:9])
    
    df_new <- data.frame(
      writer = w,
      character = 'all',
      model     = c("hotelling"),
      pval        = c(htest$pval)
    )
    
    df_all = rbind(df_all,df_new)
  }
  return(df_all)
}
character_data<-IAM_data

w.list <- sapply(unique(IAM_data$writer_id), list)

df_conjugate_all_list = list()
for (w in w.list){
  df_conjugate_all_list[[w]] <- same_source_def(IAM_data,w)
  print(w)
}

df_conjugate_all <- do.call(rbind, df_conjugate_all_list)

mean(df_conjugate_all$pval<0.01)

df_conjugate_all[df_conjugate_all$pval<0.01,]


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
    htest<- hotelling.test(x = questioned_data[,1:9], y = suspect_data[,1:9])
    
    df_new <- data.frame(
      writer1 = composition[w,1],
      writer2 = composition[w,2],
      character = 'all',
      model     = c("hotelling"),
      pval        = c(htest$pval)
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

mean(df_conjugate_all_ds$pval>0.01)

df_conjugate_all_ds[df_conjugate_all_ds$pval>0.01,]

