setwd("C:/Users/Lampis_lab/Desktop/PhD/Handwritten_Loop_characters/handwriting-evidence-evaluation-fourier-based-feature-extraction/Modelling")
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


IAM_data <- read_excel("IAM_fourier_features_dataset/DB_loop_handwriting_ls.xlsx")
IAM_data = as.data.frame(IAM_data)


#IAM_data[,2:9] = IAM_data[,2:9]/sqrt(IAM_data$area)
#IAM_data[,1] = log(IAM_data[,1])

IAM_data = cbind(scale(IAM_data[,1:9]),IAM_data[,10:ncol(IAM_data)])

writers_ids <- unique(IAM_data$writer_id)


count_ch = as.data.frame(IAM_data %>% group_by(writer_id,character) %>% 
                           count()%>% 
                           group_by(character) %>% 
                           summarise(count = sum(n>3)))

count_ch


table(IAM_data$character)




background_statistics_niw <- function(background_data){
  
  
  p=9
  nw.min = p + 2
  nw_hat = nw.min
  
  mu_hat=matrix(colMeans(do.call(rbind, lapply(unique(background_data$writer_id), function(w)
    colMeans(background_data[background_data$writer_id == w, 1:p])))), nrow = 1)
  
  S <- 0
  Sw <- 0
  n_rows <- 0
  n_writers <- 0
  for (w in unique(background_data$writer_id)){
    df_writer = background_data[(background_data$writer_id==w),]
    if (nrow(df_writer)>2){
      n_rows <- n_rows + nrow(df_writer)
      n_writers <- n_writers+1
      var_data <- unname(as.matrix(df_writer[,1:p]))
      theta_w <- matrix(colMeans(var_data), nrow = 1)
      S.this <- (t(theta_w - mu_hat) %*% (theta_w - mu_hat))
      S <- S + S.this
      Cov.this <- cov(var_data)*(nrow(df_writer)-1) 
      Sw <- Sw + Cov.this
    }
  } 
  
  B_hat = S/(n_writers - 1)
  #B_hat = cov(background_data[,1:p])
  if (!is.positive.definite(B_hat)){B_hat = as.matrix(nearPD(B_hat)$mat)}
  
  
  W_hat <- Sw/(n_rows - n_writers)
  U_hat <- W_hat*(nw_hat-p-1)
  
  
  
  log_sd_mat <- sapply(unique(background_data$writer_id), function(w) {
    df_w <- background_data[background_data$writer_id == w, 1:p]
    if (nrow(df_w) <= 3) return(rep(NA_real_, p))  # skip degenerate writers
    log(sqrt(diag(cov(as.matrix(df_w)))))
  })
  
  # Remove degenerate writers (columns with NA)
  log_sd_mat <- log_sd_mat[, colSums(is.na(log_sd_mat)) == 0]
  
  loc <- rowMeans(log_sd_mat)
  sc  <- apply(log_sd_mat, 1, sd)
  
  eta <- 1
  
  return(list(mu_hat,B_hat,nw_hat,U_hat,loc,sc,eta))
}

background_statistics_br <- function(background_data){
  
  p=9
  l = length(unique(background_data$character))
  nw.min = p + 2
  nw_hat = nw.min
  
  a_data = background_data[(background_data$character==1),]
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
    letter_data = as.matrix(unname(background_data[(
      background_data$character==l_id),1:p]))
    if (nrow(letter_data)>2){
      letter_diff = letter_data - matrix(mu_hat[col(letter_data)],ncol = p)
      beta_l = colMeans(letter_diff)
      beta_mu[l_id,] = beta_l
      S = 0
      n_rows <- 0
      for (w in unique(background_data$writer_id)){
        letter_writer = background_data[(
          background_data$character==l_id)& (background_data$writer_id==w),1:p]
        if (nrow(letter_writer)>2){
          n_rows <- n_rows + nrow(letter_writer)
          a_data_writer = background_data[(
            background_data$character==1)& (background_data$writer_id==w),1:p]
          
          mu_hat_writer=matrix(colMeans(a_data_writer),nrow = 1)
          
          letter_diff_writer = letter_writer - matrix(mu_hat_writer[col(letter_writer)],ncol = p)
          
          beta_w = matrix(colMeans(letter_diff_writer), nrow = 1)
          S.this <- (t(beta_w - beta_l) %*% (beta_w - beta_l))
          S <- S + S.this
        }
      }
      B_hat_l = S/(n_rows - 1)
      #B_hat_l = cov(letter_data)
      
      if (!is.positive.definite(B_hat_l)) { B_hat_l <- as.matrix(nearPD(B_hat_l)$mat) }
      beta_cov[,,l_id] = B_hat_l
    }
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
  
  
  log_sd_mat <- sapply(unique(background_data$writer_id), function(w) {
    df_w <- background_data[background_data$writer_id == w, 1:p]
    if (nrow(df_w) <= 3) return(rep(NA_real_, p))  # skip degenerate writers
    log(sqrt(diag(cov(as.matrix(df_w)))))
  })
  
  # Remove degenerate writers (columns with NA)
  log_sd_mat <- log_sd_mat[, colSums(is.na(log_sd_mat)) == 0]
  
  loc <- rowMeans(log_sd_mat)
  sc  <- apply(log_sd_mat, 1, sd)
  
  eta <- 1
  
  return(list(mu_hat,B_hat,beta_mu,beta_cov,nw_hat,U_hat,loc,sc,eta))
}


stan_model_niw <- stan_model(file = "Stan_models/niw.stan", model_name = "niw")
stan_model_nlkj <- stan_model(file = "Stan_models/normal_lkj_model.stan", model_name = "normal_lkj_model")
stan_model_manova_iw <- stan_model(file = "Stan_models/MANOVA_iw_model.stan", model_name = "MANOVA_iw")
stan_model_manova_lkj <- stan_model(file = "Stan_models/MANOVA_lkj_model.stan", model_name = "MANOVA_lkj")



write_xlsx(data.frame(),"Paper_experiments/same_source_results_iter.xlsx")


same_source_def <- function(character_data, w) {
  
  # ── helper: structured log to stdout (goes to outfile on workers) ──────────
  log_msg <- function(level, location, msg) {
    cat(sprintf("[%s][%s][writer=%s] %s\n",
                format(Sys.time(), "%H:%M:%S"), level, w, 
                paste0("[", location, "] ", msg)),
        file = stderr())                 # stderr also captured by outfile=
  }
  
  all_chars <- sort(unique(character_data$character))
  l         <- length(all_chars)
  df_all    <- data.frame()
  
  writer_data_all <- character_data[character_data$writer_id == w, ]
  background_data <- character_data[character_data$writer_id != w, ]
  
  for (iter_for_eval in seq_len(100)) {
    
    bf_rows <- list()
    i       <- 1
    
    # ── sample split ──────────────────────────────────────────────────────────
    sample_size <- nrow(writer_data_all)
    writer_data <- tryCatch({
      writer_data_all %>%
        add_count(character, name = "char_freq") %>%
        slice_sample(n = sample_size, weight_by = char_freq, replace = FALSE)
    }, error = function(e) {
      log_msg("ERROR", paste0("iter=", iter_for_eval, " slice_sample"), conditionMessage(e))
      return(NULL)
    })
    if (is.null(writer_data)) next
    
    questioned_data <- writer_data[seq_len(floor(sample_size / 2)), ]
    suspect_data    <- writer_data[(floor(sample_size / 2) + 1):sample_size, ]
    
    alphabet_map <- setNames(seq_along(all_chars), all_chars)
    questioned_data$character <- as.numeric(alphabet_map[questioned_data$character])
    suspect_data$character    <- as.numeric(alphabet_map[suspect_data$character])
    bg                        <- background_data
    bg$character              <- as.numeric(alphabet_map[bg$character])
    
    # ── per-character NIW models ──────────────────────────────────────────────
    for (ch in unique(bg$character)) {
      q_ch <- questioned_data[questioned_data$character == ch, ]
      s_ch <- suspect_data[suspect_data$character == ch, ]
      
      if (nrow(q_ch) <= 1 || nrow(s_ch) <= 1) {
        log_msg("WARN",
                paste0("iter=", iter_for_eval, " ch=", ch),
                paste0("skipped: nrow(Q)=", nrow(q_ch), " nrow(S)=", nrow(s_ch)))
        next
      }
      
      bg_stats_ch <- tryCatch(
        background_statistics_niw(bg[bg$character == ch, ]),
        error = function(e) {
          log_msg("ERROR", paste0("iter=", iter_for_eval, " background_statistics_niw ch=", ch),
                  conditionMessage(e))
          NULL
        })
      if (is.null(bg_stats_ch)) next
      
      niw_conj <- tryCatch(niw_conjugate(q_ch, s_ch, bg_stats_ch),
                           error = function(e) { log_msg("ERROR", paste0("niw_conjugate ch=", ch), conditionMessage(e)); NA_real_ })
      niw_val  <- tryCatch(normal_iW(q_ch, s_ch, bg_stats_ch),
                           error = function(e) { log_msg("ERROR", paste0("normal_iW ch=", ch),     conditionMessage(e)); NA_real_ })
      nlkj_val <- tryCatch(normal_lkj(q_ch, s_ch, bg_stats_ch),
                           error = function(e) { log_msg("ERROR", paste0("normal_lkj ch=", ch),    conditionMessage(e)); NA_real_ })
      
      bf_rows[[i]] <- data.frame(
        writer    = w,
        character = names(alphabet_map)[ch],
        model     = c("niw_conjugate", "niw", "nlkj"),
        BF        = c(niw_conj, niw_val, nlkj_val)
      )
      i <- i + 1
    }
    
    # ── all-character NIW ─────────────────────────────────────────────────────
    bg_stats_all <- tryCatch(
      background_statistics_niw(bg),
      error = function(e) {
        log_msg("ERROR", paste0("iter=", iter_for_eval, " background_statistics_niw ALL"), conditionMessage(e))
        NULL
      })
    
    if (!is.null(bg_stats_all)) {
      niw_conj_all <- tryCatch(niw_conjugate(questioned_data, suspect_data, bg_stats_all),
                               error = function(e) { log_msg("ERROR", "niw_conjugate ALL", conditionMessage(e)); NA_real_ })
      niw_all_val  <- tryCatch(normal_iW(questioned_data, suspect_data, bg_stats_all),
                               error = function(e) { log_msg("ERROR", "normal_iW ALL",     conditionMessage(e)); NA_real_ })
      nlkj_all     <- tryCatch(normal_lkj(questioned_data, suspect_data, bg_stats_all),
                               error = function(e) { log_msg("ERROR", "normal_lkj ALL",    conditionMessage(e)); NA_real_ })
      bf_rows[[i]] <- data.frame(writer = w, character = "all",
                                 model = c("niw_conjugate", "niw", "nlkj"),
                                 BF    = c(niw_conj_all, niw_all_val, nlkj_all))
      i <- i + 1
    }
    
    # ── all-character MANOVA ──────────────────────────────────────────────────
    bg_stats_br <- tryCatch(
      background_statistics_br(bg),
      error = function(e) {
        log_msg("ERROR", paste0("iter=", iter_for_eval, " background_statistics_br"), conditionMessage(e))
        NULL
      })
    
    if (!is.null(bg_stats_br)) {
      mv_conj <- tryCatch(MANOVA_conjugate(questioned_data, suspect_data, bg_stats_br),
                          error = function(e) { log_msg("ERROR", "MANOVA_conjugate", conditionMessage(e)); NA_real_ })
      mv_iw   <- tryCatch(MANOVA_iw(questioned_data, suspect_data, bg_stats_br),
                          error = function(e) { log_msg("ERROR", "MANOVA_iw",        conditionMessage(e)); NA_real_ })
      mv_lkj  <- tryCatch(MANOVA_LKJ(questioned_data, suspect_data, bg_stats_br),
                          error = function(e) { log_msg("ERROR", "MANOVA_LKJ",       conditionMessage(e)); NA_real_ })
      bf_rows[[i]] <- data.frame(writer = w, character = "all",
                                 model = c("manova_conjugate", "manova_iw", "manova_lkj"),
                                 BF    = c(mv_conj, mv_iw, mv_lkj))
    }
    
    df_all <- rbind(df_all, bind_rows(bf_rows))
    cat(sprintf("[%s][INFO][writer=%s] completed iter %d/%d\n",
                format(Sys.time(), "%H:%M:%S"), w, iter_for_eval, 100))
  }
  
  # ── write results (in-function, so errors here don't lose df_all) ──────────
  tryCatch({
    ssr_i <- read_excel("Paper_experiments/same_source_results_iter.xlsx")
    write_xlsx(rbind(ssr_i, df_all), "Paper_experiments/same_source_results_iter.xlsx")
  }, error = function(e) {
    log_msg("ERROR", "write_xlsx", conditionMessage(e))
    # fallback: save as RDS so results are never lost
    saveRDS(df_all, file = paste0("Paper_experiments/fallback_writer_", w, ".rds"))
  })
  
  return(df_all)
}

df_resutls<- same_source_def(IAM_data,'37')
df_resutls

detectCores()
cl <- makeCluster(5,
                  outfile="C:/Users/Lampis_lab/Desktop/PhD/Handwritten_Loop_characters/handwriting-evidence-evaluation-fourier-based-feature-extraction/Modelling/Paper_experiments/log.txt")

clusterEvalQ(cl, {
  library(dplyr)
  library(Matrix)
})

clusterExport(cl,
              list("background_statistics_br","background_statistics_niw","abind",
                   "is.positive.definite","nearPD","fitdistr","lmvgamma", "inv",
                   "stan_model_niw", "stan_model_nlkj",
                   "stan_model_manova_iw","stan_model_manova_lkj",
                   "assess_BF","sampling","extract","bridge_sampler",
                   "marginal_likelihood_niw_conjugate", "niw_conjugate",
                   "marginal_likelihood_manova_conjugate","MANOVA_conjugate",
                   "normal_iW","normal_lkj","MANOVA_iw","MANOVA_LKJ",
                   "read_excel","write_xlsx"),
              envir=globalenv())

w.list <- sapply(unique(IAM_data$writer_id), list)

system.time({saves = parLapply(cl, w.list,
                               same_source_def,
                               character_data = IAM_data)})

stopCluster(cl)

df_all <- do.call("rbind", saves)

write_xlsx(df_all,"Paper_experiments/same_source_results_ls.xlsx")


ssr <- read_excel("Paper_experiments/same_source_results_ls.xlsx")







#manova_lkj<- ssr[(ssr$model=='manova_lkj'),]
#manova_lkj[order(manova_lkj$BF),][1:10,1]
#"110" "109" "114" "89"  "315" "92"  "94"  "113" "16"  "61"

ssr = as.data.frame(ssr)
unique(ssr$writer)

indx <- apply(ssr, 2, function(x) any(is.na(x) | is.infinite(x)))
colnames(ssr)[indx]
ssr[sapply(ssr, is.infinite)] <- NA
ssr[is.na(ssr)] = 0

ssr_binary = ssr
mean(ssr_binary$BF<0)
as.data.frame(ssr %>% group_by(model, character) %>% summarise(FN = mean(BF<0)))


#plot

ssr$Prior_approach <- ifelse((ssr$model=="niw_conjugate" | ssr$model=="manova_conjugate"),"(1) NIW Conjugate",
                             ifelse((ssr$model=="niw" | ssr$model=="manova_iw"), "(2) NIW Hierarchical", 
                                    "(3) Normal-LogNormal-LKJ"))

ssr$Prior_approach <- factor(ssr$Prior_approach, 
                             levels =c("(1) NIW Conjugate",
                                       "(2) NIW Hierarchical",
                                       "(3) Normal-LogNormal-LKJ"))

library(stringr)
split_data = str_split_fixed(ssr$model, "_", 2)
ssr$model = split_data[,1]

ssr$model <- ifelse(ssr$model=='manova',"MANOVA","Normal") 

ssr$model <- paste0(ssr$model, ' ', ssr$character)
unique(ssr$model)

ssr$model <- factor(ssr$model, levels =c("Normal a", "Normal b", "Normal d", "Normal e",
                                         "Normal g", "Normal o", "Normal p", "Normal all", 
                                         "MANOVA all"))



ssr$BF = as.numeric(ssr$BF)

library(ggplot2)
library(latex2exp)
plot = ggplot(ssr,
              aes(x = Prior_approach, y = BF, fill = Prior_approach)) +
  geom_boxplot() +
  facet_wrap(~model,ncol = 9) +
  scale_y_continuous(name = TeX(r"(\textbf{LogBF})"), limits = c(-50, 200)) +
  scale_x_discrete(labels = c("(1)","(2)","(3)"), name = "Models")+
  #labs(title="Logarithmic Bayes Factors for \n Different Source Comparisons") + 
  #theme(plot.title = element_text(hjust = 0.5))+
  scale_fill_brewer(palette="Set2")+
  geom_hline(yintercept = 0, color = 'brown',lty='dashed')+ 
  labs(fill = "Prior approach") +
  theme(
    #plot.title = element_text(hjust = 0.5),
    strip.text = element_text(size = 12, face = "bold"),
    axis.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 15, face = "bold"),
    legend.position = "bottom"
  ) +
  guides(
    fill = guide_legend(byrow = TRUE),
    colour = guide_legend(override.aes = list(size = 5))
  )

#jpeg("Paper_experiments/plots/ds_boxplot.jpg",width=3920, height=2000, res=300)
plot
#dev.off()
