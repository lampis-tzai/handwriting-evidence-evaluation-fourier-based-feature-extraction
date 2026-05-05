setwd("C:/Users/lampis/Desktop/PhD/Handwritten_Loop_characters/handwriting-evidence-evaluation-fourier-based-feature-extraction/Modelling")
library(readxl)
library(dplyr)
library(MASS)
library(matrixcalc)
library(Matrix)
library(abind)
library(matlib)


ssr <- read_excel("Paper_experiments/same_source_results_ls.xlsx")
ssr = as.data.frame(ssr)
dsr <- read_excel("Paper_experiments/Extra/different_source_results_less_than-50.xlsx")
dsr = as.data.frame(dsr)


# --- Data extraction ---
lr_same <- ssr[ssr$model == 'manova_lkj', 'BF']
lr_diff <- dsr[dsr$model == 'manova_lkj', 'BF']

# --- Grid on log_e(BF) scale ---
llr_grid <- seq(-200, 200, length.out = 500)

cum_same <- sapply(llr_grid, function(t) mean(lr_same >= t, na.rm = TRUE))
cum_diff <- sapply(llr_grid, function(t) mean(lr_diff >= t, na.rm = TRUE))

# --- Base R plot ---
plot(llr_grid, cum_diff, type = "l", col = "steelblue", lwd = 2,
     xlab = expression(log[e](BF)),
     ylab = "Fraction of trials",
     main = "Tippett plot")
lines(llr_grid, cum_same, col = "firebrick", lwd = 2)
abline(v = 0, lty = 2, col = "grey50")
legend("right",
       legend = c("Different source", "Same source"),
       col    = c("steelblue", "firebrick"),
       lty = 1, lwd = 2, bty = "n")
graphics::grid()   # explicit namespace avoids conflict with 'llr_grid'

# --- ggplot2 version ---
library(ggplot2)
ssr<-ssr[ssr$character=='all',]
dsr<-dsr[dsr$character=='all',]

models <- unique(ssr$model)

llr_grid <- seq(-150, 150, length.out = 500)

# --- Build long-format data frame for all models ---
df <- do.call(rbind, lapply(models, function(m) {
  
  lr_same <- ssr[ssr$model == m, 'BF']
  lr_diff <- dsr[dsr$model == m, 'BF']
  
  cum_same <- sapply(llr_grid, function(t) mean(lr_same >= t, na.rm = TRUE))
  cum_diff <- sapply(llr_grid, function(t) mean(lr_diff >= t, na.rm = TRUE))
  
  data.frame(
    llr    = rep(llr_grid, times = 2),
    cum    = c(cum_same, cum_diff),
    source = rep(c("Same source", "Different source"), each = length(llr_grid)),
    model  = m
  )
}))

# --- Faceted Tippett plot ---
ggplot(df, aes(x = llr, y = cum, color = source)) +
  geom_line(linewidth = 0.8) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  scale_color_manual(values = c("Same source"      = "firebrick",
                                "Different source" = "steelblue")) +
  facet_wrap(~ model, ncol = 3) +
  labs(x     = expression(log[e](BF)),
       y     = "Fraction of trials",
       title = "Tippett plots by model",
       color = NULL) +
  theme_bw() +
  theme(legend.position  = "bottom",
        strip.background = element_rect(fill = "grey92"),
        strip.text       = element_text(face = "bold"))






library(dplyr)
library(ggplot2)
library(latex2exp)
library(stringr)


ssr <- read_excel("Paper_experiments/same_source_results_ls.xlsx")
ssr = as.data.frame(ssr)
dsr <- read_excel("Paper_experiments/Extra/different_source_results_less_than-50.xlsx")
dsr = as.data.frame(dsr)


## --- Reusable recoding function ---
recode_prior_and_model <- function(df) {
  df <- df %>%
    mutate(
      Prior_approach = ifelse(
        model %in% c("niw_conjugate", "manova_conjugate"), "(1) NIW Conjugate",
        ifelse(
          model %in% c("niw", "manova_iw"), "(2) NIW Hierarchical",
          "(3) Normal-LogNormal-LKJ"
        )
      ),
      Prior_approach = factor(
        Prior_approach,
        levels = c(
          "(1) NIW Conjugate",
          "(2) NIW Hierarchical",
          "(3) Normal-LogNormal-LKJ"
        )
      )
    )
  
  split_data       <- str_split_fixed(df$model, "_", 2)
  df$model         <- split_data[, 1]
  df$model         <- ifelse(df$model == "manova", "MANOVA", "Normal")
  df$model         <- paste0(df$model, " ", df$character)
  df$model         <- factor(df$model,
                             levels = c("Normal a", "Normal b", "Normal d", "Normal e",
                                        "Normal g", "Normal o", "Normal p", "Normal all",
                                        "MANOVA all"))
  df$BF <- as.numeric(df$BF)
  df
}

ssr <- recode_prior_and_model(ssr)
dsr <- recode_prior_and_model(dsr)


## --- Build long-format Tippett data frame ---
## Loop over every (model, Prior_approach) combination
llr_grid <- seq(-150, 150, length.out = 500)

combos <- ssr %>%distinct(model, Prior_approach)

df <- do.call(rbind, lapply(seq_len(nrow(combos)), function(i) {
  m <- combos$model[i]
  p <- combos$Prior_approach[i]
  
  lr_same <- ssr[ssr$model == m & ssr$Prior_approach == p, "BF"]
  lr_diff <- dsr[dsr$model == m & dsr$Prior_approach == p, "BF"]
  
  cum_same <- sapply(llr_grid, function(t) mean(lr_same >= t, na.rm = TRUE))
  cum_diff <- sapply(llr_grid, function(t) mean(lr_diff >= t, na.rm = TRUE))
  
  data.frame(
    llr            = rep(llr_grid, times = 2),
    cum            = c(cum_same, cum_diff),
    source         = rep(c("Same source", "Different source"), each = length(llr_grid)),
    model          = m,
    Prior_approach = p
  )
}))

## --- Plot: prior approach = color, source = linetype, model = facet ---
ggplot(df, aes(x = llr, y = cum,
               color    = Prior_approach,
               linetype = source)) +
  geom_line(linewidth = 0.75) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.4) +
  facet_wrap(~ model, ncol = 9) +
  scale_color_brewer(palette = "Set2") +
  scale_linetype_manual(values = c("Same source"      = "solid",
                                   "Different source" = "dashed")) +
  scale_x_continuous(limits = c(-150, 150)) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(
    x        = TeX(r"($\log_e(\text{BF})$)"),
    y        = "Fraction of trials",
    color    = "Prior approach",
    linetype = "Source"
  ) +
  theme_bw() +
  theme(
    strip.text       = element_text(size = 10, face = "bold"),
    strip.background = element_rect(fill = "grey92"),
    axis.title       = element_text(size = 11, face = "bold"),
    legend.text      = element_text(size = 10),
    legend.title     = element_text(size = 12, face = "bold"),
    legend.position  = "bottom"
  ) +
  guides(
    color    = guide_legend(nrow = 1, byrow = TRUE),
    linetype = guide_legend(nrow = 1)
  )

