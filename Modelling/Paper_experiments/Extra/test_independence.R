setwd("C:/Users/ltzai/Desktop/PhD/Handwritten_Loop_characters/handwriting-evidence-evaluation-fourier-based-feature-extraction/Modelling")
library(readxl)
library(dplyr)

IAM_data <- read_excel("IAM_fourier_features_dataset/DB_loop_handwriting.xlsx")
IAM_data <- IAM_data[IAM_data$writer_id=='0',]

library(PerformanceAnalytics)
chart.Correlation(IAM_data[,2:9], histogram = TRUE, method = "pearson")

#IAM_data[,2:9] = IAM_data[,2:9]/sqrt(IAM_data$area)

# Compute amplitudes and phases BEFORE any scaling
Ampl1 <- IAM_data$a1^2 + IAM_data$b1^2
Ampl2 <- IAM_data$a2^2 + IAM_data$b2^2
Ampl3 <- IAM_data$a3^2 + IAM_data$b3^2
Ampl4 <- IAM_data$a4^2 + IAM_data$b4^2

Phase1 <- atan2(IAM_data$b1, IAM_data$a1)
Phase2 <- atan2(IAM_data$b2, IAM_data$a2)
Phase3 <- atan2(IAM_data$b3, IAM_data$a3)
Phase4 <- atan2(IAM_data$b4, IAM_data$a4)

fourier <- cbind(Ampl1, Phase1, Ampl2, Phase2,
                 Ampl3, Phase3, Ampl4, Phase4)

library(PerformanceAnalytics)
chart.Correlation(fourier, histogram = TRUE, method = "pearson")

library(circular)
rayleigh.test(circular(Phase1))
rayleigh.test(circular(Phase2))
rayleigh.test(circular(Phase3))
rayleigh.test(circular(Phase4))

library(goftest)
ad.test(Ampl1, "pchisq", df = 2)
ad.test(Ampl2, "pchisq", df = 2)
ad.test(Ampl3, "pchisq", df = 2)
ad.test(Ampl4, "pchisq", df = 2)
