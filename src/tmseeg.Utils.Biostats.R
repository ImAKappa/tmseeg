# Name: Alex Kappen
# SciNet username: tmp_akappen
# Description:
#   Utilities for the four analyses for the assignment

analysis_1_polynomial <- function(eeg) {
  cat("Analysis 1 [Polynomial]: Finding polynomial approximations to TEP\n")

  # Compute class means
  class_means <- mean_trials_by_class(tmseeg)

  # Select a channel of interest (coi)
  coi <- 17
  # Select condition of interest
  condition <- 1 # Relax

  coi_data <- class_means$data[coi,,condition]
  time <- class_means$time
  df <- data.frame(time, mean_coi_activity)

  plot(df$time, df$mean_coi_activity)

  # Fit a nth-degree polymials to the data
  order <- 1:10
  models <- lapply(order, function(n) lm(mean_coi_activity ~ stats::poly(time, n), data = df))

  # Compare models using AIC
  aic <- sapply(models, AIC)

  df_results <- data.frame(order, aic)
  print(df_results)
}

analysis_2_cov_diff_test <- function(A, B) {
  cat("Analysis 2 [Hypothesis testing]: Are the covariance matrices the same between the Relax and Contract conditions?\n")

  library(corpcor)

  # Estimate the covariance matrices for each trial
}

analysis_2_gmfp_feature_selection < function(eeg) {
  cat("Analysis 3 [Feature selection]: Which time points are useful for differentiating Relax vs Contract\n")
}

analysis_3_classification <- function(eeg) {
  cat("Analysis 4 [Classification]: \n")
}