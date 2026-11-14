# Name: Alex Kappen
# SciNet username: tmp_akappen
# Description:
#   Utilities for the four analyses for the assignment

analysis_1_polynomial_fitting <- function(eeg) {
  cat("Analysis 1 [Polynomial fitting]: Finding polynomial approximations to TEP\n")

  # Focus only on TMS-Evoked Potential (the brain response to the single magnetic pulse)
  eeg <- filter_time(eeg, c(0.025, 0.250))

  # Compute class means
  class_means <- mean_trials_by_class(eeg)

  # Select a channel of interest (coi)
  coi <- 17
  # Select condition of interest
  condition <- 1 # Relax
  cat("Channel of Interest: ", coi, ", Condition of interest: ", condition, "\n", sep = "")

  coi_data <- class_means$data[coi,,condition]
  time <- class_means$time
  df <- data.frame(time, coi_data)

  # For plotting
  # plot(df$time, df$coi_data)

  # Fit a nth-degree polymials to the data
  # NOTE: For this specific example, any order higher than 27 will cause an error,
  #    and I'm not entirely sure I understand why. I get this error:
  #    
  #     Error in `stats::poly()`:
  #     ! 'degree' must be less than number of unique points
  order <- 1:27
  models <- lapply(order, function(n) lm(coi_data ~ stats::poly(time, n), data = df))

  # Compare models using AIC
  aic <- sapply(models, AIC)

  df_results <- data.frame(order, aic)
  print(df_results)

  # For plotting
  # for (mi in 1:length(models)) {
  #   pi <- predict(models[[mi]], data.frame(time = time))
  #   lines(time, pi, lwd = 2, col = mi)
  # }

  best_order <- df_results$order[which.min(df_results$aic)]
  cat("The best-fitting order is a ", best_order, "-degree polynomial\n", sep = "")
  # NOTE: After running this a few times, with a few different channels,
  #    I find it interesting that even such a flexible polynomial doesn't quite fit the EEG signal well
  #    Also, the AIC decreases monotically with polynomial order, which is also interesting
}

analysis_2_cov_diff_test <- function(eeg) {
  cat("Analysis 2 [Hypothesis testing]: Are the covariance matrices the same between the Relax and Contract conditions?\n")

  classes <- list(relax = 1, contract = 2)

  # Compute covariance matrices
  covs_relaxed <- cov_by_trial(filter_trialinfo(eeg, classes$relax))
  covs_contract <- cov_by_trial(filter_trialinfo(eeg, classes$contract))

  # Perform permutation test to assess if the covariance of each class comes from the null distribution
}

analysis_3_gmfp_feature_selection <- function(eeg) {
  cat("Analysis 3 [Feature selection]: Which time points are useful for differentiating Relax vs Contract\n")
}

analysis_4_classification <- function(eeg) {
  cat("Analysis 4 [Classification]: \n")
}