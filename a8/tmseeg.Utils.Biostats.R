# Name: Alex Kappen
# SciNet username: tmp_akappen
# Description:
#   Utilities for the four analyses for the assignment

analysis_1_polynomial_fitting <- function(eeg) {
  cat("=== Analysis 1: Polynomial approximations to TEP ===\n")

  # Focus only on TMS-Evoked Potential (the brain response to the single magnetic pulse)
  eeg <- filter_time(eeg, c(0.025, 0.250))

  # Compute class means
  class_means <- mean_trials_by_class(eeg)

  # In EEG analysis, its common for researchers to focus on one or a few channels of interest
  # Although, I think mutlivarate approaches (shown in later analysis functions) are usually more appropriate

  # Select a channel of interest (coi)
  coi <- 17
  # Select condition of interest
  condition <- 1 # Relax
  cat("Channel of Interest: ", coi, ", Condition of interest: ", condition, "\n", sep = "")

  coi_data <- class_means$data[coi,,condition]
  time <- class_means$time
  df <- data.frame(time, coi_data)

  # Fit a nth-degree polymials to the data
  # NOTE: For this specific channel, any order higher than 27 will cause an error,
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
  # plot(df$time, df$coi_data)


  # For plotting
  # for (mi in 1:length(models)) {
  #   pi <- predict(models[[mi]], data.frame(time = time))
  #   lines(time, pi, lwd = 2, col = mi)
  # }

  best_order <- df_results$order[which.min(df_results$aic)]
  cat("The best-fitting order is a ", best_order, "-degree polynomial\n", sep = "")
  # NOTE: After running this a few times, with a few different channels,
  #    I find it interesting that even such a flexible polynomial doesn't quite fit the TEP well
  #    Also, the AIC decreases monotically with polynomial order, which is also interesting
}

analysis_2_estimating_covariance <- function(eeg) {
  cat("=== Analysis 2: Estimating the covariance ===\n")

  # Compute covariance matrices for trials of each condition
  classes <- list(relax = 1, contract = 3)
  relax <- trial_cov(filter_trialinfo(eeg, classes$relax)$data)
  contract <- trial_cov(filter_trialinfo(eeg, classes$contract)$data)

  cat("Observed difference of mean sample covariance between 'relax' vs 'contract' conditions:\n")
  
  # Euclidean metric
  cat("\tEuclidean distance:", cov_dist_euclidean(relax, contract), "\n")
  # Riemannian metric
  geometries <- c("affine_invariant", "log_euclidean")
  for (g in geometries) {
    cat("\tRiemannian distance (", g, "): ", cov_dist_riemannian(relax, contract, geometry = g), "\n", sep = "")
  }
}

library(future.apply)

analysis_3_cov_diff_perm_test <- function(eeg, nperm = 1000) {
  cat("== Analysis 3: Permutation test for group differences in covariance ===\n")

  # Pre-compute covariance matrices to improve performance
  covs <- trial_cov(eeg$data)

  # Null hypothesis: Both groups come from the same distribution
  #    Therefore, the class labels are interchangeable
  test_stat <- function(idx) {
    group_A <- covs[,,idx, drop = FALSE]
    group_B <- covs[,,-idx, drop = FALSE]
    # Keep track of both distance metrics to compare how choice of metric affects conclusions
    return(array(c(
      cov_dist_euclidean(group_A, group_B), 
      # Log-Euclidean is an approximation to the Affine Invariant geometry which is must faster to compute
      cov_dist_riemannian(group_A, group_B, geometry = "log_euclidean"))
    ))
  }

  # NOTE FOR TA: 
  # I know we didn't cover parallelization in the course, but my covariance calculation functions are too slow!
  # With a single-thread, 10-20 samples was already taking almost a minute ...
  # I usually parallelize my code for my research, though I am more used to MATLAB and Python 
  Sys.setenv(OMP_NUM_THREADS = 1, OPENBLAS_NUM_THREADS = 1)
  n_cores <- max(1L, parallel::detectCores() - 1L)

  cat("Generating null distribution ...")
  plan(multisession, workers = n_cores)

  null_dist <- future_replicate(nperm, {
    perm_idx <- sample.int(length(eeg$trialseq), size = as.integer(length(eeg$trialseq) / 2), replace = FALSE)
    test_stat(perm_idx)
  }, future.seed = TRUE)

  plan(sequential)
  cat(" done\n")

  null_dist <- data.frame(euclidean = null_dist[1,], riemannian = null_dist[2,])
  cat("Null distribution via permutation test:\n")
  print(summary(null_dist))

  # Compute observed statistics
  classes <- list(relax = 1, contract = 3)
  relax <- trial_cov(filter_trialinfo(eeg, classes$relax)$data)
  contract <- trial_cov(filter_trialinfo(eeg, classes$contract)$data)

  cat("\n")
  cat("Observed distance:\n")
  obs_stat <- list(
    euclidean = cov_dist_euclidean(relax, contract), 
    riemannian = cov_dist_riemannian(relax, contract, geometry = "log_euclidean"))

  # Compute p-value (distances have a lower bound at 0, so I used a one-sided test)
  for (metric in c("euclidean", "riemannian")) {
    p_value <- (sum(null_dist[metric] >= obs_stat[metric]) + 1) / (nperm + 1)
    # NOTE FOR TA: I know we didn't cover `sprintf` in class, but I am familiar with it through C, MATLAB, and Python format strings.
    #    Using `cat` alone to format the string was a bit cumbersome
    cat(sprintf("\t%s: %.4f (p = %.5f)\n", metric, obs_stat[metric], p_value))
  }
}

library(caret)
library(recipes)

analysis_4_classification_with_covariance <- function(eeg) {
  cat("Analysis 4: Classification using Tangent-Space Projections of Covariance matrices\n")

  X <- trial_cov(eeg$data)
  y <- factor(eeg$trialinfo)

  # Define resamples
  folds <- createFolds(y, k = 10, list = TRUE, returnTrain = FALSE)

  oof_preds <- factor(rep(NA, length(y)), levels = levels(y))

  fit_ctrl <- trainControl(method = "none")

  for (i in seq_along(folds)) {
    test_idx <- folds[[i]]

    X_train <- X[, , -test_idx, drop = FALSE]
    X_test <- X[, , test_idx, drop = FALSE]

    y_train <- y[-test_idx]

    ref_mean <- mean_cov_riemannian(X_train, geometry = "log_euclidean")
    X_train <- map_cov_to_tangent_space(X_train, ref_mean)
    X_test <- map_cov_to_tangent_space(X_test, ref_mean)

    # To keep things simple, I won't search for the optimal hyperparameter, but just demonstrate how to feed covariance matrices features into an SVM 
    svm_fit <- train(
      x = X_train, y = y_train, method = "svmLinear", trControl = fit_ctrl,
      tuneGrid = data.frame(C = 1)
    )

    oof_preds[test_idx] <- predict(svm_fit, newdata = X_test)
  }

  conf_mat <- confusionMatrix(oof_preds, y)
  print(conf_mat)
}