# Name: Alex Kappen
# SciNet username: tmp_akappen
# Description:
#   Utilities for working with the TMS-Evoked EEG for tmseeg.Analysis.R

mean_trials_by_class <- function(eeg) {
  stopifnot(
    "Expected valid eeg data structure" = is_valid_eeg(eeg)
  )
  class_means <- apply(eeg$data, MARGIN = c(1, 2), FUN = function(x) tapply(x, eeg$trialinfo, mean))
  class_means <- aperm(class_means, c(2, 3, 1))

  eeg$data <- class_means
  eeg$trialinfo <- unique(eeg$trialinfo)
  eeg$trialseq <- 1:length(unique(eeg$trialinfo))
  return(eeg)
}

mean_class_diff <- function(eeg) {
  class_diff <- mean_trials_by_class(eeg)
  class_diff$data <- class_diff$data[,,1] - class_diff$data[,,2]
  class_diff$trialinfo <- NA
  class_diff$trialseq <- NA
  return(class_diff)
}

#' Global Mean Field Power
#' The standard deviation of amplitude over all channels, for every time point 
gmfp <- function(eeg) {
  stopifnot(
    "Expected valid eeg data structure" = is_valid_eeg(eeg)
  )
  K <- length(eeg$label)
  result <- apply(eeg$data, MARGIN = c(1, 3), FUN = sd)
  return(result)
}

library("corpcor")
library("CovTools")
library("Riemann")

cov_by_trial <- function(eeg) {
  n_trials <- dim(eeg$data)[3]
  cov_list <- vector("list", n_trials)

  for (i in 1:n_trials) {
    # Reshapeing
    trial_i <- eeg$data[,,i]
    trial_i <- t(trial_i)

    # De-meaning
    trial_i <- scale(trial_i, center = TRUE, scale = FALSE)

    # Shrink
    cov_mat <- corpcor::cov.shrink(trial_i, verbose = FALSE)
    cov_list[[i]] <- as.matrix(cov_mat)
  }
  return(simplify2array(cov_list))
}

riemann_perm_test <- function(covs_A, covs_B, nperm = 1000, seed = 42) {
  set.seed(seed)

  n1 <- dim(covs_A)[3]
  n2 <- dim(covs_B)[3]
  all_covs <- c(covs_A, covs_B)
  total_trials <- n1 + n2

  # Observed mean and distance
  # We calculate these with a Riemannian formula instead of Euclidean formulas
  #   because the geometry of covariance matrices is actually not quite flat!
  mean_A <- CovTools::CovMean(covs_A, method = "LERM")
  mean_B <- CovTools::CovMean(covs_B, method = "LERM")

  mean_A <- riem.mean(covs_A, method)
}