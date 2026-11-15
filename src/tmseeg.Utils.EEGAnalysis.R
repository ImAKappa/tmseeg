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

# === COVARIANCE ===

library("corpcor")
library("CovTools")
library("Riemann")

#' Computes covariance matrix for each trial
trial_cov <- function(eeg_array) {
  n_chans <- dim(eeg_array)[1]
  n_trials <- dim(eeg_array)[3]
  cov_flat <- apply(eeg_array, 3, function(trial_mat) {
    # Shrinking is a form of regularization (reduces variance in covariance estimate)
    as.matrix(corpcor::cov.shrink(t(trial_mat), verbose = FALSE))
  })
  cov_3d <- array(cov_flat, dim = c(n_chans, n_chans, n_trials))
  return(cov_3d)
}

mean_cov_euclidean <- function(covs) {
  n_chans <- dim(covs)[1]
  n_trials <- dim(covs)[3]
  cov_flat <- matrix(covs, nrow = n_chans * n_chans, ncol = n_trials)
  mean_cov <- matrix(rowMeans(cov_flat), nrow = n_chans, ncol = n_chans)
  return(mean_cov)
}

cov_dist_euclidean <- function(covs_A, covs_B) {
  return(norm(mean_cov_euclidean(covs_A) - mean_cov_euclidean(covs_B), type = "F"))
}

mean_cov_riemannian <- function(covs, geometry = "affine_invariant") {
  cov_riemannian <- wrap.spd(covs)
  mean_cov <- riem.mean(cov_riemannian, geometry = geometry)$mean
  return(mean_cov)
}

cov_dist_riemannian <- function(covs_A, covs_B, geometry = "affine_invariant") {
  means <- list(
    mean_cov_riemannian(covs_A, geometry = geometry),
    mean_cov_riemannian(covs_B, geometry = geometry)
  )
  dist_mat <- riem.pdist(wrap.spd(means), geometry = geometry)
  # The resulting distance matrix is a symmetric 2 x 2 matrix with zeros along the diagonal
  dist <- dist_mat[1, 2]
  return(dist)
}
