# Name: Alex Kappen
# SciNet username: tmp_akappen
# Description:
#   Utilities for working with covariance matrices

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

map_cov_to_tangent_space <- function(covs, ref_mean) {
  C = dim(covs)[1]
  N <- dim(covs)[3]

  # Some linear algebra, if we want to take the inverse square root of a matrix, we do so using a symmetric eigendecomposition
  e_ref <- eigen(ref_mean, symmetric = TRUE)
  inv_sqrt_ref <- e_ref$vectors %*% ((1 / sqrt(e_ref$values)) * t(e_ref$vectors))

  upper_mask <- upper.tri(matrix(0, C, C), diag = TRUE)
  off_diag_mask <- upper.tri(matrix(0, C, C), diag = FALSE)

  features <- matrix(0, nrow = N, ncol = sum(upper_mask))
  colnames(features) <- paste0("TS_", seq_len(ncol(features)))

  for (i in seq_len(N)) {
    # Whiten by trial matrix
    whitened <- inv_sqrt_ref %*% covs[,,i] %*% inv_sqrt_ref

    # Enforce symmetry
    whitened <- 0.5 * (whitened + t(whitened))

    # Matrix logarithm
    e_w <- eigen(whitened, symmetric = TRUE)
    log_mat <- e_w$vectors %*% (log(e_w$values) * t(e_w$vectors))

    # Weight off-diagonal elements by sqrt(2) to preserve exact Riemannian distance
    log_mat[off_diag_mask] <- log_mat[off_diag_mask] * sqrt(2)

    # Extract upper triangle features
    features[i,] <- log_mat[upper_mask]
  }
  return(features)
}

isometry_check <- function(covs) {

  # 1. Compute reference Fréchet mean on TRAINING trials only
  train_spd <- wrap.spd(covs)
  ref_mean <- riem.mean(train_spd, geometry = "log_euclidean")$mean

  # 2. Project both training and test sets using the training reference mean
  X_train <- map_cov_to_tangent_space(covs, ref_mean)

  # 3. Verify that Euclidean norm in tangent space equals exact AIRM distance to ref_mean
  ts_norm <- sqrt(sum(X_train[1, ]^2))
  pair_spd <- wrap.spd(array(c(covs[, , 1], ref_mean), dim = c(61, 61, 2)))
  airm_dist <- as.numeric(riem.pdist(pair_spd, geometry = "affine_invariant")[1, 2])

  cat(sprintf("Tangent Space L2 norm: %.6f | Exact AIRM distance: %.6f\n", ts_norm, airm_dist))

}
