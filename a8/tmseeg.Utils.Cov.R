# Name: Alex Kappen
# SciNet username: tmp.akappen
# Description:
#   Utilities for working with covariance matrices

library("corpcor")
library("Riemann")

trial.cov <- function(eeg.array) {
  stopifnot(
    "Expected numeric array" = is.numeric(eeg.array),
    "Epected 3D input" = length(dim(eeg.array)) == 3
  )

  # Get shape
  n.chans <- dim(eeg.array)[1]
  n.trials <- dim(eeg.array)[3]

  # Estimate covariance matrices for each trial
  cov.flat <- apply(eeg.array, 3, function(trial.mat) {
    # Shrinking is a form of regularization (reduces variance in covariance estimate)
    as.matrix(corpcor::cov.shrink(t(trial.mat), verbose = FALSE))
  })
  cov.3d <- array(cov.flat, dim = c(n.chans, n.chans, n.trials))
  return(cov.3d)
}

mean.cov.euclidean <- function(covs) {
  stopifnot(
    "Expected numeric array" = is.numeric(covs),
    "Epected 3D input" = length(dim(covs)) == 3
  )

  # Get shape
  n.chans <- dim(covs)[1]
  n.trials <- dim(covs)[3]
  # Compute covariance matrix for each trial
  cov.flat <- matrix(covs, nrow = n.chans * n.chans, ncol = n.trials)
  # Compute Euclidean mean
  mean.cov <- matrix(rowMeans(cov.flat), nrow = n.chans, ncol = n.chans)
  return(mean.cov)
}

cov.dist.euclidean <- function(covs.A, covs.B) {
  stopifnot(
    "Expected numeric array" = is.numeric(covs.A) & is.numeric(covs.B),
    "Epected 3D input" = length(dim(covs.A)) == 3 & length(dim(covs.B)) == 3,
    "Expected both covariance matrix data to have the same dimension" = dim(covs.A) == dim(covs.B)
  )

  # Compute Euclidean (a.k.a Frobenius norm) of the difference between two matrices
  # While subtraction is not commutative, the norm is commutative (||A - B|| == ||B - A||)
  return(norm(mean.cov.euclidean(covs.A) - mean.cov.euclidean(covs.B), type = "F"))
}

mean.cov.riemannian <- function(covs, geometry = "affine.invariant") {
  stopifnot(
    "Expected numeric array" = is.numeric(covs),
    "Epected 3D input" = length(dim(covs)) == 3
  )

  # Riemannian geometry is curvy(ish)
  # Covariance matrices are typically symmetric positive definite (SPD) matrices
  #   and SPD matrices have a particular geometry to them that can't be captured by Euclidean metrics
  # Check out https://www.youtube.com/watch?v=ZQv5FPohsm8
  cov.riemannian <- wrap.spd(covs)
  mean.cov <- riem.mean(cov.riemannian, geometry = geometry)$mean
  return(mean.cov)
}

cov.dist.riemannian <- function(covs.A, covs.B, geometry = "affine.invariant") {
  stopifnot(
    "Expected numeric array" = is.numeric(covs.A) & is.numeric(covs.B),
    "Epected 3D input" = length(dim(covs.A)) == 3 & length(dim(covs.B)) == 3,
    "Expected both covariance matrix data to have the same dimension" = dim(covs.A) == dim(covs.B)
  )

  # Riemannian distance metrics take into account the geometry of the Riemannian space in which the data lives
  # This space is not flat like a Euclidean space because covariance matrices (usually) are symmetric positive definite (SPD)

  # The geometry is affine invariant, although 
  #
  means <- list(
    mean.cov.riemannian(covs.A, geometry = geometry),
    mean.cov.riemannian(covs.B, geometry = geometry)
  )
  dist.mat <- riem.pdist(wrap.spd(means), geometry = geometry)
  # The resulting distance matrix is a symmetric 2x2 matrix with zeros along the diagonal
  #  so I index element at one of the off-diagonals ([1, 2])
  dist <- dist.mat[1, 2]
  return(dist)
}

map.cov.to.tangent.space <- function(covs, ref.mean) {
  stopifnot(
    "Expected numeric array" = is.numeric(covs) & is.numeric(ref.mean),
    "Epected 3D input for covs" = length(dim(covs)) == 3,
    "Expected 2D input for ref.mean" = length(dim(ref.mean)) == 2
  )

  # To learn more about how or why to perform these specific linear algebra operations
  #    refer to https://hal.science/hal-00681328
  #    and to https://www.seas.ucla.edu/~vandenbe/133B/lectures/symm.pdf

  # Get shape
  C = dim(covs)[1]
  N <- dim(covs)[3]

  # Use a symmetric eigendecomposition to take the inverse square root of a matrix
  e.ref <- eigen(ref.mean, symmetric = TRUE)
  inv.sqrt.ref <- e.ref$vectors %*% ((1 / sqrt(e.ref$values)) * t(e.ref$vectors))

  # Define masks to select different parts of the covariance matrices
  upper.mask <- upper.tri(matrix(0, C, C), diag = TRUE)
  off.diag.mask <- upper.tri(matrix(0, C, C), diag = FALSE)

  # Initialize a matrix buffer for features
  features <- matrix(0, nrow = N, ncol = sum(upper.mask))
  colnames(features) <- paste0("TS.", seq_len(ncol(features)))

  for (i in seq_len(N)) {
    # Whiten each trial matrix,
    #     which decorrelates, normalizes, and centers the data relative to the reference [mean] covariance
    whitened <- inv.sqrt.ref %*% covs[,,i] %*% inv.sqrt.ref

    # Enforce symmetry by taking the average of the whitened and transposed whitened matrices
    whitened <- 0.5 * (whitened + t(whitened))

    # Matrix logarithm, also via symmetric eigendecomposition
    e.w <- eigen(whitened, symmetric = TRUE)
    log.mat <- e.w$vectors %*% (log(e.w$values) * t(e.w$vectors))

    # Weight off-diagonal elements by sqrt(2) to preserve exact Riemannian distance
    log.mat[off.diag.mask] <- log.mat[off.diag.mask] * sqrt(2)

    # Extract features from upper triangle of transformed covariance matrices
    features[i,] <- log.mat[upper.mask]
  }
  return(features)
}
