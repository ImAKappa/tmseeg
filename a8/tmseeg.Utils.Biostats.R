# Name: Alex Kappen
# SciNet username: tmp.akappen
# Description:
#   Utilities for the four analyses for the assignment

analysis.1.polynomial.fitting <- function(eeg) {
  cat("=== Analysis 1: Polynomial approximation of TMS-Evoked Potential ===\n")

  stopifnot(
    "Expected valid eeg data structure" = is.valid.eeg(eeg)
  )

  # Focus only on TMS-Evoked Potential (the change in EEG amplitude in response to the single-pulse of TMS)
  eeg <- slice.time(eeg, c(0.025, 0.250))
  class.means <- timelock.avg.by.class(eeg)

  # Select a channel of interest (coi)
  #    I picked this channel to be consistent with the tutorial from which I sourced the data
  #    https://www.fieldtriptoolbox.org/tutorial/tms/tms-eeg/
  #    In EEG analysis, its common for researchers to focus on one or a few channels of interest
  #    Although, I think mutlivarate approaches (shown in later analysis functions) are usually more appropriate
  coi <- 17
  cat("Channel of Interest: ", coi, "\n")

  # Create dataframe of voltage values over time
  classes <- c("relax", "contract")
  df <- list(
    relax = data.frame(class.means$time, class.means$data[coi,,1]),
    contract = data.frame(class.means$time, class.means$data[coi,,2])
  )
  colnames(df$relax) <- c("time", "voltage")
  colnames(df$contract) <- c("time", "voltage")

  for (k in classes) {
    # Fit a nth-degree polymials to the data
    # NOTE: For this specific channel, any order higher than 27 will cause an error,
    #    and I'm not entirely sure I understand why
    order <- 1:27
    models <- lapply(order, function(n) lm(voltage ~ stats::poly(time, n), data = df[[k]]))

    # Compare models using Akaike Information Criterion (AIC)
    aic <- sapply(models, AIC)
    df.results <- data.frame(order, aic)
    best.order <- df.results$order[which.min(df.results$aic)]
    cat("For condition '", k, "', the TEP is best fit by a ", best.order, "-degree polynomial\n", sep = "")
  }

  # NOTE: After running this a few times, with a few different channels,
  #    I find it interesting that even such a flexible polynomial doesn't quite fit the TEP well
  #    Also, the AIC decreases monotically with polynomial order, which is also interesting
}

analysis.2.estimating.covariance <- function(eeg) {
  cat("=== Analysis 2: Estimate distance between class mean covariance matrices ===\n")

  stopifnot(
    "Expected valid eeg data structure" = is.valid.eeg(eeg)
  )

  # Compute covariance matrices for trials of each class
  classes <- list(relax = 1, contract = 3)
  relax <- trial.cov(slice.trialinfo(eeg, classes$relax)$data)
  contract <- trial.cov(slice.trialinfo(eeg, classes$contract)$data)

  cat("Observed difference of sample mean covariance between 'relax' and 'contract' conditions:\n")
  
  # Euclidean metric
  cat("\tEuclidean distance:", cov.dist.euclidean(relax, contract), "\n")
  # Riemannian metric
  geometries <- c("affine_invariant", "log_euclidean")
  for (g in geometries) {
    cat("\tRiemannian distance (", g, "): ", cov.dist.riemannian(relax, contract, geometry = g), "\n", sep = "")
  }
}

library("future.apply")

analysis.3.cov.diff.perm.test <- function(eeg, nperm = 1000) {
  cat("== Analysis 3: Permutation test for group differences in covariance ===\n")

  stopifnot(
    "Expected valid eeg data structure" = is.valid.eeg(eeg)
  )

  # Pre-compute covariance matrices to improve performance
  covs <- trial.cov(eeg$data)

  # Null hypothesis: Both groups come from the same distribution
  #    Therefore, the class labels are interchangeable
  test.stat <- function(idx) {
    group.A <- covs[,,idx, drop = FALSE]
    group.B <- covs[,,-idx, drop = FALSE]
    # Keep track of both distance metrics to compare how choice of metric affects conclusions
    return(array(c(
      cov.dist.euclidean(group.A, group.B), 
      # Log-Euclidean is an approximation to the Affine Invariant geometry which is must faster to compute
      cov.dist.riemannian(group.A, group.B, geometry = "log_euclidean"))
    ))
  }

  # NOTE FOR TA: 
  # I know we didn't cover parallelization in the course, but my covariance calculation functions are too slow!
  # With a single-thread, 10-20 samples was already taking almost a minute ...
  # I usually parallelize my code for my research, though I am more used to MATLAB and Python 
  Sys.setenv(OMP.NUM.THREADS = 1, OPENBLAS.NUM.THREADS = 1)
  n.cores <- max(1L, parallel::detectCores() - 1L)

  cat("Generating null distribution ...")
  plan(multisession, workers = n.cores)

  null.dist <- future_replicate(nperm, {
    perm.idx <- sample.int(length(eeg$trialseq), size = as.integer(length(eeg$trialseq) / 2), replace = FALSE)
    test.stat(perm.idx)
  }, future.seed = TRUE)

  plan(sequential)
  cat(" done\n")

  null.dist <- data.frame(euclidean = null.dist[1,], riemannian = null.dist[2,])
  cat("Null distribution via permutation test:\n")
  print(summary(null.dist))

  # Compute observed statistics
  classes <- list(relax = 1, contract = 3)
  relax <- trial.cov(slice.trialinfo(eeg, classes$relax)$data)
  contract <- trial.cov(slice.trialinfo(eeg, classes$contract)$data)

  cat("\n")
  cat("Observed distance:\n")
  obs.stat <- list(
    euclidean = cov.dist.euclidean(relax, contract), 
    riemannian = cov.dist.riemannian(relax, contract, geometry = "log_euclidean"))

  # Compute p-value (distances have a lower bound at 0, so I used a one-sided test)
  for (metric in c("euclidean", "riemannian")) {
    p.value <- (sum(null.dist[metric] >= obs.stat[metric]) + 1) / (nperm + 1)
    # NOTE FOR TA: I know we didn't cover `sprintf` in class, 
    #    but I am familiar with it through C, MATLAB, and Python format strings.
    #    Using `cat` alone to format the string was a bit cumbersome
    cat(sprintf("\t%s: %.4f (p = %.5f)\n", metric, obs.stat[metric], p.value))
  }
}

library("caret")

analysis.4.classification.with.covariance <- function(eeg) {
  stopifnot(
    "Expected valid eeg data structure" = is.valid.eeg(eeg)
  )

  cat("Analysis 4: Classification using Tangent-Space Projections of Covariance matrices\n")

  X <- trial.cov(eeg$data)
  y <- factor(eeg$trialinfo)

  # Define resamples
  folds <- createFolds(y, k = 10, list = TRUE, returnTrain = FALSE)

  # OOF Preds = Out-Of-Fold predictions 
  oof.preds <- factor(rep(NA, length(y)), levels = levels(y))

  fit.ctrl <- trainControl(method = "none")

  for (i in seq_along(folds)) {
    # Generate train and test splits
    test.idx <- folds[[i]]

    X.train <- X[, , -test.idx, drop = FALSE]
    X.test <- X[, , test.idx, drop = FALSE]

    y.train <- y[-test.idx]

    # Compute mean covariance.
    # The training-derived reference should be used for train and test to prevent data leakage
    ref.mean <- mean.cov.riemannian(X.train, geometry = "log_euclidean")
    X.train <- map.cov.to.tangent.space(X.train, ref.mean)
    X.test <- map.cov.to.tangent.space(X.test, ref.mean)

    # NOTE TO TA: 
    # To keep the code shorter, I don't search for the optimal hyperparameter
    # but instead just demonstrate how to feed covariance matrices features into an SVM
    # For an actual application, I would use nested cross-validation to tune the C parameter 
    svm.fit <- train(
      x = X.train, y = y.train, method = "svmLinear", trControl = fit.ctrl,
      tuneGrid = data.frame(C = 1)
    )

    # Save test predictions
    oof.preds[test.idx] <- predict(svm.fit, newdata = X.test)
  }

  # Summarize classification performance
  conf.mat <- confusionMatrix(oof.preds, y)
  print(conf.mat)
}