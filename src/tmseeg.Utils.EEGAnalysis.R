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
  return(cov_list)
}