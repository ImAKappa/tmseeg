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