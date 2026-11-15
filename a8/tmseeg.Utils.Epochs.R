# Name: Alex Kappen
# SciNet username: tmp.akappen
# Description:
#   Utilities for working with the epoched TMS-Evoked EEG data

# === VALIDATION ===

is.valid.eeg <- function(eeg) {
  # Check list names
  required.names <- c("label", "data", "time", "trialinfo", "trialseq")
  if (!all(required.names %in% names(eeg))) {
    cat("Missing names: ", paste(setdiff(required.names, names(eeg)), collapse = ", "), "\n")
    return(FALSE)
  }
  
  # --- Check field properties
  # Data
  if (length(dim(eeg$data)) != 3) {
    cat("Expected data to have 3 dimensions (channel x time x trial), got", dim(eeg$data), "\n")
    return(FALSE)
  }

  return(TRUE)
}

# === SLICING ===

slice.time <- function(eeg, window) {
  stopifnot(
    "Expected valid eeg data structure" = is.valid.eeg(eeg),
    "Expected window to be an interval" = length(window) == 2,
    "Expected window to be monotonic" = window[2] > window[1],
    "Expected window to be numeric" = is.numeric(window)
  )
  
  # Slice data
  time.mask <- -0.1 <= eeg$time & eeg$time <= 0.6
  eeg$data <- eeg$data[,time.mask,]
  
  # Update labels to match data
  eeg$time <- eeg$time[time.mask]
  return(eeg)
}

slice.trialseq <- function(eeg, i) {
  stopifnot(
    "Expected valid eeg data structure" = is.valid.eeg(eeg),
    "Expected i to be numeric" = is.numeric(i),
    "Expected i to be greater than 1" = all(i >= 1)
  )

  # Slice data
  trial.mask <- eeg$trialseq == i
  eeg$data <- eeg$data[,,trial.mask]
  
  # Update labels to match data
  # TODO(2026-11-11) What if no i in trialseq? Throw error
  eeg$trialseq <- eeg$trialseq[trial.mask]
  eeg$trialinfo <- eeg$trialinfo[trial.mask]
  return(eeg)
}

slice.trialinfo <- function(eeg, info) {
  stopifnot(
    "Expected valid eeg data structure" = is.valid.eeg(eeg),
    "Expected info to be numeric" = is.numeric(info)
  )

  # Slice data
  trial.mask <- eeg$trialinfo == info
  eeg$data <- eeg$data[,,trial.mask]
  
  # Update labels to match data
  # TODO(2026-11-11) What if no info in trialinfo? Throw error
  eeg$trialinfo <- eeg$trialinfo[trial.mask]
  eeg$trialseq <- eeg$trialseq[trial.mask]
  return(eeg)
}

# == EPOCH SUMMARY ===

timelock.avg.by.class <- function(eeg) {
  stopifnot(
    "Expected valid eeg data structure" = is.valid.eeg(eeg)
  )
  class.means <- apply(eeg$data, MARGIN = c(1, 2), FUN = function(x) tapply(x, eeg$trialinfo, mean))
  class.means <- aperm(class.means, c(2, 3, 1))

  eeg$data <- class.means
  eeg$trialinfo <- unique(eeg$trialinfo)
  eeg$trialseq <- 1:length(unique(eeg$trialinfo))
  return(eeg)
}