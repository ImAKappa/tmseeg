# Name: Alex Kappen
# SciNet username: tmp_akappen
# Description:
#   Utilities for working with the TMS-Evoked EEG for tmseeg.Analysis.R

# For opening MAT-Files. Under the hood, MAT-Files (v7.3) are really just HDF5 files
library("hdf5r")

load_mat <- function(filepath) {
  stopifnot(
    "Expected character vector" = is.character(filepath),
    "Expected filepath to end with .mat" = endsWith(filepath, ".mat"),
    "File not found" = file.exists(filepath)
  )

  # MATLAB MAT-Files v7.3 are really just HDF5 files
  file_h5 <- hdf5r::H5File$new(filepath, mode = "r")

  # Repackage relevant data into R data structure
  data = list(
    label = file_h5[["label"]]$read(), 
    data = file_h5[["data"]]$read(),
    time = file_h5[["time"]]$read(),
    trialinfo = file_h5[["trialinfo"]]$read(),
    trialseq = file_h5[["trialseq"]]$read()
  )

  # Clean up any open file handles
  file_h5$close_all()

  return(data)
}

is_valid_eeg <- function(eeg) {
  # Check list names
  required_names <- c("label", "data", "time", "trialinfo", "trialseq")
  if (!all(required_names %in% names(eeg))) {
    cat("Missing names: ", paste(setdiff(required_names, names(eeg)), collapse = ", "), "\n")
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

filter_time <- function(eeg, window) {
  stopifnot(
    "Expected valid eeg data structure" = is_valid_eeg(eeg),
    "Expected window to be an interval" = length(window) == 2,
    "Expected window to be monotonic" = window[2] > window[1],
    "Expected window to be numeric" = is.numeric(window)
  )
  
  # Slice data
  time_mask <- -0.1 <= eeg$time & eeg$time <= 0.6
  eeg$data <- eeg$data[,time_mask,]
  
  # Update labels to match data
  eeg$time <- eeg$time[time_mask]
  return(eeg)
}

filter_trialseq <- function(eeg, i) {
  stopifnot(
    "Expected valid eeg data structure" = is_valid_eeg(eeg),
    "Expected i to be numeric" = is.numeric(i),
    "Expected i to be greater than 1" = all(i >= 1)
  )

  # Slice data
  trial_mask <- eeg$trialseq == i
  eeg$data <- eeg$data[,,trial_mask]
  
  # Update labels to match data
  # TODO(2026-11-11) What if no i in trialseq? Throw error
  eeg$trialseq <- eeg$trialseq[trial_mask]
  eeg$trialinfo <- eeg$trialinfo[trial_mask]
  return(eeg)
}

filter_trialinfo <- function(eeg, info) {
  stopifnot(
    "Expected valid eeg data structure" = is_valid_eeg(eeg),
    "Expected info to be numeric" = is.numeric(info)
  )

  # Slice data
  trial_mask <- eeg$trialinfo == info
  eeg$data <- eeg$data[,,trial_mask]
  
  # Update labels to match data
  # TODO(2026-11-11) What if no info in trialinfo? Throw error
  eeg$trialinfo <- eeg$trialinfo[trial_mask]
  eeg$trialseq <- eeg$trialseq[trial_mask]
  return(eeg)
}