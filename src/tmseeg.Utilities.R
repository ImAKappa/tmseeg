# Name: Alex Kappen
# SciNet username: tmp_akappen
# Description: 
#   Utilities for tmseeg.Analysis.R

library("hdf5r") # For opening MAT-Files

parse_args <- function() {
  help <- "Usage: Rscript tmseeg.Analysis.R [MAT-FILE PATH]"

  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) != 1) {
    cat(help, "\n")
    quit()
  }

  filepath <- args[1]
  return(filepath)
}

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
    trial = file_h5[["trial"]]$read(),
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
  required_names <- c("label", "trial", "time", "trialinfo", "trialseq")
  if (!all(required_names %in% names(eeg))) {
    cat("Missing names: ", paste(setdiff(required_names, names(eeg)), collapse = ", "), "\n")
    return(FALSE)
  }
  
  # --- Check field properties
  # Trial
  if (length(dim(eeg$trial)) != 3) {
    cat("Expected 3 dimensions (channel x time x trial), got", dim(eeg$trial), "\n")
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
  time_mask <- -0.1 < eeg$time & eeg$time < 0.6
  eeg$trial <- eeg$trial[,time_mask,]
  
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
  eeg$trial <- eeg$trial[,,trial_mask]
  
  # Update labels to match data
  eeg$trialinfo <- eeg$trialseq[trial_mask]
  return(eeg)
}

filter_trialinfo <- function(eeg, info) {
  stopifnot(
    "Expected valid eeg data structure" = is_valid_eeg(eeg),
    "Expected info to be numeric" = is.numeric(info),
  )

  # Slice data
  trial_mask <- eeg$trialinfo == info
  eeg$trial <- eeg$trial[,,trial_mask]
  
  # Update labels to match data
  eeg$trialinfo <- eeg$trialseq[trial_mask]
  return(eeg)
}