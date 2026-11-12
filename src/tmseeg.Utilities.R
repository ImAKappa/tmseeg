# Name: Alex Kappen
# SciNet username: tmp_akappen
# Description: 
#   Utilities for tmseeg.Analysis.R

library("hdf5r") # For opening MAT-Files

### ARGUMENT PARSING

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

### EEG DATA STRUCTURE

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
  time_mask <- -0.1 < eeg$time & eeg$time < 0.6
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

### EEG PREPROCESSING

library("gsignal")

downsample <- function(eeg, fs_new) {
  stopifnot(
    "Expected valid eeg data structure" = is_valid_eeg(eeg)
  )

  # Compute decimation factor (assuming constant sampling rate)
  fs_old <- 1 / (eeg$time[2] - eeg$time[1])
  factor <- fs_old / fs_new

  stopifnot(
    "Old sampling frequency (fs_old) must be integer multiple of new sampling frequency (fs_new)" = abs(factor - round(factor)) < 1e-8 & factor > 1
  )
  factor <- as.integer(round(factor))
  cat("Downsampling from ", fs_old, " Hz to ", fs_new, " Hz (factor=", factor, ")\n", sep = "")

  # Reshape data to resample the timepoints for every channel & trial at the same time
  # dims: [channels (1) x time (2) x trials (3)]
  dims <- dim(eeg$data)
  mat_2d <- matrix(aperm(eeg$data, c(2, 1, 3)), nrow = dims[2])
  mat_resampled <- gsignal::resample(mat_2d, p = 1, q = factor)
  arr_3d <- array(mat_resampled, dim = c(nrow(mat_resampled), dims[1], dims[3]))
  eeg$data <- aperm(arr_3d, c(2, 1, 3))
  
  # Update time labels to match downsampled signals
  eeg$time <- seq(eeg$time[1], by = 1/fs_new, length.out = dim(eeg$data)[2])
  
  return(eeg)
}

### Analysis

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

### VISUALIZATION
library("latex2exp")

plot_trial <- function(eeg) {
  cmap <- hcl.colors(length(eeg$label), palette = "Viridis", alpha = 0.4)
  matplot(x = eeg$time, y = t(eeg$data), type = "l", lty = 1,
    xlab = "Time (s)", ylab = latex2exp::TeX("Amplitude ($\\mu$V)"),
    col = cmap)
}