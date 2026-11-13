# Name: Alex Kappen
# SciNet username: tmp_akappen
# Description:
#   Utilities for digital signal processing (DSP) for tmseeg.Analysis.R

# Popular signal processing library in R
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

correct_baseline <- function(eeg, window) {
  stopifnot(
    "Expected valid eeg data structure" = is_valid_eeg(eeg),
    "Window should be an interval" = length(window) == 2,
    "Window should be monotonic" = window[2] > window[1]
  )

  baseline_mask <- window[1] <= eeg$time & eeg$time <= window[2]

  # Calculate the baseline mean of all channels and trials
  # [channel x trial]
  baseline_mean <- apply(eeg$data[,baseline_mask,], MARGIN = c(1, 3), mean)
  
  # Permute dimensions to allow R to subtract the correct dimensions and recylce over time dimension
  # [channel x trial x time]
  data <- aperm(eeg$data, c(1, 3, 2))
  # Subtract mean from data
  demeaned_data <- data - as.vector(baseline_mean)

  # Permute data back to [channel x time x trial]
  eeg$data <- aperm(demeaned_data, c(1, 3, 2))
  return(eeg)
}