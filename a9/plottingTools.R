# Name: Alex Kappen
# SciNet username: tmp_akappen
# Description: 
#   Utilities file for plotting TMS-EEG data
#
#   The data comes from this tutorial https://www.fieldtriptoolbox.org/tutorial/tms/tms-eeg/
#      Specifically: https://download.fieldtriptoolbox.org/tutorial/tms/sp/data_tms_segmented.mat
#      However, I preprocessed that data in MATLAB first to produce a more R-friendly file

load.data <- function(filepath) {
  stopifnot(
    "Expected RDS file" = ends.With(filepath),
    "Data file not found" = file.exists(filepath)
  )
  
  eeg <- readRDS(filepath)
  return(eeg)
}

source("tmseeg.Utils.Epochs.R")
source("tmseeg.Utils.Cov.R")

plot.avg.tep.by.class <- function(eeg) {
  # Get classes
  classes <- unique(eeg$trialinfo)

  # Calculate class means
  class_means <- timelock.avg.by.class(eeg)

  # Plot each class
  par(mfrow = c(2, 1))

  for (i in 1:length(classes)) {
    trials <- slice.trialinfo(class_means, classes[i])
    n_trials <- length(eeg$trialinfo[eeg$trialinfo == classes[i]])
    plot.trial(trials, ylim = c(-7, 7))
    title(paste("Mean '", class_names[i], "' (n=", n_trials, ")", sep = ""))
  }
  
  par(mfrow = c(1, 1))
}

plot.mean.cov.matrice.by.class <- function(eeg) {
  
}