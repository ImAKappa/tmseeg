# Name: Alex Kappen
# SciNet username: tmp_akappen
# Description: 
#   Classifies hand contraction vs relaxation from Transcranial Magnetic Stimulation (TMS)-Evoked Electroencephalography (EEG) data
#   TMS = Magnetic brain stimulation treatment
#   EEG = recording of voltage fluctuations on the scalp from brain and muscle electrical activity
#   
#   The data comes from this tutorial https://www.fieldtriptoolbox.org/tutorial/tms/tms-eeg/
#      Specifically: https://download.fieldtriptoolbox.org/tutorial/tms/sp/data_tms_segmented.mat
#      However, I preprocessed that data in MATLAB first to produce a more R-friendly file 

suppressWarnings(
  suppressPackageStartupMessages({
    source("a8/tmseeg.Utils.Cli.R")
    source("a8/tmseeg.Utils.Epochs.R")
    source("a8/tmseeg.Utils.Cov.R")
    source("a8/tmseeg.Utils.Biostats.R")
  })
)

# Load data
filepath <- parse.args()
eeg <- readRDS(filepath)
cat("Loaded TMS-EEG dataset with fields:", paste(names(eeg), collapse = ", "), "\n\n")

# Analyses
analysis.1.polynomial.fitting(eeg)
cat("\n")

analysis.2.estimating.covariance(eeg)
cat("\n")

analysis.3.cov.diff.perm.test(eeg)
cat("\n")

analysis.4.classification.with.covariance(eeg)
