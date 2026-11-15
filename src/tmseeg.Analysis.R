# Name: Alex Kappen
# SciNet username: tmp_akappen
# Description: 
#   Classifies hand contraction vs relaxation from TMS-EEG data
#   Based on the tutorial: https://www.fieldtriptoolbox.org/tutorial/tms/tms-eeg/

library("tidyverse")
source("src/tmseeg.Utils.R")

# Load data
# filepath <- parse_args()
filepath <- file.path("data", "tmseeg.rds")
eeg <- readRDS(filepath)
cat("Loaded TMS-EEG dataset with fields:", paste(names(eeg), collapse = ", "), "\n")

# = Analyses
analysis_1_polynomial_fitting(eeg)
cat("\n")

analysis_2_estimating_covariance(eeg)
cat("\n")

analysis_3_cov_diff_perm_test(eeg)


