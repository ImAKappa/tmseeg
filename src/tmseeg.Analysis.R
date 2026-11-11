# Name: Alex Kappen
# SciNet username: tmp_akappen
# Description: 
#   Classifies hand contraction vs relaxation from TMS-EEG data
#   Based on the tutorial: https://www.fieldtriptoolbox.org/tutorial/tms/tms-eeg/

library("tidyverse")
source("src/tmseeg.Utilities.R")

# = Set paths
fig_dir <- "figures"

# = Parse CLI arguments
# filepath <- parse_args()
filepath <- "data/data_tms_clean_for_R.mat"
print(filepath)
tmseeg <- load_mat(filepath) 
cat("Loaded dataset with fields:", paste(names(tmseeg), collapse = ","), "\n")

# = Visualization

# Downsample to reduce visual clutter
tmseeg <- downsample(tmseeg, 250)

# Select window
tmseeg <- filter_time(tmseeg, c(-0.1, 0.6))

# == Single trial

trial_i <- 1
png(filename = file.path(fig_dir, "trial.png"), width = 800, height = 600, res = 100)
plot_trial(filter_trialseq(tmseeg, trial_i))
title(paste("Trial", trial_i))
dev.off()

# == Mean trial per class

class_means <- mean_trials_by_class(tmseeg)
print(dim(class_means$data))

# Plot 'eigen-time-series'
