# Name: Alex Kappen
# SciNet username: tmp_akappen
# Description: 
#   Classifies hand contraction vs relaxation from TMS-EEG data
#   Based on the tutorial: https://www.fieldtriptoolbox.org/tutorial/tms/tms-eeg/

library("tidyverse")
source("src/tmseeg.Utils.R")

# = Set paths
fig_dir <- "figures"

# = Parse CLI arguments
# filepath <- parse_args()
filepath <- "data/data_tms_clean_for_R.mat"

# = Load data
tmseeg <- load_mat(filepath)
cat("Loaded TMS-EEG dataset with fields:", paste(names(tmseeg), collapse = ", "), "\n")

# = Visualize Raw
class_names <- c("Relax", "Contract")

png(filename = file.path(fig_dir, "class_means_raw.png"), width = 800, height = 600, res = 100)
plot_class_means(tmseeg, class_names)
dev.off()
cat("Plotted raw data\n")

# = Preprocess (I)
tmseeg <- correct_baseline(tmseeg, c(-0.1, -0.001))
tmseeg <- filter_lowpass(tmseeg, 35)
cat("Applied low pass filter\n")

# = Visualize Preprocessed (I)
png(filename = file.path(fig_dir, "class_means_preproc.png"), width = 800, height = 600, res = 100)
plot_class_means(tmseeg, class_names)
dev.off()
cat("Plotted preprocessed data\n")

# = Preprocessing (II)
# Reduce dataset size
tmseeg <- downsample(tmseeg, 250)
cat("Downsampled data\n")
# Select time window of interest
tmseeg <- filter_time(tmseeg, c(-0.1, 0.6))
cat("Filtered time\n")

# = Visualized Preprocessed (II)
png(filename = file.path(fig_dir, "class_means_250Hz.png"), width = 800, height = 600, res = 100)
plot_class_means(tmseeg, class_names)
dev.off()
cat("Plotted preprocessed data\n")

# == Single trial
trial_i <- 1
png(filename = file.path(fig_dir, "trial.png"), width = 800, height = 600, res = 100)
plot_trial(filter_trialseq(tmseeg, trial_i))
title(paste("Trial", trial_i))
dev.off()

# Plot 'eigen-time-series'

# = Classification


# == Feature engineering

# --- Global Mean Field Amplitude ... verify that I get a plot like https://www.fieldtriptoolbox.org/assets/img/tutorial/tms-eeg/figure23.png

# == Calculate cross-validation