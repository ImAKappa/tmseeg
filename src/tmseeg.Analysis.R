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
print(filepath)

# = Load data

tmseeg <- load_mat(filepath) 
cat("Loaded dataset with fields:", paste(names(tmseeg), collapse = ","), "\n")

# = Visualize Raw

classes <- list(relax = 1, contract = 3)

class_means <- mean_trials_by_class(tmseeg)
print(dim(class_means$data))

png(filename = file.path(fig_dir, "relax_avg_raw.png"), width = 800, height = 600, res = 100)
plot_trial(filter_trialinfo(class_means, classes$relax))
title(paste("Mean of 'Relax' condition (n=", length(tmseeg$trialinfo[tmseeg$trialinfo == classes$relax]), ")"))
dev.off()

png(filename = file.path(fig_dir, "contract_avg_raw.png"), width = 800, height = 600, res = 100)
plot_trial(filter_trialinfo(class_means, classes$contract))
title(paste("Mean of 'Contract' condition (n=", length(tmseeg$trialinfo[tmseeg$trialinfo == classes$contract]), ")"))
dev.off()

# = Preprocess (I)
# TODO(2026-11-11) Baseline correction (de-meaning) based on [-0.05, -0.001]
tmseeg <- correct_baseline(tmseeg, c(-0.05, -0.001))
# TODO(2026-11-12) Lowpass filter of 35Hz


# = Visualize Preprocessed

classes <- list(relax = 1, contract = 3)

class_means <- mean_trials_by_class(tmseeg)
print(dim(class_means$data))

png(filename = file.path(fig_dir, "relax_avg_preproc.png"), width = 800, height = 600, res = 100)
plot_trial(filter_trialinfo(class_means, classes$relax))
title(paste("Mean of 'Relax' condition (n=", length(tmseeg$trialinfo[tmseeg$trialinfo == classes$relax]), ")"))
dev.off()

png(filename = file.path(fig_dir, "contract_avg_preproc.png"), width = 800, height = 600, res = 100)
plot_trial(filter_trialinfo(class_means, classes$contract))
title(paste("Mean of 'Contract' condition (n=", length(tmseeg$trialinfo[tmseeg$trialinfo == classes$contract]), ")"))
dev.off()

# Preprocessing (I)
tmseeg <- downsample(tmseeg, 250)

# Select window
tmseeg <- filter_time(tmseeg, c(-0.1, 0.6))

# == Single trial
trial_i <- 1
png(filename = file.path(fig_dir, "trial.png"), width = 800, height = 600, res = 100)
plot_trial(filter_trialseq(tmseeg, trial_i))
title(paste("Trial", trial_i))
dev.off()

# Average of each class
# png(filename = file.path(fig_dir, "class_avg.png"), width = 800, height = 600, res = 100)
# contract_avg <- filter_trialinfo(class_means, classes$contract)
# contract_avg <- apply(contract_avg$data, MARGIN = 3, FUN = mean)
# dev.off()

# Plot 'eigen-time-series'

# = Classification

# == Feature engineering

# --- Global Mean Field Amplitude ... verify that I get a plot like https://www.fieldtriptoolbox.org/assets/img/tutorial/tms-eeg/figure23.png

# == Calculate cross-validation