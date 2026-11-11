# Name: Alex Kappen
# SciNet username: tmp_akappen
# Description: 
#   Classifies hand contraction vs relaxation from TMS-EEG data

library("tidyverse")
source("src/tmseeg.Utilities.R")

# --- Set paths
fig_dir <- "figures"

# --- Parse CLI arguments
filepath <- parse_args()
print(filepath)
tmseeg <- load_mat(filepath)

# --- Visualize a single trial

# Select window
tmseeg <- filter_time(tmseeg, c(-0.1, 0.6))
print(dim(tmseeg$trial))

# Select trial
tmseeg <- filter_trialseq(tmseeg, 1)
print(dim(tmseeg$trial))



# Plot trial
fig <- png(filename = file.path(fig_dir, "trials.png"), width = 800, height = 600, res = 100)
matplot(t(tmseeg$trial), )

# Plot 'eigen-time-series'
