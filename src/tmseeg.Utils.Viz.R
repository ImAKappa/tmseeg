# Name: Alex Kappen
# SciNet username: tmp_akappen
# Description:
#   Utilities for visualizing data for tmseeg.Analysis.R

# For rendering math symbols in plots
library("latex2exp")

plot_trial <- function(eeg) {
  cmap <- hcl.colors(length(eeg$label), palette = "Viridis", alpha = 0.4)
  matplot(x = eeg$time, y = t(eeg$data), type = "l", lty = 1,
    xlab = "Time (s)", ylab = latex2exp::TeX("Amplitude ($\\mu$V)"),
    col = cmap)
}