# Name: Alex Kappen
# SciNet username: tmp_akappen
# Description:
#   Utilities for visualizing data for tmseeg.Analysis.R

# For rendering math symbols in plots
library("latex2exp")

plot.class.means <- function(eeg, class_names) {
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

plot.trial <- function(eeg, ylim = NULL) {
  cmap <- hcl.colors(length(eeg$label), palette = "Viridis", alpha = 0.4)
  matplot(x = eeg$time, y = t(eeg$data), type = "l", lty = 1,
    xlab = "Time (s)", ylab = latex2exp::TeX("Amplitude ($\\mu$V)"),
    ylim = ylim,
    col = cmap)
}