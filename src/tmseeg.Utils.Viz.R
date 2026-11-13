# Name: Alex Kappen
# SciNet username: tmp_akappen
# Description:
#   Utilities for visualizing data for tmseeg.Analysis.R

# For rendering math symbols in plots
library("latex2exp")

plot_class_means <- function(eeg, class_names) {
  # Get classes
  classes <- unique(eeg$trialinfo)

  # Calculate class means
  class_means <- mean_trials_by_class(tmseeg)

  # Plot each class
  par(mfrow = c(2, 1))

  for (i in 1:length(classes)) {
    trials <- filter_trialinfo(class_means, classes[i])
    n_trials <- length(eeg$trialinfo[eeg$trialinfo == classes[i]])
    plot_trial(trials, ylim = c(-7, 7))
    title(paste("Mean '", class_names[i], "' (n=", n_trials, ")", sep = ""))
  }
  
  par(mfrow = c(1, 1))
}

plot_class_differences <- function(eeg, class_names) {
  # Get classes
  classes <- unique(eeg$trialinfo)

  # Calculate class means
  class_means <- mean_trials_by_class(tmseeg)

  class_diff <- class_means
  class_diff$data <- class_means$data[,,1] - class_means$data[,,2]
  plot_trial(class_diff)
}

plot_trial <- function(eeg, ylim = NULL) {
  cmap <- hcl.colors(length(eeg$label), palette = "Viridis", alpha = 0.4)
  matplot(x = eeg$time, y = t(eeg$data), type = "l", lty = 1,
    xlab = "Time (s)", ylab = latex2exp::TeX("Amplitude ($\\mu$V)"),
    ylim = ylim,
    col = cmap)
}