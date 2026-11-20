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
    "Expected RDS file" = endsWith(filepath, ".rds"),
    "Data file not found" = file.exists(filepath)
  )
  eeg <- readRDS(filepath)
  return(eeg)
}

library(tidyr)
library(dplyr)

trial.to.df <- function(eeg, trial.number) {
  stopifnot(
    "Expected valid eeg" = is.valid.eeg(eeg),
    "Expected trial number to be numeric" = is.numeric(trial.number),
    "Expected trial number > 0" = trial.number > 0
  )

  trial <- eeg$data[,,trial.number]
  trial <- t(trial)
  colnames(trial) <- eeg$label

  trial.class <- eeg$trialinfo[trial.number]

  df <- as.data.frame(trial) %>% 
    mutate(time_ms = eeg$time) %>% 
    pivot_longer(cols = -time_ms, names_to="channel", values_to="voltage_uV") %>%
    mutate(channel = factor(channel, levels = unique(channel))) %>%
    mutate(channel_number = as.numeric(channel)) %>%
    mutate(trial_class = trial.class)
  return(df)
}

source("a9/tmseeg.Utils.Epochs.R")

library(ggplot2)

plot.trial <- function(eeg, trial.number) {
  stopifnot(
    "Expected valid eeg" = is.valid.eeg(eeg),
    "Expected trial number to be numeric" = is.numeric(trial.number),
    "Expected trial number > 0" = trial.number > 0
  )

  # Conver the trial matrices to a dataframe to feed into ggplot2
  trial.df <- trial.to.df(eeg, trial.number)

  # Amplitude range of EEG voltages for this dataset
  # I determined this during my exploratory data analysis
  amplitude.range <- c(-6, 6)

  g <- ggplot(trial.df, aes(x = time_ms, y = voltage_uV, group = channel, color = channel_number)) +
    geom_line(alpha = 0.75, linewidth = 0.5) +
    geom_vline(xintercept = 0.0, linetype = "dashed", color = "darkgrey") +
    annotate(
      "text",
      x = 0.0, y = max(amplitude.range) - 1,
      label = "TMS pulse", angle = 90, vjust = -0.5, color = "darkgrey"
    ) +
    scale_x_continuous(n.breaks = 10) +
    coord_cartesian(ylim = amplitude.range) + # TMS-Evoked potential amplitudes in this dataset within [-6, 6] uV
    labs(x = "time (ms)", y = "voltage (uV)", color = "channel", title = unique(trial.df$trial_class)) +
    scale_color_viridis_c(
      breaks = c(1, 32, length(eeg$label)),
      labels = c("Posterior", "Central", "Anterior")
    ) +
    guides(
      color = guide_colorbar(
        theme = theme(
          legend.title = element_text(hjust = 0.5, vjust = 2.0),
        )
      )
    ) +
    theme_minimal()
  
  return(g)
}

library(patchwork)

plot.avg.tep.by.class <- function(eeg) {
  stopifnot(
    "Expected valid eeg" = is.valid.eeg(eeg),
  )
  # Calculate class means  
  class.means <- timelock.avg.by.class(eeg)

  # Plot the mean trial for each condition
  relax <- plot.trial(class.means, 1)
  contract <- plot.trial(class.means, 2)

  # Combine the plots into a larger plot with a shared colorbar
  combined <- (relax / contract) +
    plot_layout(guides = "collect")

  # Add number of samples per class to each plot title
  counts <- table(eeg$trialinfo)
  for (i in seq_along(combined)) {
    current.title <- combined[[i]]$labels$title

    combined[[i]] <- combined[[i]] + 
      labs(title = paste0("mean '", current.title, "' condition (n =", as.numeric(counts[i]), ")"))
  }

  return(combined)
}

source("a9/tmseeg.Utils.Cov.R")

cov.to.df <- function(eeg, cov.mat) {
  tbl <- as.table(cov.mat)
  colnames(tbl) <- eeg$label
  rownames(tbl) <- eeg$label

  df <- as.data.frame(tbl) %>%
    rename(channel_X = Var1, channel_Y = Var2, covariance = Freq) %>%
    mutate(channel_X = factor(channel_X, levels = eeg$label)) %>%
    mutate(channel_Y = factor(channel_Y, levels = rev(eeg$label)))
  return(df)
}

plot.cov.mat <- function(eeg, cov.mat, max.abs = NULL) {

  cov.df <- cov.to.df(eeg, cov.mat)
  if (is.null(max.abs)) {
    max.abs <- max(abs(cov.mat))
  }

  # Get lowest channel number, highest channel number
  # We don't label the reference channel (always placed last in the list of channels)
  end.channels <- c(eeg$label[1], eeg$label[length(eeg$label) - 1])

  g.cov <- ggplot(cov.df, aes(x = channel_X, y = channel_Y, fill = covariance)) +
    geom_tile(color = "white", linewidth = 0.3) +
    scale_x_discrete(breaks = end.channels) +
    scale_y_discrete(breaks = end.channels) +
    scale_fill_gradient2(
      low = "blue", mid = "white", high = "red", 
      midpoint = 0, limits = c(-max.abs, max.abs),
      name = "Covariance\n(µV²)"
    ) +
    coord_fixed() +
    labs(x = "channel", y = "channel") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid = element_blank(),
      legend.title = element_text(hjust = 0.5, vjust = 2.0)
    )
  return(g.cov)
}

plot.mean.cov.matrice.by.class <- function(eeg) {
  stopifnot(
    "Expected valid eeg" = is.valid.eeg(eeg)
  )

  # Estimate covariance matrices for each class
  cov.relax <- cov.from.trials(slice.trialinfo(eeg, "relax")$data)
  cov.contract <- cov.from.trials(slice.trialinfo(eeg, "contract")$data)

  # Calculate mean Riemannian covariance matrices per class
  cov.mean.relax <- cov.mean.riemannian(cov.relax)
  cov.mean.contract <- cov.mean.riemannian(cov.contract)

  # Compute the maximum value across all covariances to properly set the colour scale
  # But only take the 98th percentile to avoid outliers from super noisy channels!
  max.abs <- as.numeric(max(sapply(list(cov.mean.relax, cov.mean.contract), function(C) quantile(abs(C), probs = 0.98))))

  # Generate covariance matrix heatmaps
  g.relax <- plot.cov.mat(eeg, cov.mean.relax, max.abs) + labs(title = "relax")
  g.contract <- plot.cov.mat(eeg, cov.mean.contract, max.abs) + labs(title = "contract")
  
  # Compute maximum as above
  max.abs <- as.numeric(quantile(abs(cov.mean.relax - cov.mean.contract), probs = 0.98, na.rm = TRUE))
  g.diff <- plot.cov.mat(eeg, cov.mean.relax - cov.mean.contract, max.abs) + labs(title = "relax - contract")

  # Combine subplots into one figure
  combined <- (g.relax + g.contract + g.diff)
  return(combined)
}