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
    source("a8/tmseeg.Utils.Viz.R")
    source("a8/tmseeg.Utils.Cov.R")
    source("a8/tmseeg.Utils.Biostats.R")
  })
)

# Load data
filepath <- parse.args()
eeg <- readRDS(filepath)
cat("Loaded TMS-EEG dataset with fields:", paste(names(eeg), collapse = ", "), "\n\n")

# Visualize data
png(filename = file.path("figures", "avg_tep_by_class.png"), width = 800, height = 600, res = 100)
plot.class.means(eeg, c("Relax", "Contract"))
invisible(dev.off())
cat("Plotted preprocessed data\n")

# Analyses
# analysis.1.polynomial.fitting(eeg)
# cat("\n")

# analysis.2.estimating.covariance(eeg)
# cat("\n")

# analysis.3.cov.diff.perm.test(eeg)
# cat("\n")

# analysis.4.classification.with.covariance(eeg)
