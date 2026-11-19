# Name: Alex Kappen
# SciNet username: tmp_akappen
# Description: 
#   Main driver script fr plotting TMS-EEG data
#   
#   The data comes from this tutorial https://www.fieldtriptoolbox.org/tutorial/tms/tms-eeg/
#      Specifically: https://download.fieldtriptoolbox.org/tutorial/tms/sp/data_tms_segmented.mat
#      However, I preprocessed that data in MATLAB first to produce a more R-friendly file 

source("plottingTools.R")

parse.args <- function() {
  help <- "Usage: Rscript generatePlots.R [DATA FILE PATH] [PLOT INDEX]\n"
  
  args <- commandArgs(trailingOnly = TRUE)

  if (length(args) != 2) {
    cat(help)
    quit()
  }

  data <- args[1]
  plot.idx <- args[2]
  if (!(plot.idx %in% c("plot1", "plot2"))) {
    cat("Expected [PLOT INDEX] to be 'plot1' or 'plot2'\n")
    cat(help)
    quit()
  }
  return(list(data=data, plot.idx=plot.idx))
}

args <- parse.args()

eeg <- load.data(args$data)

if (args$plot.idx == "plot1") {
  plot.avg.tep.by.class(eeg)
} else if (args$plot.idx == "plot2") {
  plot.mean.cov.matrice.by.class(eeg)
}