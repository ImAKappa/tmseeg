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
}

args <- parse.args()

# TODO: Check if file exists? Or just crash anyways