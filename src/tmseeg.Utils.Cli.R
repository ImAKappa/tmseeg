# Name: Alex Kappen
# SciNet username: tmp_akappen
# Description: 
#   Utilities for CLI for tmseeg.Analysis.R

parse_args <- function() {
  help <- "Usage: Rscript tmseeg.Analysis.R [MAT-FILE PATH]"

  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) != 1) {
    cat(help, "\n")
    quit()
  }

  filepath <- args[1]
  return(filepath)
}