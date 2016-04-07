## Created -/-/2015 by Daniel Beck
## Last modified 4/6/2016

## This function reformats DMR table to be acceptable to the slidingWindowCluster.R script.
## When a new clustering method is implimented, this function can be eliminated.

formatForCluster <- function(dmrList) {
  cn <- colnames(dmrList)
  cn <- gsub(pattern = "chr", replacement = "Chromosome", cn)
  cn <- gsub(pattern = "start", replacement = "cSTART", cn)
  cn <- gsub(pattern = "stop", replacement = "cSTOP", cn)
  cn <- gsub(pattern = "numSigWin", replacement = "nProbes", cn)
  colnames(dmrList) <- cn
  return(dmrList)
}
