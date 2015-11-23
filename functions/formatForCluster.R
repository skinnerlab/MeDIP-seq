
formatForCluster <- function(dmrList) {
  cn <- colnames(dmrList)
  cn <- gsub(pattern="chr", replacement="Chromosome", cn)
  cn <- gsub(pattern="start", replacement="cSTART", cn)
  cn <- gsub(pattern="stop", replacement="cSTOP", cn)
  cn <- gsub(pattern="numSigWin", replacement="nProbes", cn)
  colnames(dmrList) <- cn
  return(dmrList)
}
