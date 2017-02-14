## This function extracts all windows that meet CpG number thresholds. It also includes
## options for requiring the same threshold in neighboring windows.

identifyZero <- function(numCpG, chr=NULL, neighborhood=10, max=0) {

  counts <- swCount(numCpG=numCpG, chr=chr, neighborhood=neighborhood)

  final <- numeric(length(numCpG))
  final[which(counts <= max)] <- 1
  return(final)
}
