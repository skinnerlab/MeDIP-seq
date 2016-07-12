## Created 10/14/2015 by Daniel Beck
## Last modified 4/6/2016

## This function takes a list of base pairs (separated by chromosome) and a DMR table and 
## outputs a table of DMR that overlap with those base pairs.

bpToDmr <- function(dmrTable, bpList) {
  if (is.null(bpList)) return(NULL)
  chrList <- split.data.frame(dmrTable, f = dmrTable$chr)
  chrList <- chrList[which(sapply(chrList, nrow) > 0)]
  commonChrs <- names(which(table(c(unique(dmrTable$chr), names(bpList))) == 2))
  chrList <- chrList[match(commonChrs, names(chrList))]
  bpList <- bpList[match(commonChrs, names(bpList))]
  resultList <- list()
  # Loop over all chromosomes
  if (length(chrList)==1) {
    resultID <- chrList[[1]]$ID[which((chrList[[1]]$start <= i) & (chrList[[1]]$stop >= i))]
    resultList[[1]] <- chrList[[1]][match(unique(unlist(resultID)), chrList[[1]]$ID), ]
  } else {
    for (chr in 1:length(chrList)) {
      resultID <- lapply(bpList[[chr]], function(i) {
        chrList[[chr]]$ID[which((chrList[[chr]]$start <= i) & (chrList[[chr]]$stop >= i))]
        })
      resultList[[chr]] <- chrList[[chr]][match(unique(unlist(resultID)), chrList[[chr]]$ID), ]
    }
  }
  resultTable <- do.call(rbind, resultList)
  return(resultTable)
}
