## Created -/-/2015 by Daniel Beck
## Last modified 4/6/2016

## This function takes in two tables and merges them using a specified column.

addToTable <- function(x, toAdd, idColumn1 = 3, idColumn2 = 1, discardColumns = NULL) {
  # There are comma separated entrezIDs in the DAVID table
  toAddList <- split.data.frame(toAdd, f = toAdd[, idColumn2])
  nTAL <- unlist(strsplit(names(toAddList), split = ", "))
  toAddList <- rep(toAddList, sapply(strsplit(names(toAddList), split = ", "), length))
  names(toAddList) <- nTAL
  toAddList <- lapply(toAddList, function(i) {
    apply(i, 2, function(j) paste(j, collapse = "_"))
  })
  alignedTAL <- toAddList[match(x[, idColumn1], nTAL, nomatch = NA)]
  alignedTAL[which(sapply(alignedTAL, length) == 0)] <- "NA"
  alignedTAL <- do.call(rbind, alignedTAL)
  if (length(discardColumns) > 0) {
    alignedTAL <- alignedTAL[, -c(discardColumns)]
  }
  combined <- cbind(x, alignedTAL)
  return(combined)
}
