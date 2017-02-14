## Created -/-/2015 by Daniel Beck
## Last modified 4/6/2016

## This function tries to identify a custom category for a gene using homologs

identifyCategory <- function(query, homologs, categories, annotationTable=NULL) {
  cg <- homologs$V1[match(categories$symbol, homologs$V4)]
  queryHID <- homologs$V1[match(query, homologs$V4, nomatch=0)]
  cgTab <- categories[which(cg == queryHID),]
  atTab <- annotationTable[which(annotationTable$homologNumber == queryHID),]
  maxrow <- max(nrow(cgTab), nrow(atTab))
  if (nrow(cgTab) < maxrow) {
    cgTab[(nrow(cgTab) + 1):maxrow, ] <- c(NA, NA)
  }
  if (nrow(atTab) < maxrow) {
    atTab[(nrow(atTab) + 1):maxrow, ] <- c(NA, NA, NA, NA, NA, NA)
  }
  rTab <- cbind(cgTab, atTab)
  
  if (length(rTab) < 1) return(NA)
  return(rTab)
}
