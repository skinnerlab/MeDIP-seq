## Created -/-/2015 by Daniel Beck
## Last modified 4/6/2016

## This function matches a gene to an annotationTable using homologs

customAnnotation <- function(query, homologs, annotationTable) {
  queryHID <- homologs$V1[match(query, homologs$V4, nomatch = 0)]
  present <- match(queryHID, annotationTable$homologNumber, nomatch = 0)
  if (length(present) < 1) present <- 0
  if (present) {
    return(annotationTable[present, ])
  }
  return(NA)
}
