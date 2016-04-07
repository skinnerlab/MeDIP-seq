## Created -/-/2015 by Daniel Beck
## Last modified 4/6/2016

## This function tries to identify a custom category for a gene using homologs

identifyCategory <- function(query, homologs, categories) {
  cg <- homologs$V1[match(categories$symbol, homologs$V4)]
  queryHID <- homologs$V1[match(query, homologs$V4, nomatch=0)]
  present <- match(queryHID, cg, nomatch=0)
  if (length(present)<1) present <- 0
  if (present) {
    return(categories[present,])
  }
  return(NA)
}
