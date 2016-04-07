## Created -/-/2015 by Daniel Beck
## Last modified 4/6/2016

## This function excludes unwanted chromosomes from the DMR table.

excludeChr <- function(dmrList, exclude) {
  dmrList <- dmrList[match(dmrList$chr, exclude, nomatch =0 ) == 0, ]
  return(dmrList)
}
