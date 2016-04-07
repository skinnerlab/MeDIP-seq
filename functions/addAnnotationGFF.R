## Created -/-/2015 by Daniel Beck
## Last modified 4/6/2016

## This function adds annotation information to a DMR table. The annotation information must
## be in GFF format. There are several GFF types, currently this function only works for one
## of them.

addAnnotationGFF <- function(dmrList, gff, chrPrefix = "", maxDMR = 1000) {
  if (is.null(nrow(dmrList)))
    return(dmrList)
  if (nrow(dmrList) < maxDMR && nrow(dmrList) > 0) {
    # save original chromosome names
    ochr <- dmrList$chr
    # add prefix to chr name if necessary
    dmrList$chr <- paste(chrPrefix, dmrList$chr, sep = "")
    # convert dmrList to GRanges object
    dmrList$start <- as.numeric(dmrList$start)
    dmrList$stop <- as.numeric(dmrList$stop)
    gdmrList <- makeGRangesFromDataFrame(dmrList)
    # find overlaps
    overlaps <- findOverlaps(gdmrList, gff)
    # add annotation to dmrList
    dmrList <- cbind(dmrList, annotation = NA)
    dmrList$annotation[queryHits(overlaps)] <- as.character(gff$group[subjectHits(overlaps)])
    # put original chromosome names back
    dmrList$chr <- ochr
  }
  return(dmrList)
}
