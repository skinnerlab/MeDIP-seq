## Created 10/14/2015 by Daniel Beck
## Last modified 4/6/2016

## This function takes as many DMR tables as necessary and returns a list of base pairs that
## overlap in all DMR tables. 

overlappingBP <- function(dmrList) {
  # Split DMR tables by chromosome
  sdmrList <- lapply(dmrList, 
                     function(i) {
                       temp <- split.data.frame(i, f = i$chr)
                       i <- temp[which(sapply(temp, nrow) > 0)]
                     })
  # Find chromosomes shared by all dmrLists
  chrCounts <- table(unlist(lapply(sdmrList, names)))
  commonChr <- names(chrCounts)[chrCounts == length(sdmrList)]
  if (length(commonChr)==0) return(NULL)
  # Loop over common chromosomes and extract all overlapping BP
  overlaps <- list()  # Overlaping BP are stored in this list
  for (chr in 1:length(commonChr)) {
    smallList <- lapply(sdmrList, function(i) i[[which(names(i) == commonChr[chr])]])
    overInt <- cbind(smallList[[1]]$start, smallList[[1]]$stop)
    for (i in 2:length(smallList)) {
      if (!is.null(overInt)){
        overInt <- pairOverlap(start1 = overInt[,1], stop1 = overInt[,2], 
                               start2 = smallList[[i]]$start, stop2 = smallList[[i]]$stop)
      }
    }
    # Expand BP ranges into vectors
    if (!is.null(overInt)) {
      overlaps[[chr]] <- unlist(apply(overInt, 1, function(i) i[1]:i[2]))
    } else {
      overlaps[[chr]] <- NA
    }
  }
  names(overlaps) <- commonChr
  #overlaps<-lapply(overlaps, function(i) if (is.na(i)) i<-NULL)
  # Return list of overlapping basepairs
  return(overlaps)
}
