## Created -/-/2015 by Daniel Beck
## Last modified 4/6/2016

## This function calculates CpG density from a DMR table. The table must include "chr", 
## "start", "stop", and "length"

calcCpGdensity <- function(dmrList, maxDMR = 1000) {
  if (is.null(nrow(dmrList)))
    return(dmrList)
  if (nrow(dmrList) < maxDMR && nrow(dmrList) > 0) {
    # If single site, include entire CpG
    singleSite <- which(dmrList$start==dmrList$stop)
    dmrList$stop[singleSite] <- dmrList$stop[singleSite]+1
    dmrList$start[singleSite] <- dmrList$start[singleSite]-1
    
    cpgNum <- apply(dmrList, 1, function(i) {
      dinucleotideFrequency(subseq(
        eval(parse(text = referenceName))[[match(i["chr"],
                                                 seqnames(eval(parse(text = referenceName))))]],
        start = as.numeric(i["start"]), 
        end = as.numeric(i["stop"])))["CG"]
      })
    cpgDensity <- 100 * cpgNum / dmrList$length
    dmrList$cpgNum <- cpgNum
    dmrList$cpgDensity <- cpgDensity
    return(dmrList)
  } else {
    return(dmrList)
  }
}
