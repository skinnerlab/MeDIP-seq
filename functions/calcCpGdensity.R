## Created -/-/2015 by Daniel Beck
## Last modified 4/6/2016

## This function calculates CpG density from a DMR table. The table must include "chr", 
## "start", "stop", and "length"

calcCpGdensity <- function(dmrList, maxDMR = 1000) {
  if (is.null(nrow(dmrList)))
    return(dmrList)
  if (nrow(dmrList) < maxDMR && nrow(dmrList) > 0) {
    cpgNum <- apply(dmrList, 1, function(i) {
      dinucleotideFrequency(subseq(
        eval(parse(text = referenceName))[[match(i["chr"],
                                                 seqnames(eval(parse(text = referenceName))))]],
        start = as.numeric(i["start"]), 
        end = as.numeric(i["stop"])))["CG"]
      })
    cpgDensity <- 100 * cpgNum / dmrList$length
    return(cbind(dmrList, cpgNum, cpgDensity))
  } else {
    return(dmrList)
  }
}
