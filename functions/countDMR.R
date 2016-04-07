## Created -/-/2015 by Daniel Beck
## Last modified 4/6/2016

## This function counts the number of DMR on a set of chromosomes.

countDMR <- function(dmrList, chr) {
     numDMR <- sum(match(dmrList$chr, chr, nomatch = 0) > 0)
     return(numDMR)
}
