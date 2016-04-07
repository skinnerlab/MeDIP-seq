## Created 9/18/2015 by Daniel Beck
## Last modified 4/6/2016

## This function generates a DMR length histogram. It requires a DMR table with a "length" column.

plotDMRlength <- function(dmrList, xlab = "DMR length (kbp)", ylab = "Number of DMR", 
                          main = "", axis.fixed = FALSE, ...) {
  if (is.null(nrow(dmrList))) return()
  if (is.null(dmrList$length)) return()
  if (nrow(dmrList) > 0){
    dmrList$length <- dmrList$length / 1000
    # generate plot
    if ((length(which(dmrList$length > 10)) > 0) | (axis.fixed)) {
      # all DMR length > 10000 get lumped together into one group
      dmrList$length[dmrList$length > 10] <- 10.1
      hist(dmrList$length, breaks = 0:11, xaxt = "n", xlab = xlab, ylab = ylab, 
           xaxp = c(0, 10, 10), main = main, ...)
      axis(1, at = c(0:11), labels = c(0:10, ">10"))
    } else {
      hist(dmrList$length, xlab = xlab, ylab = ylab, main = main, ...)
    }
  } else {
    return()
  }        
}


