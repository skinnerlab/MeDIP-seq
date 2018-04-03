## Created 8/13/2015 by Daniel Beck

## This function generates a CpG density plot. It requires a dataframe with DMR in rows
## with CpG density listed in a "cpgDensity" column.

# This version of the function generates a standard histogram with CpG density > 10 lumped
# into a single category.
plotCpGdensity<-function(dmrList, xlab = "Number of CpG sites per 100bp", 
                         ylab = "Number of DMR", 
                         main = "", ...) {
  if (is.null(nrow(dmrList))) return()
  if (is.null(dmrList$cpgDensity)) return()
  if (nrow(dmrList) > 0){
    
    # all CpG density > 10.5 get lumped together into one group
    dmrList$cpgDensity[dmrList$cpgDensity > 9] <- 9.1
    
    # generate plot
    hist(dmrList$cpgDensity, breaks = 0:10, xaxt = "n", xlab = xlab, ylab = ylab, 
         xaxp = c(0, 9, 9), main = main, ...)
    axis(1, at = c(0:5*2), labels = c(0:4*2, "\u226510"), ...)
    axis(1, at = c(1:5*2-1), labels = c(1:5*2-1), ...)
    
  } else {
    return()
  }        
}


# This version of the CpG density histogram has modified intervals, per Mike's request. He 
# changed his mind, but I'm keeping it here just in case.
#plotCpGdensity<-function(dmrList, xlab = "Number of CpG sites per 100bp", 
#                         ylab = "Number of differential DNA methylation regions", ...){
#  if (is.null(dmrList$cpgDensity)) return()
#  if (is.null(nrow(dmrList))) return()
#  if (nrow(dmrList) > 0){
# 
#  # start at 0 unless there are no CpGs in the 0 group
#  firstBreak <- 0
#  if (length(table(dmrList$cpgDensity <= 0.5)) < 2) firstBreak <- 1
# 
#  # all CpG density > 10.5 get lumped together into one group
#  dmrList$cpgDensity[dmrList$cpgDensity > 10.5] <- 11
# 
#  # generate plot
#  hist(dmrList$cpgDensity, breaks = firstBreak:12 - 0.5, xaxt = "n", 
#       xlab = xlab, ylab = ylab, xaxp = c(firstBreak, 10, 10 - firstBreak), ...)
#  axis(1, at = c(firstBreak:11), labels = c(firstBreak:10, ">10"))
#  } else {
#    return()
#  }        
#}










