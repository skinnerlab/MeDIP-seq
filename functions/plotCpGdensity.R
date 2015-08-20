## Daniel Beck
## Created 8/13/2015
## Modified 

## This function generates a CpG density plot. It requires a dataframe with DMR in rows with CpG density listed in a "cpgDensity" column.

# The CpG density plot is just a histogram. However, Mike would like custom intervals.
# plotCpGdensity<-function(dmrList, xlab="Number of CpG sites per 100bp", ylab="Number of differential DNA methylation regions", ...){
#      if (is.null(dmrList$cpgDensity)) return()
#      if (is.null(nrow(dmrList))) return()
#      if (nrow(dmrList)>0){
# 
#         # start at 0 unless there are no CpGs in the 0 group
#         firstBreak<-0
#         if (length(table(dmrList$cpgDensity<=0.5))<2) firstBreak<-1
# 
#         # all CpG density > 10.5 get lumped together into one group
#         dmrList$cpgDensity[dmrList$cpgDensity>10.5]<-11
# 
#         # generate plot
#         hist(dmrList$cpgDensity, breaks=firstBreak:12-0.5, xaxt="n", 
#              xlab=xlab, ylab=ylab, xaxp=c(firstBreak,10,10-firstBreak), ...)
#         axis(1, at=c(firstBreak:11), labels=c(firstBreak:10, ">10"))
#         
#      } else {
#           return()
#      }        
# }


# Reverting to original CpG density histogram. Mike changed his mind.
plotCpGdensity<-function(dmrList, xlab="Number of CpG sites per 100bp", ylab="Number of differential DNA methylation regions", main="", ...){
     if (is.null(nrow(dmrList))) return()
     if (is.null(dmrList$cpgDensity)) return()
     if (nrow(dmrList)>0){
          
          # all CpG density > 10.5 get lumped together into one group
          dmrList$cpgDensity[dmrList$cpgDensity>10]<-11
          
          # generate plot
          hist(dmrList$cpgDensity, breaks=0:12, xaxt="n", 
               xlab=xlab, ylab=ylab, xaxp=c(0,10,10), main=main, ...)
          axis(1, at=c(0:11), labels=c(0:10, ">10"))
     } else {
          return()
     }        
}
