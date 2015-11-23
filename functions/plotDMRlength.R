## Daniel Beck
## Created 9/18/2015
## Modified 

## This function generates a DMR length histogram. It requires a dataframe with DMR in rows

plotDMRlength<-function(dmrList, xlab="DMR length (kbp)", ylab="Number of DMR", main="", axis.fixed=FALSE, ...){
     if (is.null(nrow(dmrList))) return()
     if (is.null(dmrList$length)) return()
     if (nrow(dmrList)>0){
          dmrList$length<-dmrList$length/1000
          # generate plot
          if ((length(which(dmrList$length>10))>0)|(axis.fixed)){
               # all DMR length > 10000 get lumped together into one group
               dmrList$length[dmrList$length>10]<-10.1
               hist(dmrList$length, breaks=0:11, xaxt="n",xlab=xlab, ylab=ylab, xaxp=c(0,10,10), main=main, ...)
               axis(1, at=c(0:11), labels=c(0:10, ">10"))
          } else {
               hist(dmrList$length, xlab=xlab, ylab=ylab, main=main, ...)
          }
     } else {
          return()
     }        
}

## General form of this plot (histogram with maximum)
histMax<-function(x, max = 15, ...){
  if (length(x) > 0) {
    if (length(which(x > max)) > 0){
      # all values > max get changed to max
      x[x > max]<- max + 0.01
      hist(x, breaks=0:(max + 1), xaxt="n", xaxp=c(0, max, max), ...)
      axis(1, at=c(0:(max + 1)), labels=c(0:max, paste(">", max, sep="")))
    } else {
      hist(x, ...)
    }
  } else {
    return()
  }        
}
