

# This function makes sure the window stop site is actually on the chromosome. If it isn't, the stop site is changed to the last position on the chromosome and the length is adjusted accordingly.
modifyStop <- function(dmrList, refGenome, maxDMR){
     if (is.null(nrow(dmrList))) return(dmrList)
     if (nrow(dmrList)<maxDMR && nrow(dmrList)>0){
          for (rn in 1:nrow(dmrList)){
               i<-dmrList[rn,]
               i$start<-as.numeric(i$start); i$stop<-as.numeric(i$stop)
               i$stop<-min(i$stop, eval(parse(text=paste("length(refGenome$", "\"", i$chr, "\"", ")", sep=""))))
               i$length<-i$stop-as.numeric(i$start)+1
               dmrList[rn,]<-i
          }
     }
     return(dmrList)
}
