
# function to identify and extract DMRs. DMRs are defined by several different constants. "pValueIdentify" is the p-value threshold for calling a DMR. All DMRs include at least one window with a p-value less than pValueIdentify. "pValueExtend" defines a p-value threshold for windows adjacent to the significant window(s) that should be included in the DMR. The "distance" parameter defines how far away the windows can be from each other and still be considered adjacent. The default distance=1 forces the windows to be exactly next to each other.

# I've modified the function to loop over a vector of p values for pValueIdentify and return lists
extractDMR<-function(x, distance=1, pValueIdentify=c(0.05), pValueExtend=0.1, mtcExtend=FALSE, mtcIdentify=FALSE){
     # first extract all windows that meet the pValueExtend threshold
     if (mtcExtend){
          subsetX<-x[which(x$edgeR.adj.p.value<pValueExtend),]
     } else {
          subsetX<-x[which(x$edgeR.p.value<pValueExtend),]
     }
     if (nrow(subsetX)<1) return(list(methList=NULL, methListEtc=NULL))
  
     # merge together adjacent windows
     potentialDMRtable<-MEDIPS.mergeFrames(subsetX, distance=distance)
     methList<-list()
     methListEtc<-list()
     
     # determine which potential DMRs include windows meeting the pValueIdentify threshold
     for (pV in 1:length(pValueIdentify)){
          if (mtcIdentify){
               methList[[pV]]<-x[which(x$edgeR.adj.p.value<pValueIdentify[pV]),]
          } else {
               methList[[pV]]<-x[which(x$edgeR.p.value<pValueIdentify[pV]),]
          }
          
          if (nrow(methList[[pV]])<1){ 
               methList[[pV]]<-NA; methListEtc[[pV]]<-NA
          } else {
               dmrIDs<-MEDIPS.selectROIs(results=methList[[pV]], rois=potentialDMRtable)
               # identify potentialDMRtable index for DMRs
               dmrInd<-which(!is.na(match(potentialDMRtable$ID, dmrIDs$ROI)))
               # extract DMRs
               methList[[pV]]<-potentialDMRtable[dmrInd,]
               # keep track of all p-value information
               methListEtc[[pV]]<-addToMergedResults(allWindows=subsetX, mergedWindows=methList[[pV]])
               # identify windows meeting the pValueIdentify threshold
               if (mtcIdentify){
                    pValueFlags<-lapply(strsplit(as.character(as.data.frame(rbind(methListEtc[[pV]]), stringsAsFactors=F)$edgeR.adj.p.value), split=";"), function(i) i<-as.numeric(i)<pValueIdentify[pV])
               } else {
                    pValueFlags<-lapply(strsplit(as.character(as.data.frame(rbind(methListEtc[[pV]]), stringsAsFactors=F)$edgeR.p.value), split=";"), function(i) i<-as.numeric(i)<pValueIdentify[pV])
               }
               # count number of significant windows in each DMR
               numSigWin<-sapply(pValueFlags, sum)
               # calculate length. add it and numSigWin to dmrList
               methList[[pV]]$start<-as.numeric(methList[[pV]]$start); methList[[pV]]$stop<-as.numeric(methList[[pV]]$stop)
               methList[[pV]]<-cbind(methList[[pV]], length=methList[[pV]]$stop-methList[[pV]]$start+1, numSigWin)
          }
     }
     methList<-lapply(methList, function(i){ if (is.null(nrow(i))){i<-NULL}; i})
     methListEtc<-lapply(methListEtc, function(i){ if (is.null(nrow(i))){i<-NULL}; i})
     
     # return both dmrList and dmrListEtc
     return(list(methList, methListEtc))
}
