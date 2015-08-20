
# function to identify and extract DMRs. DMRs are defined by several different constants. "pValueIdentify" is the p-value threshold for calling a DMR. All DMRs include at least one window with a p-value less than pValueIdentify. "pValueExtend" defines a p-value threshold for windows adjacent to the significant window(s) that should be included in the DMR. The "distance" parameter defines how far away the windows can be from each other and still be considered adjacent. The default distance=1 forces the windows to be exactly next to each other.
extractDMR<-function(x, distance=1, pValueIdentify=0.05, pValueExtend=0.1, mtc=FALSE){
     # first extract all windows that meet the pValueExtend threshold
     if (mtc){
          subsetX<-x[which(x$edgeR.adj.p.value<pValueExtend),]
     } else {
          subsetX<-x[which(x$edgeR.p.value<pValueExtend),]
     }
     # merge together adjacent windows
     potentialDMRtable<-MEDIPS.mergeFrames(subsetX, distance=distance)
     # determine which potential DMRs include windows meeting the pValueIdentify threshold
     if (mtc){
          strictX<-x[which(x$edgeR.adj.p.value<pValueIdentify),]
     } else {
          strictX<-x[which(x$edgeR.p.value<pValueIdentify),]
     }
     
     dmrIDs<-MEDIPS.selectROIs(results=strictX, rois=potentialDMRtable)
     # identify potentialDMRtable index for DMRs
     dmrInd<-which(!is.na(match(potentialDMRtable$ID, dmrIDs$ROI)))
     # extract DMRs
     dmrList<-potentialDMRtable[dmrInd,]
     # keep track of all p-value information
     dmrListEtc<-addToMergedResults(allWindows=subsetX, mergedWindows=dmrList)
     # identify windows meeting the pValueIdentify threshold
     if (mtc){
          pValueFlags<-lapply(strsplit(as.data.frame(dmrListEtc, stringsAsFactors=F)$edgeR.adj.p.value, split=";"), function(i) i<-as.numeric(i)<pValueIdentify)
     } else {
          pValueFlags<-lapply(strsplit(as.data.frame(dmrListEtc, stringsAsFactors=F)$edgeR.p.value, split=";"), function(i) i<-as.numeric(i)<pValueIdentify)
     }
     # count number of significant windows in each DMR
     numSigWin<-sapply(pValueFlags, sum)
     # calculate length. add it and numSigWin to dmrList
     dmrList$start<-as.numeric(dmrList$start); dmrList$stop<-as.numeric(dmrList$stop)
     dmrList<-cbind(dmrList, length=dmrList$stop-dmrList$start+1, numSigWin)
     # return both dmrList and dmrListEtc
     return(list(dmrList, dmrListEtc))
}
