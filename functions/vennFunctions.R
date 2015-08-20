## Daniel Beck
## Created 8/13/2015
## Modified 

## This script contains functions for comparing DMR lists.

## Venn diagram wrappers for using GRanges objects
vennBPfour<-function(dmrList, names, ...){
     totalBP<-sapply(dmrList, function(i) i<-sum(i$stop-i$start+1))
     
     n12<-sum(as.logical(unlist(calcOverlapsBP(dmrList[c(1,2)]))))
     n13<-sum(as.logical(unlist(calcOverlapsBP(dmrList[c(1,3)]))))
     n14<-sum(as.logical(unlist(calcOverlapsBP(dmrList[c(1,4)]))))
     n23<-sum(as.logical(unlist(calcOverlapsBP(dmrList[c(2,3)]))))
     n24<-sum(as.logical(unlist(calcOverlapsBP(dmrList[c(2,4)]))))
     n34<-sum(as.logical(unlist(calcOverlapsBP(dmrList[c(3,4)]))))
     n123<-sum(as.logical(unlist(calcOverlapsBP(dmrList[c(1,2,3)]))))
     n124<-sum(as.logical(unlist(calcOverlapsBP(dmrList[c(1,2,4)]))))
     n134<-sum(as.logical(unlist(calcOverlapsBP(dmrList[c(1,3,4)]))))
     n234<-sum(as.logical(unlist(calcOverlapsBP(dmrList[c(2,3,4)]))))
     n1234<-sum(as.logical(unlist(calcOverlapsBP(dmrList[c(1,2,3,4)]))))
     
     plot.new()
     draw.quad.venn(area1=totalBP[1], area2=totalBP[2], area3=totalBP[3], area4=totalBP[4], n12=n12, n13=n13, n14=n14, n23=n23, n24=n24, n34=n34, n123=n123, n124=n124, n134=n134, n234=n234, n1234=n1234, category = names, lwd = rep(2, 4), lty = rep("solid", 4), col = rep("black", 4), fill = c("skyblue", "pink1", "mediumorchid", "orange"), ...)
}

vennBPthree<-function(dmrList, names, ...){
     totalBP<-sapply(dmrList, function(i) i<-sum(i$stop-i$start+1))
     n12<-sum(as.logical(unlist(calcOverlapsBP(dmrList[c(1,2)]))))
     n13<-sum(as.logical(unlist(calcOverlapsBP(dmrList[c(1,3)]))))
     n23<-sum(as.logical(unlist(calcOverlapsBP(dmrList[c(2,3)]))))
     n123<-sum(as.logical(unlist(calcOverlapsBP(dmrList[c(1,2,3)]))))
  
     plot.new()
     draw.triple.venn(area1=totalBP[1], area2=totalBP[2], area3=totalBP[3], n12=n12, n13=n13, n23=n23, n123=n123, category=names, col=rep("black", 3), fill=c("skyblue", "pink1", "orange"), ...)
}


vennBPtwo<-function(dmrList, names, ...){
     totalBP<-sapply(dmrList, function(i) i<-sum(i$stop-i$start+1))
     n12<-sum(as.logical(unlist(calcOverlapsBP(dmrList[c(1,2)]))))
     
     plot.new()
     draw.pairwise.venn(area1=totalBP[1], area2=totalBP[2], cross.area=n12, category=names, col=rep("black", 2), fill=c("skyblue", "pink1"), ...)
}



calcOverlaps<-function(a,b){
     # first count number of non-overlapping DMRs
     uniqueA<-sum(!countOverlaps(a,b))
     uniqueB<-sum(!countOverlaps(b,a))
     # next find all overlapping DMRs
     allAB<-findOverlaps(a,b)
     qh<-queryHits(allAB)
     sh<-subjectHits(allAB)
     # find complex overlaps (overlapping DMRs in A that do not uniquely overlap a single DMR in B)
     shc<-!is.na(match(sh, names(table(sh)[table(sh)>1])))
     qhc<-!is.na(match(qh, names(table(qh)[table(qh)>1])))
     complexAB<-allAB[shc|qhc]
     # subtract complex overlaps from all overlapping to find simple overlaps
     simpleAB<-length(allAB)-length(complexAB)
     # calculate number of DMR involved in complex overlaps for each sample
     ncA<-length(unique(queryHits(complexAB)))
     ncB<-length(unique(subjectHits(complexAB)))
     results<-c(uniqueA, ncA, simpleAB, ncB, uniqueB)
     names(results)<-c("uniqueA", "complexA", "simpleAB", "complexB", "uniqueB")
     return(results)
}

vennBarSet<-function(gdmrList, names, col, scale=T, ...){
     comparisons<-matrix(nrow=length(gdmrList)^2, ncol=2)
     comparisons[,1]<-c(sapply(1:length(gdmrList), rep, length(gdmrList)))
     comparisons[,2]<-rep(1:length(gdmrList), length(gdmrList))
     comparisons<-comparisons[!(comparisons[,1]==comparisons[,2]),]
     overlaps<-apply(comparisons, 1, function(i){
          calcOverlaps(gdmrList[[i[1]]], gdmrList[[i[2]]])
     })
     countOverlaps<-overlaps
     if (scale){
          overlaps<-apply(overlaps, 2, function(i) i<-i/sum(i))
     }
     #colMat<-colorMatrix<-apply(comparisons, 1, function(i){
     colMat<-colorRampPalette(c(col[1], col[2]))(5)
     #})
     par(mar=c(3.1,4.1,2.1,4.1))
     a<-barplot(overlaps, beside=F, horiz=T, col=(colMat), xpd=F, space=c(0,rep(c(rep(0,length(gdmrList)-2), 0.5),length(gdmrList)-1)),...)
     mtext(side=2, at=a, c(names[comparisons[,1]]), las=1, line=0.1)
     mtext(side=4, at=a, c(names[comparisons[,2]]), las=1, line=0.1)
     text(apply(overlaps, 2, function(i) i<-i[1]/2), a, countOverlaps[1,], cex=0.7)
     text(apply(overlaps, 2, function(i) i<-i[1]+i[2]/2), a, countOverlaps[2,], cex=0.7)
     text(apply(overlaps, 2, function(i) i<-i[1]+i[2]+i[3]/2), a, countOverlaps[3,],cex=0.7)
     text(apply(overlaps, 2, function(i) i<-i[1]+i[2]+i[3]+i[4]/2), a, countOverlaps[4,],cex=0.7)
     text(apply(overlaps, 2, function(i) i<-i[1]+i[2]+i[3]+i[4]+i[5]/2), a, countOverlaps[5,], cex=0.7)
}

# This function takes a list of DMR tables and calculates the number of BP that overlap in all tables
calcOverlapsBP<-function(dmrList){
     # split DMR tables by chromosome
     sdmrList<-lapply(dmrList, function(i) split.data.frame(i, f=i$chr))
     # find chromosomes shared by all dmrLists
     commonChr<-names(sdmrList[[1]])
     for (i in 2:length(dmrList)){
          commonChr<-commonChr[!is.na(match(commonChr, names(sdmrList[[i]])))]
     }
     
     overlapBPlist<-list()
     for (chr in 1:length(commonChr)){
          chrList<-lapply(sdmrList, function(i) i<-i[[which(names(i)==commonChr[chr])]])

          pdmrBP<-lapply(chrList, function(i) lapply(1:nrow(i), function(j) j<-i$start[j]:i$stop[j]))
          dmrBP<-lapply(pdmrBP, function(i) unlist(i))
          
          overlapBP<-dmrBP[[2]][match(dmrBP[[1]], dmrBP[[2]], nomatch=0)]
          if (length(dmrBP)>2){
               for (i in 3:length(dmrBP)){
                    overlapBP<-overlapBP[match(dmrBP[[i]], overlapBP, nomatch=0)]
               }
          }
          overlapBPlist[[chr]]<-overlapBP
     }
     names(overlapBPlist)<-commonChr

     return(overlapBPlist)
}

# this function is based on draw.pairwise.venn() but is customized for DMR overlaps
vennDMRthree<-function(dmrList, names, ...){
     # split DMR tables by chromosome
     sdmrList<-lapply(dmrList, function(i) split.data.frame(i, f=i$chr))
     # find chromosomes shared by all dmrLists
     commonChr<-names(sdmrList[[1]])
     for (i in 2:length(dmrList)){
          commonChr<-commonChr[!is.na(match(commonChr, names(sdmrList[[i]])))]
     }

     uniqueAbp<-NULL; uniqueBbp<-NULL; uniqueCbp<-NULL
     uniqueAdmr<-NULL; uniqueBdmr<-NULL; uniqueCdmr<-NULL    
     overlapABbp<-NULL; overlapACbp<-NULL; overlapBCbp<-NULL; overlapABCbp<-NULL;
     overlapABdmr<-NULL; overlapACdmr<-NULL; overlapBCdmr<-NULL; overlapABCdmr<-NULL;
     
     # loop over all chromosomes
     for (chr in 1:length(commonChr)){
          # extract all dmr of given chromosome
          chrList<-lapply(sdmrList, function(i) i<-i[[which(names(i)==commonChr[chr])]])
          # expand ranges to individual bp vector
          pdmrBP<-lapply(chrList, function(i) lapply(1:nrow(i), function(j) j<-i$start[j]:i$stop[j]))
          dmrBP<-lapply(pdmrBP, function(i) unlist(i))
          # lump all bp together then table them to determine which ones overlap
          bps<-table(unlist(dmrBP))
          bpsOne<-as.numeric(names(bps[which(bps==1)]))
          bpsTwo<-as.numeric(names(bps[which(bps==2)]))
          bpsThree<-as.numeric(names(bps[which(bps==3)]))
          
          # extract double overlaps for each analysis
          aX<-dmrBP[[1]][match(bpsTwo, dmrBP[[1]], nomatch=0)]
          bX<-dmrBP[[2]][match(bpsTwo, dmrBP[[2]], nomatch=0)]
          cX<-dmrBP[[3]][match(bpsTwo, dmrBP[[3]], nomatch=0)]
          
          # determine double overlap values for venn diagram
          uBPa<-dmrBP[[1]][match(bpsOne, dmrBP[[1]], nomatch=0)]
          uBPb<-dmrBP[[2]][match(bpsOne, dmrBP[[2]], nomatch=0)]
          uBPc<-dmrBP[[3]][match(bpsOne, dmrBP[[3]], nomatch=0)]
          overlapBPab<-as.numeric(names(table(c(aX,bX))[which(table(c(aX,bX))==2)]))
          overlapBPac<-as.numeric(names(table(c(aX,cX))[which(table(c(aX,cX))==2)]))
          overlapBPbc<-as.numeric(names(table(c(bX,cX))[which(table(c(bX,cX))==2)]))
          overlapBPabc<-bpsThree

          uDMRa<-sum(unlist(lapply(pdmrBP, function(i) i<-sapply(i, function(j) length(which(!is.na(match(j,uBPa))))>0))))
          uDMRb<-sum(unlist(lapply(pdmrBP, function(i) i<-sapply(i, function(j) length(which(!is.na(match(j,uBPb))))>0))))
          uDMRc<-sum(unlist(lapply(pdmrBP, function(i) i<-sapply(i, function(j) length(which(!is.na(match(j,uBPc))))>0))))
          overlapDMRab<-sum(unlist(lapply(pdmrBP, function(i) i<-sapply(i, function(j) length(which(!is.na(match(j,overlapBPab))))>0))))
          overlapDMRac<-sum(unlist(lapply(pdmrBP, function(i) i<-sapply(i, function(j) length(which(!is.na(match(j,overlapBPac))))>0))))
          overlapDMRbc<-sum(unlist(lapply(pdmrBP, function(i) i<-sapply(i, function(j) length(which(!is.na(match(j,overlapBPbc))))>0))))
          overlapDMRabc<-sum(unlist(lapply(pdmrBP, function(i) i<-sapply(i, function(j) length(which(!is.na(match(j,overlapBPabc))))>0))))
          
          # store results in arrays
          uniqueAbp<-c(uniqueAbp, length(uBPa))
          uniqueBbp<-c(uniqueBbp, length(uBPb))
          uniqueCbp<-c(uniqueCbp, length(uBPc))
          overlapABbp<-c(overlapABbp, length(overlapBPab))
          overlapACbp<-c(overlapACbp, length(overlapBPac))
          overlapBCbp<-c(overlapBCbp, length(overlapBPbc))
          overlapABCbp<-c(overlapABCbp, length(overlapBPabc))
          
          uniqueAdmr<-c(uniqueAdmr, uDMRa)
          uniqueBdmr<-c(uniqueBdmr, uDMRb)
          uniqueCdmr<-c(uniqueCdmr, uDMRc)
          overlapABdmr<-c(overlapABdmr, overlapDMRab)
          overlapACdmr<-c(overlapACdmr, overlapDMRac)
          overlapBCdmr<-c(overlapBCdmr, overlapDMRbc)
          overlapABCdmr<-c(overlapABCdmr, overlapDMRabc)
          
     }

     area1=sum(uniqueAdmr)+sum(overlapABdmr)+sum(overlapACdmr)+sum(overlapABCdmr)
     area2=sum(uniqueBdmr)+sum(overlapABdmr)+sum(overlapBCdmr)+sum(overlapABCdmr)
     area3=sum(uniqueCdmr)+sum(overlapACdmr)+sum(overlapBCdmr)+sum(overlapABCdmr)
     n12=sum(overlapABdmr)+sum(overlapABCdmr)
     n13=sum(overlapACdmr)+sum(overlapABCdmr)
     n23=sum(overlapBCdmr)+sum(overlapABCdmr)
     n123=sum(overlapABCdmr)
     
     plot.new()
     draw.triple.venn(area1=area1, area2=area2, area3=area3, n12=n12, n13=n13, n23=n23, n123=n123, category=names, col=rep("black", 3), fill=c("skyblue", "pink1", "orange"), ...)
     
}


# this function is based on draw.pairwise.venn() but is customized for DMR overlaps. This variation allows only a single grouping for each DMR. DMRs are added to the categories in the order of 3 overlaps, 2 overlaps, unique. 
vennDMRthree2<-function(dmrList, names, scaled=F, ...){
     # split DMR tables by chromosome
     sdmrList<-lapply(dmrList, function(i) split.data.frame(i, f=i$chr))
     # find chromosomes shared by all dmrLists
     commonChr<-names(sdmrList[[1]])
     for (i in 2:length(dmrList)){
          commonChr<-commonChr[!is.na(match(commonChr, names(sdmrList[[i]])))]
     }
     
     uniqueAbp<-NULL; uniqueBbp<-NULL; uniqueCbp<-NULL
     uniqueAdmr<-NULL; uniqueBdmr<-NULL; uniqueCdmr<-NULL    
     overlapABbp<-NULL; overlapACbp<-NULL; overlapBCbp<-NULL; overlapABCbp<-NULL;
     overlapABdmr<-NULL; overlapACdmr<-NULL; overlapBCdmr<-NULL; overlapABCdmr<-NULL;
     
     # loop over all chromosomes
     for (chr in 1:length(commonChr)){
          # extract all dmr of given chromosome
          chrList<-lapply(sdmrList, function(i) i<-i[[which(names(i)==commonChr[chr])]])
          # expand ranges to individual bp vector
          pdmrBP<-lapply(chrList, function(i) lapply(1:nrow(i), function(j) j<-i$start[j]:i$stop[j]))
          dmrBP<-lapply(pdmrBP, function(i) unlist(i))
          # lump all bp together then table them to determine which ones overlap
          bps<-table(unlist(dmrBP))
          bpsOne<-as.numeric(names(bps[which(bps==1)]))
          bpsTwo<-as.numeric(names(bps[which(bps==2)]))
          bpsThree<-as.numeric(names(bps[which(bps==3)]))
          
          # extract double overlaps for each analysis
          aX<-dmrBP[[1]][match(bpsTwo, dmrBP[[1]], nomatch=0)]
          bX<-dmrBP[[2]][match(bpsTwo, dmrBP[[2]], nomatch=0)]
          cX<-dmrBP[[3]][match(bpsTwo, dmrBP[[3]], nomatch=0)]
          
          # determine double overlap values for venn diagram
          uBPa<-dmrBP[[1]][match(bpsOne, dmrBP[[1]], nomatch=0)]
          uBPb<-dmrBP[[2]][match(bpsOne, dmrBP[[2]], nomatch=0)]
          uBPc<-dmrBP[[3]][match(bpsOne, dmrBP[[3]], nomatch=0)]
          overlapBPab<-as.numeric(names(table(c(aX,bX))[which(table(c(aX,bX))==2)]))
          overlapBPac<-as.numeric(names(table(c(aX,cX))[which(table(c(aX,cX))==2)]))
          overlapBPbc<-as.numeric(names(table(c(bX,cX))[which(table(c(bX,cX))==2)]))
          overlapBPabc<-bpsThree
          
          uDMRa<-lapply(pdmrBP, function(i) i<-sapply(i, function(j) length(which(!is.na(match(j,uBPa))))>0))
          uDMRb<-lapply(pdmrBP, function(i) i<-sapply(i, function(j) length(which(!is.na(match(j,uBPb))))>0))
          uDMRc<-lapply(pdmrBP, function(i) i<-sapply(i, function(j) length(which(!is.na(match(j,uBPc))))>0))
          overlapDMRab<-lapply(pdmrBP, function(i) i<-sapply(i, function(j) length(which(!is.na(match(j,overlapBPab))))>0))
          overlapDMRac<-lapply(pdmrBP, function(i) i<-sapply(i, function(j) length(which(!is.na(match(j,overlapBPac))))>0))
          overlapDMRbc<-lapply(pdmrBP, function(i) i<-sapply(i, function(j) length(which(!is.na(match(j,overlapBPbc))))>0))
          overlapDMRabc<-lapply(pdmrBP, function(i) i<-sapply(i, function(j) length(which(!is.na(match(j,overlapBPabc))))>0))
          for (i in 1:length(chrList)){
               overlapDMRab[[i]][overlapDMRab[[i]]&overlapDMRabc[[i]]]<-FALSE
               overlapDMRac[[i]][overlapDMRac[[i]]&overlapDMRabc[[i]]]<-FALSE
               overlapDMRbc[[i]][overlapDMRbc[[i]]&overlapDMRabc[[i]]]<-FALSE
               uDMRa[[i]][uDMRa[[i]]&(overlapDMRab[[i]]|overlapDMRac[[i]]|overlapDMRabc[[i]])]<-FALSE
               uDMRb[[i]][uDMRb[[i]]&(overlapDMRab[[i]]|overlapDMRbc[[i]]|overlapDMRabc[[i]])]<-FALSE
               uDMRc[[i]][uDMRc[[i]]&(overlapDMRac[[i]]|overlapDMRbc[[i]]|overlapDMRabc[[i]])]<-FALSE
          }
          uDMRa<-sum(unlist(uDMRa))
          uDMRb<-sum(unlist(uDMRb))
          uDMRc<-sum(unlist(uDMRc))
          overlapDMRab<-sum(unlist(overlapDMRab))
          overlapDMRac<-sum(unlist(overlapDMRac))
          overlapDMRbc<-sum(unlist(overlapDMRbc))
          overlapDMRabc<-sum(unlist(overlapDMRabc))
          
          # store results in arrays
          uniqueAbp<-c(uniqueAbp, length(uBPa))
          uniqueBbp<-c(uniqueBbp, length(uBPb))
          uniqueCbp<-c(uniqueCbp, length(uBPc))
          overlapABbp<-c(overlapABbp, length(overlapBPab))
          overlapACbp<-c(overlapACbp, length(overlapBPac))
          overlapBCbp<-c(overlapBCbp, length(overlapBPbc))
          overlapABCbp<-c(overlapABCbp, length(overlapBPabc))
          
          uniqueAdmr<-c(uniqueAdmr, uDMRa)
          uniqueBdmr<-c(uniqueBdmr, uDMRb)
          uniqueCdmr<-c(uniqueCdmr, uDMRc)
          overlapABdmr<-c(overlapABdmr, overlapDMRab)
          overlapACdmr<-c(overlapACdmr, overlapDMRac)
          overlapBCdmr<-c(overlapBCdmr, overlapDMRbc)
          overlapABCdmr<-c(overlapABCdmr, overlapDMRabc)
          
          
     }
     SuniqueAdmr<-sum(uniqueAdmr)
     SuniqueBdmr<-sum(uniqueBdmr)
     SuniqueCdmr<-sum(uniqueCdmr)
     SoverlapABdmr<-sum(overlapABdmr)
     SoverlapACdmr<-sum(overlapACdmr)
     SoverlapBCdmr<-sum(overlapBCdmr)
     SoverlapABCdmr<-sum(overlapABCdmr)
     
     if (scaled){
          SoverlapABdmr<-SoverlapABdmr/2
          SoverlapACdmr<-SoverlapACdmr/2
          SoverlapBCdmr<-SoverlapBCdmr/2
          SoverlapABCdmr<-SoverlapABCdmr/3
     }
     area1=SuniqueAdmr+SoverlapABdmr+SoverlapACdmr+SoverlapABCdmr
     area2=SuniqueBdmr+SoverlapABdmr+SoverlapBCdmr+SoverlapABCdmr
     area3=SuniqueCdmr+SoverlapACdmr+SoverlapBCdmr+SoverlapABCdmr
     n12=SoverlapABdmr+SoverlapABCdmr
     n13=SoverlapACdmr+SoverlapABCdmr
     n23=SoverlapBCdmr+SoverlapABCdmr
     n123=SoverlapABCdmr
     

     plot.new()
     draw.triple.venn(area1=area1, area2=area2, area3=area3, n12=n12, n13=n13, n23=n23, n123=n123, category=names, col=rep("black", 3), fill=c("skyblue", "pink1", "orange"), ...)
     
}


plotFourVenns<-function(dmrList, names, filename){
     venn1<-vennBPthree(dmrList, names=names)
     venn2<-vennDMRthree(dmrList, names=names)
     venn3<-vennDMRthree2(dmrList, names=names, scaled=F)
     venn4<-vennDMRthree2(dmrList, names=names, scaled=T)
     dev.off()
     plot.new()
     pdf(filename)
     gl <- grid.layout(nrow=2, ncol=2)
     # grid.show.layout(gl)

     # setup viewports
     vp.1 <- viewport(layout.pos.col=1, layout.pos.row=1) 
     vp.2 <- viewport(layout.pos.col=2, layout.pos.row=1) 
     vp.3 <- viewport(layout.pos.col=1, layout.pos.row=2) 
     vp.4 <- viewport(layout.pos.col=2, layout.pos.row=2) 
     
     # init layout
     pushViewport(viewport(layout=gl))
     # access the first position
     pushViewport(vp.1)
     # start new base graphics in first viewport
     par(new=TRUE, fig=gridFIG())
     grid.draw(venn1)
     # done with the first viewport
     popViewport()
     # move to the next viewport
     pushViewport(vp.2)
     grid.draw(venn2)
     popViewport()
     # move to the next viewport
     pushViewport(vp.3)
     grid.draw(venn3)
     popViewport()
     # move to the next viewport
     pushViewport(vp.4)
     grid.draw(venn4)
     # done with this viewport
     popViewport()
     dev.off()

}