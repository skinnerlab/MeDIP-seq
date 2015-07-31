## Daniel Beck
## Created 5/28/2015
## Modified 
##	5/29/2015 Finished writing functions. Sucessful test using steelhead RBC data.
##	6/4/2015  Fixed error caused by numbers being stored as characters
##   7/6/2015  Added CpG density plot function. Changed name of file to customFunctions.R. Added addToMergedResults function to keep track of all DMR information.
##   7/7/2015  Fixed addToMergedResults function that failed in the case of empty or single line dataframes
##   7/13/2015 Added file to git repository "medipPipeline"
##   7/14/2015 Fixed errors that killed code when passed lists to plots were empty. Now returns().

# This file contains code with additional functions. These functions include those necessary for the chromosome plot and the CpG plot.


# Modified from chromPlot.R
# siteTable: table holding genomic regions to be marked. This table should have each row being a genomic region to be marked on the plot and include the following columns   
   # chr: The chromosome name
   # start: Start position on chromosome for mark
   # stop: Stop position on chromosome for mark
# chrLengths: vector of chromosome lengths, labels should be chromosome names

plotChromosomes<-function(siteTable, chrLengths, ymar=4, xmar=5, cex.axis=1, markerWidth=0, main="", chrExclude="", ...){
     if (is.null(nrow(siteTable))) return()
     if(nrow(siteTable)>0){

          chrNames <- unique(siteTable$chr)
          chrLengths <- chrLengths[match(chrNames, names(chrLengths))]
          if (!is.na(match(chrExclude, chrNames)[1])){
               rc<-match(chrExclude, chrNames)
               chrNames<-chrNames[-rc]
               chrLengths<-chrLengths[-rc]
          }
     	glen=0.4
          if (markerWidth==0){
               markerWidth<-max(seqlengths(eval(parse(text=referenceName))))/100
          }
     	# generate empty plot
     	par(mar = c(xmar,ymar,4,2) + 0.1)
     	plot(c(1, max(chrLengths)), c(1 - glen, length(chrNames) + glen), type = "n", xlab="", ylab="", axes = FALSE, las=2, main=main)
     	# add axes
     	axis(2, c(1:length(chrNames)), chrNames, las = 2, cex.axis=cex.axis)
     	axis(1, seq(0,max(chrLengths),max(chrLengths)/20),labels=paste(signif(seq(0,max(chrLengths),max(chrLengths)/20)/1e6,digits=3),"Mb"), las = 2, cex.axis=cex.axis)

     	# plot lines and points for each chromosome
     	for (i in chrNames) { 
     		rows <- siteTable[which(siteTable[,"chr"]==i),]
     		siteStart <- as.numeric(rows[,"start"])
     		siteStop <- as.numeric(rows[,"stop"])
     		plotPoly(chrName=i, siteStart=siteStart, siteStop=siteStop, chrLengths=chrLengths, markerWidth=markerWidth, ...)
     	}
     } else {
          return()
     }
}

# Modified from chromPlot.R
# Used by plotChromosomes. The plotPoly function plots lines representing chromosomes and adds site markers
plotPoly<-function(chrName, siteStart, siteStop, chrLengths, fg, markerWidth=500000, glen=0.4, color.marker="red", color.chr="blue", color.border="black"){

	scaledMax<-cScale(max(chrLengths), chrLengths, method="relative", chrName=chrName)

	# Count the number of site markers
	nSites <- length(siteStart)
    
	# Scale the marker locations
	chrNumber <- match(chrName, names(chrLengths))

	scaledStart <- siteStart*scaledMax
	scaledStop <- siteStop*scaledMax
	scaledMidpoint <- rowMeans(cbind(scaledStart,scaledStop))

        ## Determine the direction of the Y plot (+ or -)
        ypos <- rep(chrNumber, nSites)
	ytop <- ifelse(scaledStart>0, ypos+glen, ypos-glen)

    	polys <- c()
        for ( i in 1:length(scaledMidpoint) ){
		polys<- rbind(polys,
				cbind(abs(scaledMidpoint[i]),ypos),
				cbind(abs(scaledMidpoint[i])-markerWidth*scaledMax, ytop),
				cbind(abs(scaledMidpoint[i])+markerWidth*scaledMax, ytop),
				cbind(1,NA))}

        ## Plot site markers
        polygon(polys, col=color.marker,border=color.border)

	## Plot chromosome lines
	lines(c(1,chrLengths[chrNumber]*scaledMax),c(chrNumber,chrNumber),col=color.chr)
}


# From chromPlot.R
# Used to calculate scaling on the geneplots
# Uses the vector of chromosome lengths and returns a vector
# of scales.
# Inputs:
   # points - the number of points to scale the chromosomes too
   # cLengths - a vector of chromosome lengths.
cScale <- function(points, cLengths, method=c("max","relative"), chrName) {
   method <- match.arg(method)
   if (method == "max") {
      cScales <- points / cLengths[chrName];
   }
   else {
      cScales <- points / max(cLengths)
   }
   return(cScales);
}

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

# This function collects all information about the original DMR windows and puts it in a dataframe with rows identical to the merged DMR list. This allows for verification of p-value and count information.
addToMergedResults<-function(allWindows, mergedWindows){
     # if data is empty, return NA. If data has one row, return that row.
     if (nrow(allWindows)<1) return(allWindows)
     if (nrow(allWindows)==1) return(allWindows)
     
	allWithIDs<-MEDIPS.selectROIs(results=allWindows, rois=mergedWindows)
	
	allCombined<-apply(allWithIDs, 2, function(i){
	     sapply(split(i, f=allWithIDs$ROI), paste, collapse=";")
	})
	
	# If it is a single row, allCombined is collapsed into a vector. This causes an error. So only combine if there is more than one merged window
	if (nrow(mergedWindows)>1){
	     allCombined<-allCombined[match(mergedWindows$ID, rownames(allCombined)),]
	}
	return(allCombined)
}

# This function calculates CpG density from a DMR list. The list must include chr, start, stop, and length
calcCpGdensity<-function(dmrList, maxDMR=1000){
     if (is.null(nrow(dmrList))) return(dmrList)
     if (nrow(dmrList)<maxDMR && nrow(dmrList)>0){
            cpgNum<-apply(dmrList, 1, function(i){
                 dinucleotideFrequency(subseq(eval(parse(text=referenceName))[[match(i["chr"], seqnames(eval(parse(text=referenceName))))]], start=as.numeric(i["start"]), end=as.numeric(i["stop"])))["CG"]
                 })
            cpgDensity<-100*cpgNum/dmrList$length
            return(cbind(dmrList, cpgNum, cpgDensity))
       } else {
            return(dmrList)
       }                        
}
# This function adds annotation information 
addAnnotationGFF<-function(dmrList, gff, chrPrefix="", maxDMR=1000){
     if (is.null(nrow(dmrList))) return(dmrList)
     if(nrow(dmrList)<maxDMR && nrow(dmrList)>0){
          # save original chromosome names
          ochr<-dmrList$chr
          # add prefix to chr name if necessary
          dmrList$chr<-paste(chrPrefix, dmrList$chr, sep="")
          # convert dmrList to GRanges object
          dmrList$start<-as.numeric(dmrList$start)
          dmrList$stop<-as.numeric(dmrList$stop)
          gdmrList<-makeGRangesFromDataFrame(dmrList)
          # find overlaps
          overlaps<-findOverlaps(gdmrList, gff)
          # add annotation to dmrList
          dmrList<-cbind(dmrList, annotation=NA)
          dmrList$annotation[queryHits(overlaps)]<-as.character(gff$group[subjectHits(overlaps)])
          # put original chromosome names back
          dmrList$chr<-ochr
     }
     return(dmrList)
}

addAnnotationBiomart<-function(dmrList, annotationObject, chrPrefix="", maxDMR=1000){
     if (is.null(nrow(dmrList))) return(dmrList)
     if(nrow(dmrList)<maxDMR && nrow(dmrList)>0){
          # save original chromosome names
          ochr<-dmrList$chr
          # add prefix to chr name if necessary
          dmrList$chr<-paste(chrPrefix, dmrList$chr, sep="")
          dmrList<-setAnnotation(regions=dmrList, annotation=annotationObject)
          # put original chromosome names back
          dmrList$chr<-ochr
     }
     return(dmrList)
}

## This function is taken from the MEDIPS library (MEDIPS.setAnnotation). It causes an error when there are no regions that can be annotated. I've modified the function to not error in this case.
setAnnotation<-function (regions, annotation, cnv = F) {
     tmp.regions = GRanges(seqnames = regions[, 1], ranges = IRanges(start = as.numeric(regions[, 2]), end = as.numeric(regions[, 3])))
     ans = NULL
     if (is.data.frame(annotation)) 
          annotation = list(annotation = annotation)
     for (n in names(annotation)) {
          anno = annotation[[n]]
          anno.data = GRanges(anno[, 2], ranges = IRanges(start = anno[, 3], end = anno[, 4]), ID = anno[, 1])
          overlapsM = IRanges::as.matrix(findOverlaps(tmp.regions, anno.data))
          splitL = split(overlapsM[, 2], overlapsM[, 1])
          maxEle = max(c(0, unlist(lapply(splitL, length))))
          if (maxEle == 0) {
               message("no \"", n, "\"annotation overlap the provided regions")
          } else {
               tmp.ans = matrix(ncol = maxEle, nrow = length(tmp.regions))
               colnames(tmp.ans) = paste(c(1:maxEle), "_", names(anno)[1], sep = "")
               j = rep(names(splitL), sapply(splitL, length))
               k = unlist(lapply(splitL, function(x) {1:length(x)}))
               tmp.ans[matrix(c(as.integer(j), as.integer(as.vector(k))), ncol = 2)] = as.character(values(anno.data)[overlapsM[, 2], 1])
               ans = cbind(ans, tmp.ans)
          }
     }
     if (cnv) {
          ans = data.frame(regions, CNV.log2.ratio = as.numeric(ans), stringsAsFactors = F)
     } else {
          if (!is.null(nrow(ans))){      # I added this test
               ans = cbind(regions, ans, stringsAsFactors = F)
          } else {
               ans <- regions           # If the test fails return the unannotated dmrList
          }
     }
     return(ans)
}

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

## Venn diagram wrappers for using GRanges objects
vennCustomFour<-function(a, b, c, d, names, ...){
        ab<-subsetByOverlaps(a,b); ba<-subsetByOverlaps(b,a)
        ac<-subsetByOverlaps(a,c); ca<-subsetByOverlaps(c,a)
        ad<-subsetByOverlaps(a,d); da<-subsetByOverlaps(d,a)
        bc<-subsetByOverlaps(b,c); cb<-subsetByOverlaps(c,b)
        bd<-subsetByOverlaps(b,d); db<-subsetByOverlaps(d,b)
        cd<-subsetByOverlaps(c,d); dc<-subsetByOverlaps(d,c)
        abc<-subsetByOverlaps(ab,c); bac<-subsetByOverlaps(ba,c); cab<-subsetByOverlaps(c,ab)
        cba<-subsetByOverlaps(c,ba)
        abd<-subsetByOverlaps(ab,d); bad<-subsetByOverlaps(ba,d); dab<-subsetByOverlaps(d,ab)
        dba<-subsetByOverlaps(d,ba)
        acd<-subsetByOverlaps(ac,d); cad<-subsetByOverlaps(ca,d); dac<-subsetByOverlaps(d,ac)
        dca<-subsetByOverlaps(d,ca)
        bcd<-subsetByOverlaps(bc,d); cbd<-subsetByOverlaps(cb,d); dbc<-subsetByOverlaps(d,bc)
        dcb<-subsetByOverlaps(d,cb)
        abcd<-subsetByOverlaps(abc, d); cabd<-subsetByOverlaps(cab, d)
        cbad<-subsetByOverlaps(cba, d); bacd<-subsetByOverlaps(bac, d)
        
        dbca<-subsetByOverlaps(dbc, a); dcba<-subsetByOverlaps(dcb, a)
        cbda<-subsetByOverlaps(cbd, a); bcda<-subsetByOverlaps(bcd, a)
        
        acdb<-subsetByOverlaps(acd, b); cadb<-subsetByOverlaps(cad, b)
        dacb<-subsetByOverlaps(dac, b); dcab<-subsetByOverlaps(dca, b)
        
        abdc<-subsetByOverlaps(abd, c); dabc<-subsetByOverlaps(dab, c)
        dbac<-subsetByOverlaps(dba, c); badc<-subsetByOverlaps(bad, c)

        n1234<-min(length(abcd),length(cabd),length(cbad),length(bacd),
                   length(dbca),length(dcba),length(cbda),length(bcda),
                   length(acdb),length(cadb),length(dacb),length(dcab),
                   length(abdc),length(dabc),length(dbac),length(badc))
     
        n123<-min(length(abc),length(bac),length(cab),length(cba))
        n124<-min(length(abd),length(bad),length(dab),length(dba))
        n134<-min(length(acd),length(cad),length(dac),length(dca))
        n234<-min(length(bcd),length(cbd),length(dbc),length(dcb))
        n12<-min(length(ab),length(ba))
        n13<-min(length(ac),length(ca))
        n14<-min(length(ad),length(da))
        n23<-min(length(bc),length(cb))
        n24<-min(length(bd),length(db))
        n34<-min(length(cd),length(dc))


        area1<-length(a)
        area2<-length(b)
        area3<-length(c)
        area4<-length(d)

        plot.new()
        draw.quad.venn(area1=area1, area2=area2, area3=area3, area4=area4, n12=n12, n13=n13, n14=n14, n23=n23, n24=n24, n34=n34, n123=n123, n124=n124, n134=n134, n234=n234, n1234=n1234, category = names, lwd = rep(2, 4), lty = rep("solid", 4), col = rep("black", 4), fill = c("skyblue", "pink1", "mediumorchid", "orange"), ...)
}

vennCustomThree<-function(a, b, c, names, ...){
        ab<-subsetByOverlaps(a,b); ba<-subsetByOverlaps(b,a)
        ac<-subsetByOverlaps(a,c); ca<-subsetByOverlaps(c,a)
        bc<-subsetByOverlaps(b,c); cb<-subsetByOverlaps(c,b)
        abc<-subsetByOverlaps(ab,c); bac<-subsetByOverlaps(ba,c); cab<-subsetByOverlaps(c,ab)
        cba<-subsetByOverlaps(c,ba)

        n123<-max(length(abc),length(bac),length(cab),length(cba))
        n12<-max(length(ab),length(ba))
        n13<-max(length(ac),length(ca))
        n23<-max(length(bc),length(cb))

        plot.new()
        draw.triple.venn(area1=length(a), area2=length(b), area3=length(c), n12=n12, n13=n13, n23=n23, n123=n123, category=names, col=rep("black", 3), fill=c("skyblue", "pink1", "orange"), ...)
}


vennCustomTwo<-function(a, b, names, ...){
     ab<-subsetByOverlaps(a,b); ba<-subsetByOverlaps(b,a)
     n12<-max(length(ab),length(ba))
     
     plot.new()
     draw.pairwise.venn(area1=length(a), area2=length(b), cross.area=n12, category=names, col=rep("black", 2), fill=c("skyblue", "pink1"), ...)
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


extendDMR<-function(dmrList, windowSize=1, pValueCutoff=0.9){
     dmrList
     
     
}