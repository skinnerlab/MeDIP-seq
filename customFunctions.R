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

plotChromosomes<-function(siteTable, chrLengths, ymar=4, xmar=5, cex.axis=1, main="", ...){
     if (is.null(nrow(siteTable))) return()
     if(nrow(siteTable)>0){

     	chrNames <- unique(siteTable$chr)
     	chrLengths <- chrLengths[match(chrNames, names(chrLengths))]

     	glen=0.4

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
     		plotPoly(chrName=i, siteStart=siteStart, siteStop=siteStop, chrLengths=chrLengths, ...)
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
plotCpGdensity<-function(dmrList, xlab="Number of CpG sites per 100bp", ylab="Number of differential DNA methylation regions", ...){
     if (is.null(dmrList$cpgDensity)) return()
     if (is.null(nrow(dmrList))) return()
     if (nrow(dmrList)>0){

        # start at 0 unless there are no CpGs in the 0 group
        firstBreak<-0
        if (length(table(dmrList$cpgDensity<=0.5))<2) firstBreak<-1

        # all CpG density > 10.5 get lumped together into one group
        dmrList$cpgDensity[dmrList$cpgDensity>10.5]<-11

        # generate plot
        hist(dmrList$cpgDensity, breaks=firstBreak:12-0.5, xaxt="n", 
             xlab=xlab, ylab=ylab, xaxp=c(firstBreak,10,10-firstBreak), ...)
        axis(1, at=c(firstBreak:11), labels=c(firstBreak:10, ">10"))
        
     } else {
          return()
     }        
}

# This function collects all information about the original DMR windows and puts it in a dataframe with rows identical to the merged DMR list. This allows for verification of p-value and count information.
addToMergedResults<-function(allWindows, mergedWindows){
     # if data is empty, return NA. If data has one row, return that row.
     if (nrow(allWindows)<1) return(NA)
     if (nrow(allWindows)==1) return(allWindows)
     
	allWithIDs<-MEDIPS.selectROIs(results=allWindows, rois=mergedWindows)
	allSplit<-split.data.frame(allWithIDs, f=allWithIDs$ROI)

	allCombined<-apply(allWithIDs, 2, function(i){
		sapply(split(i, f=allWithIDs$ROI), paste, collapse=";")
	})
	allCombined<-allCombined[match(mergedWindows$ID, rownames(allCombined)),]
	return(allCombined)
}

# This function calculates CpG density from a DMR list. The list must include chr, start, stop, and length
calcCpGdensity<-function(dmrList, maxDMR){
     if (is.null(nrow(dmrList))) return(NA)
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
addAnnotationGFF<-function(dmrList, gff, maxDMR=1000){
     if (is.null(nrow(dmrList))) return(NA)
     if(nrow(dmrList)<maxDMR && nrow(dmrList)>0){
          # convert dmrList to GRanges object
          dmrList$start<-as.numeric(dmrList$start)
          dmrList$stop<-as.numeric(dmrList$stop)
          gdmrList<-makeGRangesFromDataFrame(dmrList)
          # find overlaps
          overlaps<-findOverlaps(gdmrList, gff)
          # add annotation to dmrList
          dmrList<-cbind(dmrList, annotation=NA)
          dmrList$annotation[queryHits(overlaps)]<-as.character(gff$group[subjectHits(overlaps)])
          return(dmrList)
     } else {
          return(dmrList)
     }
}

addAnnotationBiomart<-function(dmrList, annotationObject, maxDMR=1000){
     if (is.null(nrow(dmrList))) return(NA)
     if(nrow(dmrList)<maxDMR && nrow(dmrList)>0){
          dmrList<-setAnnotation(regions=dmrList, annotation=annotationObject)
     } else {
          return(dmrList)
     }
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
               return(NULL) ### I added this line. It also skips the "if (cnv)". I don't know if this is desired.
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
          ans = cbind(regions, ans, stringsAsFactors = F)
     }
     return(ans)
}