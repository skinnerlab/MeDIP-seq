# Modified from chromPlot.R
# siteTable: table holding genomic regions to be marked. This table should have each row being a genomic region to be marked on the plot and include the following columns   
# chr: The chromosome name
# start: Start position on chromosome for mark
# stop: Stop position on chromosome for mark
# chrLengths: vector of chromosome lengths, labels should be chromosome names

plotChromosomes<-function(siteTable=NULL, clusters=NULL, chrLengths, ymar=4, 
                          xmar=5, cex.axis=1, markerWidth=0, main="", ...) {
     # Check for null DMR table
     if (is.null(nrow(siteTable))) return()
     if(nrow(siteTable)>0){
          
          chrNames <- names(chrLengths)
          
          glen=0.4
          if (markerWidth==0){
               markerWidth<-max(chrLengths)/100
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
#               if (nrow(rows) < 1)
               siteStart <- as.numeric(rows[,"start"])
               siteStop <- as.numeric(rows[,"stop"])
               if(!is.null(clusters)){
                    if (!(nrow(clusters)==0)){
                         clusterRows <- clusters[which(clusters[,"chr"]==i),]
                         clusterStart <- as.numeric(clusterRows[,"start"])
                         clusterStop <- as.numeric(clusterRows[,"stop"])
                         plotRect(chrName=i, siteStart=clusterStart, siteStop=clusterStop, chrLengths=chrLengths)
                    }
               }
               plotPoly(chrName=i, siteStart=siteStart, siteStop=siteStop, chrLengths=chrLengths, markerWidth=markerWidth, ...)
          }
     } else {
          return()
     }
}

plotRect<-function(chrName=i, siteStart=siteStart, siteStop=siteStop, chrLengths=chrLengths, glen=0.4){
     
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
     ytop <- ypos-glen
     rect(scaledStart, ypos, scaledStop, ytop, col="black")
     
}


# Modified from chromPlot.R
# Used by plotChromosomes. The plotPoly function plots lines representing chromosomes and adds site markers
plotPoly<-function(chrName, siteStart, siteStop, chrLengths, markerWidth=500000, glen=0.4, color.marker="red", color.chr="blue", color.border="black"){
     
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
     if (nSites > 0){
          for ( i in 1:length(scaledMidpoint) ){
               polys<- rbind(polys,
                             cbind(abs(scaledMidpoint[i]),ypos),
                             cbind(abs(scaledMidpoint[i])-markerWidth*scaledMax, ytop),
                             cbind(abs(scaledMidpoint[i])+markerWidth*scaledMax, ytop),
                             cbind(1,NA))
          }
     }
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

