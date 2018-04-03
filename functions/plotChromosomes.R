## Created -/-/- by Md Haque
## Last modified 4/6/2016

## Haque's original function was named chromPlot.R. I have modified this function to more
## easily fit into the MEDIP pipeline. This function takes a DMR table and outputs the 
## chromosome plot. It requires the DMR table to have "chr", "start", and "stop" column
## names. The variable chrLengths must hold the lengths of the chromosomes to be plotted.
## The names of the chrLengths entries must correspond to the chromosome names.

## Several other functions are included in this file. They are used by the plotChromosomes
## function. I've left them here, as they are not used anywhere else in the MEDIP pipeline.

# Main chromosome plot function
plotChromosomes <- function(siteTable=NULL, clusters=NULL, chrLengths, ymar=4, 
                            xmar=6, cex.axis=1, markerWidth=0, main="", cex.lab=1, ...) {
  # Check for null DMR table
  if (is.null(nrow(siteTable))) return()
  if(nrow(siteTable) > 0){
    chrNames <- names(chrLengths)
    glen = 0.4
    if (markerWidth == 0) {
      markerWidth <- max(chrLengths) / 100
    }
    # generate empty plot
    par(mar = c(xmar, ymar, 4, 2) + 0.1)
    plot(c(1, max(chrLengths)), c(1 - glen, length(chrNames) + glen), type = "n", 
         xlab = "", ylab = "Chromosome", axes = FALSE, las = 2, main = main, cex.lab=cex.lab)
    # add axes
    axis(2, c(1:length(chrNames)), chrNames, las = 2, cex.axis = cex.axis)
    axis(1, seq(0,max(chrLengths), max(chrLengths)/20), 
         labels = paste(signif(seq(0, max(chrLengths), 
                                   max(chrLengths) / 20) / 1e6, digits=3)), 
         las = 2, cex.axis = cex.axis)
    mtext("Chromosome length (Mb)", side=1, line=4, cex=cex.lab)

    # plot lines and points for each chromosome
    for (i in chrNames) { 
      rows <- siteTable[which(siteTable[, "chr"] == i),]
      siteStart <- as.numeric(rows[, "start"])
      siteStop <- as.numeric(rows[, "stop"])
      if (!is.null(clusters)) {
        if (!(nrow(clusters) == 0)) {
          clusterRows <- clusters[which(clusters[, "chr"] == i), ]
          clusterStart <- as.numeric(clusterRows[, "start"])
          clusterStop <- as.numeric(clusterRows[, "stop"])
          plotRect(chrName = i, siteStart = clusterStart, siteStop = clusterStop, 
                   chrLengths = chrLengths)
        }
      }
      plotPoly(chrName = i, siteStart = siteStart, siteStop = siteStop, 
               chrLengths = chrLengths, markerWidth = markerWidth, ...)
    }
  } else {
    return()
  }
}

plotRect <- function(chrName = i, siteStart = siteStart, siteStop = siteStop, 
                     chrLengths = chrLengths, glen = 0.4) {
  scaledMax <- cScale(max(chrLengths), chrLengths, method = "relative", chrName = chrName)
  
  # Count the number of site markers
  nSites <- length(siteStart)
  
  # Scale the marker locations
  chrNumber <- match(chrName, names(chrLengths))
  
  scaledStart <- siteStart*scaledMax
  scaledStop <- siteStop*scaledMax
  scaledMidpoint <- rowMeans(cbind(scaledStart, scaledStop))
  
  ## Determine the direction of the Y plot (+ or -)
  ypos <- rep(chrNumber, nSites)
  ytop <- ypos - glen
  rect(scaledStart, ypos, scaledStop, ytop, col = "black")
}

# Used by plotChromosomes. This function plots lines representing chromosomes and adds site markers
plotPoly <- function(chrName, siteStart, siteStop, chrLengths, markerWidth = 500000, glen = 0.4, 
                     color.marker = "red", color.chr = "blue", color.border = "black") {
  
  scaledMax<-cScale(max(chrLengths), chrLengths, method="relative", chrName=chrName)
  
  # Count the number of site markers
  nSites <- length(siteStart)
  
  # Scale the marker locations
  chrNumber <- match(chrName, names(chrLengths))
  
  scaledStart <- siteStart * scaledMax
  scaledStop <- siteStop * scaledMax
  scaledMidpoint <- rowMeans(cbind(scaledStart, scaledStop))
  
  ## Determine the direction of the Y plot (+ or -)
  ypos <- rep(chrNumber, nSites)
  ytop <- ifelse(scaledStart > 0, ypos + glen, ypos - glen)
  
  polys <- c()
  if (nSites > 0) {
    for (i in 1:length(scaledMidpoint)) {
      polys <- rbind(polys, 
                     cbind(abs(scaledMidpoint[i]), ypos),
                     cbind(abs(scaledMidpoint[i]) - markerWidth * scaledMax, ytop),
                     cbind(abs(scaledMidpoint[i]) + markerWidth * scaledMax, ytop),
                     cbind(1, NA))
    }
  }
  ## Plot site markers
  polygon(polys, col = color.marker, border = color.border)
  
  ## Plot chromosome lines
  lines(c(1, chrLengths[chrNumber] * scaledMax), c(chrNumber, chrNumber), col = color.chr)
}


# This function is used to calculate scaling on the geneplots. It uses the vector of chromosome 
# lengths and returns a vector of scales. "points" should be the number of points to scale the 
# chromosomes to. "cLengths" should be a vector of chromosome lengths.
cScale <- function(points, cLengths, method = c("max", "relative"), chrName) {
  method <- match.arg(method)
  if (method == "max") {
    cScales <- points / cLengths[chrName]
  } else {
    cScales <- points / max(cLengths)
  }
  return(cScales)
}




# Experimental function to add different samples to same plot in different colors.
plotChromosomesMulti <- function(siteTable=NULL, clusters=NULL, chrLengths, ymar=4, tableNames=NULL,
                            xmar=5, cex.axis=1, markerWidth=0, main="", colors=NULL, oAng=NULL, ...) {
  # Check for null DMR table
  if (!is.list(siteTable)) {warning("siteTable must be a list"); return()}
  if (is.null(nrow(siteTable[[1]]))) {warning("first siteTable is empty"); return()}
  chrNames <- names(chrLengths)
  glen = 0.4
  if (markerWidth == 0) {
    markerWidth <- max(chrLengths) / 200
  }
  if (is.null(oAng)) oAng <- seq(from=-45, to=45, length.out=length(siteTable))
  # generate empty plot
  par(mar = c(xmar, ymar, 4, 2) + 0.1)
  plot(c(1, max(chrLengths)), c(1 - glen, length(chrNames) + glen), type = "n", 
       xlab = "", ylab = "", axes = FALSE, las = 2, main = main)
  # add axes
  axis(2, c(1:length(chrNames)), chrNames, las = 2, cex.axis = cex.axis)
  axis(1, seq(0,max(chrLengths), max(chrLengths)/20), 
       labels = paste(signif(seq(0, max(chrLengths), 
                                 max(chrLengths) / 20) / 1e6, digits=3), "Mb"), 
       las = 2, cex.axis = cex.axis)
  
  # plot lines and points for each chromosome
  if (is.null(colors)) colors <- 1:length(siteTable)
  for (j in 1:length(siteTable)) {
    st <- siteTable[[j]]
    for (i in chrNames) { 
      rows <- st[which(st[, "chr"] == i),]
      siteStart <- as.numeric(rows[, "start"])
      siteStop <- as.numeric(rows[, "stop"])
      if (!is.null(clusters)) {
        if (!(nrow(clusters) == 0)) {
          clusterRows <- clusters[which(clusters[, "chr"] == i), ]
          clusterStart <- as.numeric(clusterRows[, "start"])
          clusterStop <- as.numeric(clusterRows[, "stop"])
          plotRect(chrName = i, siteStart = clusterStart, siteStop = clusterStop, 
                   chrLengths = chrLengths)
        }
      }
      plotLine(chrName = i, siteStart = siteStart, siteStop = siteStop, color.marker = colors[j],
               chrLengths = chrLengths, markerWidth = markerWidth, oAngle = oAng[j], ...)
    }
  }
  legend("topright", legend=tableNames, fill=colors, bty="n")

}


plotLine <- function(chrName, siteStart, siteStop, chrLengths, markerWidth = 500000,
                     color.marker = "red", color.chr = "blue", oAngle=0) {

  # Count the number of site markers
  nSites <- length(siteStart)
  
  chrNumber <- match(chrName, names(chrLengths))
  
  # x position
  midpoint <- rowMeans(cbind(siteStart, siteStop))
  
  # top and bottom y positions
  xscale <- max(chrLengths)/(length(chrLengths)*1.2)
  ypos <- rep(chrNumber, nSites)
  yang.adj <- 0.4 * cos(pi * oAngle / 180)
  xang.adj <- 0.4 * sin(pi * oAngle / 180)
  ytop <- ypos + yang.adj
  xtop <- abs(midpoint) + (xscale*xang.adj)
  ybot <- ypos - yang.adj
  xbot <- abs(midpoint) - (xscale*xang.adj)

  ## Plot site markers
  if (nSites > 0) {
    for (i in 1:length(midpoint)) {
      segments(xbot[i], ybot[i], xtop[i], ytop[i], col=color.marker, lwd=2)

    }
  }

  ## Plot chromosome lines
  lines(c(1, chrLengths[chrNumber]), c(chrNumber, chrNumber), col = color.chr)
}