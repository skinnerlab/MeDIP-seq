
# This function generates the depth plot for the bam region of interest. 
plotDepth <- function(bams, bamNames = NULL, chr, start, end, lg = "y", seqSum, col, 
                      xlab = "Chromosome position", ylab = "Read depth", legendPos = "topright",
                      bsgenomePackageName, pVx = NULL, pVy = NULL, addLegend = T, ...) {
  
  if (is.null(bamNames)) bamNames <- names(bams)
  if (is.null(bamNames)) bamNames <- paste("S", 1:length(bams), sep="")
  
  id <- seqSum$ID[match(chr, seqSum$SN)]
  
  ref <- getSeq(eval(parse(text = bsgenomePackageName)), start = start, 
                end = end, names = as.character(chr))
  cgPos <- as.numeric(gregexpr("CG", ref)[[1]])
  
  ad <- lapply(bams, function(i) alignDepth(bamRange(i, c(id, start, end))))
  maxDepth <- max(c(sapply(ad, function(i) max(i@depth)), 3))
  
  x <- ad[[1]]@pos
  y <- ad[[1]]@depth
  yS <- 0
  if (lg == "y") {
    y[which(y == 0)] <- 1
    yS <- 1
  }
  par(mar = c(5, 4, 4, 5) + 0.1)
  plot(1, type = "n", xlim = c(start, end), ylim = c(yS, maxDepth), 
       las = 1, log = lg, xlab = xlab, ylab = ylab, ...)
  
  abline(v=cgPos+start, col="lightgray")
  
  for (i in 1:length(bams)) {
    x <- ad[[i]]@pos
    y <- ad[[i]]@depth
    if (lg == "y") y[which(y == 0)] <- 1
    points(x, y, col = col[[i]], type = "l", bty = "n")
  }
  if (!is.null(pVx) & !is.null(pVy)) {
    par(new=TRUE)
    plot(pVx, pVy, col = "green",  pch = 15, yaxt = "n", xlab = "", ylab = "", 
         ylim = c(0, 0.01), xlim = c(start, end))
    axis(4)
    mtext("edgeR p-value", side = 4, line = 3)
  }
  if (addLegend)
    legend(legendPos, col = col, bamNames, lty = 1, lwd = 3)
   #legend("topright", col=c(col, "lightgray"), c(bamNames, "CpG"), 
   #       pch=c(rep("", length(bamNames)), "|"), lty = c(rep(1, length(bamNames)), NA))
  
}
