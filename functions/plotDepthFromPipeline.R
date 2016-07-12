# This function is a wrapper for plotDepth. It must be run with methLists and methListEtc
# objects in the environment.
plotDepthFromPipeline <- function(pV, analysis, multipleWindow, dmrID = 1, ...) {

  ## Extract DMR table of interest
  if (multipleWindow) { 
    dmrTable <- methList2p[[pV]]
  } else { 
    dmrTable <- methList[[pV]]
  }

  # Set up bam file objects
  bamList <- lapply(paste(dataDirectory, sbamFileName, sep=""), bamReader)
  for (i in 1:length(bamList)) {
    loadIndex(bamList[[i]], paste(dataDirectory, sbamFileName[i], ".bai", sep = ""))
  }
  seqSum <- getRefData(bamList[[1]])


  ####################
  ####################
  # plot depth

  pVdata <- as.data.frame(methListEtc[[pV]][match(dmrTable$ID, 
                                                  row.names(methListEtc[[pV]]), 
                                                  nomatch = 0), ], stringsAsFactors=F)
  pVdata <- apply(pVdata, 1, function(i) {
    data.frame("edgeR.p.value" = unlist(strsplit(i["edgeR.p.value"], split = ";")),
               "start" = unlist(strsplit(i["start"], split = ";")),
               "chr" = unlist(strsplit(i["chr"], split=";")),
               stringsAsFactors=F)
  })

  col <- c("cadetblue1", "cadetblue3", "cadetblue4", "brown1", "brown3", "brown4")
  i <- dmrID
  plotDepth(bams = bamList, chr = dmrTable$chr[i], start = dmrTable$start[i], 
            end = dmrTable$start[i] + dmrTable$length[i], 
            lg = "y", seqSum = seqSum, col = col, 
            main = paste("Chromosome ", dmrTable$chr[i], sep = ""), 
            bsgenomePackageName = bsgenomePackageName, addLegend=T,
            bamNames = seqFiles$sampleName, pVx = as.numeric(pVdata[[i]]$start) + (ws / 2), 
            pVy = as.numeric(pVdata[[i]]$edgeR.p.value), ...)

}

