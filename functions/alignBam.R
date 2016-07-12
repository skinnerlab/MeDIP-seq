# Extract fasta alignment from BAM file
alignBam <- function(analysis, multipleWindow = TRUE, dmrNum, pV) {
  
  ## Extract DMR table of interest
  if (multipleWindow) { 
    dmrTable <- methList2p[[pV]]
  } else { 
    dmrTable <- methList[[pV]]
  }
  
  bamList <- lapply(paste(dataDirectory, sbamFileName, sep=""), bamReader)
  for (i in 1:length(bamList)) {
    loadIndex(bamList[[i]], paste(dataDirectory, sbamFileName[i], ".bai", sep = ""))
  }
  seqSum <- getRefData(bamList[[1]])
  
  chr <- dmrTable$chr[dmrNum]
  start <- dmrTable$start[dmrNum]
  end <- dmrTable$start[dmrNum] + dmrTable$length[dmrNum]
  
  subBam <- list()
  for (i in 1:length(sbamFileName)) {
    subfilename <- paste(dataDirectory, "tempSub_", sbamFileName[i], sep="")
    command <- paste("samtools view -b ", dataDirectory, sbamFileName[i],
                     " ", chr, ":", start, "-", end,
                     " > ", subfilename, sep="")
    system(command)
    subBam[[i]] <- scanBam(subfilename)
  }
  
  aSeqs <- lapply(1:length(subBam), function(i) {
    alignSeqs(seqs = subBam[[i]][[1]]$seq, cigars = subBam[[i]][[1]]$cigar,
              startPos = subBam[[i]][[1]]$pos, sampleName = seqFiles$sampleName[[i]],
              refStart = start, refEnd = end)
  })
  
  combined <- do.call(c, sapply(aSeqs, function(i) i[[1]]))
  cpos <- do.call(c, lapply(subBam, function(i) i[[1]]$pos))
  combined <- combined[order(cpos)]
  return(combined)
}
