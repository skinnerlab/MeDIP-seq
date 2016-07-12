
# This function takes a DNAStringSet formatted set of sequence reads, a corresponding set of
# cigar strings, and a vector of sequence start positions. It returns the aligned sequences
# and a number corresponding to the number of bases deleted from the reads that were not present
# in the reference. Since getGSeq() isn't working, this doesn't work well either.
alignSeqs <- function(seqs, cigars, startPos, sampleName = "", refStart = NULL, refEnd = NULL) {
  
  minPos <- min(c(startPos, refStart))
  maxPos <- max(c(startPos, refEnd))
  
  startGap <- sapply(startPos - minPos, function(i) paste(rep("-", i), collapse = ""))
  endGap <- sapply(maxPos - startPos, function(i) paste(rep("-", i), collapse = ""))
  
  processed <- sapply(1:length(seqs), function(i) {
    list(getGSeq(seq = seqs[i], cigar = cigars[i]))
  })
  cSeqs <- sapply(processed, function(i) i[[1]])
  nDel <- sapply(processed, function(i) i[[2]])
  alignedSeqs <- DNAStringSet(paste(startGap, cSeqs, endGap, sep=""))
  names(alignedSeqs) <- paste(sampleName, "_" , 1:length(alignedSeqs), sep="")
  alignedSeqs <- subseq(alignedSeqs, start = refStart - minPos +1, width=refEnd-refStart) # trim seqs
  return(list(alignedSeqs = alignedSeqs, nDel = sum(nDel)))
}