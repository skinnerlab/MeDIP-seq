## Created 2/10/2016 by Daniel Beck
## Last modified 5/5/2016

## This function converts DMR to underlying reference sequences.

dmr2seq <- function(dmrTable, bsgenomeName, extendEdges=0){
     # First extend edges by set amount (we often want DMR in region of genes, not directly 
     # overlapping). For overlapping, set extendEdges to 0
     if (extendEdges > 0) {
          seqLengths <- seqlengths(eval(as.symbol(bsgenomeName)))
          # extend edges to maximum required
          dmrTable$start <- dmrTable$start - 10000
          dmrTable$stop <- dmrTable$stop + 10000
          # adjust edges that run off the chromosome
          dmrTable$start <- sapply(dmrTable$start, function(i) max(i, 1))
          maxLen <- seqLengths[match(dmrTable$chr, names(seqLengths))]
          dmrTable$stop <- apply(cbind(dmrTable$stop, maxLen), 1, min)
     }

     chrs <- match(dmrTable$chr, seqnames(eval(as.symbol(bsgenomeName))))

     seqs <- sapply(1:length(chrs), function(i) {
          chr <- eval(as.symbol(bsgenomeName))[[chrs[i]]]
          seq <- subseq(chr, start=dmrTable$start[i], end=dmrTable$stop[i])
          seq})
     ## To avoid a mysterious bug that eventually results in a segmentation fault, 
     ## the loop below replaces the more simple code commented here. I don't know why
     ## the simple code fails, but it does in a mysterious way.
     # seqs <- DNAStringSet(seqs, use.names=T)
     a <- DNAStringSet(rep("N", length(seqs)))
     for (i in 1:length(seqs)) {
          a[[i]] <- seqs[[i]]
     }
     seqs <- a
     names(seqs) <- dmrTable$ID

     return(seqs)
}


