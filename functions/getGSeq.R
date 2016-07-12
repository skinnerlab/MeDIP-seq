
# This function takes a sequence (DNAStringSet) and a cigar string (from BAM file) and returns
# a properly gapped sequence. (It adds gaps according to the cigar string). This fuction
# doesn't actually seem to work. Something is off with the interpretion of the cigar string.
getGSeq <- function(seq, cigar) {
  nDel <- 0
  sc <- unlist(strsplit(cigar, split = ""))
  letterPos <- as.logical(match(sc, LETTERS, nomatch = 0))
  letters <- sc[letterPos]
  sc[letterPos] <- "_"
  numbers <- as.numeric(unlist(strsplit(paste(sc, collapse=""), split="_")))
  
  fragments <- list()
  for (i in 1:length(letters)) {
    if (letters[i] == "M") fragments[[i]] <- subseq(seq, 1, width = numbers[i])
    if (letters[i] == "D") fragments[[i]] <- ""; paste(rep("-", numbers[i]), collapse = "")
    if (letters[i] == "I") { 
      fragments[[i]] <- ""
      nDel <- nDel + numbers[i]
    }
  }
  gappedSeq <- paste(sapply(fragments, as.character), collapse = "")
  return(list(gappedSeq, nDel))    
}
