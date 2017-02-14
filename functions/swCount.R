## Function to use a sliding window to count things in windows
swCount <- function(numCpG, chr, neighborhood) {
  # If no chromosome information, set to a single chromosome
  if (is.null(chr)) chr <- rep(1, length(numCpG))
  lnumCpG <- split(x=numCpG, f=chr)
  # Create vector to store results
  countResults <- numeric(length(numCpG))

  # Add results in original order
  for (i in 1:length(lnumCpG)){
    count <- numeric(length(lnumCpG[[i]]))
    # The begining and end of the chromosome is done separatly to account for edge effects.
    if (length(count) > 2*neighborhood) {
      # Count begining of chromosome
      count[1:neighborhood] <- sapply(1:neighborhood, function(j) sum(lnumCpG[[i]][1:(neighborhood+j)]))
      # Count middle of chromosome
      count[(neighborhood+1):(length(count)-neighborhood)] <- sapply(1:(length(lnumCpG[[i]])-(2*neighborhood)), function(j) sum(lnumCpG[[i]][j:(j+(2*neighborhood))]))
      # Count end of chromosome
      count[(length(count)-neighborhood):length(count)] <- sapply((length(count)-neighborhood):length(count), function(j) sum(lnumCpG[[i]][(j-neighborhood):length(count)]))
    } else {
      warning("Chromosome too short for current identifyZero function implimentation.")
    }
    countResults[which(chr==names(lnumCpG)[i])] <- count
  }
  return(countResults)
}
