## Created 10/14/2015 by Daniel Beck
## Last modified 4/6/2016

## This function takes start and stop positions and returns overlaps. Positions must be 
## equilivant (on the same chromosome).

pairOverlap <- function(start1, stop1, start2, stop2) {
  if (length(start1) != length(stop1)) stop("start1 and stop1 must have identical length")
  if (length(start2) != length(stop2)) stop("start2 and stop2 must have identical length")
  overlaps <- list()  # list to hold overlaping regions
  index <- 1  # next index in overlaps list that is currently empty
  for (pos1 in 1:length(start1)) {
    for (pos2 in 1:length(start2)) {
      if ((start1[pos1] < stop2[pos2]) & (stop1[pos1] > start2[pos2])) {
        overlaps[[index]] <- c(max(start1[pos1], start2[pos2]), min(stop1[pos1], stop2[pos2]))
        index <- index + 1
      }
    }
  }
  overlaps <- do.call(rbind, overlaps)
  if (!is.null(overlaps)) {
    colnames(overlaps) <- c("start", "stop")
  }
  return(overlaps)
}

