## Created 10/14/2015 by Daniel Beck
## Last modified 4/6/2016

## This function takes two matrices or data.frames and removes lines in matrix 1 that are also 
## in matrix 2. It returns matrix 1 without the lines of matrix 2.

removeDuplicates <- function(mat1, mat2) {
  # If mat1 is NULL or NA, return empty matrix (so that nrow still works and returns 0)
  if (is.null(mat1)) return(cbind(list.index = 1, dmr.index = 1)[0, ])
  if (is.na(mat1)) return(cbind(list.index = 1, dmr.index = 1)[0, ])
  # If mat2 is NULL or NA, return mat1
  if (is.null(mat2)) return(mat1)
  if (is.na(mat2)) return(mat1)
  
  mat1c <- apply(mat1, 1, paste, collapse = "_")
  mat2c <- apply(mat2, 1, paste, collapse = "_")
  
  overlaps<-match(mat2c, mat1c, nomatch = 0)
  # If there are no overlaps, the subset returns an empty matrix.
  # This line prevents this by removing a non-existant line (out of bounds).
  # This method seems strange and may be fragile.
  if ((unique(overlaps)[1] == 0) & (length(unique(overlaps)) == 1)) overlaps <- nrow(mat1) + 10
  uniqueSet <- mat1[-c(overlaps), ]
  if (length((1:length(mat1c))[-c(overlaps)]) < 2) uniqueSet <- as.data.frame(t(uniqueSet))
  
  return(uniqueSet)
}


