## Created 9/18/2015 by Daniel Beck
## Last modified 4/6/2016

## This function generates a histogram with a maximum value. All values above this maximum
## are lumped into a >max category. This function should be able to replace the CpG density
## and DMR length histograms. However, this replacement is not yet complete.

histMax <- function(x, max = 15, ...) {
  if (length(x) > 0) {
    if (length(which(x > max)) > 0) {
      # all values > max get changed to max
      x[x > max] <- max + 0.01
      hist(x, breaks = 0:(max + 1), xaxt = "n", xaxp = c(0, max, max), ...)
      axis(1, at = c(0:(max + 1)), labels = c(0:max, paste(">", max, sep = "")))
    } else {
      hist(x, ...)
    }
  } else {
    return()
  }        
}
