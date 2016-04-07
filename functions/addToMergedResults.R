## Created -/-/2015 by Daniel Beck
## Last modified 4/6/2016

## This function collects all information about the original DMR windows and puts 
## it in a dataframe with rows identical to the merged DMR list. This allows for 
## verification of p-value and count information.

addToMergedResults <- function(allWindows, mergedWindows) {
  # if data is empty, return NA. If data has one row, return that row.
  if (nrow(allWindows) <= 1) {
    return(allWindows)
  }
  allWithIDs <- MEDIPS.selectROIs(results = allWindows, rois = mergedWindows)
  allCombined <- apply(allWithIDs, 2, function(i) {
    sapply(split(i, f = allWithIDs$ROI), paste, collapse = ";")
  })
  # If it is a single row, allCombined is collapsed into a vector. 
  # This causes an error. So only combine if there is more than one merged window
  if (nrow(mergedWindows) > 1) {
    allCombined <- allCombined[match(mergedWindows$ID, rownames(allCombined)),]
  }
  allCombined <- rbind(allCombined)
  return(allCombined)
}
