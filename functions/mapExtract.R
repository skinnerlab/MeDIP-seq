## Created -/-/2015 by Daniel Beck
## Last modified 4/6/2016


## This function extracts mapping information from the prepareData.R log file. It is used for 
## the report generation.

mapExtract <- function(cleaned = TRUE, sn=NA) {
  pctLine <- 15
  if (cleaned) pctLine <- 19
  prepareData <- readLines("prepareData.Rout")
  sampleIndex <- sapply(sn, function(i) {
    grep(pattern = paste("Running Bowtie2 on ", i, '\"', sep = ""), x = prepareData)
    })
  
  pct <- sapply(strsplit(prepareData[sampleIndex + pctLine], split = " "), function(i) i[1])
  rn <- sapply(strsplit(prepareData[sampleIndex + 1], split = " "), function(i) i[1])
  outTable <- as.data.frame(rbind(rn, pct))
  colnames(outTable) <- sn
  rownames(outTable) <- c("read number", "overall alignment rate")
  return(outTable)
}