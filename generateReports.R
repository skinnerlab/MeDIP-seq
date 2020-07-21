## Created 8/10/2015 by Daniel Beck
## Last modified 11/9/2016

## This script automatically generates all relevant reports for selected comparisons and
## p-values. The number of analyses, p-values, MTC p-values, and calcVar flags should be
## the same.

source("dataNames.R")
source("customFunctions.R")
library("rmarkdown")

standardReport <- c("analysis.name")
standardPvalues <- c("p.value")

cpgMaxV <- rep(NA, length(standardReport))
lenMaxV <- rep(NA, length(standardReport))
topNV <- rep(NA, length(standardReport))
calcVar <- rep(FALSE, length(standardReport))

## For generating standard reports
for (i in 1:length(standardReport)){
  analysisName <- standardReport[i]
  reportPvalue <- standardPvalues[i]
  cVar <- calcVar[i]
  cpgMax <- cpgMaxV[i]
  lenMax <- lenMaxV[i]
  topN <- topNV[i]
  mrs <- minRowSum[which(comparisonNames==standardReport[i])]
  
  reportFileName <- paste(resultsDirectory, standardReport, "/", gsub(" ", "_", projectName),
                          "_", analysisName, "_", standardPvalues[i], 
                          "_report2.pdf", sep = "")
  save(analysisName, reportPvalue, cVar, cpgMax, lenMax, topN, mrs,
       file = paste(codeDirectory, "/reportValues.Rdata", sep = ""))
  render(input = "medipReport.Rmd", output_file = reportFileName[i])
}



