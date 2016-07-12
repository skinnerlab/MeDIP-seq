## Created 8/10/2015 by Daniel Beck
## Last modified 7/12/2016

## This script automatically generates all relevant reports for selected comparisons and
## p-values. The number of analyses, p-values, MTC p-values, and calcVar flags should be
## the same.

source("dataNames.R")
source("customFunctions.R")
library("rmarkdown")

standardReport <- c("allF1", "allF2", "allF3")
apoReport <- c("apoF1", "apoF2", "apoF3")
standardPvalues <- rep(1e-06, 3)
apoPvalues <- rep(1e-03, 3)
standardMTCpValues <- rep(0.05, 3)
apoMTCpValues <- rep(0.1, 3)

calcVar <- rep(TRUE, 3)

## For generating standard reports
for (i in 1:length(standardReport)){
  analysisName <- standardReport[i]
  reportPvalue <- standardPvalues[i] 
  MTCreportPvalue <- standardMTCpValues[i]
  cVar <- calcVar[i]
  reportFileName <- paste(resultsDirectory, standardReport, "/", gsub(" ", "_", projectName), 
                          "_", analysisName, "_", standardPvalues[i], "_", standardMTCpValues[i], 
                          "_report.pdf", sep = "")
  save(analysisName, MTCreportPvalue, reportPvalue, cVar, 
       file = paste(codeDirectory, "/reportValues.Rdata", sep = ""))
  render(input = "medipReport.Rmd", output_file = reportFileName[i])
}

## For generating intersection (APO) reports
if (!is.null(apoReport)) {
  for (i in 1:length(apoReport)) {
    analysisName <- apoReport[i]
    reportPvalue <- apoPvalues[i] 
    MTCreportPvalue <- apoMTCpValues[i]
    reportFileName <- paste(resultsDirectory, apoReport, "/", gsub(" ", "_", projectName), 
                            "_", analysisName, "_", apoPvalues[i], "_", apoMTCpValues[i], 
                            "_intersection_report.pdf", sep = "")
    save(analysisName, MTCreportPvalue, reportPvalue,
         file = paste(codeDirectory, "/reportValues.Rdata", sep = ""))
    render(input = "medipReport_APO.Rmd", output_file = reportFileName[i])
 
  } 
}



