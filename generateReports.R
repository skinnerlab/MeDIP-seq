## Created 8/10/2015 by Daniel Beck
## Last modified 4/6/2016

## This script automatically generates all relevant reports for selected comparisons and
## p-values. The number of analyses, p-values, MTC p-values, and calcVar flags should be
## the same.

source("dataNames.R")
source("customFunctions.R")
library("rmarkdown")

rAnalysis <- c("all")
rpValues <- c(1e-04)
rMTCpValues <- c(0.1)
calcVar <- c(TRUE)

for (i in 1:length(rAnalysis)){
  analysisName <- rAnalysis[i]
  reportPvalue <- rpValues[i] 
  MTCreportPvalue <- rMTCpValues[i]
  cVar <- calcVar[i]
  fnReport <- paste(resultsDirectory, rAnalysis, "/", gsub(" ", "_", projectName), 
                    "_", analysisName, "_", rpValues[i], "_", rMTCpValues[i], 
                    "_report.pdf", sep = "")
  save(analysisName, MTCreportPvalue, reportPvalue, 
       file = paste(codeDirectory, "/reportValues.Rdata", sep = ""))
  render(input = "medipReport.Rmd", output_file = fnReport[i])
}
     
     
