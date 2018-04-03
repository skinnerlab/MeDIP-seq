## Created 8/10/2015 by Daniel Beck
## Last modified 11/9/2016

## This script automatically generates all relevant reports for selected comparisons and
## p-values. The number of analyses, p-values, and option flags should be the same.

source("dataNames.R")
source("customFunctions.R")
library("rmarkdown")

report.analyses <- c("DDT.Caput", "Vin.Caput")
report.pvalues <- rep(1e-5, 2)
report.filenames <- paste(resultsDirectory, report.analyses, "/report_", gsub(" ", "_", projectName), "_", 
                         report.analyses, "_", report.pvalues, ".pdf", sep="")

cpgMaxV <- rep(NA, length(report.analyses))    # Y-axis maximum for CpG density histogram (NA for auto).
lenMaxV <- rep(NA, length(report.analyses))    # Y-axis maximum for DMR length histogram (NA for auto).
topNV <- rep(NA, length(report.analyses))      # Generate figures using top N DMR (NA for all DMR).

## For generating reports
for (i in 1:length(report.analyses)){
  analysisName <- report.analyses[i]
  reportPvalue <- report.pvalues[i] 
  cpgMax <- cpgMaxV[i]
  lenMax <- lenMaxV[i]
  topN <- topNV[i]
  save(analysisName, reportPvalue, cpgMax, lenMax, topN,
       file = paste(codeDirectory, "/reportValues.Rdata", sep = ""))
  render(input = "medipReport.Rmd", output_file = reportFileName[i])
}





