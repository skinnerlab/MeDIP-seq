## Daniel Beck
## Created 8/10/2015
## Modified

## This script automatically generates all relevant reports for selected comparisons and p-values

source("dataNames.R")
source("customFunctions.R")
library("rmarkdown")

rAnalysis<-c("all", "all")
rpValues<-c(1e-05, 1e-06)
rMTCpValues<-c(0.05, 0.05)

for (i in 1:length(rAnalysis)){
     analysisName<-rAnalysis[i]; reportPvalue<-rpValues[i]; MTCreportPvalue<-rMTCpValues[i]; 
     fnReport<-paste(resultsDirectory, rAnalysis,"/", gsub(" ", "_", projectName), "_", analysisName, "_", rpValues[i], "_", rMTCpValues[i], "_report.pdf", sep="")
     fnSummary<-paste(resultsDirectory, rAnalysis,"/", gsub(" ", "_", projectName), "_", analysisName, "_", rpValues[i], "_", rMTCpValues[i], "_summary.pdf", sep="")
     
     save(analysisName, MTCreportPvalue, reportPvalue, file=paste(codeDirectory, "/reportValues.Rdata", sep=""))
     render(input="medipReport.Rmd", output_file=fnReport[i])
     render(input="briefSummary.Rmd", output_file=fnSummary[i])
}
     
     