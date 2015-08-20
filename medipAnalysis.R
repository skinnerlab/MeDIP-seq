## Daniel Beck
## Created 5/27/2015
## Modified
##	6/2/2015  Removed the CS object to prevent CpG normalization. This is probably unnecessary, but rerunning analysis to be sure. Also commented out BSgenome creation. This only needs to be done once.
##	6/22/2015 Added back the CS object to ensure its availability. Modified code to analyze finch dataset
##	6/29/2015 Split file processing (up to BAM file) and BSgenome package creation to separate script.
##   6/30/2015 Continued development of new pipeline structure. Added CpG density calculation. Fixed CpG length calculation (was previously 1bp short)
##   7/1/2015  Continued development
##   7/2/2015  Tested on steelheadTrial dataset
##   7/6/2015  Made CpG density calculation generic to all analyses
##   7/7/2015  Removed CpG density calculation to function in customFunctions.R. Added ability to limit calculation to smaller DMR numbers. Fixed several NA/NULL logical tests.
##   7/9/2015  Added annotation for gff and biomart based data. Added DMR number summary output as dmrNumber.csv in results folder. Moved pValues to dataNames.R for use in medipReport.Rmd.
##   7/10/2015 Added CS normalization option.
##   7/13/2015 Added file to git repository "medipPipeline"
##   8/5/2015  Split all steps after medips analysis to medipProcessing.R file

## This code is intended to perform the MeDIP analysis. It is intended to work as
## part of a pipeline including customFunctions.R and prepareData.R. This script
## requires the presence of an R object file containing several object names.

# Load necessary packages
library(MEDIPS)
source("dataNames.R")
source("customFunctions.R")
library(bsgenomePackageName, character.only=T, lib.loc=genomeDirectory)
library(biomaRt)
library(rmarkdown)

###############################
## Read files and preprocess ##
###############################

processedBamFiles <-
  lapply(
    X = paste(dataDirectory, sbamFileName, sep = ""), FUN = MEDIPS.createSet, BSgenome =
      bsgenomePackageName, extend = extend, shift = shift, uniq = uniq, window_size =
      ws, chr.select = chr.select
  )

###################
## Identify DMRs ##
###################

# Normalization if necessary
if (CScalc){
        CS<-MEDIPS.couplingVector(pattern = "CG", refObj = processedBamFiles[[1]])
} else {
        CS<-NULL
}

mset1<-processedBamFiles[comparison[[analysis]]$mset1]
mset2<-processedBamFiles[comparison[[analysis]]$mset2]
if (length(mset1)==0) mset1<-NULL
if (length(mset2)==0) mset2<-NULL

# do all comparisons
for (analysis in 1:length(comparison)){  
  methResults <- MEDIPS.meth(
  	MSet1=mset1,
  	MSet2=mset2,
  	p.adj=p.adj,
  	diff.method=diff.method,
  	MeDIP=MeDIP,
  	CNV=CNV,
  	CSet=CS,
  	minRowSum=minRowSum
  )
  
  system(paste("mkdir ", resultsDirectory, comparisonNames[analysis], sep=""))
  save(methResults, file=paste(resultsDirectory, comparisonNames[analysis], "/methResults.RData", sep=""))

  # clean up unnecessary objects
  rm(methResults)
}


