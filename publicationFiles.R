## Daniel Beck
## Created 9/15/2015
## Modified

## This code is intended to produce all publication figures and tables. It repeats many of the report figures. This file should be modified for each project so that includes the exact figures produced for the final manuscript.

source("dataNames.R")
source("customFunctions.R")
library(bsgenomePackageName, character.only=T, lib.loc=genomeDirectory)
library(biomaRt)

#########################
##### Figure values #####
#########################

# These values need to be modified. They define which analysis and p-value to use to generate result figures.
analysisName<-"all"
reportPvalue<-1e-04
MTC=FALSE

#########################
#########################
#########################


# select the desired comparison
an<-which(comparisonNames==analysisName)

# list position for these p-values
pvc<-which(pValues==reportPvalue)
MTCpvc<-which(pValues==MTCreportPvalue)

# load results
load(paste(resultsDirectory, comparisonNames[an], "/methLists.RData", sep=""))
load(paste(resultsDirectory, comparisonNames[an], "/methResults.RData", sep=""))