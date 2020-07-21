## Created 5/27/2015 by Daniel Beck
## Last modified 2/22/2017

## This code is performs the MeDIP analysis. It is the second step in the analysis
## pipeline. The prepareData.R script should typically be run first. The dataNames.R
## configuration file is also used for this script.

## This latest version only loads the samples needed for the comparisons selected. It 
## does not perform any QC. It can be used when additional comparisons are added
## to the project without the addition of new samples.


# Load relevant libraries, custom functions, and the configuration script. Note, the
# order is important, as dataNames.R holds information about the bsgenomePackageName.
library(MEDIPS)
source("dataNames.R")
source("customFunctions.R")
library(bsgenomePackageName, character.only = T, lib.loc = genomeDirectory)


#############################
## Read files and reformat ##
#############################
## This step reads in the sorted BAM sample files and converts them to a matrix of
## genomic windows with associated coverage (read counts). These options are defined
## in the dataNames.R configuration file.

select.ind <- unique(unlist(sapply(comparison, function(i) c(i$mset1, i$mset2))))
select.ind <- select.ind[which(select.ind > 0)]

select.files <- paste(dataDirectory, sbamFileName, sep = "")[select.ind]

processedBamFiles <- vector("list", length(sbamFileName))
for (i in 1:length(select.ind)) {
  processedBamFiles[select.ind[i]] <- MEDIPS.createSet(file=select.files[i],
                                                       BSgenome = bsgenomePackageName, 
                                                       extend = extend, 
                                                       shift = shift, 
                                                       uniq = uniq, 
                                                       window_size = ws, 
                                                       chr.select = chr.select)
}


###################
## Normalization ##
###################
## This step normalizes the data by calculating local CpG density. This is set to false
## by default in our pipline. This was taken from Haque's analysis and hasn't been 
## sufficiently explored.

if (CScalc) {
  CS <- MEDIPS.couplingVector(pattern = "CG", refObj = processedBamFiles[[select.ind[1]]])
} else {
  CS <- NULL
}


###################
## Identify DMRs ##
###################
## This step performs the actual DMR analysis. Each genomic window is assigned a 
## probability of being a DMR. All samples in mset1 are compared with all samples 
## in mset2. This code loops over all comparisons specified in the dataNames.R 
## configuration script.

for (analysis in 1:length(comparison)) {  
  mset1<-processedBamFiles[comparison[[analysis]]$mset1]
  mset2<-processedBamFiles[comparison[[analysis]]$mset2]
  # I'm not sure when this would occur. I think it was an early test.
  if (length(mset1)==0) {
    mset1 <- NULL
  }
  if (length(mset2)==0) {
    mset2 <- NULL
  }
  # Perform analysis
  methResults <- MEDIPS.meth(MSet1 = mset1, 
                             MSet2 = mset2,
                             p.adj = p.adj,
                             diff.method = diff.method,
                             MeDIP = MeDIP,
                             CNV = CNV,
                             CSet = CS,
                             minRowSum = minRowSum[analysis])
  # Save results to a comparison specific folder in the results directory
  system(paste("mkdir ", resultsDirectory, comparisonNames[analysis], sep=""))
  save(methResults, file=paste(resultsDirectory, 
                               comparisonNames[analysis], 
                               "/methResults.RData", sep=""))
  # Clean up unnecessary objects
  rm(methResults)
}
