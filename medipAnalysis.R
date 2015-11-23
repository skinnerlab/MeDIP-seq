## Daniel Beck
## Created 5/27/2015
## Modified 10/8/2015

## This code is intended to perform the MeDIP analysis. It is intended to work as
## part of a pipeline including customFunctions.R and prepareData.R. This script
## requires the presence of an R object file containing several object names.

# Load necessary packages
library(MEDIPS)
source("dataNames.R")
source("customFunctions.R")
library(bsgenomePackageName, character.only = T, lib.loc = genomeDirectory)
library(biomaRt)
library(rmarkdown)

###############################
## Read files and preprocess ##
###############################

processedBamFiles <- lapply(X = paste(dataDirectory, sbamFileName, sep = ""), 
                            FUN = MEDIPS.createSet, 
                            BSgenome = bsgenomePackageName, 
                            extend = extend, 
                            shift = shift, 
                            uniq = uniq, 
                            window_size = ws, 
                            chr.select = chr.select)


###################
## Identify DMRs ##
###################

# Normalization if necessary
if (CScalc) {
  CS <- MEDIPS.couplingVector(pattern = "CG", refObj = processedBamFiles[[1]])
} else {
  CS <- NULL
}



# do all comparisons
for (analysis in 1:length(comparison)) {  
  mset1<-processedBamFiles[comparison[[analysis]]$mset1]
  mset2<-processedBamFiles[comparison[[analysis]]$mset2]
  if (length(mset1)==0) {
    mset1 <- NULL
  }
  if (length(mset2)==0) {
    mset2 <- NULL
  }
  methResults <- MEDIPS.meth(MSet1 = mset1, 
                             MSet2 = mset2,
                             p.adj = p.adj,
                             diff.method = diff.method,
                             MeDIP = MeDIP,
                             CNV = CNV,
                             CSet = CS,
                             minRowSum = minRowSum)
  system(paste("mkdir ", resultsDirectory, comparisonNames[analysis], sep=""))
  save(methResults, file=paste(resultsDirectory, 
                               comparisonNames[analysis], 
                               "/methResults.RData", sep=""))
  # clean up unnecessary objects
  rm(methResults)
}

# Next calculate saturation for each of the bam files (for quality control)
satList<-lapply(paste(dataDirectory, sbamFileName, sep = ""), 
                function(i) {
                  MEDIPS.saturation(file = i, BSgenome = bsgenomePackageName, uniq = uniq, 
                                    extend = extend, shift = shift, window_size = ws, 
                                    chr.select = chr.select, nit = 10, nrit = 1, empty_bins = TRUE,
                                    rank = FALSE)
                })

enrichList<-lapply(paste(dataDirectory, sbamFileName, sep = ""),
                   function(i) {
                     MEDIPS.CpGenrich(file = i, BSgenome = bsgenomePackageName, uniq = uniq,
                                      extend = extend, shift = shift, chr.select = chr.select)
                   })

coverList <- lapply(paste(dataDirectory, sbamFileName, sep=""),
                    function(i) {
                      MEDIPS.seqCoverage(file =i, pattern = "CG", 
                                         BSgenome = bsgenomePackageName, 
                                         chr.select = chr.select, extend = extend, 
                                         shift = shift, uniq = uniq)
                    })


corMatrix = MEDIPS.correlation(MSets = processedBamFiles , plot = F, method = "pearson")


save(enrichList, satList, coverList, corMatrix, file = paste(resultsDirectory, "/qcLists.RData", sep=""))

#a<-paste("Coverage\n", cr$numberReadsWO, " of ", cr$numberReads, " reads (", format(cr$numberReadsWO/cr$numberReads*100, digits=3), "%) do not cover a pattern", sep="")

