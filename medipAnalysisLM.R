## Created 5/27/2015 by Daniel Beck
## Last modified 7/20/2016

## This code is performs the MeDIP analysis. It is the second step in the analysis
## pipeline. The prepareData.R script should typically be run first. The dataNames.R
## configuration file is also used for this script. This version of the script reads
## in the bam files for each comparison, lowering memory use and slowing down the 
## analysis considerably in some cases. This should only be used when memory use must
## be kept low.

# Load relevant libraries, custom functions, and the configuration script. Note, the
# order is important, as dataNames.R holds information about the bsgenomePackageName.
library(MEDIPS)
source("dataNames.R")
source("customFunctions.R")
library(bsgenomePackageName, character.only = T, lib.loc = genomeDirectory)


###################
## Identify DMRs ##
###################
## This step performs the actual DMR analysis. Each genomic window is assigned a 
## probability of being a DMR. All samples in mset1 are compared with all samples 
## in mset2. This code loops over all comparisons specified in the dataNames.R 
## configuration script.

processedBamFiles1 <- list()
processedBamFiles2 <- list()
for (analysis in 1:length(comparison)) {
  
  #############################
  ## Read files and reformat ##
  #############################
  ## This step reads in the sorted BAM sample files and converts them to a matrix of
  ## genomic windows with associated coverage (read counts). These options are defined
  ## in the dataNames.R configuration file.
  
  read.mset1 <- TRUE
  read.mset2 <- TRUE
  if (analysis > 1) {
    read.mset1 <- !Reduce("&", comparison[[analysis]]$mset1 == comparison[[analysis-1]]$mset1)
    read.mset2 <- !Reduce("&", comparison[[analysis]]$mset2 == comparison[[analysis-1]]$mset2)
  }
  if (read.mset1) {
    processedBamFiles1 <- lapply(X = paste(dataDirectory, sbamFileName[comparison[[analysis]]$mset1], sep = ""), 
                                 FUN = MEDIPS.createSet, 
                                 BSgenome = bsgenomePackageName, 
                                 extend = extend, 
                                 shift = shift, 
                                 uniq = uniq, 
                                 window_size = ws, 
                                 chr.select = chr.select)
  }
  if (read.mset2) {
    processedBamFiles2 <- lapply(X = paste(dataDirectory, sbamFileName[comparison[[analysis]]$mset2], sep = ""), 
                                 FUN = MEDIPS.createSet, 
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
  ## sufficiently explored. This only needs to be done once.
  if (analysis == 1) {
    if (CScalc) {
      CS <- MEDIPS.couplingVector(pattern = "CG", refObj = processedBamFiles1[[1]])
    } else {
      CS <- NULL
    }
  }
  
  # I'm not sure when this would occur. I think it was an early test.
  if (length(processedBamFiles1)==0) {
    processedBamFiles1 <- NULL
  }
  if (length(processedBamFiles2)==0) {
    processedBamFiles2 <- NULL
  }
  # Perform analysis
  methResults <- MEDIPS.meth(MSet1 = processedBamFiles1, 
                             MSet2 = processedBamFiles2,
                             p.adj = p.adj,
                             diff.method = diff.method,
                             MeDIP = MeDIP,
                             CNV = CNV,
                             CSet = CS,
                             minRowSum = minRowSum)
  # Save results to a comparison specific folder in the results directory
  system(paste("mkdir ", resultsDirectory, comparisonNames[analysis], sep=""))
  save(methResults, file=paste(resultsDirectory, 
                               comparisonNames[analysis], 
                               "/methResults.RData", sep=""))
  # Clean up unnecessary objects
  rm(methResults)
}


#####################
## Quality control ##
#####################
## This section calculates various quality control metrics and visualizations. While 
## several are included here, many of them take considerable time. I've commented out
## the more computationally intensive ones. These can be uncommented if more extensive
## QC is required.

## Next calculate saturation for each of the bam files (for quality control)
# satList<-lapply(paste(dataDirectory, sbamFileName, sep = ""), 
#                function(i) {
#                  MEDIPS.saturation(file = i, BSgenome = bsgenomePackageName, uniq = uniq, 
#                                    extend = extend, shift = shift, window_size = ws, 
#                                    chr.select = chr.select, nit = 10, nrit = 1, 
#                                    empty_bins = TRUE, rank = FALSE)
#                })

## Next calculate CpG enrichment
# enrichList<-lapply(paste(dataDirectory, sbamFileName, sep = ""),
#                   function(i) {
#                     MEDIPS.CpGenrich(file = i, BSgenome = bsgenomePackageName, uniq = uniq,
#                                      extend = extend, shift = shift, chr.select = chr.select)
#                   })

# This looks at coverage levels for the reference genome CpGs.
coverList <- lapply(paste(dataDirectory, sbamFileName, sep=""),
                    function(i) {
                      MEDIPS.seqCoverage(file =i, pattern = "CG", 
                                         BSgenome = bsgenomePackageName, 
                                         chr.select = chr.select, extend = extend, 
                                         shift = shift, uniq = uniq)
                    })

# This measures the correlation in read depth between samples
# corMatrix = MEDIPS.correlation(MSets = processedBamFiles , plot = F, method = "pearson")

# The QC results are saved to a RData file in the results directory. This line should be changed
# if the commented analyses above are used (to include satList and/or enrichList).
save(coverList, file = paste(resultsDirectory, "/qcLists.RData", sep=""))


