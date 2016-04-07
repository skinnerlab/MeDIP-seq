## Created 5/27/2015 by Daniel Beck
## Last modified 4/5/2016

## This code is performs the MeDIP analysis. It is the second step in the analysis
## pipeline. The prepareData.R script should typically be run first. The dataNames.R
## configuration file is also used for this script.

# Load relevant libraries, custom functions, and the configuration script. Note, the
# order is important, as dataNames.R holds information about the bsgenomePackageName.
library(MEDIPS)
source("dataNames.R")
source("customFunctions.R")
library(bsgenomePackageName, character.only = T, lib.loc = genomeDirectory)
library(biomaRt)


#############################
## Read files and reformat ##
#############################
## This step reads in the sorted BAM sample files and converts them to a matrix of
## genomic windows with associated coverage (read counts). These options are defined
## in the dataNames.R configuration file.

processedBamFiles <- lapply(X = paste(dataDirectory, sbamFileName, sep = ""), 
                            FUN = MEDIPS.createSet, 
                            BSgenome = bsgenomePackageName, 
                            extend = extend, 
                            shift = shift, 
                            uniq = uniq, 
                            window_size = ws, 
                            chr.select = chr.select)


###################
## Normalization ##
###################
## This step normalizes the data by calculating local CpG density. This is set to false
## by default in our pipline. This was taken from Haque's analysis and hasn't been 
## sufficiently explored.

if (CScalc) {
  CS <- MEDIPS.couplingVector(pattern = "CG", refObj = processedBamFiles[[1]])
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
corMatrix = MEDIPS.correlation(MSets = processedBamFiles , plot = F, method = "pearson")

# The QC results are saved to a RData file in the results directory. This line should be changed
# if the commented analyses above are used (to include satList and/or enrichList).
save(coverList, corMatrix, file = paste(resultsDirectory, "/qcLists.RData", sep=""))


