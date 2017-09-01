## Created 5/27/2015 by Daniel Beck
## Last modified 9/2/2016

## This is a modification of the medipAnalysis.R script. It does not look for
## windows with differential coverage, it just counts the number of reads in
## each sample for every genomic window. These results can be used to look at
## read depth.

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

processedBamFiles <- lapply(X = paste(dataDirectory, sbamFileName, sep = ""), 
                            FUN = MEDIPS.createSet, 
                            BSgenome = bsgenomePackageName, 
                            extend = extend, 
                            shift = shift, 
                            uniq = uniq, 
                            window_size = ws, 
                            chr.select = chr.select)


###################
## 

  
## Create list of offsets per sample
for (analysis in 1:length(comparisonNames)) {
  
  # Combine counts for every sample into a data.frame
  countDF <- do.call(cbind, lapply(processedBamFiles, genome_count))
  colnames(countDF) <- sapply(processedBamFiles, sample_name)
  row.names(countDF) <- paste("w", 1:nrow(countDF), sep="")
  
  # Base table (taken from MEDIPS code)
  window_size = window_size(processedBamFiles[[1]])
  no_chr_windows = ceiling(chr_lengths(processedBamFiles[[1]])/window_size(processedBamFiles[[1]]))
  supersize_chr = cumsum(no_chr_windows)	
  GRanges.genome = MEDIPS.GenomicCoordinates(supersize_chr, no_chr_windows, chr_names(processedBamFiles[[1]]), chr_lengths(processedBamFiles[[1]]), window_size(processedBamFiles[[1]]))
  base = data.frame(chr=as.vector(seqnames(GRanges.genome)), start=start(GRanges.genome), stop=end(GRanges.genome), stringsAsFactors=F)
  
  counts = NULL
  rpkm = NULL
  ##Add counts
  for(i in 1:length(processedBamFiles)){
    cat(paste("Calculating RPKM ", i, " in MSet1...\n", sep=""))
    counts = cbind(counts, MSet1=genome_count(processedBamFiles[[i]]))
    rpkm = cbind(rpkm, round(((genome_count(processedBamFiles[[i]])*10^9)/(window_size*number_regions(processedBamFiles[[i]]))), digits=2))
  }
  
  cn <- seqFiles$sampleName
  colnames(counts) <- paste(cn, ".counts", sep="")
  colnames(rpkm) <- paste(cn, ".rpkm", sep="")
  base <- cbind(base, counts, rpkm)
  
  # Save results to a comparison specific folder in the results directory
  system(paste("mkdir ", resultsDirectory, comparisonNames[analysis], sep=""))
  save(base,file=paste(resultsDirectory, 
                       comparisonNames[analysis], 
                       "/readDepths.RData", sep=""))
  
}
