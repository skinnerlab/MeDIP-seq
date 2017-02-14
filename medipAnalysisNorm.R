## Created 5/27/2015 by Daniel Beck
## Last modified 4/5/2016

## This code is performs the MeDIP analysis. It is the second step in the analysis
## pipeline. The prepareData.R script should typically be run first. The dataNames.R
## configuration file is also used for this script.

# Load relevant libraries, custom functions, and the configuration script. Note, the
# order is important, as dataNames.R holds information about the bsgenomePackageName.
library(MEDIPS)
library(RUVSeq)
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



# Identify control windows

  ## This was done separately. It should be incorporated into this pipeline eventually.
  ## This will require CpG counts genome wide, and perhaps the addition of a specific
  ## script and datastructure to avoid recalculating this for ever genome.
  # load("/projects/ratVinclozolin/results/allF1/methResults2.RData")
  # numCpG <- mrwCpG$cpgNum
  # chr <- mrwCpG$chr
  # cwin <- identifyZero(chr=chr, numCpG=numCpG, neighborhood=10, max=0)
  # save(cwin, file="zeroWin_100_10_0.Rdata")

load("zeroWin_100_10_0.Rdata")



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
  
  ## Create list of offsets per sample
  
  # Combine counts for every sample into a data.frame
  countDF <- do.call(cbind, lapply(c(mset1, mset2), genome_count))
  colnames(countDF) <- sapply(c(mset1, mset2), sample_name)
  row.names(countDF) <- paste("w", 1:nrow(countDF), sep="")
  
  # Base table (taken from MEDIPS code)
  window_size = window_size(mset1[[1]])
  no_chr_windows = ceiling(chr_lengths(mset1[[1]])/window_size(mset1[[1]]))
  supersize_chr = cumsum(no_chr_windows)	
  GRanges.genome = MEDIPS.GenomicCoordinates(supersize_chr, no_chr_windows, chr_names(mset1[[1]]), chr_lengths(mset1[[1]]), window_size(mset1[[1]]))
  base = data.frame(chr=as.vector(seqnames(GRanges.genome)), start=start(GRanges.genome), stop=end(GRanges.genome), stringsAsFactors=F)
  
  counts.medip = NULL
  rpkm.medip = NULL
  ##Add counts
  if(!is.null(mset1)){
    for(i in 1:length(mset1)){
      cat(paste("Preprocessing MEDIPS SET ", i, " in MSet1...\n", sep=""))
      counts.medip = cbind(counts.medip, MSet1=genome_count(mset1[[i]]))
      rpkm.medip = cbind(rpkm.medip, round(((genome_count(mset1[[i]])*10^9)/(window_size*number_regions(mset1[[i]]))), digits=2))
    }
  }	
  if(!is.null(mset2)){
    for(i in 1:length(mset2)){
      cat(paste("Preprocessing MEDIPS SET ", i, " in MSet2...\n", sep=""))
      counts.medip = cbind(counts.medip, MSet2=genome_count(mset2[[i]]))
      rpkm.medip = cbind(rpkm.medip, round(((genome_count(mset2[[i]])*10^9)/(window_size*number_regions(mset2[[i]]))), digits=2))

    }
  }	
  cn <- c(seqFiles$sampleName[comparison[[analysis]]$mset1], seqFiles$sampleName[comparison[[analysis]]$mset2])
  colnames(counts.medip) <- paste(cn, ".counts", sep="")
  colnames(rpkm.medip) <- paste(cn, ".rpkm", sep="")
  base <- cbind(base, counts.medip, rpkm.medip)
  
  # Filter rows of both the count table and the control window vector
  filter <- which(rowSums(countDF) >= minRowSum)
  fcountDF <- countDF[filter,]
  fcwin <- cwin[filter]
  
  genes <- rownames(fcountDF)
  zeros <- rownames(fcountDF[which(fcwin==1),])
  
  # Create SeqExpressionSet object holding genomic window read counts (modified from RUVSeq documents)
  x <- as.factor(c(seqFiles$ctFlag[comparison[[analysis]]$mset1], 
                   seqFiles$ctFlag[comparison[[analysis]]$mset2]))
  set <- newSeqExpressionSet(as.matrix(fcountDF), phenoData = data.frame(x, row.names=colnames(fcountDF)))

  # Upper quantile normalization (for library size differences)
  set <- betweenLaneNormalization(set, which="upper")

  # Calculate offsets using RUVseq (RUVg) with zeros as background
  set1 <- RUVg(set, zeros, k=1)

  # plotRLE(set1, outline=F, ylim=c(-4,4), col=colors[x])
  # plotPCA(set1, col=colors[x], cex=1.2)

  # Convert offsets and counts to DGEList for analysis with edgeR
  design <- model.matrix(~x + W_1, data=pData(set1))
  y <- DGEList(counts=counts(set1), group=x)
  
  # Perform analysis
  y <- calcNormFactors(y, method="upperquartile")
  y <- estimateGLMCommonDisp(y, design)
  y <- estimateGLMTagwiseDisp(y, design)
  
  fit <- glmFit(y, design)
  lrt <- glmLRT(fit, coef=2)
  
  mr <- topTags(lrt, n=nrow(lrt$table), sort.by="none")

  # Reformat results
  methResults <- data.frame(base, "logFC"=NA, "logCPM"=NA, "LR"=NA, "edgeR.p.value"=NA, "edgeR.adj.p.value"=NA)
  methResults[filter, (ncol(methResults)-4):ncol(methResults)] <- mr$table

  # Save results to a comparison specific folder in the results directory
  system(paste("mkdir ", resultsDirectory, comparisonNames[analysis], sep=""))
  save(methResults, set1, file=paste(resultsDirectory, 
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


