## Created 6/29/2015 by Daniel Beck
## Last modified 5/31/2016

## This is the configuration script for the MeDIP-seq analysis pipeline. It holds
## most analysis parameters and options. It also holds sample/filename information
## and defines which analyses are performed. It also performs some preliminary 
## checks to ensure the analysis can procede. 


################################
## General project attributes ##
################################
## This section sets the folder structure. The default is to have a single project
## folder with three main subfolders for code, data, and results. An additional 
## genome directory holds the reference genome and is located within the data folder.
## It should be fine to modify this structure, however, extensive testing of any
## alternatives has not been done.

projectName <- "project name"

# Project directory (main folder)
projectDirectory <- "/projects/projectDirectory/"
# Data directory (holds data files)
dataDirectory <- paste(projectDirectory, "data/", sep="")
# Code directory (holds all code used for the MeDIP-seq analysis)
codeDirectory <- paste(projectDirectory, "code/", sep="")
# Results directory (folder for all results)
resultsDirectory <- paste(projectDirectory, "results/", sep="")
# Genome directory (holds reference genome and BSgenome package)
genomeDirectory <- paste(dataDirectory, "genome/", sep="")


####################
## Analysis flags ##
####################
## These are TRUE/FALSE flags for whether a particular part of the analysis should be
## performed. If these are set to false, the resulting output files should already be 
## present in the appropriate directories.

# Download the reference genome from ftp/http site?
downloadGenome <- TRUE
# Generate FastQC reports for raw datafiles?
generateFastQC <- TRUE
# Clean raw reads using default parameters?
cleanReads <- TRUE
# Use cleaned reads for mapping? (This allows the read cleaning step to be separated
# from the mapping step.)
useCleanReads <- TRUE
# Build Bowtie2 index?
buildBTindex <- TRUE
# Build BSgenome library?
buildBSgenome <- TRUE
# Map reads with Bowtie2?
mapReads <- TRUE
# Convert SAM files to sorted BAM files?
convertToBam <- TRUE
# Generate index files for BAM files?
generateIndex <- TRUE
# Download annotation files?
downloadAnnotation <- TRUE


######################
## Reference genome ##
######################
## These variables store the source of the reference genome files and associated annotation.
## Even if they aren't used, these variables should be specified to ensure the source of 
## the files is recorded and easily determinable for all projects.

# Reference genome source
genomeSourceFiles <- "ftp or http address"
# Annotation source (if an annotation file is used)
annotationSource <- "http or ftp site/name.gff.gz"


################################
## BSgenome package variables ##
################################
## These options are used for the BSgenome package creation. The seqnames parameter
## is often the most tedious (especially for large numbers of scaffolds). This will
## need to be improved in the future. It can be easier to manually create the 
## BSgenome package.

# The seed file name for BSgenome package creation
seedName <- "referenceName.seed"
# The package name
bsgenomePackageName <- "BSgenome.referenceName.Source.Version"
# Name for the reference genome
referenceName <- "Gspecies"
# Species name
species <- "Genus species"
# Source of reference genome
provider <- "Source"
# Reference genome version
version <- "Version"
# Release date
release_date <- "date"
# Names of all scaffolds or chromosomes
seqnames <- 'c("seq1", "seq2", "etc")'
# Any prefix needed to regenerate filename from seqnames
seqfiles_prefix <- ""
# Any suffix needed to regenerate filename from seqnames
seqfiles_suffix <- ".fa"


###############################
## Raw file / sample pairing ##
###############################
## This section associates files with samples. The files are assumed to be located in the
## dataDirectory folder. The full path will be based on that directory.

# The seqFiles data frame holds all sample information. The ctFlag is an identifier for the 
# treatment group, but any attribute of the sample can be used.
seqFiles <- data.frame(
  sampleName = c("S1", "S2"), 
  p1FileName = c("S1_1.fq", "S2_1.fq"),
  p2FileName = c("S1_2.fq", "S2_2.fq"),
  ctFlag = c("C", "T"),
  stringsAsFactors = F
)

# The Trimmomatic cleaned files include four files for each sample, two paired files and 
# two unpaired files. By default, the file names are generated from the raw file names in 
# seqFiles (above). If the cleaning was done manually, or if the clean data filenames are
# different than the default below, cleanFileNames needs to be modified.
cleanFileNames <- data.frame(
  cp1FileName = paste("pair_", seqFiles$p1FileName, sep=""),
  cp2FileName = paste("pair_", seqFiles$p2FileName, sep=""),
  cs1FileName = paste("sngl_", seqFiles$p1FileName, sep=""),
  cs2FileName = paste("sngl_", seqFiles$p2FileName, sep=""),
  stringsAsFactors = F
)

# Files are automatically given names when converted to SAM and sorted BAM formats. All
# files are kept. These can be changed if necessary.
samFileName <- paste(seqFiles$sampleName, ".sam", sep="")
bamFileName <- paste(seqFiles$sampleName, ".bam", sep="")
sbamFileName <- paste("s", bamFileName, sep="")


##################################
## Bowtie2 and SAMtools options ##
##################################
## These options determine how Bowtie2 maps the sample reads to the reference. Any
## Bowtie2 options can be specified using the otherMappingOptions variable.

# This defines the number of processing threads used by Bowtie2 for the mapping. It
# is also used for Trimmomatic read cleanning and trimming if necessary.
numThreads <- 10 
# All other parameters can be added here
otherMappingOptions <- "--no-unal"
# The index name can be specified, if different from referenceName due to manual creation.
indexName <- referenceName
# Multithreaded samtools may use excessive memory. Used for view and sort (conversion of
# SAM to BAM and sorted BAM.
numSamThreads <- 10


#######################
## MEDIPS parameters ##
#######################
## These parameters modify the MEDIPS analysis. The defaults shown here have not been fully
## explored. It isn't clear what the ideal values might be. Most of these are inherited from
## Haque's code. See MEDIPS documentation for details.

uniq <- 1
extend <- 50
shift <- 0
ws = 100
chr.select <- NULL
p.adj <- "fdr"
diff.method <- "edgeR"
MeDIP <- FALSE
CNV <- FALSE
CScalc <- FALSE
minRowSum <- 1

# This vector holds all raw p-value thresholds to use for the analyses. 
pValues <- c(1e-03, 1e-04, 1e-05, 1e-06, 1e-07)
# This vector holds all multiple testing adjusted p-value thresholds to use for the analyses. 
MTCpValues <- c(0.3, 0.2, 0.1, 0.05)

####################
## DMR parameters ##
####################
## These parameters define how DMR edges are defined. These values are completely arbitrarily
## chosen. It isn't clear how these values should be chosen to accuratly reflect a biologically
## significant feature.

# This p-value threshold defines DMR boundaries
dmrBoundPvalue <- 0.1
# Adjacency distance (=1 when windows must be exactly adjacent). This determines how far apart
# significant windows can be and remain in the same DMR.
adjDist <- 1000

# The maxDMRnum variable gives a maximum number of DMRs on which to calculate CpG density and other
# information. This speeds up the pipeline. However, this will need to be increased if the p-value 
# of interest has more DMR than this number.
maxDMRnum <- 5000

###########################
## Annotation parameters ##
###########################
## These variables hold information on the annotation. Three types of annotation are possible.
## The annotation is still being incorporated into the pipeline.

# The annotation can be from Biomart ("biomart"), from a GFF file ("gff"), or from BLAST results
# ("blast").
annotationType <- "gff"
# If annotation is from a GFF file, include the filename.
annotationGFF <- ""
# If the chromosomes in the annotation file have an additional prefix, specify here.
chrPrefix <- ""

# If annotation is from Biomart, include biomaRt dataset information.
biomartHost <- ""
biomartDataset <- ""


############################
## Comparisons / analyses ##
############################
## This section defines which comparisons should be made. All analyses will iterate through
## each of these.

# Names for the comparisons
comparisonNames <- 
  c("all",
    "pair_1-2", "pair_1-3", "pair_2-3",
    "pair_4-5", "pair_4-6", "pair_5-6",
    "pair_1-4", "pair_1-5", "pair_1-6",
    "pair_2-4", "pair_2-5", "pair_2-6",
    "pair_3-4", "pair_3-5", "pair_3-6"
  )

## This vector should be a subset of comparisonNames and will be used to perform an APO 
## analysis using the identifyApoDmr.R script.
pair.analysis.names <- c("pair_1-4", "pair_2-5", "pair_3-6")

# Which samples are being compared. The pairs flag was included to allow for a pairwise
# type analysis. It isn't currently functional, but some code requires it.
comparison <- list()
comparison[[1]] <- data.frame(mset1=c(1:3), mset2=c(4:6), pairs=F)
comparison[[2]] <- data.frame(mset1=c(1), mset2=c(2), pairs=T)
comparison[[3]] <- data.frame(mset1=c(1), mset2=c(3), pairs=T)
comparison[[4]] <- data.frame(mset1=c(2), mset2=c(3), pairs=T)
comparison[[5]] <- data.frame(mset1=c(4), mset2=c(5), pairs=T)
comparison[[6]] <- data.frame(mset1=c(4), mset2=c(6), pairs=T)
comparison[[7]] <- data.frame(mset1=c(5), mset2=c(6), pairs=T)
comparison[[8]] <- data.frame(mset1=c(1), mset2=c(4), pairs=T)
comparison[[9]] <- data.frame(mset1=c(1), mset2=c(5), pairs=T)
comparison[[10]] <- data.frame(mset1=c(1), mset2=c(6), pairs=T)
comparison[[11]] <- data.frame(mset1=c(2), mset2=c(4), pairs=T)
comparison[[12]] <- data.frame(mset1=c(2), mset2=c(5), pairs=T)
comparison[[13]] <- data.frame(mset1=c(2), mset2=c(6), pairs=T)
comparison[[14]] <- data.frame(mset1=c(3), mset2=c(4), pairs=T)
comparison[[15]] <- data.frame(mset1=c(3), mset2=c(5), pairs=T)
comparison[[16]] <- data.frame(mset1=c(3), mset2=c(6), pairs=T)


#################
## Data checks ##
#################
## These are very basic checks to detect problems in this configuration file. This section
## should be extended each time an unexpected error occurs in this script.

if (length(comparison)!=length(unique(comparisonNames))) { 
  warning("The number of comparison names does not match the number of comparisons.\n",
          "This could be caused by duplicate comparison names.")
}
if (CScalc & !MeDIP) {
  warning("The CSset coupling set is being calculated (CScalc = TRUE) 
          but not used (MeDIP = FALSE).")
}
if (!CScalc & MeDIP) {
  warning("The CSset coupling set is not calculated (CScalc = FALSE) 
          but is needed (MeDIP = TRUE).")
}
