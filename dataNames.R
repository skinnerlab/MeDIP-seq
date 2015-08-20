## Daniel Beck
## Created 6/29/2015
## Modified
##	6/30/2015 Continued development. Ready for testing.
##	7/2/2015	Modified for the steelheadTrial dataset
##   7/6/2015  Code failed due to combined genome fasta file. I fixed this manually and am rerunning this script without the parts that completed sucessfully. Modified for use on steelheadSperm project. Added flags for analysis customization
##   7/9/2015  Added multithreading variable for samtools. Added annotation and result director variables for updated medipAnalysis.R script.
##   7/10/2015 Added CS normalization option with data checks.
##   7/13/2015 Removed project specific information and added file to git repository "medipPipeline"

## This file contains file names, file locations, and analysis options/constants
## 


################################
## General project attributes ##
################################
projectName<-"project name"

projectDirectory<-"/projects/projectDirectory/"
dataDirectory<-paste(projectDirectory, "data/", sep="")
codeDirectory<-paste(projectDirectory, "code/", sep="")
resultsDirectory<-paste(projectDirectory, "results/", sep="")
genomeDirectory<-paste(dataDirectory, "genome/", sep="")


####################
## Analysis flags ##
####################
# If these are set to false, the resulting output files should already be present in the appropriate directories.

# download genome from ftp/http site
downloadGenome<-FALSE
# generate fastqc report for raw datafiles
generateFastQC<-FALSE
# build bowtie2 index
buildBTindex<-FALSE
# build BSgenome library
buildBSgenome=FALSE
# map reads with bowtie2
mapReads=FALSE
# convert SAM files to sorted BAM files
convertToBam=FALSE
# download annotation files
downloadAnnotation=FALSE

genomeSourceFiles<-"ftp or http address"

################################
## BSgenome package variables ##
################################

# seed file name for BSgenome package creation
seedName<-"referenceName.seed"
bsgenomePackageName<-"BSgenome.referenceName.Source.Version"
referenceName<-"Gspecies"
species<-"Genus species"
provider<-"Source"
version<-"Version"
release_date<-"date"
seqnames<-'c("seq1", "seq2", "etc")'
seqfiles_prefix<-""
seqfiles_suffix<-".fa"

###############################
## Raw file / sample pairing ##
###############################

# Sequence file names associated with sample names. 
seqFiles<-data.frame(
  sampleName=c("S1", "S2"), 
  p1FileName=c("S1_1.fq", "S2_1.fq"),
  p2FileName=c("S1_2.fq", "S2_2.fq"),
  ctFlag=c("C", "T"),
  stringsAsFactors=F
)

samFileName<-paste(seqFiles$sampleName, ".sam", sep="")
bamFileName<-paste(seqFiles$sampleName, ".bam", sep="")
sbamFileName<-paste("s", bamFileName, sep="")

########################
## Mapping parameters ##
########################
numThreads<-25
otherMappingOptions<-"--no-unal"
indexName<-referenceName

#########################
## Samtools parameters ##
#########################

# Multithreaded samtools may use excessive memory. Used for view and sort.
numSamThreads<-25

#######################
## MEDIPS parameters ##
#######################

uniq <- TRUE
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

# This vector holds all p-value thresholds to use for the analyses
pValues<-c(5e-02, 1e-02, 1e-03, 1e-04, 1e-05, 1e-06, 1e-07)

# allow for multiple comparisons. loop over comparison list in medipAnalysis.R
comparisonNames<-c(
     "all",
     "pair_1-2", "pair_1-3", "pair_2-3",
     "pair_4-5", "pair_4-6", "pair_5-6",
     "pair_1-4", "pair_1-5", "pair_1-6",
     "pair_2-4", "pair_2-5", "pair_2-6",
     "pair_3-4", "pair_3-5", "pair_3-6"
)

comparison<-list()
## all comparisons
comparison[[1]]<-data.frame(mset1=c(1:3), mset2=c(4:6), pairs=F)
comparison[[2]]<-data.frame(mset1=c(1), mset2=c(2), pairs=T)
comparison[[3]]<-data.frame(mset1=c(1), mset2=c(3), pairs=T)
comparison[[4]]<-data.frame(mset1=c(2), mset2=c(3), pairs=T)
comparison[[5]]<-data.frame(mset1=c(4), mset2=c(5), pairs=T)
comparison[[6]]<-data.frame(mset1=c(4), mset2=c(6), pairs=T)
comparison[[7]]<-data.frame(mset1=c(5), mset2=c(6), pairs=T)
comparison[[8]]<-data.frame(mset1=c(1), mset2=c(4), pairs=T)
comparison[[9]]<-data.frame(mset1=c(1), mset2=c(5), pairs=T)
comparison[[10]]<-data.frame(mset1=c(1), mset2=c(6), pairs=T)
comparison[[11]]<-data.frame(mset1=c(2), mset2=c(4), pairs=T)
comparison[[12]]<-data.frame(mset1=c(2), mset2=c(5), pairs=T)
comparison[[13]]<-data.frame(mset1=c(2), mset2=c(6), pairs=T)
comparison[[14]]<-data.frame(mset1=c(3), mset2=c(4), pairs=T)
comparison[[15]]<-data.frame(mset1=c(3), mset2=c(5), pairs=T)
comparison[[16]]<-data.frame(mset1=c(3), mset2=c(6), pairs=T)


####################
## DMR parameters ##
####################

# p-value threshold for defining DMR boundaries
dmrBoundPvalue<-0.1
# Adjacency distance (=1 when windows must be exactly adjacent)
adjDist<-1

# The maxDMRnum variable gives a maximum number of DMRs on which to calculate CpG density and other information. This speeds up the pipeline. However, this will need to be increased if the p-value of interest has more DMR than this number.
maxDMRnum<-5000

###########################
## Annotation parameters ##
###########################
# annotationType can be from Biomart ("biomart") or GFF file ("gff")
annotationType<-"gff"
# if annotation is from Biomart, include biomaRt dataset information, otherwise use gff file name
# if GFF file, include file name
annotationSource<-"http or ftp site/name.gff.gz"
annotationGFF<-"name.gff"
chrPrefix<-""

# if annotation is from Biomart, include host and dataset
biomartHost<-""
biomartDataset<-""

#################
## Data checks ##
#################
if (length(comparison)!=length(unique(comparisonNames))){ 
  warning("The number of comparison names does not match the number of comparisons.\n",
          "This could be caused by duplicate comparison names.")
}
if (CScalc & !MeDIP){
  warning("The CSset coupling set is being calculated (CScalc = TRUE) but not used (MeDIP = FALSE).")
}
if (!CScalc & MeDIP){
  warning("The CSset coupling set is not calculated (CScalc = FALSE) but is needed (MeDIP = TRUE).")
}