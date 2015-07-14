## Daniel Beck
## Created 6/29/2015
## Modified
##   6/30/2015 Continued initial development. Script ready for testing.
##   7/1/2015  Fixed misnamed variable error.
##   7/6/2015  Added flags for analysis customization
##   7/8/2015  Added annotation gff file download if necessary
##   7/9/2015  Added multithreading for samtools view and sort

## This script prepares the raw data for analysis by medipAnalysis.R. It may include raw file cleaning, bowtie2-build index generation, mapping, BSgenome creation, etc.
## This script requires a dataset specific R script that defines filenames, locations, etc.

# load relevant libraries
library(BSgenome)

# source dataset specific code
source("dataNames.R")

#####################
## Download Genome ##
#####################
if (downloadGenome) {
  system(paste("wget -nH -t inf -P ", genomeDirectory, " ", genomeSourceFiles, sep=""))
  system(paste("gzip -d ", genomeDirectory, "*.gz", sep=""))
}

if (downloadAnnotation) {
     system(paste("wget -nH -t inf -P ", genomeDirectory, " ", annotationSource, sep=""))
     system(paste("gzip -d ", genomeDirectory, "*.gz", sep=""))
}


############################
## Generate FastQC report ##
############################
if (generateFastQC){
  for (i in 1:length(seqFiles$sampleName)){
    system(paste("fastqc ", dataDirectory, seqFiles$p1FileName[i], sep=""))
    system(paste("fastqc ", dataDirectory, seqFiles$p2FileName[i], sep=""))
  }
}

###############
## Map reads ##
###############
setwd(genomeDirectory)
if (buildBTindex){
# build bowtie2 index
     system(paste("bowtie2-build -f ", paste(paste(genomeDirectory, seqfiles_prefix, eval(parse(text=seqnames)), ".fa",sep=""), collapse=","), " ", indexName, sep=""))
}
# allign illumina reads to reference using Bowtie2. Convert resulting SAM files to sorted BAM files
for (i in 1:length(seqFiles$sampleName)){
     if (mapReads){
          print(paste("Running Bowtie2 on", seqFiles$sampleName[i], sep=" "))
          system(paste("bowtie2 -x ", indexName, " ", otherMappingOptions, " -p ",numThreads, " -1 ", dataDirectory, seqFiles$p1FileName[i], " -2 ", dataDirectory, seqFiles$p2FileName[i], " -S ", dataDirectory, samFileName[i],sep=""))
     }
     if (convertToBam){
  # convert to sorted BAM format using samtools
          system(paste("samtools view -@ ", numSamThreads, " -bS ", dataDirectory, samFileName[i], " > ", dataDirectory, bamFileName[i], sep=""))
          system(paste("samtools sort -@ ", numSamThreads, " ", dataDirectory, bamFileName[i], " ", dataDirectory, paste("s",seqFiles$sampleName[i], sep=""), sep=""))
     }
}

####################
## Forge BSgenome ##
####################

# Build BSgenome, install it locally, then load it. This is a convoluted way of doing things, but hopefully it will work on the server. It is also generic, requireing only a specific file structure. This structure should likely be done for every project for consistency anyway. If this becomes a problem, the paths can be abstracted to variables.
if (buildBSgenome){
  # There must be a better way of creating a package than writing to file then reading from file.
  seedFile<-c(
    paste("Package: ", bsgenomePackageName, sep=""),
    paste("Title: ", species, " genome", sep=""),
    paste("Description: ", species, " genome downloaded from: ", genomeSourceFiles, sep=""),
    paste("Version: 1.0.0", sep=""),
    paste("organism: ", species, sep=""),
    paste("common_name: ", species, sep=""),
    paste("organism_biocview: ", referenceName, sep=""),
    paste("provider: ", provider, sep=""),
    paste("provider_version: ", version, sep=""),
    paste("release_date: ", release_date, sep=""),
    paste("release_name: ", species, " genome", sep=""),
    paste("source_url: ", genomeSourceFiles, sep=""),
    paste("BSgenomeObjname: ", referenceName, sep=""),
    paste("seqnames: ", seqnames, sep=""),
    paste("SrcDataFiles: ", genomeSourceFiles, sep=""),
    paste("seqs_srcdir: ", genomeDirectory, sep=""),
    paste("seqfiles_prefix: ", seqfiles_prefix, sep=""),
    paste("seqfiles_suffix: ", seqfiles_suffix, sep="")

  )
  writeLines(seedFile, paste(genomeDirectory, seedName, sep=""))
  forgeBSgenomeDataPkg(paste(genomeDirectory, seedName, sep=""), destdir=genomeDirectory)
  system(paste("R CMD build ", genomeDirectory, bsgenomePackageName, sep=""))
  system(paste("export R_LIBS=", genomeDirectory, sep=""))
  system(paste("R CMD INSTALL -l ", genomeDirectory, " ", bsgenomePackageName, "_1.0.0.tar.gz", sep=""))
}
