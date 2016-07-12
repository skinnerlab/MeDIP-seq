## Created 6/29/2015 by Daniel Beck
## Last modified 4/5/2016

## This script prepares the raw data for analysis by medipAnalysis.R. It may include 
## raw file cleaning, bowtie2 index generation, mapping, and BSgenome creation. This 
## script requires a configuration file (dataNames.R) that defines filenames, locations,
## and other parameter values.

# Load relevant libraries and configuration scripts
library(BSgenome)
source("dataNames.R")


#####################
## Download Genome ##
#####################
## This section allows for genome and/or annotation download from http or ftp sites.
## There is currently no check for download success or file integrety. This will need to 
## be checked manually if there are any concerns.

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
## This section generates the FastQC reports for each sample file. Seperate reports are 
## generated for each paired end. This step is currently performed prior to data cleaning
## and filtering. This may need to be moved to after cleaning depending on question.

if (generateFastQC) {
  for (i in 1:length(seqFiles$sampleName)) {
    system(paste("fastqc ", dataDirectory, seqFiles$p1FileName[i], sep=""))
    system(paste("fastqc ", dataDirectory, seqFiles$p2FileName[i], sep=""))
  }
}


########################
## Clean/filter reads ##
########################
## This section runs Trimmomatic to clean and filter raw sequence reads. The adapter file is
## curently hardcoded. This may need to be changed in the future to allow different adapters
## to be specified in dataNames.R. Trimmomatic options are also hardcoded and are based on 
## options used by the Arizona sequencing lab.

for (i in 1:length(seqFiles$sampleName)) {
  if (cleanReads) {
    # run Trimmomatic
    system(paste("java -jar /apps/Trimmomatic-0.36/trimmomatic-0.36.jar ", 
                 "PE ", "-threads ", numThreads, " -phred33 ", 
                 dataDirectory, seqFiles$p1FileName[i], " ", 
                 dataDirectory, seqFiles$p2FileName[i], " ", 
                 dataDirectory, cleanFileNames$cp1FileName[i], " ", 
                 dataDirectory, cleanFileNames$cs1FileName[i], " ", 
                 dataDirectory, cleanFileNames$cp2FileName[i], " ", 
                 dataDirectory, cleanFileNames$cs2FileName[i], " ",
                 "ILLUMINACLIP:/apps/Trimmomatic-0.36/adapters/TruSeq3-PE-2.fa:2:30:11:1:true ",
                 "TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:20", sep=""))
  }
}

## The command above is hard coded for several Trimmomatic parameters and options. These are
## summarized below. The specific paramter values are copied directly from the Arizon lab
## procedures.

# java -jar trimmomatic-0.36.jar              -- Trimmomatic JAR file
#  PE                                         -- Paired end input files
#  -threads                                   -- Number of computation threads to use
#  -phred33                                   -- Specifies the base quality encoding
#  Sample_R1.fq.gz                            -- Sample R1 file name (input)
#  Sample_R2.fq.gz                            -- Sample R2 file name (input)
#  pair_Sample_R1.fq.gz                       -- Cleaned R1 pairs file name (output)
#  sngl_Sample_R1.fq.gz                       -- Cleaned R1 singles file name (output)
#  pair_Sample_R2.fq.gz                       -- Cleaned R2 pairs file name (output)
#  sngl_Sample_R2.fq.gz                       -- Cleaned R2 singles file name (output)
#  ILLUMINACLIP                               -- Remove adapters
#      :TruSeq3-PE-2.fa                       -- Fasta file with adapters
#      :2                                     -- Maximum seed mismatches
#      :30                                    -- Palindrome clip threshold
#      :11                                    -- Simple clip threshold
#      :1                                     -- Minimum adapter length
#      :true                                  -- Keep both reads
#  TRAILING:3                                 -- Minimum quality required to keep base
#  SLIDINGWINDOW:4:15                         -- windowSize:requiredQuality
#  MINLEN:20                                  -- Minimum length of reads to keep


###############
## Map reads ##
###############
## This section maps the reads to the reference genome using Bowtie2. Two steps are required.
## First, the Bowtie2 index is built from the reference genome. Second, the sample reads are 
## mapped to the reference. This step also converts the mapped read files from SAM to a sorted
## BAM format using Samtools. All mapping options are set in dataNames.R. Similarily, each of
## the steps performed here are conditional based on flags set in dataNames.R.

setwd(genomeDirectory)
if (buildBTindex) {
# build bowtie2 index
     system(paste("bowtie2-build -f ", 
                  paste(paste(genomeDirectory, 
                              seqfiles_prefix, 
                              eval(parse(text=seqnames)), 
                              ".fa",
                              sep=""), 
                        collapse=","), 
                  " ", 
                  indexName, 
                  sep=""))
}

# map illumina reads to reference using Bowtie2. 
# convert resulting SAM files to sorted BAM files
for (i in 1:length(seqFiles$sampleName)) {
  if (mapReads) {
    print(paste("Running Bowtie2 on", seqFiles$sampleName[i], sep=" "))
    # Run bowtie2 mapping. If cleaning reads, use cleaned names. Otherwise, use original files.
    if (useCleanReads) {
      system(paste("bowtie2 -x ", indexName, 
                   " ", otherMappingOptions, 
                   " -p ", numThreads, 
                   " -1 ", dataDirectory, cleanFileNames$cp1FileName[i], 
                   " -2 ", dataDirectory, cleanFileNames$cp2FileName[i], 
                   " -U ", dataDirectory, cleanFileNames$cs1FileName[i], ",", 
                           dataDirectory, cleanFileNames$cs2FileName[i],
                   " -S ", dataDirectory, samFileName[i],
                   sep=""))
    } else {
      # Use original file names if cleaned read files are not used
      system(paste("bowtie2 -x ", indexName, 
                   " ", otherMappingOptions, 
                   " -p ", numThreads, 
                   " -1 ", dataDirectory, seqFiles$p1FileName[i], 
                   " -2 ", dataDirectory, seqFiles$p2FileName[i], 
                   " -S ", dataDirectory, samFileName[i],
                   sep=""))
    }
  }
  
  if (convertToBam) {
    # convert to sorted BAM format using samtools
    system(paste("samtools view -@ ", numSamThreads, 
                 " -bS ", dataDirectory, samFileName[i], 
                 " > ", dataDirectory, bamFileName[i], sep=""))
    system(paste("samtools sort -@ ", numSamThreads, 
                 " ", dataDirectory, bamFileName[i], 
                 " ", dataDirectory, paste("s",seqFiles$sampleName[i], sep=""), 
                 sep=""))
  }

  if (generateIndex) {
    system(paste("samtools index ", 
                 dataDirectory, sbamFileName[i], " ", 
                 dataDirectory, paste(sbamFileName[i], ".bai", sep=""), sep=""))
  }
}

####################
## Forge BSgenome ##
####################
## This section builds the BSgenome package, installs it locally, and loads it. This is a 
## slightly convoluted way of doing things, but is generic, requiring only a specific file 
## structure. This structure should likely be done for every project for consistency anyway. 
## If this becomes a problem, the paths can be abstracted to variables. Nearly all options
## to the BSgenome package creation can be specified in the dataNames.R configuration file.

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
  system(paste("R CMD INSTALL -l ", genomeDirectory, 
               " ", bsgenomePackageName, "_1.0.0.tar.gz", sep=""))
}
