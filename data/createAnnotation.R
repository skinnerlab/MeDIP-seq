## Created 2/25/2016 by Daniel Beck
## Modified 2/26/2016

## This script builds a table with custom annotation information. 
## It will be used by the medipPipeline to annotate DMR.

####
## Load data from various sources
####

# lab category -> gene name (This is a curated lab-specific list)
symCat<-read.csv("symbolCategory_v2.csv", stringsAsFactors=F)

# gene -> homolog (This is the homoloGene database from NCBI)
homo <- read.table("homologene_v1.data", sep="\t", stringsAsFactors=F)

# panther category -> lab category (This was manually created by Mike)
pantherCat<-read.csv("pantherCategories.csv", stringsAsFactors=F)

# gene name -> panther category (This is downloaded from the DAVID database)
davidAll<-read.table("davidAll.txt", sep="\t", head=T, quote="", stringsAsFactors=F, fill=T)


####
## Generate vectors with needed information in order
####

# convert panter_mf_all to lab categories
davidAllAdd <- sapply(davidAll$PANTHER_MF_ALL, function(i) {
  pcats <- unlist(strsplit(i, split=","))
  lcats <- pantherCat$Classification[match(pcats, pantherCat$panther_mf_all)]
  i <- paste(unique(lcats), collapse=";")
})
davidAll <- cbind(davidAll, davidAllAdd)

# add lab categories to homolog table
chomo <- cbind(homo, symCat[match(homo$V4, symCat$symbol, nomatch=NA),])
firstOGS <- sapply(strsplit(davidAll$OFFICIAL_GENE_SYMBOL, split=","), function(i) i[1])
chomo <- cbind(chomo, davidAll[match(homo$V4, firstOGS, nomatch=NA),])

# rows of created table will be entries in this list
shomo <- split.data.frame(chomo, f=chomo$V1)

# extract all needed information for annotationTable
homologNumber <- sapply(shomo, function(i) i$V1[1])
humanSymbol <- sapply(shomo, function(i) unique(unlist(strsplit(i$OFFICIAL_GENE_SYMBOL, split=","))))
humanSymbol <- sapply(sapply(humanSymbol, function(i) i[which(!is.na(i))]), paste, collapse=";")
humanSummary <- sapply(shomo, function(i) unique(i$ENTREZ_GENE_SUMMARY))
humanSummary <- sapply(sapply(humanSummary, function(i) i[which(!is.na(i))]), paste, collapse=";")
humanPanther <- sapply(shomo, function(i) unique(unlist(strsplit(i$PANTHER_MF_ALL, split=","))))
humanPanther <- sapply(sapply(humanPanther, function(i) i[which(!is.na(i))]), paste, collapse=";")
pantherLabCat <- sapply(shomo, function(i) unique(unlist(strsplit(as.character(i$davidAllAdd), split=","))))
pantherLabCat <- sapply(sapply(pantherLabCat, function(i) i[which(!is.na(i))]), paste, collapse=";")
labCat <- sapply(shomo, function(i) unique(unlist(strsplit(as.character(i$category), split=","))))
labCat <- sapply(sapply(labCat, function(i) i[which(!is.na(i))]), paste, collapse=";")

####
## Combine vectors into annotationTable
####

annotationTable <- cbind(homologNumber, humanSymbol, labCat, pantherLabCat, humanPanther, humanSummary)

write.table(annotationTable, file="annotationTable_v2.csv", quote=F, row.names=F, sep="\t")
