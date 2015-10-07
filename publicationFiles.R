## Daniel Beck
## Created 9/15/2015
## Modified

## This code is intended to produce all publication figures and tables. It repeats many of the report figures. This file should be modified for each project so that includes the exact figures produced for the final manuscript.


#####################
##### Load data #####
#####################

resultsFiles<-c(
     "/projects/finch/results/allForRBC_CDRSvsEG/",
     "/projects/finch/results/allForSperm_CDRSvsEG/",
     "/projects/finch/results/allFulRBC_CDRSvsEG/",
     "/projects/finch/results/allFulSperm_CDRSvsEG/"
)

dataNamesFiles<-c(
     "/projects/finch/code/dataNames.R",
     "/projects/finch/code/dataNames.R",
     "/projects/finch/code/dataNames.R",
     "/projects/finch/code/dataNames.R"
)

analysisNames<-c("fortis RBC", "fortis sperm", "fuliginosa RBC", "fuliginosa sperm")
selectedPvalue<-c(1e-3, 1e-3, 1e-3, 1e-3)

dataList<-list()
resultsList<-list()
dmrList<-list()
annList<-list()
pvc<-numeric(length(analysisNames))
for (i in 1:length(resultsFiles)){
     dataList[[i]]<-new.env()
     source(dataNamesFiles[[i]], local=dataList[[i]])
     resultsList[[i]]<-local({load(paste(resultsFiles[[i]], "methLists.RData", sep="")); environment()})
     pvc[i]<-which(dataList[[i]]$pValues==selectedPvalue[[i]])
     if (length(resultsList[[i]]$annMat)>0) annList[[i]]<-resultsList[[i]]$annMat[[pvc[i]]]
     if (length(resultsList[[i]]$methList2p)>0) dmrList[[i]]<-resultsList[[i]]$methList2p[[pvc[i]]]
     
     library(dataList[[i]]$bsgenomePackageName, character.only=T, lib.loc=dataList[[i]]$genomeDirectory)
}


## Load necessary libraries
library(GenomicRanges)
library(VennDiagram)

## load current version of customFunctions.R
setwd("/projects/medipPipeline")
source("customFunctions.R")

setwd("/projects/finch/results/publicationFiles")

#########################
##### Basic figures #####
#########################

# This section should include all standard figures produced for every project
## Chromosome plot
pdf(paste("chromosomePlot_", analysisNames[1], ".pdf", sep=""))
     methList<-excludeChr(dmrList[[1]], exclude="")
     clusterObject<-slidingWindowCluster(formatForCluster(methList), windowLength=2000000, incrementLength=50000)
     colnames(clusterObject)<-c("chr", "start", "stop", "minP")
     chrLengths<-seqlengths(eval(parse(text=dataList[[1]]$referenceName)))
     plotChromosomes(siteTable=methList, chrLengths=chrLengths, ymar=5, clusters=clusterObject)
dev.off()
pdf(paste("chromosomePlot_", analysisNames[2], ".pdf", sep=""))
     methList<-excludeChr(dmrList[[2]], exclude="")
     clusterObject<-slidingWindowCluster(formatForCluster(methList), windowLength=2000000, incrementLength=50000)
     colnames(clusterObject)<-c("chr", "start", "stop", "minP")
     chrLengths<-seqlengths(eval(parse(text=dataList[[2]]$referenceName)))
     plotChromosomes(siteTable=methList, chrLengths=chrLengths, ymar=5, clusters=clusterObject)
dev.off()
pdf(paste("chromosomePlot_", analysisNames[3], ".pdf", sep=""))
     methList<-excludeChr(dmrList[[3]], exclude="")
     clusterObject<-slidingWindowCluster(formatForCluster(methList), windowLength=2000000, incrementLength=50000)
     colnames(clusterObject)<-c("chr", "start", "stop", "minP")
     chrLengths<-seqlengths(eval(parse(text=dataList[[3]]$referenceName)))
     plotChromosomes(siteTable=methList, chrLengths=chrLengths, ymar=5, clusters=clusterObject)
dev.off()
pdf(paste("chromosomePlot_", analysisNames[4], ".pdf", sep=""))
     methList<-excludeChr(dmrList[[4]], exclude="")
     clusterObject<-slidingWindowCluster(formatForCluster(methList), windowLength=2000000, incrementLength=50000)
     colnames(clusterObject)<-c("chr", "start", "stop", "minP")
     chrLengths<-seqlengths(eval(parse(text=dataList[[4]]$referenceName)))
     plotChromosomes(siteTable=methList, chrLengths=chrLengths, ymar=5, clusters=clusterObject)
dev.off()

## CpG density histogram
pdf(paste("cpgDensity_", analysisNames[1], ".pdf", sep=""))
     par(mar=c(5, 5, 4, 2) + 0.1)
     plotCpGdensity(dmrList[[1]], main="", xlab="Number of CpG sites per 100bp", ylab="Number of differential DNA \n methylation regions")
dev.off()
pdf(paste("cpgDensity_", analysisNames[2], ".pdf", sep=""))
     par(mar=c(5, 5, 4, 2) + 0.1)
     plotCpGdensity(dmrList[[2]], main="", xlab="Number of CpG sites per 100bp", ylab="Number of differential DNA \n methylation regions")
dev.off()
pdf(paste("cpgDensity_", analysisNames[3], ".pdf", sep=""))
     par(mar=c(5, 5, 4, 2) + 0.1)
     plotCpGdensity(dmrList[[3]], main="", xlab="Number of CpG sites per 100bp", ylab="Number of differential DNA \n methylation regions")
dev.off()
pdf(paste("cpgDensity_", analysisNames[4], ".pdf", sep=""))
     par(mar=c(5, 5, 4, 2) + 0.1)
     plotCpGdensity(dmrList[[4]], main="", xlab="Number of CpG sites per 100bp", ylab="Number of differential DNA \n methylation regions")
dev.off()

## DMR length histogram
pdf(paste("dmrLength_", analysisNames[1], ".pdf", sep=""))
     plotDMRlength(dmrList[[1]])
dev.off()
pdf(paste("dmrLength_", analysisNames[2], ".pdf", sep=""))
     plotDMRlength(dmrList[[2]])
dev.off()
pdf(paste("dmrLength_", analysisNames[3], ".pdf", sep=""))
     plotDMRlength(dmrList[[3]])
dev.off()
pdf(paste("dmrLength_", analysisNames[4], ".pdf", sep=""))
     plotDMRlength(dmrList[[4]])
dev.off()
########################
##### Basic tables #####
########################

## DMR table
write.csv(dmrList[[1]][,-4], file=paste("dmrTable_", analysisNames[1], "_", "2p", "_", selectedPvalue[1], ".csv", sep=""), quote=F, row.names=F)

write.csv(dmrList[[2]][,-4], file=paste("dmrTable_", analysisNames[2], "_", "2p", "_", selectedPvalue[2], ".csv", sep=""), quote=F, row.names=F)

write.csv(dmrList[[3]][,-4], file=paste("dmrTable_", analysisNames[3], "_", "2p", "_", selectedPvalue[3], ".csv", sep=""), quote=F, row.names=F)

write.csv(dmrList[[4]][,-4], file=paste("dmrTable_", analysisNames[4], "_", "2p", "_", selectedPvalue[4], ".csv", sep=""), quote=F, row.names=F)
## Annotation tables

write.csv(annList[[1]], file=paste("annTable_", analysisNames[1], "_", "2p", "_", selectedPvalue[1], ".csv", sep=""), quote=F, row.names=F)

write.csv(annList[[2]], file=paste("annTable_", analysisNames[2], "_", "2p", "_", selectedPvalue[2], ".csv", sep=""), quote=F, row.names=F)

write.table(c(paste(unique(annList[[1]]$entrezgene), collapse=" "), paste(unique(annList[[1]]$external_gene_name), collapse=" ")), row.names=F, col.names=F, quote=F, file=paste("listForKegg_", analysisNames[2], "_", "2p", "_", selectedPvalue[2], ".txt"))
##############################
##### Comparison figures #####
##############################

## Venn diagram
pdf("vennDiagram_ForFulSpermRBC.pdf")
     vennDMRfour(dmrList=dmrList, names=analysisNames, scaled=TRUE, ext.text=F, cat.pos=0, cex=2, cat.cex=2)
dev.off()

#################
##### Other #####
#################

# This section includes nonstandard figures and tables



