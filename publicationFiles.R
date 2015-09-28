## Daniel Beck
## Created 9/15/2015
## Modified

## This code is intended to produce all publication figures and tables. It repeats many of the report figures. This file should be modified for each project so that includes the exact figures produced for the final manuscript.


#####################
##### Load data #####
#####################

resultsFiles<-c(
     "/projects/steelheadSperm/results/all/",
     "/projects/steelheadRBC/results/all/"
)

dataNamesFiles<-c(
     "/projects/steelheadSperm/code/dataNames.R",
     "/projects/steelheadRBC/code/dataNames.R"
)

analysisNames<-c("steelheadSperm", "steelheadRBC")
selectedPvalue<-c(1e-4, 1e-4)

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

setwd("/projects/steelheadSperm/results/publicationFiles")

#########################
##### Basic figures #####
#########################

# This section should include all standard figures produced for every project
## Chromosome plot
pdf(paste("chromosomePlot_", analysisNames[1], ".pdf", sep=""))
     methList<-excludeChr(dmrList[[1]], exclude="chrUn")
     clusterObject<-slidingWindowCluster(formatForCluster(methList), windowLength=2000000, incrementLength=50000)
     colnames(clusterObject)<-c("chr", "start", "stop", "minP")
     chrLengths<-seqlengths(eval(parse(text=dataList[[1]]$referenceName)))
     plotChromosomes(siteTable=methList, chrLengths=chrLengths, ymar=5, clusters=clusterObject)
dev.off()

pdf(paste("chromosomePlot_", analysisNames[2], ".pdf", sep=""))
     methList<-excludeChr(dmrList[[2]], exclude="chrUn")
     clusterObject<-slidingWindowCluster(formatForCluster(methList), windowLength=2000000, incrementLength=50000)
     colnames(clusterObject)<-c("chr", "start", "stop", "minP")
     chrLengths<-seqlengths(eval(parse(text=dataList[[2]]$referenceName)))
     plotChromosomes(siteTable=methList, chrLengths=chrLengths, ymar=5, clusters=clusterObject)
dev.off()

## CpG density histogram
par(mar=c(5, 5, 4, 2) + 0.1)
plotCpGdensity(methList2p[[pvc]], main="", xlab="Number of CpG sites per 100bp", ylab="Number of differential DNA \n methylation regions")

## DMR length histogram
plotDMRlength(methList2p[[pvc]])

########################
##### Basic tables #####
########################

## DMR table
methList2p[[]]
## Annotation table
annMat # modify to get 2p

##############################
##### Comparison figures #####
##############################


## Venn diagram
pdf("vennDiagram_RBCvsSperm.pdf")
     vennDMRtwo(dmrList=dmrList, names=c("Sperm", "RBC"), scaled=TRUE, ext.text=F, cat.pos=0)
dev.off()

#################
##### Other #####
#################

# This section includes nonstandard figures and tables



