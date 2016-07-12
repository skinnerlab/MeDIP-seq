
## Created 6/29/2016 by Daniel Beck
## Last modified 6/29/2016

## The goal of this script is to develop a simple PCA plot generator from the raw
## genomic window count data. This will be used for the incorporation of a normalization
## step into the MeDIP-seq analysis.
library(vegan)

load("/bigDisk/epigenesys/chemo1/results/all/methResults.RData")
mr1 <- methResults[,4:9]
load("/bigDisk/epigenesys/chemo2/results/all/methResults.RData")
mr2 <- methResults[,4:9]
load("/bigDisk/epigenesys/chemo3/results/all/methResults.RData")
mr3 <- methResults[,4:9]
load("/bigDisk/epigenesys/chemo4/results/all/methResults.RData")
mr4 <- methResults[,4:9]

allMr <- cbind(mr1, mr2, mr3, mr4)
ntemp <- c(colnames(mr1), colnames(mr2), colnames(mr3), colnames(mr4))
ntemp <- sapply(strsplit(ntemp, split="\\."), function(i) i[1])
colnames(allMr) <- ntemp

cpsObj <- capscale(allMr~1)
cpsObj()$CA$u

pcaObj <- princomp(allMr)

plot(pcaObj$loadings[,1], pcaObj$loadings[,2], xlab = "PC1", ylab="PC2")
(pcaObj$sdev^2) / sum(pcaObj$sdev^2) 

setwd("/bigDisk/epigenesys/chemo1/code")
source("dataNames.R")

library(RUVSeq)
a<-plotPCA(as.matrix(methResults[,4:9]))

## 35.3, 20.19