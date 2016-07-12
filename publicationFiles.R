## Created 9/15/2015 by Daniel Beck
## Last modified 4/6/2016

## This code is intended to produce all publication figures and tables. It repeats many 
## of the report figures. This file should be modified for each project so that includes 
## the exact figures produced for the final manuscript.


#####################
##### Load data #####
#####################

resultsFiles <- c(
  "/projects/human/results/all/"
)

dataNamesFiles <- c(
  "/projects/human/code/dataNames.R"
)

analysisNames <- c("human")
selectedPvalue <- c(1e-4)

dataList <- list()
resultsList <- list()
methList <- list()
dmrList <- list()
annList <- list()
annList2p <- list()
pvc <- numeric(length(analysisNames))
for (i in 1:length(resultsFiles)) {
     dataList[[i]] <- new.env()
     source(dataNamesFiles[[i]], local = dataList[[i]])
     resultsList[[i]] <- local({load(paste(resultsFiles[[i]], 
                                           "methLists.RData", 
                                           sep="")); environment()})
     pvc[i] <- which(dataList[[i]]$pValues == selectedPvalue[i])
     if (length(resultsList[[i]]$annMat) > 0) annList[[i]] <- resultsList[[i]]$annMat[[pvc[i]]]
     if (length(resultsList[[i]]$methList2p) > 0) {
          dmrList[[i]] <- resultsList[[i]]$methList2p[[pvc[i]]]
          dmrList[[i]] <- cbind(dmrName = paste("DMR",
                                                dmrList[[i]]$chr, ":",
                                                dmrList[[i]]$start, sep = ""), dmrList[[i]])
     }
     if (length(resultsList[[i]]$methList)>0) {
       methList[[i]] <- resultsList[[i]]$methList[[pvc[i]]]
     }
     annList2p[[i]] <- annList[[i]][as.logical(match(as.character(annList[[i]]$ID), 
                                                     as.character(dmrList[[i]]$ID), 
                                                     nomatch=0)), ]
     annList2p[[i]]$description <- gsub(pattern = ",", 
                                        replacement = " -", 
                                        annList2p[[i]]$description)
     namesList <- paste("DMR", dmrList[[1]]$chr, ":", dmrList[[1]]$start, sep = "")
     annList2p[[i]] <- cbind(dmrName = namesList[match(annList2p[[i]]$ID, 
                                                       dmrList[[i]]$ID)], annList2p[[i]])
     library(dataList[[i]]$bsgenomePackageName, character.only = T, 
             lib.loc=dataList[[i]]$genomeDirectory)
}

## Load necessary libraries
library(GenomicRanges)
library(VennDiagram)

## load current version of customFunctions.R
setwd("/projects/medipPipeline")
source("customFunctions.R")

setwd("/projects/human/results/publicationFiles")

#########################
##### Basic figures #####
#########################

# This section should include all standard figures produced for every project
## Chromosome plot
pdf(paste("chromosomePlot_1p_", analysisNames[1], "_", selectedPvalue[1], ".pdf", sep = ""))
smethList <- excludeChr(methList[[1]], exclude = "")
clusterObject <- slidingWindowCluster(formatForCluster(smethList), 
                                      windowLength = 2000000, 
                                      incrementLength = 50000)
colnames(clusterObject) <- c("chr", "start", "stop", "minP")
chrLengths <- seqlengths(eval(parse(text = dataList[[1]]$referenceName)))
plotChromosomes(siteTable = smethList, chrLengths = chrLengths, 
                ymar = 5, clusters = clusterObject)
dev.off()

pdf(paste("chromosomePlot_2p_", analysisNames[1], "_", selectedPvalue[1], ".pdf", sep = ""))
smethList <- excludeChr(dmrList[[1]], exclude = "")
clusterObject <- slidingWindowCluster(formatForCluster(smethList), 
                                      windowLength = 2000000, 
                                      incrementLength = 50000)
colnames(clusterObject) <- c("chr", "start", "stop", "minP")
chrLengths <- seqlengths(eval(parse(text = dataList[[1]]$referenceName)))
plotChromosomes(siteTable = smethList, chrLengths = chrLengths, 
                ymar = 5, clusters = clusterObject)
dev.off()

## CpG density histogram
pdf(paste("cpgDensity_2p_", analysisNames[1], "_", selectedPvalue[1], ".pdf", sep = ""))
par(mar = c(5, 5, 4, 2) + 0.1)
plotCpGdensity(dmrList[[1]], main = "", xlab = "Number of CpG sites per 100bp", 
               ylab = "Number of differential DNA \n methylation regions")
dev.off()


## DMR length histogram
pdf(paste("dmrLength_2p_", analysisNames[1], "_", selectedPvalue[1], ".pdf", sep = ""))
plotDMRlength(dmrList[[1]])
dev.off()

########################
##### Basic tables #####
########################

# Alignment percentage table
tmp <- getwd()
setwd("/projects/human/code")
seqFiles <- dataList[[1]]$seqFiles
comparison <- dataList[[1]]$comparison
an <- 1
alignPct <- mapExtract()
setwd(tmp)
write.csv(alignPct, file = paste("alignPct_", analysisNames[1], ".csv", sep = ""), quote = F)

# DMR number table
write.csv(resultsList[[1]]$dmrNumberTable, 
          file = paste("dmrNumberTable_", analysisNames[1], ".csv", sep = ""), 
          quote = F)

# DMR number of significant windows table
nWin <- table(methList[[1]]$numSigWin)
nWinTable <- as.data.frame(matrix(ncol = length(nWin), nrow = 2, c(names(nWin), nWin), byrow = T))
row.names(nWinTable) <- c("Number of significant windows", "Number of DMR")
colnames(nWinTable) <- rep(" ", ncol(nWinTable))
write.csv(nWinTable, file = paste("nWinTable_", analysisNames[[1]], "_", 
                                  selectedPvalue[1], ".csv", sep = ""), quote = F) 

## DMR table
write.csv(dmrList[[1]][,-4], file = paste("dmrTable_", analysisNames[1], "_", "2p", 
                                          "_", selectedPvalue[1], ".csv", sep = ""), 
          quote = F, row.names = F)

## Annotation table

write.csv(annList2p[[1]], file = paste("annTable_", analysisNames[1], "_", "2p", 
                                       "_", selectedPvalue[1], ".csv", sep = ""), 
          quote = F, row.names = F)

write.table(c(paste(unique(annList2p[[1]]$entrezgene), collapse = " "), 
              paste(unique(annList2p[[1]]$external_gene_name), collapse = " ")), 
            row.names = F, col.names = F, quote = F, 
            file = paste("listForKegg_", analysisNames[1], "_", "2p", "_", 
                         selectedPvalue[1], ".txt", sep = ""))

goTable <- sort(apply(t(apply(table(as.character(annList2p[[1]]$name_1006), 
                                    as.character(annList2p[[1]]$ID)), 1, 
                              as.numeric)), 1, sum), decreasing = T)
goTable <- goTable[-which(names(goTable) == "")]
goTable <- cbind(names(goTable), goTable)
colnames(goTable) <- c("GO category", "DMR number")
write.csv(goTable, quote = F, row.names = F, 
          file = paste("goCatNum_", analysisNames[1], "_", "2p", "_", 
                       selectedPvalue[1], ".txt", sep = ""))

##############################
##### Comparison figures #####
##############################

## Venn diagram
#pdf("vennDiagram_RBCvsSperm.pdf", colormodel="cmyk")
#  vennDMR(dmr.list = dmrList, group.names = c("Sperm", "RBC"), 
#             ext.text = F, cat.pos = 0, cex = 2, cat.cex = 2)
#dev.off()

#################
##### Other #####
#################

# This section includes nonstandard figures and tables



