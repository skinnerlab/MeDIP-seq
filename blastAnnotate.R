## Daniel Beck
## Created 2/10/2016
## Modified 

## The purpose of this script is to annotate DMR in a denovo manner.
## I want to take the DMR sequence all the way to a gene name and category.



# temporary loading of results for test
library(MEDIPS)
source("dataNames.R")
source("customFunctions.R")
library(bsgenomePackageName, character.only = T, lib.loc = genomeDirectory)
library(biomaRt)
library(ShortRead)
library(rentrez)

load("/projects/human/results/all/methLists.RData")
dmrTable <- methList2p[[2]]

# After some exploration, I've decided to use NCBI's blast+ executables to query the database.
# The results can then be imported into R. The following steps will be necessary.
#		1. Write DMR sequences to fasta file
#		2. Blast sequences
#		3. Read results into R
#		4. Filter results
#		5. Associate with genes
#		6. Find gene category
#		7. Generate necessary tables

## Write DMR sequences to fasta file

dmrSequences <- dmr2seq(dmrTable=dmrTable, bsgenomeName=bsgenomePackageName, extendEdges=10000)

writeFasta(dmrSequences, file = "dmrTable.fa")

# for the human, I ran this manually on the command line (using nohup) as it has a long runtime and several errors

system('blastn -query dmrTable.fa -out dmrTable.br -db /projects/blast/nt/nt -outfmt 6 -best_hit_overhang 0.25 -max_target_seqs 20 -num_threads 20')
br <- read.table("dmrTable2.br", sep="\t", stringsAsFactor=F)
colnames(br) <- c("query", "subject", "pctId", "alignLength", "mismatches", "gapOpenings", "q.start", "q.end", "s.start", "s.end", "eValue","bitScore")

#br <- split.data.frame(br, f=br$query)

##### Filter br here
# filteredBr <- lapply(splitBr, function(i) i[which(i$eValue==min(i$eValue)),])

# extract id number from subject
ids <- sapply(strsplit(br$subject, split="\\|"), function(i) i[2])
uids <- unique(ids)  #2078
blastSummary <- sapply(uids, function(i) {
     ind <- which(uids==i); if (ind/100 == floor(ind/100)) print(ind)
     i<-entrez_summary(db="nuccore", id=i, always_return_list=TRUE)
     i
     })

blastTitle <- sapply(blastSummary, function(i) i$title)
blastOrganism <- sapply(blastSummary, function(i) i$organism)
blastCaption <- sapply(blastSummary, function(i) i$caption)
blastExtra <- sapply(blastSummary, function(i) i$extra)
blastSubname <- sapply(blastSummary, function(i) i$subname)


geneLinks <- lapply(uids, function(i) {
     ind <- which(uids==i); if (ind/100 == floor(ind/100)) print(ind)
     i<-entrez_link(dbfrom="nuccore",db="gene", id=i)
     i
     })

linkedGenes <- sapply(geneLinks, function(i) i$links[1])
linkedGenes <- sapply(linkedGenes, function(i) {if(is.null(i)) i<-NA; i})
# ## I've changed this to a loop due to a server error. Hopefully this will
# ## prevent a single error from killing the entire query.
linkedSummary <- list()
for (i in 1:length(linkedGenes)) {
     if (!is.na(linkedGenes[[i]])) {
          linkedSummary[[i]] <- try(entrez_summary(db="gene", id=linkedGenes[[i]]))
     } else {
          linkedSummary[[i]] <- NA
     }
     if (i/100 == floor(i/100)) print(i)
     }
# 
# save.image(file="tempSave.RData")
# 
# # extract uid and name (1st and 2nd entry)
# 
homoloGene <- read.table("/projects/medipPipeline/data/homologene_v1.data", sep="\t", stringsAsFactors=F)
annotationTable <- read.table("/projects/medipPipeline/data/annotationTable_v1.txt", stringsAsFactors=F, quote="", sep="\t", header=T)
# 
# ## Replace errors
linkedSummary[c(which(sapply(linkedSummary, typeof)=="character"))] <- NA
# 
linkedGS <- list()
for (i in 1:length(linkedSummary)) {
     if (typeof(linkedSummary[[i]][[1]])=="list") {
          lst <- list()
          for (j in 1:length(linkedSummary[[i]])){
               lst[[j]] <- linkedSummary[[i]][[j]]$name
          }
          linkedGS[[i]] <- unlist(lst)
     } else if (typeof(linkedSummary[[i]][[1]])=="character") {
          linkedGS[[i]] <- linkedSummary[[i]]$name
     } else {
          linkedGS[[i]] <- NA
     }
}

# linked gene titles
linkedGT <- list()
for (i in 1:length(linkedSummary)) {
     if (typeof(linkedSummary[[i]][[1]])=="list") {
          lst <- list()
          for (j in 1:length(linkedSummary[[i]])){
               lst[[j]] <- linkedSummary[[i]][[j]]$description
          }
          linkedGT[[i]] <- unlist(lst)
     } else if (typeof(linkedSummary[[i]][[1]])=="character") {
          linkedGT[[i]] <- linkedSummary[[i]]$description
     } else {
          linkedGT[[i]] <- NA
     }
}

# extract annotation information for linked genes
annolinkedGS <- list()
for (i in 1:length(linkedGS)) {
     if (!is.na(linkedGS[[i]][1])) {
          a <- lapply(linkedGS[[i]], function(j) customAnnotation(query=j, homologs=homoloGene, annotationTable=annotationTable))
          annolinkedGS[[i]] <- do.call(rbind, a)
     } else {
          annolinkedGS[[i]] <- annotationTable[0,]
     }
     if (i/1000 == floor(i/1000)) print(i/length(linkedGS))
}
annolinkedGS <- lapply(annolinkedGS, function(i) { if(nrow(i)<1) i[1,]<-NA; i})

# number of linked genes for every uid
nl <- sapply(annolinkedGS, nrow)
#nl <- sapply(nl, function(i) {if(i==0) i<-1; i})


# The next step is to combine lab annotation and homolog information with linked gene title/description
cBlast <- list()

for (i in 1:length(annolinkedGS)) {
     cBlast[[i]] <- data.frame("uid"=rep(uids[i], nl[i]), "blastTitle"=rep(blastTitle[i], nl[i]), "blastExtra"=rep(blastExtra[i], nl[i]), "linkedGeneSymbol"=linkedGS[[i]], "linkedGeneTitle"=linkedGT[[i]], annolinkedGS[[i]], stringsAsFactors=F)
}
cBlast <- lapply(cBlast, function(i) { if (ncol(i)>7) i<-i[,-c(which(colnames(i)=="humanSummary"))]; i })

# Combind annotations into data.frame
finalAnn <- do.call(rbind, cBlast)


## finalAnn can now be used as a database to associate with dmrID using the uid
brAnn <- finalAnn[match(ids, finalAnn$uid),]

# original human annotations
directAnnotation <- dmrTable$annotation[match(br$query, dmrTable$ID)]



dmrAnn <- data.frame("dmrID"=br$query, directAnnotation, "eValue"=br$eValue, "bitScore"=br$bitScore, brAnn, stringsAsFactors=F)
dmrAnn <- dmrAnn[which(!is.na(dmrAnn$linkedGeneSymbol)),]
dmrAnnList <- split.data.frame(dmrAnn, f=dmrAnn$dmrID)

#fdmrAnnList[[i]][,c(1, 2, 4, 5, 8, 10, 11, 12, 13)]; i<-i+1

# filter out low eValue blast hits (This should be done with the blastn command if this threshold works for other datasets)
fdmrAnnList <- lapply(dmrAnnList, function(i) i[which(i$eValue < 1e-100),])
# keep only best hit for each target sequence
fdmrAnnList <- lapply(fdmrAnnList, function(i) {
     split <- split.data.frame(i, f=i$uid)
     split <- lapply(split, function (j) {
          j <- j[which(j$bitScore == max(j$bitScore))[1],]
          j
     })
     do.call(rbind, split)
})


allCombined <- do.call(rbind, fdmrAnnList)
allCombined$linkedGeneTitle <- gsub(pattern=",", replacement=";", x=allCombined$linkedGeneTitle)
write.csv(file="blastAnn_1e-04_2p.csv", allCombined[,c(1, 2, 4, 5, 8, 9, 10, 11, 12, 13)], quote=F, row.names=F)





# Filter results to include only top homolog and top non-homolog
load("tempSave.RData")
allCombined <- read.csv("blastAnn_1e-04_2p.csv")
sac <- split.data.frame(allCombined, f=allCombined$dmrID)

fsac <- lapply(sac, function(i) {
     nhp <- which(is.na(i$homologNumber))
     hp <- which(!is.na(i$homologNumber))
     
     nhResult <- i[nhp,]
     hResult <- i[hp,]
     
     nhResult <- nhResult[which(nhResult$bitScore == max(nhResult$bitScore)),]
     hResult <- hResult[which(hResult$bitScore==max(hResult$bitScore)),]
     
     rbind(nhResult, hResult)
})


allCombinedFiltered <- do.call(rbind, fsac)
write.csv(file="blastAnnFiltered_1e-04_2p.csv", allCombinedFiltered, quote=F, row.names=F)


## Verify human results by quantifying overlaps.

## One problem is the different DMRs used in the different analyses (Extended in new blast results).


## Four aspects related to false positive/fase negative.
#    1. Genes that are found identically in both
#    2. Genes found in first but not second
#    3. Genes found in second but not first
#    4. Genes found differently in both


# First generate two tables, ensembl annotations in blast results, blast results in ensembl annotations
# Number 1 above








