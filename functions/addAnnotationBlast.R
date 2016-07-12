## Created 2/10/2016 by Daniel Beck
## Last modified 5/5/2016

## This function adds annotation information to a DMR table. The annotation information is
## estimated using BLAST to associate DMR sequences with known sequences in the NCBI nr database.
## This script relies on the NCBI's blast+ executables.

# temporary loading of results for test. These should be replaced by passing to function


#dmrTable <- methList[[1]]

# a<- addAnnotationBlast(dmrTable=dmrTable, bsgenomePackageName = bsgenomePackageName, extendEdges = 0, outFile = "test", outDir = "/projects/snail/")


addAnnotationBlast <- function(dmrTable, maxDMR = 1000, bsgenomePackageName, 
                               extendEdges = 0, outFile, outDir, max_eValue = 1) {
  
  library(bsgenomePackageName, character.only = T, lib.loc = genomeDirectory)
  library(biomaRt)
  library(ShortRead)
  library(rentrez)
  
  if (is.null(nrow(dmrTable))) {
    return(list(dmrTable = NA, annMat = NA))
  }
  if (nrow(dmrTable) < maxDMR && nrow(dmrTable) > 0) {
    # Convert DMR range to reference sequence DNAStringSet
    dmrSequences <- dmr2seq(dmrTable = dmrTable, bsgenomeName = bsgenomePackageName, 
                            extendEdges = extendEdges)
    ## Write DMR sequences to fasta file
    writeFasta(dmrSequences, file = paste(outDir, outFile, ".fasta", sep=""))
    ## Run blastn against local database
    system(paste('blastn -query ', outDir, outFile, '.fasta -out ', outDir, outFile,
                 '.br -db /projects/blast/nt/nt -outfmt 6 -best_hit_overhang 0.25',
                 ' -max_target_seqs 20 -num_threads 10', sep = ""))
    br <- read.table(paste(outDir, outFile, ".br", sep = ""), sep = "\t", stringsAsFactor = F)
    
    colnames(br) <- c("query", "subject", "pctId", "alignLength", 
                      "mismatches", "gapOpenings", "q.start", "q.end", 
                      "s.start", "s.end", "eValue", "bitScore")
    
    # extract id number from subject
    ids <- sapply(strsplit(br$subject, split = "\\|"), function(i) i[2])
    uids <- unique(ids)
    blastSummary <- sapply(uids, function(i) {
      ind <- which(uids == i)
      i <- entrez_summary(db = "nuccore", id = i, always_return_list = TRUE)
      i
    })
    
    blastTitle <- sapply(blastSummary, function(i) i$title)
    blastOrganism <- sapply(blastSummary, function(i) i$organism)
    blastCaption <- sapply(blastSummary, function(i) i$caption)
    blastExtra <- sapply(blastSummary, function(i) i$extra)
    blastSubname <- sapply(blastSummary, function(i) i$subname)
    
    geneLinks <- lapply(uids, function(i) {
      ind <- which(uids == i)
      i <- entrez_link(dbfrom = "nuccore", db = "gene", id = i)
      i
    })
    linkedGenes <- sapply(geneLinks, function(i) i$links[1])
    linkedGenes <- sapply(linkedGenes, function(i) { if(is.null(i)) i <- NA; i })
    # ## I've changed this to a loop due to a server error. Hopefully this will
    # ## prevent a single error from killing the entire query.
    linkedSummary <- list()
    for (i in 1:length(linkedGenes)) {
      if (!is.na(linkedGenes[[i]])) {
        linkedSummary[[i]] <- try(entrez_summary(db = "gene", id = linkedGenes[[i]]))
      } else {
        linkedSummary[[i]] <- NA
      }
    }
    homoloGene <- read.table("/projects/medipPipeline/data/homologene_v1.data", 
                             sep = "\t", stringsAsFactors = F)
    annotationTable <- read.table("/projects/medipPipeline/data/annotationTable_v1.txt", 
                                  stringsAsFactors = F, quote = "", sep = "\t", header = T)
    
    # ## Replace errors
    linkedSummary[c(which(sapply(linkedSummary, typeof) == "character"))] <- NA
    linkedGS <- list()
    for (i in 1:length(linkedSummary)) {
      if (typeof(linkedSummary[[i]][[1]]) == "list") {
        lst <- list()
        for (j in 1:length(linkedSummary[[i]])){
          lst[[j]] <- linkedSummary[[i]][[j]]$name
        }
        linkedGS[[i]] <- unlist(lst)
      } else if (typeof(linkedSummary[[i]][[1]]) == "character") {
        linkedGS[[i]] <- linkedSummary[[i]]$name
      } else {
        linkedGS[[i]] <- NA
      }
    }
    # linked gene titles
    linkedGT <- list()
    for (i in 1:length(linkedSummary)) {
      if (typeof(linkedSummary[[i]][[1]]) == "list") {
        lst <- list()
        for (j in 1:length(linkedSummary[[i]])){
          lst[[j]] <- linkedSummary[[i]][[j]]$description
        }
        linkedGT[[i]] <- unlist(lst)
      } else if (typeof(linkedSummary[[i]][[1]]) == "character") {
        linkedGT[[i]] <- linkedSummary[[i]]$description
      } else {
        linkedGT[[i]] <- NA
      }
    }
    
    # extract annotation information for linked genes
    annolinkedGS <- list()
    for (i in 1:length(linkedGS)) {
      if (!is.na(linkedGS[[i]][1])) {
        a <- lapply(linkedGS[[i]], function(j) {
          customAnnotation(query = j, homologs = homoloGene, 
          annotationTable = annotationTable)})
        annolinkedGS[[i]] <- do.call(rbind, a)
      } else {
        annolinkedGS[[i]] <- annotationTable[0, ]
      }
    }
    annolinkedGS <- lapply(annolinkedGS, function(i) { if(nrow(i) < 1) i[1,] <- NA; i})
    
    # number of linked genes for every uid
    nl <- sapply(annolinkedGS, nrow)

    # The next step is to combine lab annotation and homolog information with linked gene title/description
    cBlast <- list()
    
    for (i in 1:length(annolinkedGS)) {
      cBlast[[i]] <- data.frame("uid" = rep(uids[i], nl[i]), 
                                "blastTitle" = rep(blastTitle[i], nl[i]),
                                "blastExtra" = rep(blastExtra[i], nl[i]), 
                                "linkedGeneSymbol" = linkedGS[[i]], 
                                "linkedGeneTitle" = linkedGT[[i]], 
                                annolinkedGS[[i]], stringsAsFactors = F)
    }
    cBlast <- lapply(cBlast, function(i) { 
      if (ncol(i) > 7) i <- i[, -c(which(colnames(i) == "humanSummary"))]
      i})
    
    # Combind annotations into data.frame
    finalAnn <- do.call(rbind, cBlast)
    
    
    ## finalAnn can now be used as a database to associate with dmrID using the uid
    brAnn <- finalAnn[match(ids, finalAnn$uid), ]
    
    dmrAnn <- data.frame("dmrID" = br$query, "eValue" = br$eValue, "bitScore" = br$bitScore, 
                         brAnn, stringsAsFactors = FALSE)
    dmrAnn <- dmrAnn[which(!is.na(dmrAnn$linkedGeneSymbol)), ]
    dmrAnnList <- split.data.frame(dmrAnn, f = dmrAnn$dmrID)
    # Filter out low eValue blast hits (This should be done with the blastn command if this 
    # threshold works for other datasets)
    fdmrAnnList <- lapply(dmrAnnList, function(i) i[which(i$eValue < max_eValue),])
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
    allCombined$linkedGeneTitle <- gsub(pattern = ",", replacement = ";", 
                                        x = allCombined$linkedGeneTitle)
    allCombined$blastTitle <- gsub(pattern = ",", replacement = ";", 
                                   x = allCombined$blastTitle)
    
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
    
    # I need to produce two data frames, annMat and dmrTable. 
    #     annMat: dmrID, gene information needed for categorization/presentation
    #     dmrTable: dmrID, DMR information, gene symbol(s)
    
    annIds <- sapply(fsac, function(i) i$dmrID[1])
    annGene <- sapply(fsac, function(i) paste(unique(i$linkedGeneSymbol), collapse=";"))
    
    dmrTable <- data.frame(dmrTable, "geneAssociation" = annGene[match(dmrTable$ID, annIds)], 
                           stringsAsFactors = F)
    dmrName <- apply(dmrTable[match(annIds, dmrTable$ID), c("chr", "start")], 1, function(i) {
      paste("DMR", i[1], ":", i[2], collapse="")
      })
    dmrName <- gsub(dmrName, pattern=" ", replacement="")
    
    annMat <- data.frame("dmrName"=dmrName[match(allCombinedFiltered$dmrID, annIds)], 
                         allCombinedFiltered, stringsAsFactors=F)
    annMat <- annMat[,-c(which(colnames(annMat)=="eValue"))]
    annMat <- annMat[,-c(which(colnames(annMat)=="uid"))]
    annMat <- annMat[,-c(which(colnames(annMat)=="dmrID"))]
    
    dmrName2 <- paste("DMR", dmrTable$chr, ":", dmrTable$start, paste="")
    dmrName2 <- gsub(dmrName2, pattern=" ", replacement="")
    dmrTable <- data.frame("dmrName"=dmrName2, dmrTable, stringsAsFactors=F)
  }
  return(list(dmrTable = dmrTable, annMat = annMat))
}


