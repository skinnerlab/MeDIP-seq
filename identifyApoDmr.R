## Daniel Beck
## Created 5/31/2016
## Modified 5/31/2016

## This code is intended to identify arbitrary pair overlap (APO) DMR. It reads 
## in the pairwise comparison results files and identifies APO DMR using standard
## overlapping functions. The end result is a methLists_APO Rdata file stored in 
## the results directory

# Load necessary packages
library(MEDIPS)
source("dataNames.R")
source("customFunctions.R")
library(bsgenomePackageName, character.only = T, lib.loc = genomeDirectory)
library(biomaRt)

###############################
## Read files and preprocess ##
###############################
# Loop over all analyses in apo.names (from dataNames.R)
for (analysis in 1:length(apo.names)) {
  ## Load relevant analyses
  apo.pair.analysis.env <- list()  # list of environments that will hold paired analyses
  for (i in 1:length(apo.pair.analysis.names[[analysis]])) {
    apo.pair.analysis.env[[i]] <- local({load(paste(resultsDirectory, 
                                                apo.pair.analysis.names[[analysis]][i], 
                                                "/methLists.RData", sep = "")); environment()})
  }

  methList.list <- lapply(apo.pair.analysis.env, function(i) i$methList)
  methList2p.list <- lapply(apo.pair.analysis.env, function(i) i$methList2p)

  MTCmethList.list <- lapply(apo.pair.analysis.env, function(i) i$MTCmethList)
  MTCmethList2p.list <- lapply(apo.pair.analysis.env, function(i) i$MTCmethList2p)

  ###################
  ## Find overlaps ##
  ###################
  ## This creates analogs to methList, methList2p, MTCmethList, and MTCmethList2p
  methListAPO <- list()
  methList2pAPO <- list()
  for (i in 1:length(pValues)) {
    tML <- lapply(methList.list, function(j) j[[i]])
    tML2p <- lapply(methList2p.list, function(j) j[[i]])
    
    ##### WARNING: arbitrary selection that affects the results!
    methListAPO[[i]] <- bpToDmr(bpList = overlappingBP(tML), dmrTable = tML[[1]])
    methList2pAPO[[i]] <- bpToDmr(bpList = overlappingBP(tML2p), dmrTable = tML2p[[1]])
  }

  MTCmethListAPO <- list()
  MTCmethList2pAPO <- list()
  for (i in 1:length(MTCpValues)) {
    tML <- lapply(MTCmethList.list, function(j) j[[i]])
    tML2p <- lapply(MTCmethList2p.list, function(j) j[[i]])
  
    ##### WARNING: arbitrary selection that affects the results!
    MTCmethListAPO[[i]] <- bpToDmr(bpList = overlappingBP(tML), dmrTable = tML[[1]])
    MTCmethList2pAPO[[i]] <- bpToDmr(bpList = overlappingBP(tML2p), dmrTable = tML2p[[1]])
  }

  ###############
  ## Count DMR ##
  ###############
  ## This creates the dmrNumberTable and MTCdmrNumberTable objects
  dmrNumberTableAPO <- cbind(
    "p-value" = pValues,
    "allWindow" = sapply(methListAPO, nrow),
    "twoWindow" = sapply(methList2pAPO, nrow)
  )
  MTCdmrNumberTableAPO <- cbind(
    "p-value" = MTCpValues,
    "allWindow" = sapply(MTCmethListAPO, nrow),
    "twoWindow" = sapply(MTCmethList2pAPO, nrow)
  )


  ##################
  ## Save results ##
  ##################
  dest.dir <- paste(resultsDirectory, apo.names[analysis], sep = "")
  system(paste("mkdir ", dest.dir, sep = ""))
  save(methListAPO, methList2pAPO, MTCmethListAPO, MTCmethList2pAPO, 
       dmrNumberTableAPO, MTCdmrNumberTableAPO, 
       file = paste(dest.dir, "/methListsAPO.RData", sep = ""))

}




