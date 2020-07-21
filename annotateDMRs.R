## Daniel Beck
## Created 10/5/2015
## Modified 10/8/2015

## This code is intended to add annotation to DMRs identified by the medipProcessing script.

library(MEDIPS)
source("dataNames.R")
source("customFunctions.R")
library(bsgenomePackageName, character.only = T, lib.loc = genomeDirectory)
library(biomaRt)

# Moving this here so it is only called once. It failes occationally for unknown reasons. 
if (annotationType == "biomart") {
  annotationObject <- useMart(biomart = "ensembl",
                              dataset = biomartDataset, 
                              host = biomartHost)
}

for (analysis in 1:length(comparisonNames)) {
  load(paste(resultsDirectory, comparisonNames[analysis], "/methLists.RData", sep = ""))
  pvt <- ann.pvt[analysis]; pvi <- which(pValues==pvt)
  ####################
  ## Add Annotation ##
  ####################
  if (!is.na(pvt)){
    annMat <- list()
    # add Annotation information for each DMR. This is separate from the loop so that the 
    # import.gff and getAnnotation functions are only run once
    if (annotationType == "gff") {
      gff <- import.gff(paste(genomeDirectory, annotationGFF, sep = ""))
      methList[[pvi]] <- addAnnotationGFF(dmrList = methList[[pvi]], gff = gff, 
                                            maxDMR = maxDMRnum, 
                                            chrPrefix = chrPrefix)
    }
      
    if (annotationType == "biomart") {
      a <- addAnnotationBiomart(dmrList = methList[[pvi]], extension = 10000,
                                annotationObject = annotationObject, 
                                maxDMR = maxDMRnum, chrPrefix = chrPrefix)
      methList[[pvi]] <- a$dmrList
      annMat[[pvi]] <- a$annMat
    }
  
       
    ## Re-calculate twoWindow DMRs to include annotation information
         
    ##################################
    ## Multiple significant windows ##
    ##################################
    methList2p <- lapply(methList, function(i) {
      if(!is.null(i)) {
        if (!is.na(i)) {
          i <- i[which(i$numSigWin >= 2), ]
        } else {
          i <- NA
        }
      } else {
        i <- NULL
      }
    })
    MTCmethList2p <- lapply(MTCmethList, function(i) {
      if (!is.null(i)) {
        if (!is.na(i)) {
          i <- i[which(i$numSigWin >= 2),]
        } else {
          i <- NA
        }
      } else {
        i <- NULL
      }
    })
         
    save(methList, methList2p, MTCmethList, MTCmethList2p, 
         dmrNumberTable, MTCdmrNumberTable, annMat,
         file = paste(resultsDirectory, comparisonNames[analysis], "/methLists.RData", sep = ""))
  }
  print(paste("Progress: ", analysis / length(comparison), "%", sep=""))
   
}
