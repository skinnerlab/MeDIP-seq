## Daniel Beck
## Created 10/5/2015
## Modified

## This code is intended to add annotation to DMRs identified by the medipProcessing script.

library(MEDIPS)
source("dataNames.R")
source("customFunctions.R")
library(bsgenomePackageName, character.only=T, lib.loc=genomeDirectory)
library(biomaRt)


for (analysis in 1:length(comparison)){  
     load(paste(resultsDirectory, comparisonNames[analysis], "/methLists.RData", sep=""))
     
     ####################
     ## Add Annotation ##
     ####################
     annMat<-list()
     MTCannMat<-list()
     # add Annotation information for each DMR. This is separate from the loop so that the import.gff and getAnnotation functions are only run once
     if (annotationType=="gff"){
          gff<-import.gff(paste(genomeDirectory, annotationGFF, sep=""))
          methList<-lapply(methList, function(i){
               addAnnotationGFF(dmrList=i, gff=gff, maxDMR=maxDMRnum, chrPrefix=chrPrefix)
          })
          MTCmethList<-lapply(MTCmethList, function(i){
               addAnnotationGFF(dmrList=i, gff=gff, maxDMR=maxDMRnum, chrPrefix=chrPrefix)
          })
     }
     if (annotationType=="biomart"){
          annotationObject<-annMart<-useMart(biomart="ensembl", dataset=biomartDataset, host=biomartHost)
          
          for (i in 1:length(methList)){
               a<-addAnnotationBiomart(dmrList=methList[[i]], annotationObject=annotationObject, maxDMR=maxDMRnum, chrPrefix=chrPrefix)
               methList[[i]]<-a$dmrList
               annMat[[i]]<-a$annMat
          }
          for (i in 1:length(MTCmethList)){
               b<-addAnnotationBiomart(dmrList=MTCmethList[[i]], annotationObject=annotationObject, maxDMR=maxDMRnum, chrPrefix=chrPrefix)
               MTCmethList[[i]]<-b$dmrList
               MTCannMat[[i]]<-b$annMat
          }
     }
     
     ## Re-calculate twoWindow DMRs to include annotation information
     
     ##################################
     ## Multiple significant windows ##
     ##################################
     methList2p<-lapply(methList, function(i){
          if(!is.null(i)){
               if (!is.na(i)){
                    i<-i[which(i$numSigWin>=2),]
               } else {
                    i<-NA
               }
          } else {
               i<-NULL
          }
     })
     MTCmethList2p<-lapply(MTCmethList, function(i){
          if(!is.null(i)){
               if (!is.na(i)){
                    i<-i[which(i$numSigWin>=2),]
               } else {
                    i<-NA
               }
          } else {
               i<-NULL
          }
     })
     
     save(methList, methList2p, MTCmethList, MTCmethList2p, dmrNumberTable, MTCdmrNumberTable, annMat, MTCannMat, file=paste(resultsDirectory, comparisonNames[analysis], "/methLists.RData", sep=""))
     
}