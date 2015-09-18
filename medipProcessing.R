## Daniel Beck
## Created 8/5/2015
## Modified

## This code is intended to continue the MeDIP analysis started by medipAnalysis.R

library(MEDIPS)
source("dataNames.R")
source("customFunctions.R")
library(bsgenomePackageName, character.only=T, lib.loc=genomeDirectory)
library(biomaRt)


for (analysis in 1:length(comparison)){  
     
     load(paste(resultsDirectory, comparisonNames[analysis], "/methResults.RData", sep=""))

     methList<-list()
     methListEtc<-list()
     MTCmethList<-list()
     MTCmethListEtc<-list()
     
     temp<-extractDMR(x=methResults, distance=adjDist, pValueIdentify=pValues, pValueExtend=dmrBoundPvalue, mtcIdentify=FALSE, mtcExtend=FALSE)
     methList<-temp[[1]]
     methListEtc<-temp[[2]]
     rm(temp)
     for (pV in 1:length(methList)){     
          if (!is.null(methList[[pV]])){
               methList[[pV]]$minP<-sapply(strsplit(as.character(rbind(methListEtc[[pV]])[,"edgeR.p.value"]), split=";"), function(i) min(as.numeric(i)))
               # correct window bounds outside of chromosome
               methList[[pV]]<-modifyStop(dmrList=methList[[pV]], refGenome=eval(parse(text=referenceName)), maxDMR=maxDMRnum)
               # add CpG density
               methList[[pV]]<-calcCpGdensity(methList[[pV]], maxDMRnum)
          }
     }

     temp<-extractDMR(x=methResults, distance=adjDist, pValueIdentify=MTCpValues, pValueExtend=dmrBoundPvalue, mtcIdentify=TRUE, mtcExtend=FALSE)
     MTCmethList<-temp[[1]]
     MTCmethListEtc<-temp[[2]]
     rm(temp)
     for (pV in 1:length(MTCmethList)){
          if (!is.null(MTCmethList[[pV]])){
               MTCmethList[[pV]]$minP<-sapply(strsplit(as.character(rbind(MTCmethListEtc[[pV]])[,"edgeR.p.value"]), split=";"), function(i) min(as.numeric(i)))
               # correct window bounds outside of chromosome
               MTCmethList[[pV]]<-modifyStop(dmrList=MTCmethList[[pV]], refGenome=eval(parse(text=referenceName)), maxDMR=maxDMRnum)
               # add CpG density
               MTCmethList[[pV]]<-calcCpGdensity(MTCmethList[[pV]], maxDMRnum)
          }
     }
     
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
          annotationObject<-MEDIPS.getAnnotation(host=biomartHost, dataset=biomartDataset, annotation=c("GENE"), chr=chr.select)
          
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
     
     ##################################
     ## Multiple significant windows ##
     ##################################
     methList2p<-lapply(methList, function(i){
          if (!is.na(i)){
            i<-i[which(i$numSigWin>=2),]
          }
     })
     MTCmethList2p<-lapply(MTCmethList, function(i){
       if (!is.na(i)){
          i<-i[which(i$numSigWin>=2),]
       }
     })

     # create table of DMR numbers and save to file in results. This will be used to select the p-value to use for further analysis.
     
     dmrNumberTable<-cbind(
          "p-value"=pValues,
          "allWindow"=sapply(methList, nrow),
          "twoWindow"=sapply(methList2p, nrow)
     )
     MTCdmrNumberTable<-cbind(
          "p-value"=MTCpValues,
          "allWindow"=sapply(MTCmethList, nrow),
          "twoWindow"=sapply(MTCmethList2p, nrow)
     )
     write.csv(file=paste(resultsDirectory, comparisonNames[analysis], "/dmrNumber.csv", sep=""), x=dmrNumberTable, quote=F, row.names=FALSE) 
     write.csv(file=paste(resultsDirectory, comparisonNames[analysis], "/MTCdmrNumber.csv", sep=""), x=MTCdmrNumberTable, quote=F, row.names=FALSE) 
     
     save(methList, methList2p, annMat, MTCmethList, MTCmethList2p, MTCannMat, dmrNumberTable, MTCdmrNumberTable, file=paste(resultsDirectory, comparisonNames[analysis], "/methLists.RData", sep=""))
     save(methListEtc, MTCmethListEtc, file=paste(resultsDirectory, comparisonNames[analysis], "/methListEtc.RData", sep=""))
}


