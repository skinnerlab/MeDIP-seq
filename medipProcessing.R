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
     
     for (pV in 1:length(pValues)){
          temp<-extractDMR(x=methResults, distance=adjDist, pValueIdentify=pValues[pV], pValueExtend=dmrBoundPvalue, mtc=FALSE)
          mtctemp<-extractDMR(x=methResults, distance=adjDist, pValueIdentify=pValues[pV], pValueExtend=dmrBoundPvalue, mtc=TRUE)
          temp[[1]]$minP<-sapply(strsplit(temp[[2]][,"edgeR.p.value"], split=";"), function(i) min(as.numeric(i)))
          mtctemp[[1]]$minP<-sapply(strsplit(mtctemp[[2]][,"edgeR.p.value"], split=";"), function(i) min(as.numeric(i)))
          
          methList[[pV]]<-temp[[1]]
          methListEtc[[pV]]<-temp[[2]]
          MTCmethList[[pV]]<-mtctemp[[1]]
          MTCmethListEtc[[pV]]<-mtctemp[[2]]
          rm(temp); rm(mtctemp)
          
          # There appears to be a bug in the MEDIPS package that defines genomic window stop coordinates by the window size. This is a problem when the window runs off the end of the chromosome. If this last window is determined to be differentially methylated, the CpG calculation on the window errors out. I fix this problem here by setting the last window stop site to the end of the chromosome. I'm not sure if this is an ideal solution. Hopefully the MEDIPS package isn't doing anything else weird at the end of chromosomes.
          methList[[pV]]<-modifyStop(dmrList=methList[[pV]], refGenome=eval(parse(text=referenceName)), maxDMR=maxDMRnum)
          MTCmethList[[pV]]<-modifyStop(dmrList=MTCmethList[[pV]], refGenome=eval(parse(text=referenceName)), maxDMR=maxDMRnum)
          
          # calculate CpG density
          methList[[pV]]<-calcCpGdensity(methList[[pV]], maxDMRnum)
          MTCmethList[[pV]]<-calcCpGdensity(MTCmethList[[pV]], maxDMRnum)
          
     }
  
     
     
     ####################
     ## Add Annotation ##
     ####################
     
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
          methList<-lapply(methList, function(i){
               addAnnotationBiomart(dmrList=i, annotationObject=annotationObject, maxDMR=maxDMRnum, chrPrefix=chrPrefix)
          })
          MTCmethList<-lapply(MTCmethList, function(i){
               addAnnotationBiomart(dmrList=i, annotationObject=annotationObject, maxDMR=maxDMRnum, chrPrefix=chrPrefix)
          })
     }
     
     ##################################
     ## Multiple significant windows ##
     ##################################
     methList2p<-lapply(methList, function(i){
          i<-i[which(i$numSigWin>=2),]
     })
     MTCmethList2p<-lapply(MTCmethList, function(i){
          i<-i[which(i$numSigWin>=2),]
     })
     methList3p<-lapply(methList, function(i){
          i<-i[which(i$numSigWin>=3),]
     })
     MTCmethList3p<-lapply(MTCmethList, function(i){
          i<-i[which(i$numSigWin>=3),]
     })
     # create table of DMR numbers and save to file in results. This will be used to select the p-value to use for further analysis.
     dmrNumberTable<-cbind(
          "p-value"=pValues,
          "allWindow"=sapply(methList, nrow),
          "twoWindow"=sapply(methList2p, nrow),
          "threeWindow"=sapply(methList3p,nrow),
          "MTCallWindow"=sapply(MTCmethList, nrow),
          "MTCtwoWindow"=sapply(MTCmethList2p, nrow),
          "MTCthreeWindow"=sapply(MTCmethList3p, nrow)
     )
     write.csv(file=paste(resultsDirectory, comparisonNames[analysis], "/dmrNumber.csv", sep=""), x=dmrNumberTable, quote=F, row.names=FALSE) 
     save(methList, methList2p, methList3p, MTCmethList, MTCmethList2p, MTCmethList3p, dmrNumberTable, file=paste(resultsDirectory, comparisonNames[analysis], "/methLists.RData", sep=""))
     save(methListEtc, MTCmethListEtc, file=paste(resultsDirectory, comparisonNames[analysis], "/methListEtc.RData", sep=""))
}


