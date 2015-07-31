## Daniel Beck
## Created 5/27/2015
## Modified
##	6/2/2015  Removed the CS object to prevent CpG normalization. This is probably unnecessary, but rerunning analysis to be sure. Also commented out BSgenome creation. This only needs to be done once.
##	6/22/2015 Added back the CS object to ensure its availability. Modified code to analyze finch dataset
##	6/29/2015 Split file processing (up to BAM file) and BSgenome package creation to separate script.
##   6/30/2015 Continued development of new pipeline structure. Added CpG density calculation. Fixed CpG length calculation (was previously 1bp short)
##   7/1/2015  Continued development
##   7/2/2015  Tested on steelheadTrial dataset
##   7/6/2015  Made CpG density calculation generic to all analyses
##   7/7/2015  Removed CpG density calculation to function in customFunctions.R. Added ability to limit calculation to smaller DMR numbers. Fixed several NA/NULL logical tests.
##   7/9/2015  Added annotation for gff and biomart based data. Added DMR number summary output as dmrNumber.csv in results folder. Moved pValues to dataNames.R for use in medipReport.Rmd.
##   7/10/2015 Added CS normalization option.
##   7/13/2015 Added file to git repository "medipPipeline"

## This code is intended to perform the MeDIP analysis. It is intended to work as
## part of a pipeline including customFunctions.R and prepareData.R. This script
## requires the presence of an R object file containing several object names.

# Load necessary packages
library(MEDIPS)
source("dataNames.R")
source("customFunctions.R")
library(bsgenomePackageName, character.only=T, lib.loc=genomeDirectory)
library(biomaRt)
library(rmarkdown)

###############################
## Read files and preprocess ##
###############################

processedBamFiles <-
  lapply(
    X = paste(dataDirectory, sbamFileName, sep = ""), FUN = MEDIPS.createSet, BSgenome =
      bsgenomePackageName, extend = extend, shift = shift, uniq = uniq, window_size =
      ws, chr.select = chr.select
  )

###################
## Identify DMRs ##
###################

# Normalization if necessary
if (CScalc){
        CS<-MEDIPS.couplingVector(pattern = "CG", refObj = processedBamFiles[[1]])
} else {
        CS<-NULL
}

# do all comparisons
resultList<-list()
for (analysis in 1:length(comparison)){  
  methResults = MEDIPS.meth(
  	MSet1=processedBamFiles[comparison[[analysis]]$mset1],
  	MSet2=processedBamFiles[comparison[[analysis]]$mset2],
  	p.adj=p.adj,
  	diff.method=diff.method,
  	MeDIP=MeDIP,
  	CNV=CNV,
  	CSet=CS,
  	minRowSum=minRowSum
  )

  #################
  ## Filter DMRs ##
  #################

  methResultsSubset <- lapply(pValues, function (i) methResults[which(methResults[,"edgeR.p.value"] < i),])
  MTCmethResultsSubset <- lapply(pValues, function (i) methResults[which(methResults[,"edgeR.adj.p.value"] < i),])

  ## Merge adjacent windows
  methList <- lapply(methResultsSubset, function (i) {
  	if (nrow(i)>0){
  		merged<-MEDIPS.mergeFrames(frames=i, distance=1)
  		merged$length<-as.integer(as.character(merged$stop)) - as.integer(as.character(merged$start)) + 1
  	} else {
  		merged<-NA
  	}
  	merged
  })
  MTCmethList <- lapply(MTCmethResultsSubset, function (i) {
  	if (nrow(i)>0){
  		merged<-MEDIPS.mergeFrames(frames=i, distance=1)
  		merged$length<-as.integer(as.character(merged$stop)) - as.integer(as.character(merged$start)) + 1
  	} else {
  		merged<-NA
  	}
  	merged
  })
  
  # create dataframes with identical rows as the merged dataframes that include all available information about the DMRs
  methListEtc<-lapply(1:length(methResultsSubset), function(i){
       if ((nrow(methResultsSubset[[i]])<maxDMRnum)&&(nrow(methResultsSubset[[i]])>0)){
            addToMergedResults(allWindows=methResultsSubset[[i]], mergedWindows=methList[[i]])
       }else{
            NA
       }
  })
  MTCmethListEtc<-lapply(1:length(MTCmethResultsSubset), function(i){
       if ((nrow(MTCmethResultsSubset[[i]])<maxDMRnum)&&(nrow(MTCmethResultsSubset[[i]])>0)){
            addToMergedResults(allWindows=MTCmethResultsSubset[[i]], mergedWindows=MTCmethList[[i]])
       }else{
            NA
       }
  })
  
  # extract regions with multiple merged windows
  methList3p <- lapply(methList, function (i) {
          if (!is.na(i)){
                  i[which(i$length>100),]
          } else {
                  NA
          }
  })
  MTCmethList3p <- lapply(MTCmethList, function (i) {
          if (!is.na(i)){
                  i[which(i$length>100),]
          } else {
                  NA
          }
  })

  # Filter additional information matrix in same manner
  methList3pEtc <- lapply(1:length(methList), function (i) {
          if (nrow(methList3p[[i]])>0 && !is.na(methList3p[[i]])){
               if (nrow(methList3p[[i]])>0 && !is.na(methList3p[[i]])){
                    if (!is.null(nrow(methListEtc[[i]]))){
                       methListEtc[[i]][which(methList[[i]]$length>100),]
                    } else { methListEtc[[i]] }
               } else { NA }
          } else { NA }
       })
  MTCmethList3pEtc<-lapply(1:length(MTCmethList), function (i) {
          if (nrow(MTCmethList3p[[i]])>0 && !is.na(MTCmethList3p[[i]])){
               if (!is.na(MTCmethList3p[[i]])&&!is.na(MTCmethListEtc[[i]])){
                    if (!is.null(nrow(MTCmethListEtc[[i]]))){
                       MTCmethListEtc[[i]][which(MTCmethList[[i]]$length>100),]
                         ## Below is temporary fix. Do this right later!
                    } else { MTCmethListEtc[[i]] }
               } else { NA }
           } else { NA }
       })

  #########################
  ## Chromosome edge fix ##
  #########################
  
  # There appears to be a bug in the MEDIPS package that defines genomic window stop coordinates by the window size. This is a problem when the window runs off the end of the chromosome. If this last window is determined to be differentially methylated, the CpG calculation on the window errors out. I fix this problem here by setting the last window stop site to the end of the chromosome. I'm not sure if this is an ideal solution. Hopefully the MEDIPS package isn't doing anything else weird at the end of chromosomes.

   methList<-lapply(methList, function(i){
        modifyStop(dmrList=i, refGenome=eval(parse(text=referenceName)), maxDMR=maxDMRnum)
   })
   MTCmethList<-lapply(MTCmethList, function(i){
        modifyStop(dmrList=i, refGenome=eval(parse(text=referenceName)), maxDMR=maxDMRnum)
   })
   
   methList3p<-lapply(methList3p, function(i){
        modifyStop(dmrList=i, refGenome=eval(parse(text=referenceName)), maxDMR=maxDMRnum)
   })
   
   MTCmethList3p<-lapply(MTCmethList3p, function(i){
        modifyStop(dmrList=i, refGenome=eval(parse(text=referenceName)), maxDMR=maxDMRnum)
   })
  
  #############################
  ## CpG density calculation ##
  #############################

  methList<-lapply(methList, function(i){
       calcCpGdensity(i, maxDMRnum)
       })
  
  MTCmethList<-lapply(MTCmethList, function(i){
       calcCpGdensity(i, maxDMRnum)
       })

  methList3p<-lapply(methList3p, function(i){
       calcCpGdensity(i, maxDMRnum)
       })
  
  MTCmethList3p<-lapply(MTCmethList3p, function(i){
       calcCpGdensity(i, maxDMRnum)
       })
  
  ####################
  ## Add Annotation ##
  ####################
  
  # add Annotation information for each DMR
  if (annotationType=="gff"){
       gff<-import.gff(paste(genomeDirectory, annotationGFF, sep=""))
       methList<-lapply(methList, function(i){
            addAnnotationGFF(dmrList=i, gff=gff, maxDMR=maxDMRnum, chrPrefix=chrPrefix)
       })
       MTCmethList<-lapply(MTCmethList, function(i){
            addAnnotationGFF(dmrList=i, gff=gff, maxDMR=maxDMRnum, chrPrefix=chrPrefix)
       })
       methList3p<-lapply(methList3p, function(i){
            addAnnotationGFF(dmrList=i, gff=gff, maxDMR=maxDMRnum, chrPrefix=chrPrefix)
       })
       MTCmethList3p<-lapply(MTCmethList3p, function(i){
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
       methList3p<-lapply(methList3p, function(i){
            addAnnotationBiomart(dmrList=i, annotationObject=annotationObject, maxDMR=maxDMRnum, chrPrefix=chrPrefix)
       })
       MTCmethList3p<-lapply(MTCmethList3p, function(i){
            addAnnotationBiomart(dmrList=i, annotationObject=annotationObject, maxDMR=maxDMRnum, chrPrefix=chrPrefix)
       })
  }
  
  system(paste("mkdir ", resultsDirectory, comparisonNames[analysis], sep=""))
  # create table of DMR numbers and save to file in results. This will be used to select the p-value to use for further analysis.
  dmrNumberTable<-cbind(
       "p-value"=pValues,
       "singleWindows"=sapply(methList, nrow),
       "multipleWindows"=sapply(methList3p,nrow),
       "MTCsingleWindows"=sapply(MTCmethList, nrow),
       "MTCmultipleWindows"=sapply(MTCmethList3p, nrow)
  )
  write.csv(file=paste(resultsDirectory, comparisonNames[analysis], "/dmrNumber.csv", sep=""), x=dmrNumberTable, quote=F, row.names=FALSE) 
  
  save.image(paste(resultsDirectory, comparisonNames[analysis], "/results.RData", sep=""))
}



