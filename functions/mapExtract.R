# function to extract mapping information. Used in report generation files
mapExtract<-function(){
     prepareData<-readLines("prepareData.Rout")
     sn<-c(seqFiles$sampleName[comparison[[an]]$mset1],seqFiles$sampleName[comparison[[an]]$mset2])
     pct<-sapply(strsplit(prepareData[sapply(sn, function(i) grep(pattern=paste("Running Bowtie2 on ", i, '\"', sep=""), x=prepareData))+15], split=" "), function(i) i[1])
     rn<-sapply(strsplit(prepareData[sapply(sn, function(i) grep(pattern=paste("Running Bowtie2 on ", i, '\"', sep=""), x=prepareData))+1], split=" "), function(i) i[1])
     outTable<-as.data.frame(rbind(rn, pct))
     colnames(outTable)<-sn
     rownames(outTable)<-c("read number", "overall alignment rate")
     return(outTable)
}