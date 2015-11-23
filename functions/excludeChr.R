# function to exclude chromosomes from DMR list
excludeChr<-function(dmrList, exclude){
     dmrList<-dmrList[match(dmrList$chr, exclude, nomatch=0)==0,]
     return(dmrList)
}

# function to count the number of DMR of a set of chromosomes
countDMR<-function(dmrList, chr){
     numDMR<-sum(match(dmrList$chr, chr, nomatch=0)>0)
     return(numDMR)
}