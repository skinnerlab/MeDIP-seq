excludeChr<-function(dmrList, exclude){
     dmrList<-dmrList[match(dmrList$chr, exclude, nomatch=0)==0,]
     return(dmrList)
}