## tablify function taken from old phd code. It converts a list into a data.frame.
tablify<-function(lst){
  cn<-unique(names(unlist(lst)))
  tlst<-as.data.frame(matrix(ncol=length(cn), nrow=length(lst), numeric(length(cn)*length(lst))))
  colnames(tlst)<-cn
  for (i in 1:length(lst)){
    for (j in 1:length(cn)){
      index<-which(names(lst[[i]])==cn[j])
      if (length(index)) {tlst[i,j]<-lst[[i]][index]}
    }
  }
  return(tlst)
}
