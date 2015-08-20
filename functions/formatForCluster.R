


formatForCluster<-function(dmrList){
     colnames(dmrList)<-c("Chromosome", "cSTART", "cSTOP", "ID", "length", "nProbes", "minP")
     return(dmrList)
}
