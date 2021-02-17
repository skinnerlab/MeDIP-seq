## This is the code to perform expanded overlaps. 



## This function calculates the minimum p-value in the resultTable that overlaps with the DMRs in the dmrTable
dmrMinp <- function(dmrTable, resultTable){
  minp <- foreach (i=1:nrow(dmrTable)) %dopar% {
    dmr <- dmrTable[i,]
    ## TODO
    ## This section can be greatly simplified with a single !(resultTable$stop < dmr$start) & !(resultTable$start > dmr$stop). 
    ## This should speed everything up. 
    ta <- which((resultTable$chr==dmr$chr) & (resultTable$start >= dmr$start) & (resultTable$stop <= dmr$stop))
    tb <- which((resultTable$chr==dmr$chr) & (resultTable$start >= dmr$start) & (resultTable$start < dmr$stop))
    tc <- which((resultTable$chr==dmr$chr) & (resultTable$top > dmr$start) & (resultTable$stop <= dmr$stop))
    td <- which((resultTable$chr==dmr$chr) & (resultTable$start <= dmr$start) & (resultTable$stop >= dmr$stop))
    te <- which((resultTable$chr==dmr$chr) & (resultTable$start <= dmr$start) & (resultTable$stop > dmr$start))
    tf <- which((resultTable$chr==dmr$chr) & (resultTable$start < dmr$stop) & (resultTable$stop >= dmr$stop))
    
    dmr.win <- unique(c(ta, tb, tc, td, te, tf))
    
    pvs <- resultTable$edgeR.p.value[dmr.win]
    pvs <- pvs[!is.na(pvs)]
    if (length(pvs) >= 1) {
      a <- min(pvs)
    } else {
      a <- NA
    }
    a
  }
  return(unlist(minp))
}

## Calculates expanded overlap table. This function still needs to be tested and optimized.
expandedOverlap <- function(resultTable.list, dmrTable.list, ncore=1, minp=0.05){
  
  library(doMC)
  registerDoMC(ncore)
  
  list.names <- c("ldp.testes", "ldp.kidney", "ldp.multiple")
  
  lo <- matrix(nrow=length(dmrTable.list), ncol=length(resultTable.list), numeric(1))
  row.names(lo) <- list.names
  colnames(lo) <- list.names
  
  for (i in 1:length(dmrTable.list)){
    for (j in 1:length(resultTable.list)){
      print(lo) ## Print progress flag. This should be removed after optimization
      lo[i,j] <- length(which(dmrMinp(dmrTable=dmrTable.list[[i]], resultTable=resultTable.list[[j]])<=minp))
    }
  }
  lo.pct <- format(apply(lo, 2, function(i) i/sapply(dmrTable.list,nrow))*100, digits=2)
  lo.pct <- apply(lo.pct, 2, function(i) paste("(", gsub(pattern=" ", replacement="", i), "%)", sep=""))
  lo.all <- lo
  for (col in 1:ncol(lo)){
    lo.all[,col] <- paste(lo[,col], " ", lo.pct[,col], sep="")
  }

return(lo.all)

}