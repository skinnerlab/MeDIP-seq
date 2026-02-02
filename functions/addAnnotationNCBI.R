## Created -/-/2015 by Daniel Beck
## Last modified 4/6/2016

## This function adds annotation information to a DMR table.

extractRegion <- function(chr, start, stop, feature.table){
  features <- feature.table[which(feature.table$chr==chr & !(feature.table$stop < start) & !(feature.table$start > stop)),]
  return(features)
}

addAnnotationNCBI <- function(dmrList, chrPrefix = "", extension = 0, ann=NA, panther=panther, plc=NA) {
  if (is.null(nrow(dmrList))) {
    return(list(dmrList = NA, annMat = NA))
  }
  if (nrow(dmrList) > 0) {
    # extend DMR edges if necessary
    dmrList.ext <- dmrList
    dmrList.ext$start <- sapply(dmrList.ext$start-extension, function(i) max(i, 1))
    dmrList.ext$stop <- dmrList.ext$stop+extension
    
    annotation <- list()
    for (i in 1:nrow(dmrList.ext)) {
      
      annotation[[i]] <- extractRegion(chr=dmrList.ext[i,"chr"], start=dmrList.ext[i,"start"], 
                                    stop=dmrList.ext[i,"stop"], feature.table=ann)
    }
      
    atribute.list <- lapply(annotation, function(i) strsplit(i$attributes, split=";"))
    ## get gene id
    gene.idlist <- lapply(atribute.list, function(i) unique(unlist(lapply(i, function(j) gsub(j[grep(j, pattern="gene=")], pattern="gene=", replacement="")))))
    ## get NCBI GeneID (used for  KEGG pathway search)
    gene.ncbi.pre <- lapply(atribute.list, function(i) unique(unlist(lapply(i, function(j) gsub(j[grep(j, pattern="Dbxref=")], pattern="Dbxref=", replacement="")))))
    gene.ncbi.pre <- lapply(gene.ncbi.pre, function(i) i[grep(i, pattern="GeneID:")])
    gene.ncbi <- lapply(gene.ncbi.pre, function(i) {
      all <- unlist(strsplit(as.character(i), split=",")); 
      all <- all[grep(all, pattern="GeneID:")]; 
      unique(gsub(all, pattern="GeneID:", replacement=""))
    })
    ## add gene id to main table
    dmrList$annotation <- sapply(gene.idlist, function(i) paste(i, collapse=";"))
    ## get panther category
    pcat <- sapply(gene.idlist, function(i) paste(sapply(strsplit(unique(panther[match(i, panther$V3),10]), split="#"), function(j) j[1]), "", sep=""))
    ## get lab category
    lcat <- sapply(pcat, function(i) paste(plc$newLab[match(i, plc$pantherClass, nomatch=0)], collapse=";"))
    ## add lcat to main table
    dmrList$category <- lcat
    ## add ncbi ID to table
    dmrList$ncbi <- sapply(gene.ncbi, function(i) paste(i, collapse=";"))
    
    return(list(dmrList = dmrList, annMat = NA))
    
  } else {
    return(list(dmrList = NA, annMat = NA))
  }
}

