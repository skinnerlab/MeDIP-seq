## Created -/-/2015 by Daniel Beck
## Last modified 4/6/2016

## This function adds annotation information to a DMR table. The annotation information is
## pulled from a remote Biomart database.

addAnnotationBiomart <- function(dmrList, annotationObject, chrPrefix = "", maxDMR = 1000) {
  if (is.null(nrow(dmrList))) {
    return(list(dmrList = NA, annMat = NA))
  }
  if (nrow(dmrList) < maxDMR && nrow(dmrList) > 0) {
      
    annotation <- list()
    for (i in 1:nrow(dmrList)) {
      try(
           annotation[[i]] <- getBM(
                   attributes = c("external_gene_name", "entrezgene", "chromosome_name", 
                                  "start_position", "end_position", "ensembl_gene_id", 
                                  "description"),
                   filters = c("chromosome_name", "start", "end"),
                   values = list(dmrList$chr[i], dmrList$start[i], dmrList$stop[i]),
                   mart = annotationObject)
      )
      if (nrow(annotation[[i]]) > 0) {
        annotation[[i]] <- cbind(annotation[[i]], ID = dmrList$ID[i])
      }
        
    }
    # Flatten list into dataframe
    annMat <- do.call(rbind.data.frame, annotation)
    
    # Add category information from custom list      
    symCat <- read.csv("/projects/medipPipeline/data/symbolCategory_v1.csv", 
                       stringsAsFactors = F)
    homoloGene <- read.table("/projects/medipPipeline/data/homologene_v1.data", 
                             sep = "\t", stringsAsFactors = F)
    genes <- sapply(annMat$external_gene_name, function(j) { 
      identifyCategory(query = j, homologs = homoloGene, categories = symCat)
    })
    genes <- do.call(rbind, genes)
    colnames(genes) <- c("cglSymbol", "cglCategory")
    annMat <- cbind(annMat, genes)
      
    # Add gene names to dmrList
    dmrList$annotation <- sapply(annotation, function(i) { 
      paste(unique(i$external_gene_name), collapse = ";")})
  }
  return(list(dmrList = dmrList, annMat = annMat))
}

