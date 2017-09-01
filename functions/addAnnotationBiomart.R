## Created -/-/2015 by Daniel Beck
## Last modified 4/6/2016

## This function adds annotation information to a DMR table. The annotation information is
## pulled from a remote Biomart database.

addAnnotationBiomart <- function(dmrList, annotationObject, chrPrefix = "", maxDMR = 1000, extension = 0) {
  if (is.null(nrow(dmrList))) {
    return(list(dmrList = NA, annMat = NA))
  }
  if (nrow(dmrList) < maxDMR && nrow(dmrList) > 0) {
      
    annotation <- list()
    for (i in 1:nrow(dmrList)) {
      # There is a sporadic error (possibly some sort of request timeout) that causes the getBM function to fail.
      # This failure is unrelated to the actual function values and rerunning the same function generally works.
      # To accomidate this (fairly low) error rate, I've implimented a while loop to retry the getBM function 
      # when it fails. I added a try limit to avoid infinite loops.
      a <- "start"
      ca <- 1
      tryLimit <- 3
      while((ca == 1) | ((class(a) == "try-error") & ca <= tryLimit)) {
        a <- try(
             annotation[[i]] <- getBM(
                     attributes = c("external_gene_name", "entrezgene", "chromosome_name", 
                                    "start_position", "end_position", "ensembl_gene_id", 
                                    "description"),
                     filters = c("chromosome_name", "start", "end"),
                     values = list(dmrList$chr[i], dmrList$start[i] - extension, dmrList$stop[i] + extension),
                     mart = annotationObject)
        )
        ca <- ca + 1
      }
      if ((ca > tryLimit) & (class(a) == "try-error")) warning(paste("Entry", i, "of the dmrList failed too many times."))
      
      if (length(annotation) < i) {
        annotation[[i]] <- NA
      } else {
        if (nrow(annotation[[i]]) > 0) {
          annotation[[i]] <- cbind(annotation[[i]], ID = dmrList$ID[i])
        }
      }
        
    }
    # Flatten list into dataframe
    if(is.list(annotation)) {
      annMat <- do.call(rbind.data.frame, annotation)
    } else {
      annMat <- NA
      print("annMat failed")
    }

    # Add category information from custom list      
    symCat <- read.csv("/projects/medipPipeline/data/symbolCategory_v3.csv", 
                       stringsAsFactors = F)
    homoloGene <- read.table("/projects/medipPipeline/data/homologene_v1.data", 
                             sep = "\t", stringsAsFactors = F)
    annotationTable <- read.csv("/projects/medipPipeline/data/annotationTable_v3.csv", 
                                  stringsAsFactors = F)
    
    genes <- lapply(annMat$external_gene_name, function(j) { 
      g <- identifyCategory(query = j, homologs = homoloGene, categories = symCat, annotationTable = annotationTable)
      if (nrow(g) < 1) g[1,] <- NA
      if (nrow(g) > 1) {
        g <- t(as.data.frame(apply(g, 2, function(i) paste(unique(i), collapse=";"))))
      }
      g
    })
    
    if (is.list(genes)) {
      genes <- do.call(rbind, genes)
    } else {
      genes <- NA
      print("genes failed")
    }

    annMat <- cbind(annMat, genes)
    # Get rid of any commas, they are annoying when writing to csv
    annMat$description <- gsub(annMat$description, pattern=",", replacement=";")
    
    # Add gene names to dmrList
    dmrList$annotation <- sapply(annotation, function(i) { 
      if (length(i) == 1) {
        ann <- ""
      } else {
        ann <- paste(unique(i$external_gene_name), collapse = ";")
      }
      ann
      })
  }
  return(list(dmrList = dmrList, annMat = annMat))
}

