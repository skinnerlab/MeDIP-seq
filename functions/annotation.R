

# This function adds annotation information
addAnnotationGFF <- function(dmrList, gff, chrPrefix = "", maxDMR = 1000) {
  if (is.null(nrow(dmrList)))
    return(dmrList)
  if (nrow(dmrList) < maxDMR && nrow(dmrList) > 0) {
    # save original chromosome names
    ochr <- dmrList$chr
    # add prefix to chr name if necessary
    dmrList$chr <- paste(chrPrefix, dmrList$chr, sep = "")
    # convert dmrList to GRanges object
    dmrList$start <- as.numeric(dmrList$start)
    dmrList$stop <- as.numeric(dmrList$stop)
    gdmrList <- makeGRangesFromDataFrame(dmrList)
    # find overlaps
    overlaps <- findOverlaps(gdmrList, gff)
    # add annotation to dmrList
    dmrList <- cbind(dmrList, annotation = NA)
    dmrList$annotation[queryHits(overlaps)] <- as.character(gff$group[subjectHits(overlaps)])
    # put original chromosome names back
    dmrList$chr <- ochr
  }
  return(dmrList)
}

addAnnotationBiomart <- function(dmrList, annotationObject, chrPrefix = "", maxDMR = 1000) {
    if (is.null(nrow(dmrList))) {
      return(list(dmrList = NA, annMat = NA))
    }
    if (nrow(dmrList) < maxDMR && nrow(dmrList) > 0) {
      # keep original dmrList unchanged
      # dmrListOriginal<-dmrList
      # add prefix to chr name if necessary
      # dmrList$chr<-paste(chrPrefix, dmrList$chr, sep="")
      # download associations
      annotation <- list()
      for (i in 1:nrow(dmrList)) {
        annotation[[i]] <-
          getBM(attributes = c("external_gene_name", "entrezgene", "chromosome_name", 
                               "start_position", "end_position",
                               "ensembl_gene_id", "description"),
                filters = c("chromosome_name", "start", "end"),
                values = list(dmrList$chr[i], dmrList$start[i], dmrList$stop[i]),
                mart = annotationObject)
        if (nrow(annotation[[i]]) > 0) {
          annotation[[i]] <- cbind(annotation[[i]], ID = dmrList$ID[i])
        }
      }
      # Flatten list into dataframe
      annMat <- do.call(rbind.data.frame, annotation)
      # Add gene names to dmrList
      dmrList$annotation <- sapply(annotation, function(i) { 
        paste(unique(i$external_gene_name), collapse = ";")})
    }
    return(list(dmrList = dmrList, annMat = annMat))
  }
