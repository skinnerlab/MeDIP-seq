cleanAnnMat <- function(annMat) {
  annMat$catEstimate <- sapply(strsplit(annMat$pantherLabCat, split=";"), function(i) i[1])
  lc.temp <- sapply(strsplit(annMat$labCat, split=";"), function(i) i[1])
  annMat$catEstimate[!is.na(lc.temp)] <- lc.temp[!is.na(lc.temp)]
  annMat <- annMat[,c(1, 7, 13, 14, 16, 17)]
  annMat$catEstimate[which(annMat$catEstimate == "miscellaneous and unknown")] <- "Unknown"
  annMat$catEstimate[which(annMat$catEstimate == "receptors and binding proteins")] <- "Receptor"
  annMat$catEstimate[which(annMat$catEstimate == "translation and protein modification")] <- "Translation"
  annMat$catEstimate[which(annMat$catEstimate == "metabolism and transport")] <- "Metabolism"
  annMat$catEstimate[which(annMat$catEstimate == "Miscellaneous")] <- "Unknown"
  annMat$catEstimate[which(annMat$catEstimate == "NA")] <- NA  
  annMat$catEstimate[which(annMat$catEstimate == "signaling")] <- "Signaling"
  annMat$catEstimate[which(annMat$catEstimate == "development")] <- "Development"
  annMat$catEstimate[which(annMat$catEstimate == "transcription")] <- "Transcription"
  annMat$catEstimate[which(annMat$catEstimate == "cytoskeleton")] <- "Cytoskeleton"
  annMat$catEstimate[which(annMat$catEstimate == "proteolysis")] <- "Proteolysis"
  annMat$catEstimate[which(annMat$catEstimate == "growth factor")] <- "Growth factor"
  annMat$catEstimate[which(annMat$catEstimate == "immune response")] <- "Immune response"
  annMat$catEstimate[which(annMat$catEstimate == "epigenetics")] <- "Epigenetics"
  annMat$catEstimate[which(annMat$catEstimate == "apoptosis")] <- "Apoptosis"

  annMat <- annMat[match(unique(annMat$external_gene_name), annMat$external_gene_name),]     # remove duplicates
  annMat$description <- gsub(annMat$description, pattern=",", replacement=";")
  annMat$humanSummary <- gsub(annMat$humanSummary, pattern=",", replacement=";")
  return(annMat)
}
addCatDMR <- function(dmrTable, annMat) {
  dmrcat <- unlist(lapply(lapply(strsplit(dmrTable$annotation, split=";"), function(i) {
    unique(annMat$catEstimate[match(i, annMat$external_gene_name)])
  }), function(j) {
    paste(na.omit(j), collapse=";")
  })
  )
  dmrTable$category <- dmrcat
  dmrTable$dmrName <- paste("DMR", dmrTable$chr, ":", dmrTable$start, sep="")
  dmrTable <- dmrTable[,c("dmrName", "chr", "start", "length", "numSigWin", "minP", "maxLFC", "cpgNum", "cpgDensity", "annotation", "category")]
  return(dmrTable)
}
