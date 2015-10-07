
# This function adds annotation information 
addAnnotationGFF<-function(dmrList, gff, chrPrefix="", maxDMR=1000){
     if (is.null(nrow(dmrList))) return(dmrList)
     if(nrow(dmrList)<maxDMR && nrow(dmrList)>0){
          # save original chromosome names
          ochr<-dmrList$chr
          # add prefix to chr name if necessary
          dmrList$chr<-paste(chrPrefix, dmrList$chr, sep="")
          # convert dmrList to GRanges object
          dmrList$start<-as.numeric(dmrList$start)
          dmrList$stop<-as.numeric(dmrList$stop)
          gdmrList<-makeGRangesFromDataFrame(dmrList)
          # find overlaps
          overlaps<-findOverlaps(gdmrList, gff)
          # add annotation to dmrList
          dmrList<-cbind(dmrList, annotation=NA)
          dmrList$annotation[queryHits(overlaps)]<-as.character(gff$group[subjectHits(overlaps)])
          # put original chromosome names back
          dmrList$chr<-ochr
     }
     return(dmrList)
}

addAnnotationBiomart<-function(dmrList, annotationObject, chrPrefix="", maxDMR=1000){
     if (is.null(nrow(dmrList))) return(list(dmrList=NA, annMat=NA))
     if(nrow(dmrList)<maxDMR && nrow(dmrList)>0){
          # keep original dmrList unchanged
          # dmrListOriginal<-dmrList
          # add prefix to chr name if necessary
          # dmrList$chr<-paste(chrPrefix, dmrList$chr, sep="")
          # download associations
          annotation<-list()
          for (i in 1:nrow(dmrList)) {
               annotation[[i]] <- getBM(attributes=c("ensembl_gene_id", "entrezgene", "chromosome_name", 
                                                  "start_position", "end_position", "external_gene_name", 
                                                  "namespace_1003", "name_1006", "go_id"), 
                                        filters=c("chromosome_name", "start", "end"), 
                                        values=list(dmrList$chr[i], dmrList$start[i], dmrList$stop[i]), 
                                        mart=annotationObject)
               if (nrow(annotation[[i]]) > 0){
                    annotation[[i]] <- cbind(annotation[[i]], ID=dmrList$ID[i])
               }
          }
          # Flatten list into dataframe    
          annMat <- do.call(rbind.data.frame, annotation)
          # Add gene names to dmrList
          dmrList$annotation<-sapply(annotation, function(i) paste(unique(i$external_gene_name), collapse=";"))
     }
     return(list(dmrList=dmrList, annMat=annMat))
}

## This function is taken from the MEDIPS library (MEDIPS.setAnnotation). It causes an error when there are no regions that can be annotated. I've modified the function to not error in this case. THIS FUNCTION NOT NECESSARY WITH NEW ANNOTATIONS (10/6/2015)
# setAnnotation<-function (regions, annotation, cnv = F) {
#      tmp.regions = GRanges(seqnames = regions[, 1], ranges = IRanges(start = as.numeric(regions[, 2]), end = as.numeric(regions[, 3])))
#      ans = NULL
#      if (is.data.frame(annotation)) 
#           annotation = list(annotation = annotation)
#      for (n in names(annotation)) {
#           anno = annotation[[n]]
#           anno.data = GRanges(anno[, 2], ranges = IRanges(start = anno[, 3], end = anno[, 4]), ID = anno[, 1])
#           overlapsM = IRanges::as.matrix(findOverlaps(tmp.regions, anno.data))
#           splitL = split(overlapsM[, 2], overlapsM[, 1])
#           maxEle = max(c(0, unlist(lapply(splitL, length))))
#           if (maxEle == 0) {
#                message("no \"", n, "\"annotation overlap the provided regions")
#           } else {
#                tmp.ans = matrix(ncol = maxEle, nrow = length(tmp.regions))
#                colnames(tmp.ans) = paste(c(1:maxEle), "_", names(anno)[1], sep = "")
#                j = rep(names(splitL), sapply(splitL, length))
#                k = unlist(lapply(splitL, function(x) {1:length(x)}))
#                tmp.ans[matrix(c(as.integer(j), as.integer(as.vector(k))), ncol = 2)] = as.character(values(anno.data)[overlapsM[, 2], 1])
#                ans = cbind(ans, tmp.ans)
#           }
#      }
#      if (cnv) {
#           ans = data.frame(regions, CNV.log2.ratio = as.numeric(ans), stringsAsFactors = F)
#      } else {
#           if (!is.null(nrow(ans))){      # I added this test
#                ans = cbind(regions, ans, stringsAsFactors = F)
#           } else {
#                ans <- regions           # If the test fails return the unannotated dmrList
#           }
#      }
#      return(ans)
# }
