## Created 12/2/2015 by Daniel Beck
## Modified

## This script will modify the DMRs for the steelhead datasets to keep them from spanning scaffolds.
## Annotation should be done after this script to ensure correct gene calls.

## abstract to function that can be incorporated into medipProcessing
fixDMR <- function(methResults, methList, methListETC, pValues, MTC=FALSE, scaffolds, maxDMR=5000) {
     
     # add ID column to scaffolds. this is necessary for overlap functions
     scaffolds$ID <- 1:nrow(scaffolds)
     
     # error if input does not match
     if (length(pValues) != length(methList)) stop("Number of p-values does not match length of DMR table")
     modifiedML <- list()
     # loop over all p-values
     for (pV in 1:length(pValues)) {
          print(paste("Fixing DMR edges for threshold: ", pValues[pV], sep=""))
          # get data.frame of all significant windows
          if (MTC) {
               sigWin <- methResults[which(methResults$edgeR.adj.p.value < pValues[pV]),]
          } else {
               sigWin <- methResults[which(methResults$edgeR.p.value < pValues[pV]),]
          }
          if (nrow(sigWin) > 0) {
          # loop over every significant window
          dmrList <- list()
          for (i in 1:nrow(sigWin)){
     	     # find associated DMR
     	     dmri <- bpToDmr(bpList=overlappingBP(list(methList[[pV]], sigWin[i, ])), dmrTable=methList[[pV]]) 
     	     # find associated scaffold
     	     scai<- bpToDmr(bpList=overlappingBP(list(scaffolds, sigWin[i, ])), dmrTable=scaffolds)
     	     # This conditional accounts for when there is no associated scaffold. This occurs
     	     # when the significant window is between contigs. These DMRs are removed.
     	     if (nrow(scai) > 0) {
     	       dmri$start <- max(scai$start, dmri$start)
     	       dmri$stop <- min(scai$stop, dmri$stop)
     	     } else {
     	       dmri <- dmri[0,]
     	     }
     	     dmrList[[i]] <- dmri
     	     if (i/100 == floor(i/100)) print(i/nrow(sigWin))
          }
          dmrList <- do.call(rbind, dmrList)
          dmrList <- split.data.frame(dmrList, f=dmrList$ID)
          # flag for if the DMR was split
          multiple <- sapply(dmrList, function(i) length(table(i$start)))
          # reduce dmr list to unique rows
          uniqueRows <- lapply(dmrList, function(i) i[!duplicated(i),])
          # modify ID to remain unique with split DMR
          uniqueRows <- lapply(uniqueRows, function(i) {
                               if (nrow(i)>1) i$ID <- paste(i$ID, letters[1:nrow(i)], sep="_")
                               i
                               })
          # transform to data.frame
          dmrList2 <- do.call(rbind, uniqueRows)
          # recalculate length
          dmrList2$length <- dmrList2$stop-dmrList2$start+1
          # recalculate CpG density information
          dmrList2 <- calcCpGdensity(dmrList=dmrList2, maxDMR=maxDMR)
          # the first step in recalculating everything else is to associate the dmrList2 IDs with the
          # methListEtc IDs. This requires stripping off the letters added to the dmrList2 IDs

          dl2id <- sapply(strsplit(dmrList2$ID, split="_"), function(i) i[2])
          mleid <- sapply(strsplit(row.names(methListETC[[pV]]), split="_"), function(i) i[2])
          
          # mle contains the reordered Etc information
          mle <- as.data.frame(methListETC[[pV]][match(dl2id, mleid),], stringsAsFactors=F)

          for (i in 1:nrow(dmrList2)) {
               sts<-sort(as.numeric(unlist(strsplit(mle$start[i], split=";"))))
               sps<-sort(as.numeric(unlist(strsplit(mle$stop[i], split=";"))))
               # identify windows within dmrList stop and start site
               windows <- which(((sts>=dmrList2$start[i])&(sts<=dmrList2$stop[i]))|
                                ((dmrList2$start[i]>=sts)&(dmrList2$start[i]<=sps)))
               if (MTC) {
                    minP <- min(as.numeric(unlist(strsplit(mle$edgeR.adj.p.value[i], split=";")))[windows])
                    numSigWin <- sum(as.numeric(unlist(strsplit(mle$edgeR.adj.p.value[i], split=";")))[windows]<pValues[pV])
               } else {
                    minP <- min(as.numeric(unlist(strsplit(mle$edgeR.p.value[i], split=";")))[windows])
                    numSigWin <- sum(as.numeric(unlist(strsplit(mle$edgeR.p.value[i], split=";")))[windows]<pValues[pV]) 
               }
               dmrList2$minP[i] <- minP
               dmrList2$numSigWin[i] <- numSigWin
          }

          # This removes annotation, if present. I'm not sure it is necessary 
          a <- max(which(colnames(dmrList2)=="annotation"), 0)
          if (a) dmrList2 <- dmrList2[,-c(a)]
          modifiedML[[pV]] <- dmrList2
          } else {
            modifiedML[[pV]] <- NA
          }
     }
     return(modifiedML)
}

