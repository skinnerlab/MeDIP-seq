## Daniel Beck
## Created 10/14/2015

# These functions are used to determine and extract overlapping DMR / BP

# Takes as many dmr tables as necessary and returns overlapping BP
# This was the original function. It was too slow on large overlap regions.
overlappingBPo <- function(dmrList) {
  # Split DMR tables by chromosome
  sdmrList <- lapply(dmrList, 
                     function(i) {
                       temp <- split.data.frame(i, f = i$chr)
                       i <- temp[which(sapply(temp, nrow) > 0)]
                     })
  # Find chromosomes shared by all dmrLists
  chrCounts <- table(unlist(lapply(sdmrList, names)))
  commonChr <- names(chrCounts)[chrCounts == length(sdmrList)]
  
  # Loop over common chromosomes and extract all overlapping BP
  overlaps <- list()  # Overlaping BP are stored in this list
  for (chr in 1:length(commonChr)) {
    smallList <- lapply(sdmrList, function(i) i[[which(names(i)==commonChr[chr])]])
    # Expand BP ranges into vectors
    expandedBP <- lapply(smallList, 
                         function(i) lapply(1:nrow(i), 
                                            function(j) j<-i$start[j]:i$stop[j]))
    # Unlist BP so they are not separated by DMR
    dmrBP <- lapply(expandedBP, function(i) unlist(i))
    # Lump all bp together then table them to determine which ones overlap
    allBP <- table(unlist(dmrBP))  # This is the slowest step so far
    # Extract BP that are present in all DMR tables
    overlaps[[chr]] <- as.numeric(names(allBP[which(allBP == length(dmrList))]))
  }
  names(overlaps)<-commonChr
  # Return list of overlapping basepairs
  return(overlaps)
}

# Takes as many dmr tables as necessary and returns overlapping BP
overlappingBP <- function(dmrList) {
  # Split DMR tables by chromosome
  sdmrList <- lapply(dmrList, 
                     function(i) {
                       temp <- split.data.frame(i, f = i$chr)
                       i <- temp[which(sapply(temp, nrow) > 0)]
                     })
  # Find chromosomes shared by all dmrLists
  chrCounts <- table(unlist(lapply(sdmrList, names)))
  commonChr <- names(chrCounts)[chrCounts == length(sdmrList)]
  
  # Loop over common chromosomes and extract all overlapping BP
  overlaps <- list()  # Overlaping BP are stored in this list
  for (chr in 1:length(commonChr)) {
    smallList <- lapply(sdmrList, function(i) i[[which(names(i) == commonChr[chr])]])
    overInt <- cbind(smallList[[1]]$start, smallList[[1]]$stop)
    for (i in 2:length(smallList)) {
      if (!is.null(overInt)){
        overInt <- pairOverlap(start1 = overInt[,1], stop1 = overInt[,2], 
                               start2 = smallList[[i]]$start, stop2 = smallList[[i]]$stop)
      }
    }
    # Expand BP ranges into vectors
    if (!is.null(overInt)) {
      overlaps[[chr]] <- unlist(apply(overInt, 1, function(i) i[1]:i[2]))
    } else {
      overlaps[[chr]] <- NA
    }
  }
  names(overlaps) <- commonChr
  #overlaps<-lapply(overlaps, function(i) if (is.na(i)) i<-NULL)
  # Return list of overlapping basepairs
  return(overlaps)
}

# Takes start and stop positions and returns overlaps. Positions must be equilivant (on same chromosome)
pairOverlap <- function(start1, stop1, start2, stop2) {
  if (length(start1) != length(stop1)) stop("start1 and stop1 must have identical length")
  if (length(start2) != length(stop2)) stop("start2 and stop2 must have identical length")
  overlaps <- list()  # list to hold overlaping regions
  index <- 1  # next index in overlaps list that is currently empty
  for (pos1 in 1:length(start1)) {
    for (pos2 in 1:length(start2)) {
      if ((start1[pos1] < stop2[pos2]) & (stop1[pos1] > start2[pos2])) {
        overlaps[[index]] <- c(max(start1[pos1], start2[pos2]), min(stop1[pos1], stop2[pos2]))
        index <- index + 1
      }
    }
  }
  overlaps <- do.call(rbind, overlaps)
  if (!is.null(overlaps)) {
    colnames(overlaps) <- c("start", "stop")
  }
  return(overlaps)
}

# Takes a list of base pairs (separated by chromosome) and a DMR table and outputs ta table of DMR that overlap with those base pairs

bpToDmr <- function(dmrTable, bpList) {
  chrList <- split.data.frame(dmrTable, f = dmrTable$chr)
  chrList <- chrList[which(sapply(chrList, nrow) > 0)]
  commonChrs<-names(which(table(c(unique(dmrTable$chr), names(bpList)))==2))
  chrList<-chrList[match(commonChrs, names(chrList))]
  bpList<-bpList[match(commonChrs, names(bpList))]
  resultList <- list()
  # Loop over all chromosomes
  for (chr in 1:length(chrList)) {
    resultID <- lapply(bpList[[chr]], 
                       function(i) chrList[[chr]]$ID[which((chrList[[chr]]$start <= i) & 
                                                             (chrList[[chr]]$stop >= i))])
    resultList[[chr]] <- chrList[[chr]][match(unique(unlist(resultID)), chrList[[chr]]$ID), ]
  }
  resultTable <- do.call(rbind, resultList)
  return(resultTable)
}

# Takes two matrices or data.frames, removes lines in mat1 that are also in mat2. 
# Returns mat1 without lines in mat2
removeDuplicates <- function(mat1, mat2) {
  # If mat1 is NULL or NA, return empty matrix (so that nrow still works and returns 0)
  if (is.null(mat1)) return(cbind(list.index=1, dmr.index=1)[0,])
  if (is.na(mat1)) return(cbind(list.index=1, dmr.index=1)[0,])
  # If mat2 is NULL or NA, return mat1
  if (is.null(mat2)) return(mat1)
  if (is.na(mat2)) return(mat1)
  
  mat1c <- apply(mat1, 1, paste, collapse="_")
  mat2c <- apply(mat2, 1, paste, collapse="_")
  
  overlaps<-match(mat2c, mat1c, nomatch=0)
  # If there are no overlaps, the subset returns an empty matrix.
  # This line prevents this by removing a non-existant line (out of bounds).
  # This method seems strange and may be fragile.
  if ((unique(overlaps)[1] == 0) & (length(unique(overlaps)) == 1)) overlaps <- nrow(mat1) + 10
  uniqueSet <- mat1[-c(overlaps),]
  
  return(uniqueSet)
}




vennDMR <- function(dmr.list, group.names=NULL, ...) {
  if (is.null(group.names)) group.names <- names(dmr.list)  # First try to extract names from dmrList
  if (is.null(group.names)) group.names <- as.character(1:length(dmr.list))  # Then just make up names
  
  if (length(dmr.list) < 2) stop("Can't make a single group Venn. Try drawing a circle instead.")
  
  # First get list of all DMR in a table
  unique <- list()
  for (i in 1:length(dmr.list)) {
    unique[[i]] <- cbind(list.index = i, dmr.index = 1:nrow(dmr.list[[i]]))
  }
 
  # Next get list of DMRs involved in pairwise overlaps
  pairs <- list()
  pair.flag <- 1  # Counter for list index
  for (i in 1:(length(dmr.list) - 1)) {
    for (j in (i + 1):length(dmr.list)) {
      bps <- overlappingBP(dmr.list[c(i, j)])
      dmrs1 <- bpToDmr(dmrTable = dmr.list[[i]], bpList = bps)
      dmrs2 <- bpToDmr(dmrTable = dmr.list[[j]], bpList = bps)
      ids1 <- cbind(list.index=i, dmr.index=match(dmrs1$ID, dmr.list[[i]]$ID))
      ids2 <- cbind(list.index=j, dmr.index=match(dmrs2$ID, dmr.list[[j]]$ID))
      pairs[[pair.flag]] <- rbind(ids1, ids2)
      pair.flag <- pair.flag + 1  # Increment index counter
    }
  }
  pairs[which(sapply(pairs, ncol) == 1)] <- NA
  
  # If there are at least 3 groups, count triple overlaps
  if (length(dmr.list) > 2) {
    # list of DMRs involved in triple overlaps
    triples <- list()
    triple.flag <- 1  # Counter for list index
    for (i in 1:(length(dmr.list) - 2)) {
      for (j in (i + 1):(length(dmr.list) - 1)) {
        for (k in (j + 1):length(dmr.list)) {
          bps <- overlappingBP(dmr.list[c(i, j, k)])
          dmrs1 <- bpToDmr(dmrTable = dmr.list[[i]], bpList = bps)
          dmrs2 <- bpToDmr(dmrTable = dmr.list[[j]], bpList = bps)
          dmrs3 <- bpToDmr(dmrTable = dmr.list[[k]], bpList = bps)
          ids1 <- cbind(list.index=i, dmr.index=match(dmrs1$ID, dmr.list[[i]]$ID))
          ids2 <- cbind(list.index=j, dmr.index=match(dmrs2$ID, dmr.list[[j]]$ID))
          ids3 <- cbind(list.index=k, dmr.index=match(dmrs3$ID, dmr.list[[k]]$ID))
          triples[[triple.flag]] <- rbind(ids1, ids2, ids3)
          triple.flag <- triple.flag + 1  # Increment index counter
        }
      }
    }
    triples[which(sapply(triples, ncol) == 1)] <- NA
  }
  # If there are at least 4 groups, count quadruple overlaps
  if (length(dmr.list) > 3) {
    # list of DMRs involved in triple overlaps
    quads <- list()
    quad.flag <- 1  # Counter for list index
    for (i in 1:(length(dmr.list) - 3)) {
      for (j in (i + 1):(length(dmr.list) - 2)) {
        for (k in (j + 1):(length(dmr.list) - 1)) {
          for (l in (k + 1):length(dmr.list)) {
            bps <- overlappingBP(dmr.list[c(i, j, k, l)])
            dmrs1 <- bpToDmr(dmrTable = dmr.list[[i]], bpList = bps)
            dmrs2 <- bpToDmr(dmrTable = dmr.list[[j]], bpList = bps)
            dmrs3 <- bpToDmr(dmrTable = dmr.list[[k]], bpList = bps)
            dmrs4 <- bpToDmr(dmrTable = dmr.list[[l]], bpList = bps)
            ids1 <- cbind(list.index=i, dmr.index=match(dmrs1$ID, dmr.list[[i]]$ID))
            ids2 <- cbind(list.index=j, dmr.index=match(dmrs2$ID, dmr.list[[j]]$ID))
            ids3 <- cbind(list.index=k, dmr.index=match(dmrs3$ID, dmr.list[[k]]$ID))
            ids4 <- cbind(list.index=l, dmr.index=match(dmrs4$ID, dmr.list[[l]]$ID))
            quads[[quad.flag]] <- rbind(ids1, ids2, ids3, ids4)
            quad.flag <- quad.flag + 1  # Increment index counter
          }
        }
      }
    }
    quads[which(sapply(quads, ncol) == 1)] <- NA
  }
  # If there are at least 5 groups, count quintuple overlaps
  if (length(dmr.list) > 4) {
    # list of DMRs involved in triple overlaps
    quins <- list()
    quin.flag <- 1  # Counter for list index
    for (i in 1:(length(dmr.list) - 4)) {
      for (j in (i + 1):(length(dmr.list) - 3)) {
        for (k in (j + 1):(length(dmr.list) - 2)) {
          for (l in (k + 1):(length(dmr.list) - 1)) {
            for (m in (l + 1):length(dmr.list)) {
              bps <- overlappingBP(dmr.list[c(i, j, k, l, m)])
              dmrs1 <- bpToDmr(dmrTable = dmr.list[[i]], bpList = bps)
              dmrs2 <- bpToDmr(dmrTable = dmr.list[[j]], bpList = bps)
              dmrs3 <- bpToDmr(dmrTable = dmr.list[[k]], bpList = bps)
              dmrs4 <- bpToDmr(dmrTable = dmr.list[[l]], bpList = bps)
              dmrs5 <- bpToDmr(dmrTable = dmr.list[[m]], bpList = bps)
              ids1 <- cbind(list.index=i, dmr.index=match(dmrs1$ID, dmr.list[[i]]$ID))
              ids2 <- cbind(list.index=j, dmr.index=match(dmrs2$ID, dmr.list[[j]]$ID))
              ids3 <- cbind(list.index=k, dmr.index=match(dmrs3$ID, dmr.list[[k]]$ID))
              ids4 <- cbind(list.index=l, dmr.index=match(dmrs4$ID, dmr.list[[l]]$ID))
              ids5 <- cbind(list.index=m, dmr.index=match(dmrs5$ID, dmr.list[[m]]$ID))
              quins[[quin.flag]] <- rbind(ids1, ids2, ids3, ids4, ids5)
              quin.flag <- quin.flag + 1  # Increment index counter
            }
          }
        }
      }
    }
    quins[which(sapply(quins, ncol) == 1)] <- NA
  }
  # Filter DMRs so they only occur on highest level of overlap and plot Venn
  if (length(dmr.list) == 2) {
    pM <- do.call(rbind, pairs)
    if (ncol(pM)>1) {
      pM <- pM[which(!is.na(pM[ ,1])), ]  # Get rid of NA rows
      n1 <- nrow(removeDuplicates(unique[[1]], pM))
      n2 <- nrow(removeDuplicates(unique[[2]], pM))
      n12 <- nrow(pM) / 2
    } else {
      n1 <- nrow(unique[[1]])
      n2 <- nrow(unique[[2]])
      n12 <- 0
    }
    area1 <- n1 + n12; area2 <- n2 + n12;
    plot.new()
    draw.pairwise.venn(area1=area1, area2=area2, cross.area=n12, 
                       category = group.names, col=rep("black", 2), 
                       fill=c("skyblue", "orange"), ...)
  }
  
  if (length(dmr.list) == 3) {
    pM <- do.call(rbind, pairs)
    if (ncol(pM) > 1) {
      pM <- pM[which(!is.na(pM[ ,1])), ]  # Get rid of NA rows
      n1 <- nrow(removeDuplicates(unique[[1]], pM))
      n2 <- nrow(removeDuplicates(unique[[2]], pM))
      n3 <- nrow(removeDuplicates(unique[[3]], pM))
    } else {
      n1 <- nrow(unique[[1]])
      n2 <- nrow(unique[[2]])
      n3 <- nrow(unique[[3]])
    }
    tM <- do.call(rbind, triples)
    if (ncol(tM) > 1) {
      tM <- tM[which(!is.na(tM[ ,1])), ]  # Get rid of NA rows
      
      n12 <- nrow(removeDuplicates(pairs[[1]], tM)) / 2
      n13 <- nrow(removeDuplicates(pairs[[2]], tM)) / 2
      n23 <- nrow(removeDuplicates(pairs[[3]], tM)) / 2
      n123 <- nrow(tM) / 3
    } else {
      n12 <- max(nrow(pairs[[1]]), 0) / 2
      n13 <- max(nrow(pairs[[2]]), 0) / 2
      n23 <- max(nrow(pairs[[3]]), 0) / 2
      n123 <- 0
    }
    area1 <- n1 + n12 + n13 + n123; area2 <- n2 + n12 + n23 + n123; area3 <- n3 + n13 + n23 + n123
    plot.new()
    draw.triple.venn(area1 = area1, area2 = area2, area3 = area3, n12 = n12 + n123, 
                     n13 = n13 + n123, n23 = n23 + n123, n123 = n123, category = group.names, 
                     col=rep("black", 3), fill=c("skyblue", "pink1", "orange"), ...)
  }
  if (length(dmr.list) == 4) {
    pM <- do.call(rbind, pairs)
    if (ncol(pM) > 1) {
      pM <- pM[which(!is.na(pM[ ,1])), ]  # Get rid of NA rows
      n1 <- nrow(removeDuplicates(unique[[1]], pM))
      n2 <- nrow(removeDuplicates(unique[[2]], pM))
      n3 <- nrow(removeDuplicates(unique[[3]], pM))
      n4 <- nrow(removeDuplicates(unique[[4]], pM))
    } else {
      n1 <- nrow(unique[[1]])
      n2 <- nrow(unique[[2]])
      n3 <- nrow(unique[[3]])
      n4 <- nrow(unique[[4]])
    }
    tM <- do.call(rbind, triples)
    if (ncol(tM) > 1) {
      tM <- tM[which(!is.na(tM[ ,1])), ]  # Get rid of NA rows
      n12 <- nrow(removeDuplicates(pairs[[1]], tM)) / 2
      n13 <- nrow(removeDuplicates(pairs[[2]], tM)) / 2
      n14 <- nrow(removeDuplicates(pairs[[3]], tM)) / 2
      n23 <- nrow(removeDuplicates(pairs[[4]], tM)) / 2
      n24 <- nrow(removeDuplicates(pairs[[5]], tM)) / 2
      n34 <- nrow(removeDuplicates(pairs[[6]], tM)) / 2
    } else {
      n12 <- max(nrow(pairs[[1]]), 0) / 2
      n13 <- max(nrow(pairs[[2]]), 0) / 2
      n14 <- max(nrow(pairs[[3]]), 0) / 2
      n23 <- max(nrow(pairs[[4]]), 0) / 2
      n24 <- max(nrow(pairs[[5]]), 0) / 2
      n34 <- max(nrow(pairs[[6]]), 0) / 2
    }
    qM <- do.call(rbind, quads)
    if (ncol(qM) > 1) {
      qM <- qM[which(!is.na(qM[ ,1])), ]  # Get rid of NA rows
      n123 <- nrow(removeDuplicates(triples[[1]], qM)) / 3
      n124 <- nrow(removeDuplicates(triples[[2]], qM)) / 3
      n134 <- nrow(removeDuplicates(triples[[3]], qM)) / 3
      n234 <- nrow(removeDuplicates(triples[[4]], qM)) / 3
      n1234 <- nrow(qM) / 4
    } else {
      n123 <- max(nrow(triples[[1]]), 0) / 3
      n124 <- max(nrow(triples[[2]]), 0) / 3
      n134 <- max(nrow(triples[[3]]), 0) / 3
      n234 <- max(nrow(triples[[4]]), 0) / 3
      n1234 <- 0
    }
    area1 <- n1 + n12 + n13 + n14 + n123 + n124 + n134 + n1234
    area2 <- n2 + n12 + n23 + n24 + n123 + n124 + n234 + n1234
    area3 <- n3 + n13 + n23 + n34 + n123 + n134 + n234 + n1234 
    area4 <- n4 + n14 + n24 + n34 + n124 + n134 + n234 + n1234
    plot.new()
    draw.quad.venn(area1 = area1, area2 = area2, area3 = area3, area4 = area4, 
                   n12 = n12 + n123 + n124 + n1234, n13 = n13 + n123 + n134 + n1234, 
                   n14 = n14 + n124 + n134 + n1234, n23 = n23 + n123 + n234 + n1234, 
                   n24 = n24 + n124 + n234 + n1234, n34 = n34 + n134 + n234 + n1234, 
                   n123 = n123 + n1234, n124 = n124 + n1234, n134 = n134 + n1234, 
                   n234 = n234 + n1234, n1234 = n1234, category=group.names, col=rep("black", 4), 
                   fill=c("skyblue", "pink1", "mediumorchid", "orange"), ...)
  }
  if (length(dmr.list) == 5) {
    pM <- do.call(rbind, pairs)
    if (ncol(pM) > 1) {
      pM <- pM[which(!is.na(pM[ ,1])), ]  # Get rid of NA rows
      n1 <- nrow(removeDuplicates(unique[[1]], pM))
      n2 <- nrow(removeDuplicates(unique[[2]], pM))
      n3 <- nrow(removeDuplicates(unique[[3]], pM))
      n4 <- nrow(removeDuplicates(unique[[4]], pM))
      n5 <- nrow(removeDuplicates(unique[[5]], pM))
    } else {
      n1 <- nrow(unique[[1]])
      n2 <- nrow(unique[[2]])
      n3 <- nrow(unique[[3]])
      n4 <- nrow(unique[[4]])
      n5 <- nrow(unique[[5]])
    }
    tM <- do.call(rbind, triples)
    if (ncol(tM) > 1) {
      tM <- tM[which(!is.na(tM[ ,1])), ]  # Get rid of NA rows
      n12 <- nrow(removeDuplicates(pairs[[1]], tM)) / 2
      n13 <- nrow(removeDuplicates(pairs[[2]], tM)) / 2
      n14 <- nrow(removeDuplicates(pairs[[3]], tM)) / 2
      n15 <- nrow(removeDuplicates(pairs[[4]], tM)) / 2
      n23 <- nrow(removeDuplicates(pairs[[5]], tM)) / 2
      n24 <- nrow(removeDuplicates(pairs[[6]], tM)) / 2
      n25 <- nrow(removeDuplicates(pairs[[7]], tM)) / 2
      n34 <- nrow(removeDuplicates(pairs[[8]], tM)) / 2
      n35 <- nrow(removeDuplicates(pairs[[9]], tM)) / 2
      n45 <- nrow(removeDuplicates(pairs[[10]], tM)) / 2
    } else {
      n12 <- max(nrow(pairs[[1]]), 0) / 2
      n13 <- max(nrow(pairs[[2]]), 0) / 2
      n14 <- max(nrow(pairs[[3]]), 0) / 2
      n15 <- max(nrow(pairs[[4]]), 0) / 2
      n23 <- max(nrow(pairs[[5]]), 0) / 2
      n24 <- max(nrow(pairs[[6]]), 0) / 2
      n25 <- max(nrow(pairs[[7]]), 0) / 2
      n34 <- max(nrow(pairs[[8]]), 0) / 2
      n35 <- max(nrow(pairs[[9]]), 0) / 2
      n45 <- max(nrow(pairs[[10]]), 0) / 2
      
    }
    qM <- do.call(rbind, quads)
    if (ncol(qM) > 1) {
      qM <- qM[which(!is.na(qM[ ,1])), ]  # Get rid of NA rows
      n123 <- nrow(removeDuplicates(triples[[1]], qM)) / 3
      n124 <- nrow(removeDuplicates(triples[[2]], qM)) / 3
      n125 <- nrow(removeDuplicates(triples[[3]], qM)) / 3
      n134 <- nrow(removeDuplicates(triples[[4]], qM)) / 3
      n135 <- nrow(removeDuplicates(triples[[5]], qM)) / 3
      n145 <- nrow(removeDuplicates(triples[[6]], qM)) / 3
      n234 <- nrow(removeDuplicates(triples[[7]], qM)) / 3
      n235 <- nrow(removeDuplicates(triples[[8]], qM)) / 3
      n245 <- nrow(removeDuplicates(triples[[9]], qM)) / 3
      n345 <- nrow(removeDuplicates(triples[[10]], qM)) / 3
      
    } else {
      n123 <- max(nrow(triples[[1]]), 0) / 3
      n124 <- max(nrow(triples[[2]]), 0) / 3
      n125 <- max(nrow(triples[[3]]), 0) / 3
      n134 <- max(nrow(triples[[4]]), 0) / 3
      n135 <- max(nrow(triples[[5]]), 0) / 3
      n145 <- max(nrow(triples[[6]]), 0) / 3
      n234 <- max(nrow(triples[[7]]), 0) / 3
      n235 <- max(nrow(triples[[8]]), 0) / 3
      n245 <- max(nrow(triples[[9]]), 0) / 3
      n345 <- max(nrow(triples[[10]]), 0) / 3
    }
    qnM <- do.call(rbind, quins)
    if (ncol(qnM) > 1) {
      qnM <- qnM[which(!is.na(qnM[ ,1])), ]  # Get rid of NA rows
      n1234 <- nrow(removeDuplicates(quads[[1]], qnM)) / 4
      n1235 <- nrow(removeDuplicates(quads[[2]], qnM)) / 4
      n1245 <- nrow(removeDuplicates(quads[[3]], qnM)) / 4
      n1345 <- nrow(removeDuplicates(quads[[4]], qnM)) / 4
      n2345 <- nrow(removeDuplicates(quads[[5]], qnM)) / 4
      n12345 <- nrow(qnM) / 5
    } else {
      n1234 <- max(nrow(quads[[1]]), 0) / 4
      n1235 <- max(nrow(quads[[2]]), 0) / 4
      n1245 <- max(nrow(quads[[3]]), 0) / 4
      n1345 <- max(nrow(quads[[4]]), 0) / 4
      n2345 <- max(nrow(quads[[5]]), 0) / 4
      n12345 <- 0
    }
    area1 <- n1 + n12 + n13 + n14 + n15 + n123 + n124 + n125 + 
             n134 + n135 + n145 + n1234 + n1235 + n1245 + n1345 + n12345
    area2 <- n2 + n12 + n23 + n24 + n25 + n123 + n124 + n125 + 
             n234 + n235 + n245 + n1234 + n1235 + n1245 + n2345 + n12345
    area3 <- n3 + n13 + n23 + n34 + n35 + n123 + n134 + n135 + 
             n234 + n235 + n345 + n1234 + n1235 + n1345 + n2345 + n12345
    area4 <- n4 + n14 + n24 + n34 + n45 + n124 + n134 + n145 + 
             n234 + n245 + n345 + n1234 + n1245 + n1345 + n2345 + n12345
    area5 <- n5 + n15 + n25 + n35 + n45 + n125 + n135 + n145 +
             n235 + n245 + n345 + n1235 + n1245 + n1345 + n2345 + n12345
    
    plot.new()
    draw.quintuple.venn(area1 = area1, area2 = area2, area3 = area3, area4 = area4, area5 = area5,
                        n12 = n12 + n123 + n124 + n125 + n1234 + n1235 + n1245 + n12345, 
                        n13 = n13 + n123 + n134 + n135 + n1234 + n1235 + n1345 + n12345,
                        n14 = n14 + n124 + n134 + n145 + n1234 + n1245 + n1345 + n12345, 
                        n15 = n15 + n125 + n135 + n145 + n1235 + n1245 + n1345 + n12345, 
                        n23 = n23 + n123 + n234 + n235 + n1234 + n1235 + n2345 + n12345,
                        n24 = n24 + n124 + n234 + n245 + n1234 + n1245 + n2345 + n12345, 
                        n25 = n25 + n125 + n135 + n145 + n1235 + n1245 + n2345 + n12345, 
                        n34 = n34 + n134 + n234 + n345 + n1234 + n1345 + n2345 + n12345, 
                        n35 = n35 + n135 + n235 + n345 + n1235 + n1345 + n2345 + n12345, 
                        n45 = n45 + n145 + n245 + n345 + n1245 + n1345 + n2345 + n12345, 
                        n123 = n123 + n1234 + n1235 + n12345, 
                        n124 = n124 + n1234 + n1245 + n12345, 
                        n125 = n125 + n1235 + n1245 + n12345, 
                        n134 = n134 + n1234 + n1345 + n12345,  
                        n135 = n135 + n1235 + n1345 + n12345, 
                        n145 = n145 + n1245 + n1345 + n12345,
                        n234 = n234 + n1234 + n2345 + n12345, 
                        n235 = n235 + n1235 + n2345 + n12345,
                        n245 = n245 + n1245 + n2345 + n12345, 
                        n345 = n345 + n1345 + n2345 + n12345,
                        n1234 = n1234 + n12345, n1235 = n1235 + n12345, n1245 = n1245 + n12345, 
                        n1345 = n1345 + n12345, n2345 = n2345 + n12345, n12345 = n12345, 
                        category=group.names,  
                        fill=c("skyblue", "pink1", "mediumorchid", "orange", "green"), ...)
  }
  
}
