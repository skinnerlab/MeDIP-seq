## Created 10/14/2015 by Daniel Beck
## Last modified 4/6/2016

## This is the main Venn diagram function. It can produce Venn diagrams with up to five groups.
## It calls several associated functions, however, I've split them into seperate files since
## the other functions are often used for other purposes.

vennDMR <- function(dmr.list, group.names = NULL, ...) {
  # First try to extract names from dmrList
  if (is.null(group.names)) group.names <- names(dmr.list)  
  # Then just make up names if that doesn't work
  if (is.null(group.names)) group.names <- as.character(1:length(dmr.list))  
  
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
      ids1 <- cbind(list.index = i, dmr.index = match(dmrs1$ID, dmr.list[[i]]$ID))
      ids2 <- cbind(list.index = j, dmr.index = match(dmrs2$ID, dmr.list[[j]]$ID))
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
          ids1 <- cbind(list.index = i, dmr.index = match(dmrs1$ID, dmr.list[[i]]$ID))
          ids2 <- cbind(list.index = j, dmr.index = match(dmrs2$ID, dmr.list[[j]]$ID))
          ids3 <- cbind(list.index = k, dmr.index = match(dmrs3$ID, dmr.list[[k]]$ID))
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
            ids1 <- cbind(list.index = i, dmr.index = match(dmrs1$ID, dmr.list[[i]]$ID))
            ids2 <- cbind(list.index = j, dmr.index = match(dmrs2$ID, dmr.list[[j]]$ID))
            ids3 <- cbind(list.index = k, dmr.index = match(dmrs3$ID, dmr.list[[k]]$ID))
            ids4 <- cbind(list.index = l, dmr.index = match(dmrs4$ID, dmr.list[[l]]$ID))
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
              ids1 <- cbind(list.index = i, dmr.index = match(dmrs1$ID, dmr.list[[i]]$ID))
              ids2 <- cbind(list.index = j, dmr.index = match(dmrs2$ID, dmr.list[[j]]$ID))
              ids3 <- cbind(list.index = k, dmr.index = match(dmrs3$ID, dmr.list[[k]]$ID))
              ids4 <- cbind(list.index = l, dmr.index = match(dmrs4$ID, dmr.list[[l]]$ID))
              ids5 <- cbind(list.index = m, dmr.index = match(dmrs5$ID, dmr.list[[m]]$ID))
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
    draw.pairwise.venn(area1 = area1, area2 = area2, cross.area = n12, 
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
      n123 <- round(nrow(tM) / 3, digits=1)
    } else {
      n12 <- max(nrow(pairs[[1]]), 0) / 2
      n13 <- max(nrow(pairs[[2]]), 0) / 2
      n23 <- max(nrow(pairs[[3]]), 0) / 2
      n123 <- 0
    }
    area1 <- n1 + n12 + n13 + n123
    area2 <- n2 + n12 + n23 + n123
    area3 <- n3 + n13 + n23 + n123
    plot.new()
    draw.triple.venn(area1 = area1, area2 = area2, area3 = area3, n12 = n12 + n123, 
                     n13 = n13 + n123, n23 = n23 + n123, n123 = n123, category = group.names, 
                     col = rep("black", 3), fill = c("skyblue", "pink1", "orange"), ...)
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
      n123 <- round(nrow(removeDuplicates(triples[[1]], qM)) / 3, digits = 1)
      n124 <- round(nrow(removeDuplicates(triples[[2]], qM)) / 3, digits = 1)
      n134 <- round(nrow(removeDuplicates(triples[[3]], qM)) / 3, digits = 1)
      n234 <- round(nrow(removeDuplicates(triples[[4]], qM)) / 3, digits = 1)
      n1234 <- round(nrow(qM) / 4, digits = 1)
    } else {
      n123 <- round(max(nrow(triples[[1]]), 0) / 3, digits = 1)
      n124 <- round(max(nrow(triples[[2]]), 0) / 3, digits = 1)
      n134 <- round(max(nrow(triples[[3]]), 0) / 3, digits = 1)
      n234 <- round(max(nrow(triples[[4]]), 0) / 3, digits = 1)
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
                   n234 = n234 + n1234, n1234 = n1234, category = group.names, 
                   col = rep("black", 4), fill = c("skyblue", "pink1", "mediumorchid", "orange"),
                   ...)
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
      n123 <- round(nrow(removeDuplicates(triples[[1]], qM)) / 3, digits = 1)
      n124 <- round(nrow(removeDuplicates(triples[[2]], qM)) / 3, digits = 1)
      n125 <- round(nrow(removeDuplicates(triples[[3]], qM)) / 3, digits = 1)
      n134 <- round(nrow(removeDuplicates(triples[[4]], qM)) / 3, digits = 1)
      n135 <- round(nrow(removeDuplicates(triples[[5]], qM)) / 3, digits = 1)
      n145 <- round(nrow(removeDuplicates(triples[[6]], qM)) / 3, digits = 1)
      n234 <- round(nrow(removeDuplicates(triples[[7]], qM)) / 3, digits = 1)
      n235 <- round(nrow(removeDuplicates(triples[[8]], qM)) / 3, digits = 1)
      n245 <- round(nrow(removeDuplicates(triples[[9]], qM)) / 3, digits = 1)
      n345 <- round(nrow(removeDuplicates(triples[[10]], qM)) / 3, digits = 1)
      
    } else {
      n123 <- round(max(nrow(triples[[1]]), 0) / 3, digits = 1)
      n124 <- round(max(nrow(triples[[2]]), 0) / 3, digits = 1)
      n125 <- round(max(nrow(triples[[3]]), 0) / 3, digits = 1)
      n134 <- round(max(nrow(triples[[4]]), 0) / 3, digits = 1)
      n135 <- round(max(nrow(triples[[5]]), 0) / 3, digits = 1)
      n145 <- round(max(nrow(triples[[6]]), 0) / 3, digits = 1)
      n234 <- round(max(nrow(triples[[7]]), 0) / 3, digits = 1)
      n235 <- round(max(nrow(triples[[8]]), 0) / 3, digits = 1)
      n245 <- round(max(nrow(triples[[9]]), 0) / 3, digits = 1)
      n345 <- round(max(nrow(triples[[10]]), 0) / 3, digits = 1)
    }
    qnM <- do.call(rbind, quins)
    if (ncol(qnM) > 1) {
      qnM <- qnM[which(!is.na(qnM[ ,1])), ]  # Get rid of NA rows
      n1234 <- round(nrow(removeDuplicates(quads[[1]], qnM)) / 4, digits = 1)
      n1235 <- round(nrow(removeDuplicates(quads[[2]], qnM)) / 4, digits = 1)
      n1245 <- round(nrow(removeDuplicates(quads[[3]], qnM)) / 4, digits = 1)
      n1345 <- round(nrow(removeDuplicates(quads[[4]], qnM)) / 4, digits = 1)
      n2345 <- round(nrow(removeDuplicates(quads[[5]], qnM)) / 4, digits = 1)
      n12345 <- round(nrow(qnM) / 5, digits = 1)
    } else {
      n1234 <- round(max(nrow(quads[[1]]), 0) / 4, digits = 1)
      n1235 <- round(max(nrow(quads[[2]]), 0) / 4, digits = 1)
      n1245 <- round(max(nrow(quads[[3]]), 0) / 4, digits = 1)
      n1345 <- round(max(nrow(quads[[4]]), 0) / 4, digits = 1)
      n2345 <- round(max(nrow(quads[[5]]), 0) / 4, digits = 1)
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
                        category = group.names,  
                        fill = c("skyblue", "pink1", "mediumorchid", "orange", "green"), ...)
  }
  
}
