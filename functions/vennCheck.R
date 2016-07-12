## Created 6/16/2016 by Daniel Beck
## Last modified 6/16/2016

## This function generates several Venn diagrams as a test for the underlying functions



## The resulting Venn diagrams should have areas of gradually increasing volume
dmrTableList <- list()
dmrTableList[[1]] <- data.frame("chr" = rep(1, 10),
                                "start" = c(1, 11, 21, 31, 41, 51, 61, 71, 81, 91),
                                "stop" = c(1, 11, 21, 31, 41, 51, 61, 71, 81, 91) + 1,
                                "ID" = paste("dmr", 1:10, sep="")
                               )
dmrTableList[[2]] <- data.frame("chr" = rep(1, 10),
                                "start" = c(1, 11, 21, 31, 41, 51, 61, 71, 81, 91) + 2,
                                "stop" = c(1, 11, 21, 31, 41, 51, 61, 71, 81, 91) + 3,
                                "ID" = paste("dmr", 1:10, sep="")
)
dmrTableList[[3]] <- data.frame("chr" = rep(1, 10),
                                "start" = c(1, 11, 21, 31, 41, 51, 61, 71, 81, 91) + 4,
                                "stop" = c(1, 11, 21, 31, 41, 51, 61, 71, 81, 91) + 5,
                                "ID" = paste("dmr", 1:10, sep="")
)
dmrTableList[[4]] <- data.frame("chr" = rep(1, 10),
                                "start" = c(1, 11, 21, 31, 41, 51, 61, 71, 81, 91) + 6,
                                "stop" = c(1, 11, 21, 31, 41, 51, 61, 71, 81, 91) + 7,
                                "ID" = paste("dmr", 1:10, sep="")
)
dmrTableList[[5]] <- data.frame("chr" = rep(1, 10),
                                "start" = c(1, 11, 21, 31, 41, 51, 61, 71, 81, 91) + 8,
                                "stop" = c(1, 11, 21, 31, 41, 51, 61, 71, 81, 91) + 9,
                                "ID" = paste("dmr", 1:10, sep="")
)
dmrTableList[[6]] <- data.frame("chr" = rep(1, 10),
                                "start" = c(1, 11, 21+2, 31+2, 41+4, 51+4, 61+6, 71+6, 81+8, 91+8),
                                "stop" = c(1+1, 11+1, 21+3, 31+3, 41+5, 51+5, 61+7, 71+7, 81+9, 91+9),
                                "ID" = paste("dmr", 1:10, sep="")
)

## Five groups: all unique 10
vennDMR(dmr.list = dmrTableList[c(1:5)], group.names = 1:5)
## Five groups: perfectly overlapping 10
vennDMR(dmr.list = dmrTableList[c(1,1,1,1,1)], group.names = 1:5)
## Five groups: Perfect overlap between groups 1 and 2, 3 and 4, all unique for group 5
vennDMR(dmr.list = dmrTableList[c(1,1,2,2,3)], group.names = 1:5)
## Five groups: Two overlaps with group 5 and other sets
vennDMR(dmr.list = dmrTableList[c(1,2,3,4,6)], group.names = 1:5)
## Five groups: Triple overlaps
vennDMR(dmr.list = dmrTableList[c(1,2,3,6,6)], group.names = 1:5)
## Five groups: Quadruple overlaps
vennDMR(dmr.list = dmrTableList[c(1,2,6,6,6)], group.names = 1:5)
## Five groups: Two complete overlaps
vennDMR(dmr.list = dmrTableList[c(1,6,6,6,6)], group.names = 1:5)

## Four groups: all unique 10
vennDMR(dmr.list = dmrTableList[c(1:4)], group.names = 1:4)
## Four groups: perfectly overlapping 10
vennDMR(dmr.list = dmrTableList[c(1,1,1,1)], group.names = 1:4)
## Four groups: Perfect overlap between groups 1 and 2, 3 and 4
vennDMR(dmr.list = dmrTableList[c(1,1,2,2)], group.names = 1:4)
## Four groups: Two double overlaps with group 5 and 4
vennDMR(dmr.list = dmrTableList[c(1,2,3,6)], group.names = 1:4)
## Four groups: Two triple overlaps
vennDMR(dmr.list = dmrTableList[c(1,2,6,6)], group.names = 1:4)
## Four groups: Two complete overlaps
vennDMR(dmr.list = dmrTableList[c(1,6,6,6)], group.names = 1:4)

vennDMR(dmr.list = dmrTableList[c(1,2,3)], group.names = 1:3, euler.d=F, scaled=F)
vennDMR(dmr.list = dmrTableList[c(1,1,1)], group.names = 1:3, euler.d=F, scaled=F)
vennDMR(dmr.list = dmrTableList[c(1,2,6)], group.names = 1:3, euler.d=F, scaled=F)
vennDMR(dmr.list = dmrTableList[c(1,6,6)], group.names = 1:3, euler.d=F, scaled=F)




