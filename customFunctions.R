## Created 5/28/2015 by Daniel Beck
## Last modified 4/7/2016

## This script acts as an index file for all necessary functions. This is not a standard
## way of doing things, however. It may be best to create an R package in the future that
## contains the MEDIP pipeline. It could then be called with library().

## Keeps track of p-values and other information when merging DMR
source("./functions/addToMergedResults.R")
## Adds annotation to DMR table from GFF source
source("./functions/addAnnotationGFF.R")
## Adds annotation to DMR table using remote Biomart database
source("./functions/addAnnotationBiomart.R")
## Adds annotation to DMR table using BLAST
source("./functions/addAnnotationBlast.R")
## Extracts reference sequence corresponding to DMR
source("./functions/dmr2seq.R")
## Finds a custom category for a gene using homologs
source("./functions/identifyCategory.R")
## Matches gene to annotationTable using homologs
source("./functions/customAnnotation.R")
## Takes in two tables and merges them using a specified column.
source("./functions/addToTable.R")
## Haque's function to perform sliding window cluster detection
source("./functions/slidingWindowCluster.R")
## Functions used to generate chromosome plot
source("./functions/plotChromosomes.R")
## Generates CpG density histogram
source("./functions/plotCpGdensity.R")
## Calculates CpG density
source("./functions/calcCpGdensity.R")
## Modifies incorrectly labeled stop site in DMR lists
source("./functions/modifyStop.R")
## Identifies and returns DMRs
source("./functions/extractDMR.R")
## Extracts mapping % from prepareData log file
source("./functions/mapExtract.R")
## Modifies dmrList with correct column names for slidingWindowCluster function.
source("./functions/formatForCluster.R")
## Excludes selected chromosomes from the DMR table
source("./functions/excludeChr.R")
## Counts the number of DMR on a set of chromosomes
source("./functions/countDMR.R")
## Generates DMR length histogram
source("./functions/plotDMRlength.R")
## Generates a histogram with a maximum value
source("./functions/histMax.R")
## Generates a Venn diagram from a list of DMR tables
source("./functions/vennDMR.R")
## Returns list of base pairs that overlap between all input DMR tables
source("./functions/overlappingBP.R")
## Takes start and stop positions and returns overlaps
source("./functions/pairOverlap.R")
## Returns DMR that include specified base pairs
source("./functions/bpToDmr.R")
## Removes lines in matrix 1 that are also in matrix 2
source("./functions/removeDuplicates.R")
## These functions modify the VennDiagram package functions to allow easier integration
## with the medipPipeline code.
source("./functions/draw.quintuple.venn.dbmod.R")
source("./functions/draw.quad.venn.dbmod.R")
source("./functions/draw.triple.venn.dbmod.R")
## Use a sliding window to count things in genomic windows
source("./functions/swCount.R")
## Extracts all windows that meet CpG number thresholds
source("./functions/identifyZero.R")

## These are no longer used. I've kept them in this file for documentation and on the off 
## chance they become useful.
# source("./functions/oldFunctions.R")




