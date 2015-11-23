## Daniel Beck
## Created 5/28/2015
## Modified 10/8/2015

# This file contains code with additional functions. These functions include 
# those necessary for the chromosome plot and the CpG plot.

## Function to perform sliding window cluster detection
source("./functions/slidingWindowCluster.R")
## Functions used to generate chromosome plot
source("./functions/chromosomePlot.R")
## Function used to generate CpG density histogram
source("./functions/plotCpGdensity.R")
## Function to keep track of p-values and other information when merging DMR
source("./functions/addToMergedResults.R")
## Script to hold small generic functions
source("./functions/miscFunctions.R")
## Function to calculate CpG density
source("./functions/calcCpGdensity.R")
## Functions to add annotation results to DMR lists
source("./functions/annotation.R")
## Function to modify incorrectly labeled stop site in DMR lists
source("./functions/modifyStop.R")
## Venn diagram functions
source("./functions/vennFunctions.R")
source("./functions/overlap.R")
## Functions to identify and extract DMRs
source("./functions/extractDMR.R")
# Function that extracts mapping % from prepareData.Rout
source("./functions/mapExtract.R")
# Modifies dmrList with correct column names for slidingWindowCluster function.
source("./functions/formatForCluster.R")
# Excludes selected chromosomes from dmrList
source("./functions/excludeChr.R")
# Plot DMR length
source("./functions/plotDMRlength.R")

