## Daniel Beck
## Created 5/28/2015
## Modified 
##	5/29/2015 Finished writing functions. Sucessful test using steelhead RBC data.
##	6/4/2015  Fixed error caused by numbers being stored as characters
##   7/6/2015  Added CpG density plot function. Changed name of file to customFunctions.R. Added addToMergedResults function to keep track of all DMR information.
##   7/7/2015  Fixed addToMergedResults function that failed in the case of empty or single line dataframes
##   7/13/2015 Added file to git repository "medipPipeline"
##   7/14/2015 Fixed errors that killed code when passed lists to plots were empty. Now returns().
##   8/13/2015 I'm splitting all functions into separate files. This script will stay here and act as a hub for sourcing all the other functions

# This file contains code with additional functions. These functions include those necessary for the chromosome plot and the CpG plot.

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
## Functions to identify and extract DMRs
source("./functions/extractDMR.R")
# Function that extracts mapping % from prepareData.Rout
source("./functions/mapExtract.R")
# Modifies dmrList with correct column names for slidingWindowCluster function.
source("./functions/formatForCluster.R")
# Excludes selected chromosomes from dmrList
source("./functions/excludeChr.R")
