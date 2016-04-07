# medipPipeline
#### A collection of scripts for MeDIP-seq data analysis

This repository currently contains a collection of R scripts and functions for analyzing MeDIP-seq datasets. Eventually these scripts will form a pipeline that can be used on generic datasets.

## Availability

medipPipeline is available at https://github.com/danlbek/medipPipeline.

## Overview

**dataNames.R:** Configuration file for the pipeline. It contains filenames, file locations, constants, and other adjustable parameters.

**prepareData.R:** Cleans, filteres (Trimmomatic), and maps (Bowtie2) reads as well as generates read quality reports (FastQC).

**medipAnalysis.R:** Performs the basic MeDIP-seq analysis using the MEDIPS library.

**medipProcessing.R:** Identifies DMR and generates DMR tables.

**annotateDMRs.R:** Adds annotation information to DMR tables.

**customFunctions.R:** Wrapper for loading all custom functions used by the pipeline. These functions are stored in the functions folder.

**generateReports.R:** Wrapper for generating reports summarizing the MeDIP-seq results. This script renders the report using the **mediepReport.Rmd** code.

**publicationFiles.R:** Produces separate files for selected tables and plots. 

