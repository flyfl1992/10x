# 10x genomics 2,700 PBMC GO term exploration

This repository...

## Contents

`10x_PBMC_dataAnalysis_Seurat.R` and `10x_PBMC_dataAnalysis_CellRangerR.R` reflect example scripts from tutorials for Seurat and CellRangerR packages, respectively. 

## Getting set up

### Installing dependencies

Install the following dependencies:

- `CellrangerR`, follow the installation instructions at https://support.10xgenomics.com/single-cell/software/pipelines/latest/rkit . Tested with version 1.1.0.
- `seurat`, run the following in R: `devtools::install_url("https://github.com/satijalab/seurat/releases/download/v1.4.0/Seurat_1.4.0.9.tgz", binary = TRUE)`

### Obtaining the data

All data for the analyses are in the repo.

Data were obtained from https://support.10xgenomics.com/single-cell/software/pipelines/latest/rkit

### Obtaining the respository

Download or clone the repsoitory from https://github.com/zrlewis/10x


## Running the analyses


Each time you want to execute an analysis:

- `cd` to the root of the repository, eg `cd 10x`
- Run the code, eg lanch `R` and run 
  `source("10x_PBMC_addingGOterms.R")`

