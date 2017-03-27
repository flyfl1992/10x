#####################
## Attempting to add GO terms to 10X genomics PBMC data
#####################

# From non-cran sources
library(cellrangerRkit) # https://support.10xgenomics.com/single-cell/software/pipelines/latest/rkit
library(Seurat) # http://satijalab.org/seurat/

# From cran
library(tidyverse)
library(Matrix)

##############################
# Reading data into CellRangerR:
##############################

#specifying file path
pipestance_path <- getwd()

## Download data (if you haven't already):
#   (Note that this is not necessary if you have cloned the repo)
# download_sample(sample_name="pbmc3k",sample_dir=pipestance_path, 
#                 host="https://s3-us-west-2.amazonaws.com/10x.files/samples/cell/") 

#The following code loads data into CellRangerR:
gbm <- load_cellranger_matrix(pipestance_path) #gbm has 2700 cells and 32,738 genes

# #Accessing data within gbm
# exprs(gbm) # expression matrix
# fData(gbm) # data frame of genes
# pData(gbm) # data frame of cell barcodes

###########################
# reading data into Seurat:
###########################

#relative filepath for matrices
filepath <- "./outs/filtered_gene_bc_matrices/hg19/"
pbmc.data <- Read10X(filepath)

#In Seurat, the following command is used to read raw data into S4 class
#pbmc <- new("seurat", raw.data = pbmc.data)

#################
# Loading GO Terms
#################

#import list of gene IDs from 10x data (Note that this is not necessary if you have cloned the repo)
# gene_names2 <- read_delim("genes.tsv", delim = "\t",col_names = c("Ensembl","ID"))
# 
# write.csv(gene_names2[,1],file = "genewrite",row.names=FALSE) #write Ensembl IDs to file

#use Retrieve/ID mapping to convert from Ensembl to UniprotKB
# http://www.uniprot.org/uploadlists/ # paste in list of genes, or upload (may have to remove quotation marks first)
#  filter by "Reviewed / Swiss-Prot"
#   download results as a tab separated file
#    rename downloaded file to "uniprot.tab"

uniprot<-read_delim("uniprot.tab", delim = "\t")
summary(uniprot)
#rename first two columns
names(uniprot)[1:2] = c("id","isomap")
head(uniprot)
#remove duplicated Ensembl IDs
#(note: coarse way to do this, should investigate why some IDs are duplicated)
uniprot[! duplicated(uniprot$id), ] -> uniprot.noDup
head(uniprot.noDup)

##################
## Relating GO terms to PBMC Data
##################

#relational database from the genes database to the uniprot.noDup
gene_annotations <- left_join( fData(gbm), uniprot.noDup, by = "id") %>% as_tibble()
#the problem is that the CellrangerR package doesn't support addition of metadata
#  need to instead use an non-sparse matrix 

#I think that as.matrix parses the Seurat pbmc.data
denseMatrix <-as.matrix(pbmc.data)
#however,
names(denseMatrix)
# NULL

#instead tried to join this database to the raw reads parsed by Seurat
pbmc.tibble <- as_tibble(denseMatrix)
join_reads <- left_join(gene_annotations, pbmc.tibble, by = NULL) #throwing error, even though there should be a shared Gene ID column
#using as.matrix(pbmc.data) seems to elminate the gene IDs present in pbmc.data

head(pbmc.tibble)
head(pbmc.data)
head(gene_annotations)
head(denseMatrix)

#In Seurat there is a function to add metadata, perhaps could use that
#AddMetaData adds columns to object@data.info, and is a great place to stash QC stats
#mito.genes <- grep("^MT-", rownames(pbmc@data), value = T)
#percent.mito <- colSums(expm1(pbmc@data[mito.genes, ]))/colSums(expm1(pbmc@data))
pbmc <- new("seurat", raw.data = pbmc.data)
#pbmc <- AddMetaData(pbmc, percent.mito, "percent.mito")

#can potentially use the AddMetaData function from Seurat to add columns
pbmc.join <- AddMetaData(pbmc, gene_annotations)
#no errors, but
names(pbmc.join)
#NULL




