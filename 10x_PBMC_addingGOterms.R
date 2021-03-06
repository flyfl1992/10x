#####################
## Attempting to add GO terms to 10X genomics PBMC data
#####################

# From non-cran sources
library(cellrangerRkit) # https://support.10xgenomics.com/single-cell/software/pipelines/latest/rkit
library(Seurat) # http://satijalab.org/seurat/

# From cran
library(tidyverse)
library(Matrix)
library(stringr)

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

################
# Searching for a GO term and its ancestors
################

go_map <- read_tsv("go_map.tsv")

#go_map[go_map$go_id == "GO:1903506", ]

#transcription factor activity, protein binding
# can be accessed as follows:
# go_map[go_map$go_id == "GO:0000988", ]


#' Returns ENSEMBL IDs from a GO term search
#' 
#' @param Character indicating a go term, eg "GO:0003712"
#' @return A vector of ENSEMBL IDs
#' @export
id_search = function( go_term ){
  
  row_return = go_map %>%
    filter(str_detect(go_id, go_term)) 
  
  row_return = row_return[! duplicated(row_return$ensembl_gene_id), ]
  
  id_vector = as.vector(row_return$ensembl_gene_id)
  
  return(id_vector)
}

#store your search results in a variable, 
# such as: 
transcription_factor_vector <- id_search("GO:0000988")

##################
#Subset by GO_Term
##################

# Need to make sure that the ENSEMBL ids from GO_Search are actually in the gbm matrix
# `transcription_factor_vector` represents the results of the id_search() function
transcription_factor_vector2 = transcription_factor_vector[ transcription_factor_vector %in% fData(gbm)$id ]
#then subset by these GO terms
subset_by_GO_term <-gbm[c(transcription_factor_vector2),]


#can use fData, pData and exprs functions on this subset
table(exprs(subset_by_GO_term)[1,]) #summary
table(fData(subset_by_GO_term)[1:20,])

#the expression matrix would be accessed with:
expression_subset <- exprs(subset_by_GO_term)
head(expression_subset)

###
# Next step is to re-run sequencing analysis using subsetted data...

#analyzing gene expression signatures:
# Instead of using raw UMI counts for downstream differential gene analysis, 
# Filter out unexpressed genes, then normalize the UMI counts for the barcodes
# This log-transformed gene-barcode matrix can be used for analysis.
# (Raw UMI counts are not recommended for downstream analysis)

use_genes <- get_nonzero_genes(subset_by_GO_term)
gbm_bcnorm <- normalize_barcode_sums_to_median(subset_by_GO_term[use_genes,])
gbm_log <- log_gene_bc_matrix(gbm_bcnorm,base=10)
print(dim(gbm_log))

# ###########################
# # reading data into Seurat:
# ###########################
# 
# #relative filepath for matrices
# filepath <- "./outs/filtered_gene_bc_matrices/hg19/"
# pbmc.data <- Read10X(filepath)
# 
# #In Seurat, the following command is used to read raw data into S4 class
# #pbmc <- new("seurat", raw.data = pbmc.data)
# 
# #################
# # Loading GO Terms
# #################
# 
# #import list of gene IDs from 10x data (Note that this is not necessary if you have cloned the repo)
# # gene_names2 <- read_delim("genes.tsv", delim = "\t",col_names = c("Ensembl","ID"))
# # 
# # write.csv(gene_names2[,1],file = "genewrite",row.names=FALSE) #write Ensembl IDs to file
# 
# #use Retrieve/ID mapping to convert from Ensembl to UniprotKB
# # http://www.uniprot.org/uploadlists/ # paste in list of genes, or upload (may have to remove quotation marks first)
# #  filter by "Reviewed / Swiss-Prot"
# #   download results as a tab separated file
# #    rename downloaded file to "uniprot.tab"
# 
# uniprot<-read_delim("uniprot.tab", delim = "\t")
# summary(uniprot)
# #rename first two columns
# names(uniprot)[1:2] = c("id","isomap")
# head(uniprot)
# #remove duplicated Ensembl IDs
# #(note: coarse way to do this, should investigate why some IDs are duplicated)
# uniprot[! duplicated(uniprot$id), ] -> uniprot.noDup
# head( uniprot.noDup )
# 
# ##################
# ## Relating GO terms to PBMC Data
# ##################
# 
# #relational database from the genes database to the uniprot.noDup
# gene_annotations <- 
# 	left_join( fData(gbm), uniprot.noDup, by = "id") %>% 
# 	as_tibble()
# 
# #the problem is that the CellrangerR package doesn't support addition of metadata
# #  need to instead use an non-sparse matrix 
# 
# #I think that as.matrix parses the Seurat pbmc.data
# denseMatrix <- as.matrix( pbmc.data )
# #however,
# names( denseMatrix )
# # NULL
# 
# #instead tried to join this database to the raw reads parsed by Seurat
# pbmc.tibble <- as_tibble( denseMatrix )
# join_reads <- left_join( gene_annotations, pbmc.tibble, by = NULL ) #throwing error, even though there should be a shared Gene ID column
# #using as.matrix(pbmc.data) seems to elminate the gene IDs present in pbmc.data
# 
# head(pbmc.tibble)
# head(pbmc.data)
# head(gene_annotations)
# 
# #In Seurat there is a function to add metadata, perhaps could use that
# #AddMetaData adds columns to object@data.info, and is a great place to stash QC stats
# #mito.genes <- grep("^MT-", rownames(pbmc@data), value = T)
# #percent.mito <- colSums(expm1(pbmc@data[mito.genes, ]))/colSums(expm1(pbmc@data))
# pbmc <- new("seurat", raw.data = pbmc.data)
# #pbmc <- AddMetaData(pbmc, percent.mito, "percent.mito")
# 
# #can potentially use the AddMetaData function from Seurat to add columns
# pbmc.join <- AddMetaData(pbmc, gene_annotations)
# #no errors, but
# names(pbmc.join)
# #NULL




