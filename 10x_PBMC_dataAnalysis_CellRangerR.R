#10X genomics CellRangerR tutorial here:
#https://support.10xgenomics.com/single-cell/software/pipelines/latest/rkit

source("http://s3-us-west-2.amazonaws.com/10x.files/code/rkit-install-1.1.0.R")

library(cellrangerRkit)

setwd("~/Documents/Dunn_Lab/10x/")

packageVersion("cellrangerRkit")


library(devtools)

pipestance_path <- "/Users/Zack/Documents/Dunn_Lab/10x"
setwd(pipestance_path)
download_sample(sample_name="pbmc3k",sample_dir=pipestance_path, 
                host="https://s3-us-west-2.amazonaws.com/10x.files/samples/cell/") 
gbm <- load_cellranger_matrix(pipestance_path) 
analysis_results <- load_cellranger_analysis_results(pipestance_path)

exprs(gbm) # expression matrix
fData(gbm) # data frame of genes
pData(gbm) # data frame of cell barcodes


#t-SNE projection and plot the cells colored by UMI counts
# generates a t-SNE projection where each cell is colored by log10 of UMI counts.
#  Each point in the scatter plot represents a cell in the coordinates specified by the two t-SNE components. 
#  The color of each point plotted by visualize_umi_counts (Figure 1) indicates the total number of UMIs for each cell,
#  and these count values are displayed in log10 scale.
tsne_proj <- analysis_results$tsne 
visualize_umi_counts(gbm,tsne_proj[c("TSNE.1","TSNE.2")], limits=c(3,4),marker_size=0.05)

#analyzing gene expression signatures:
# Instead of using raw UMI counts for downstream differential gene analysis, 
# we recommend that you filter unexpressed genes, normalize the UMI counts 
# for each barcode, and use the log-transformed gene-barcode matrix
use_genes <- get_nonzero_genes(gbm)
gbm_bcnorm <- normalize_barcode_sums_to_median(gbm[use_genes,])
gbm_log <- log_gene_bc_matrix(gbm_bcnorm,base=10)
print(dim(gbm_log))

#Fishing with known gene markers using visualize_gene_markers
genes <- c("CD79A","NKG7","CD3D","CST3","CD8A","PF4")
tsne_proj <- analysis_results$tsne 
visualize_gene_markers(gbm_log,genes,tsne_proj[c("TSNE.1","TSNE.2")],limits=c(0,1.5))

#trying with PCA
genes <- c("CD79A","NKG7","CD3D","CST3","CD8A","PF4")
pca_proj <- analysis_results$pca 
visualize_gene_markers(gbm_log,genes,pca_proj[c("PC.1","PC.2")],limits=c(0,1.5))


#visualizing clusters at different kmeans clustering levels:

n_clu <- 2:10
km_res <- analysis_results$kmeans # load pre-computed kmeans results
clu_res <- sapply(n_clu, function(x) km_res[[paste(x,"clusters",sep="_")]]$Cluster) 
colnames(clu_res) <- sapply(n_clu, function(x) paste("kmeans",x,sep=".")) 
visualize_clusters(clu_res,tsne_proj[c("TSNE.1","TSNE.2")])

#analyzing cluster-specific genes
# DE between cell types can then be examined
example_K <- 5 # number of clusters (use "Set3" for brewer.pal below if example_K > 8) 
example_col <- rev(brewer.pal(example_K,"Set2")) # customize plotting colors
cluster_result <- analysis_results$kmeans[[paste(example_K,"clusters",sep="_")]] 
visualize_clusters(cluster_result$Cluster,tsne_proj[c("TSNE.1","TSNE.2")],colour=example_col)

#for example, comparing the mean expression between one cell class and all other classes
# prioritize_top_genes() will output genes upregulated in a particular cluster

# sort the cells by the cluster labels
cells_to_plot <- order_cell_by_clusters(gbm, cluster_result$Cluster)
# order the genes from most up-regulated to most down-regulated in each cluster 
prioritized_genes <- prioritize_top_genes(gbm, cluster_result$Cluster, "sseq", min_mean=0.5)

#now output the top genes for each cluster to a folder
output_folder <-pipestance_path 
write_cluster_specific_genes(prioritized_genes, output_folder, n_genes=10)

#making a heatmap
# create values and axis annotations for pheatmap
# genes are rows, cells are columns. Top 3 DE genes for each cluster are displayed
gbm_pheatmap(log_gene_bc_matrix(gbm), prioritized_genes, cells_to_plot, n_genes=5, colour=example_col, limits=c(-1,2))


#from your bait genes, make calls as to what cell type each cluster represents 
# and output composition of entire population
cell_composition(cluster_result$Cluster,
                 anno=c("monocytes","T cells","NK cells","megakaryocytes","B cells"))



#Subsetting data
# examples of how to subset:
# subset_by_cell <- gbm[,c("AAACATACAACCAC-1", "AAACATTGAGCTAC-1")]
# subset_by_gene_id <- gbm["ENSG00000167286",]
# subset_by_gene_symbol <- gbm[which(fData(gbm)$symbol == 'CD3D'),]
# subset_by_cell_and_gene <- gbm["ENSG00000167286", c("AAACATACAACCAC-1", "AAACATTGAGCTAC-1")]


#Subset by GO_Term
# Need to make sure that the ENSEMBL ids from GO_Search are actually in the gbm matrix
# `transcription_factor_vector` represents the results of the `GO_Search.R` results
transcription_factor_vector2 = transcription_factor_vector[ transcription_factor_vector %in% fData(gbm)$id ]
#then subset by these GO terms
subset_by_GO_term <-gbm[c(transcription_factor_vector2),]






#examining expression data over cells for a certain gene
exprs(gbm[6,400:500])
fData(gbm[1,1])

#################
# Loading GO Terms
#################

library(tidyverse)
#import list of gene IDs from 10x data
gene_names2 <- read_delim("genes.tsv", delim = "\t",col_names = c("Ensembl","ID"))

write.csv(gene_names2[,1],file = "genewrite",row.names=FALSE) #write Ensembl IDs to file

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
#(note: course way to do this, should investigate why some IDs are duplicated)
uniprot[! duplicated(uniprot$id), ] -> uniprot.noDup
head(uniprot.noDup)

#relational database from the genes database to the uniprot.noDup
join<- left_join(fData(gbm),uniprot.noDup, by = NULL)

#what now? Can I output this as a tsv, then import with pipestance path?
write_tsv(x = join, path ="joined.tsv")
# no, the load_cellranger_matrix_from_files() doesn't like the new genes.tsv file
# error message claims that its dimensions don't match the matrix file, although ...

dim(exprs(gbm))
# [1] 32738  2700
dim(fData(gbm))
# [1] 32738     2
dim(join)
# [1] 32738    15








###############
# Scratch below
###############


#to examine the code for loading in the databases:
load_cellranger_matrix_from_files #no brackets

readMM # appears to be function for reading matrix

#perhaps could use the Seurat Read10X() function to read in matrix.

#example coerce sparse matrix to data.frame, not exactly correct

summ <- summary(pbmc.data)
pbmc.dataframe <- data.frame(Origin      = rownames(mat)[summ$i],
                             Destination = colnames(mat)[summ$j],
                             Weight      = summ$x)

#I think that as.matrix parses it okay. pbmc.data from Seurat.

denseMatrix <-as.matrix(pbmc.data)
dim(denseMatrix)

denseGBM <- as.data.frame(fData(gbm))

############
# general code for elbow plots using sum of squared error plot

wss <- (nrow(mydata)-1)*sum(apply(mydata,2,var))
for (i in 2:15) wss[i] <- sum(kmeans(mydata,
                                     centers=i)$withinss)
plot(1:15, wss, type="b", xlab="Number of Clusters",
     ylab="Within groups sum of squares") 







