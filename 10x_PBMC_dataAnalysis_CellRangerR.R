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


###
# Next step is to re-run sequencing analysis using subsetted data...


# #examining expression data over cells for a certain gene
# exprs(gbm[6,400:500])
# fData(gbm[1,1])
# 


