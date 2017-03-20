###########
# Using Seurat to analyze PBMC data from 10x
###########
#from this link:
# http://satijalab.org/seurat/pbmc-tutorial.html

library(devtools)
#for Seurat install: 
#install_url("https://github.com/satijalab/seurat/releases/download/v1.4.0/Seurat_1.4.0.9.tgz", binary = TRUE)

library(Seurat)
library(tidyverse)
library(Matrix)

# Load the PBMC dataset
pbmc.data <- Read10X("~/Documents/Dunn_Lab/10x/outs/filtered_gene_bc_matrices/hg19/")

#Examine the memory savings between regular and sparse matrices
dense.size <- object.size(as.matrix(pbmc.data))
dense.size

sparse.size <- object.size(pbmc.data)
sparse.size
dense.size/sparse.size

# Initialize the Seurat object with the raw (non-normalized data)
# Note that this is slightly different than the older Seurat workflow, where log-normalized values were passed in directly.
# You can continue to pass in log-normalized values, just set do.logNormalize=F in the next step.
pbmc <- new("seurat", raw.data = pbmc.data)

# Keep all genes expressed in >= 3 cells, keep all cells with >= 200 genes
# Perform log-normalization, first scaling each cell to a total of 1e4 molecules (as in Macosko et al. Cell 2015)
pbmc <- Setup(pbmc, min.cells = 3, min.genes = 200, do.logNormalize = T, total.expr = 1e4, project = "10X_PBMC")

# Basic QC and selecting cells for further analysis
# 
# While the setup function imposes a basic minimum gene-cutoff, 
#you may want to filter out cells at this stage based on technical or biological parameters.
#Seurat allows you to easily explore QC metrics and filter cells based on any user-defined criteria.
#In the example below, we visualize gene and molecule counts, plot their relationship, 
#and exclude cells with a clear outlier number of genes detected as potential multiplets. 
#Of course this is not a guarenteed method to exclude cell doublets,
#but we include this as an example of filtering user-defined outlier cells.
#We also filter cells based on the percentage of mitochondrial genes present.

#nGene and nUMI are automatically calculated for every object by Seurat. For non-UMI data, nUMI represents the sum of the non-normalized values within a cell
# We calculate the percentage of mitochondrial genes here and store it in percent.mito using the AddMetaData. The % of UMI mapping to MT-genes is a common scRNA-seq QC metric.
mito.genes <- grep("^MT-", rownames(pbmc@data), value = T)
percent.mito <- colSums(expm1(pbmc@data[mito.genes, ]))/colSums(expm1(pbmc@data))

#AddMetaData adds columns to object@data.info, and is a great place to stash QC stats
pbmc <- AddMetaData(pbmc, percent.mito, "percent.mito")
VlnPlot(pbmc, c("nGene", "nUMI", "percent.mito"), nCol = 3)

#GenePlot is typically used to visualize gene-gene relationships, but can be used for anything calculated by the object, i.e. columns in object@data.info, PC scores etc.
#Since there is a rare subset of cells with an outlier level of high mitochondrial percentage, and also low UMI content, we filter these as well
par(mfrow = c(1, 2))
GenePlot(pbmc, "nUMI", "percent.mito")
GenePlot(pbmc, "nUMI", "nGene")

#We filter out cells that have unique gene counts over 2,500
#Note that accept.high and accept.low can be used to define a 'gate', and can filter cells not only based on nGene but on anything in the object (as in GenePlot above)
pbmc <- SubsetData(pbmc, subset.name = "nGene", accept.high = 2500)
pbmc <- SubsetData(pbmc, subset.name = "percent.mito", accept.high = 0.05)

#
###regress out unwanted sources of variation

# Seurat implements a basic version of this by constructing linear models to predict gene expression 
# based on user-defined variables. Seurat stores the z-scored residuals of these models in the scale.data slot, 
# and they are used for dimensionality reduction and clustering.

#We typically regress out cell-cell variation in gene expression driven by batch (if applicable), 
#cell alignment rate (as provided by Drop-seq tools for Drop-seq data), 
#the number of detected molecules, and mitochondrial gene expression. 

#note that this overwrites pbmc@scale.data. Therefore, if you intend to use RegressOut, you can set do.scale=F and do.center=F in the original object to save some time.
pbmc <- RegressOut(pbmc, latent.vars = c("nUMI", "percent.mito"))

par(mfrow =c(1,1))
#detection of variable genes
#MeanVarPlot(), which works by calculating the average expression and dispersion for each gene, placing these genes into bins, and then calculating a z-score for dispersion within each bin. This helps control for the relationship between variability and average expression.
pbmc <- MeanVarPlot(pbmc ,fxn.x = expMean, fxn.y = logVarDivMean, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5, do.contour = F)


