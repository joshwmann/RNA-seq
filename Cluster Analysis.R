install.packages("BiocManager")
BiocManager::install()

#RUN THROUGH CLUSTERING
#Testing
#two waysto install apckages in R for bioinformatics

#1 Load Data
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("STexampleData")

library(SpatialExperiment)
library(STexampleData)
# load object
spe <- Visium_humanDLPFC()
######OUTPUT NOT THE SAME AS ONLINE#########
# save object(s)
saveRDS(spe, file = "spe_load.rds")

# check object
spe
# number of genes (rows) and spots (columns)
dim(spe)
# names of 'assays'
assayNames(spe)
# row (gene) data
head(rowData(spe))
# column (spot) data
head(colData(spe))
# spatial coordinates
head(spatialCoords(spe))
# image data
imgData(spe)

# create data
n_genes <- 200
n_spots <- 100
counts <- matrix(0, nrow = n_genes, ncol = n_spots)
row_data <- DataFrame(
  gene_name = paste0("gene", sprintf("%03d", seq_len(n_genes)))
)
col_data <- DataFrame(
  sample_id = rep("sample01", n_spots)
)
spatial_coords <- matrix(0, nrow = n_spots, ncol = 2)
colnames(spatial_coords) <- c("x", "y")
# create SpatialExperiment object
spe <- SpatialExperiment(
  assays = list(counts = counts), 
  colData = col_data, 
  rowData = row_data, 
  spatialCoords = spatial_coords
)



#QUALITY CONTROL
#Load previously saved data
library(SpatialExperiment)
spe <- readRDS("spe_load.rds")

#Plot data
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ggspavis")
library(ggspavis)
# plot spatial coordinates (spots)
plotSpots(spe)

#Calculate QC Metrics
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("scater")
library(scater)
# subset to keep only spots over tissue
spe <- spe[, colData(spe)$in_tissue == 1]
dim(spe)
# identify mitochondrial genes
is_mito <- grepl("(^MT-)|(^mt-)", rowData(spe)$gene_name)
table(is_mito)
rowData(spe)$gene_name[is_mito]
# calculate per-spot QC metrics and store in colData
spe <- addPerCellQC(spe, subsets = list(mito = is_mito))
head(colData(spe))
# histogram of library sizes
hist(colData(spe)$sum, breaks = 20)
# plot library size vs. number of cells per spot
plotSpotQC(spe, plot_type = "scatter", 
           x_metric = "cell_count", y_metric = "sum", 
           y_threshold = 600)
# select QC threshold for library size
qc_lib_size <- colData(spe)$sum < 600
table(qc_lib_size)
colData(spe)$qc_lib_size <- qc_lib_size
# check spatial pattern of discarded spots
plotSpotQC(spe, plot_type = "spot", 
           annotate = "qc_lib_size")
# check spatial pattern of discarded spots if threshold is too high
qc_lib_size_2000 <- colData(spe)$sum < 2000
colData(spe)$qc_lib_size_2000 <- qc_lib_size_2000
plotSpotQC(spe, plot_type = "spot", 
           annotate = "qc_lib_size_2000")
# plot reference (manually annotated) layers
plotSpots(spe, annotate = "ground_truth", 
          pal = "libd_layer_colors")
# histogram of numbers of expressed genes
hist(colData(spe)$detected, breaks = 20)
# plot number of expressed genes vs. number of cells per spot
plotSpotQC(spe, plot_type = "scatter", 
           x_metric = "cell_count", y_metric = "detected", 
           y_threshold = 400)
# select QC threshold for number of expressed genes
qc_detected <- colData(spe)$detected < 400
table(qc_detected)
colData(spe)$qc_detected <- qc_detected
# check spatial pattern of discarded spots
plotSpotQC(spe, plot_type = "spot", 
           annotate = "qc_detected")
# check spatial pattern of discarded spots if threshold is too high
qc_detected_1000 <- colData(spe)$detected < 1000
colData(spe)$qc_detected_1000 <- qc_detected_1000
plotSpotQC(spe, plot_type = "spot", 
           annotate = "qc_detected_1000")
# histogram of mitochondrial read proportions
hist(colData(spe)$subsets_mito_percent, breaks = 20)
# plot mitochondrial read proportion vs. number of cells per spot
plotSpotQC(spe, plot_type = "scatter", 
           x_metric = "cell_count", y_metric = "subsets_mito_percent", 
           y_threshold = 28)
# select QC threshold for mitochondrial read proportion
qc_mito <- colData(spe)$subsets_mito_percent > 28
table(qc_mito)
colData(spe)$qc_mito <- qc_mito
# check spatial pattern of discarded spots
plotSpotQC(spe, plot_type = "spot", 
           annotate = "qc_mito")
# check spatial pattern of discarded spots if threshold is too high
qc_mito_25 <- colData(spe)$subsets_mito_percent > 25
colData(spe)$qc_mito_25 <- qc_mito_25
plotSpotQC(spe, plot_type = "spot", 
           annotate = "qc_mito_25")
# histogram of cell counts
hist(colData(spe)$cell_count, breaks = 20)
# distribution of cells per spot
tbl_cells_per_spot <- table(colData(spe)$cell_count)
# plot number of expressed genes vs. number of cells per spot
plotSpotQC(spe, plot_type = "scatter", 
           x_metric = "cell_count", y_metric = "detected", 
           x_threshold = 10)
# select QC threshold for number of cells per spot
qc_cell_count <- colData(spe)$cell_count > 10
table(qc_cell_count)
colData(spe)$qc_cell_count <- qc_cell_count
# check spatial pattern of discarded spots
plotSpotQC(spe, plot_type = "spot", 
           annotate = "qc_cell_count")
# number of discarded spots for each metric
apply(cbind(qc_lib_size, qc_detected, qc_mito, qc_cell_count), 2, sum)
# combined set of discarded spots
discard <- qc_lib_size | qc_detected | qc_mito | qc_cell_count
table(discard)
# store in object
colData(spe)$discard <- discard
# check spatial pattern of combined set of discarded spots
plotSpotQC(spe, plot_type = "spot", 
           annotate = "discard")
# remove combined set of low-quality spots
spe <- spe[, !colData(spe)$discard]
dim(spe)
# distribution of cells per spot
tbl_cells_per_spot[1:13]
# as proportions
prop_cells_per_spot <- round(tbl_cells_per_spot / sum(tbl_cells_per_spot), 2)
prop_cells_per_spot[1:13]
# save object(s)
saveRDS(spe, file = "spe_qc.rds")
#ERROR IN SAVING


#NORMALIZATION
library(SpatialExperiment)
spe <- readRDS("spe_qc.rds")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("scran")
library(scran)
# calculate library size factors
spe <- computeLibraryFactors(spe)
summary(sizeFactors(spe))
hist(sizeFactors(spe), breaks = 20)
# calculate logcounts and store in object
spe <- logNormCounts(spe)
# check
assayNames(spe)
dim(counts(spe))
dim(logcounts(spe))
# save object(s)
saveRDS(spe, file = "spe_logcounts.rds")
#ERROR WITH SAVING


#FEATURE SELECTION
library(SpatialExperiment)
spe <- readRDS("spe_logcounts.rds")
# identify mitochondrial genes
is_mito <- grepl("(^MT-)|(^mt-)", rowData(spe)$gene_name)
table(is_mito)
# remove mitochondrial genes
spe <- spe[!is_mito, ]
dim(spe)

library(scran)
# fit mean-variance relationship
dec <- modelGeneVar(spe)
# visualize mean-variance relationship
fit <- metadata(dec)
plot(fit$mean, fit$var, 
     xlab = "mean of log-expression", ylab = "variance of log-expression")
curve(fit$trend(x), col = "dodgerblue", add = TRUE, lwd = 2)
# select top HVGs
top_hvgs <- getTopHVGs(dec, prop = 0.1)
length(top_hvgs)
# save object(s)
saveRDS(spe, file = "spe_hvgs.rds")
saveRDS(top_hvgs, file = "top_hvgs.rds")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("nnSVG")
library(nnSVG)
# subsample spots for faster runtime in this example
# note: skip this step in full analysis
n <- 100
set.seed(123)
ix <- sample(seq_len(n), n)
spe_nnSVG <- spe[, ix]

# filter low-expressed and mitochondrial genes
# using stringent filtering for faster runtime in this example
# note: use default filtering in full analysis
spe_nnSVG <- filter_genes(
  spe_nnSVG, filter_genes_ncounts = 10, filter_genes_pcspots = 3
)
# re-calculate logcounts after filtering
spe_nnSVG <- logNormCounts(spe_nnSVG)
# run nnSVG
set.seed(123)
spe_nnSVG <- nnSVG(spe_nnSVG)
# investigate results
# show results
head(rowData(spe_nnSVG), 3)
# number of significant SVGs
table(rowData(spe_nnSVG)$padj <= 0.05)
# show results for top n SVGs
rowData(spe_nnSVG)[order(rowData(spe_nnSVG)$rank)[1:6], ]
# identify top-ranked SVG
rowData(spe_nnSVG)$gene_name[which(rowData(spe_nnSVG)$rank == 1)]



#DIMENSIONALITY REDUCTION
library(SpatialExperiment)
spe <- readRDS("spe_hvgs.rds")
top_hvgs <- readRDS("top_hvgs.rds")
library(scater)
# compute PCA
set.seed(123)
spe <- runPCA(spe, subset_row = top_hvgs)
reducedDimNames(spe)
dim(reducedDim(spe, "PCA"))
# compute UMAP on top 50 PCs
set.seed(123)
spe <- runUMAP(spe, dimred = "PCA")
reducedDimNames(spe)
dim(reducedDim(spe, "UMAP"))
# update column names for easier plotting
colnames(reducedDim(spe, "UMAP")) <- paste0("UMAP", 1:2)
# save object(s)
saveRDS(spe, file = "spe_reduceddims.rds")
library(ggspavis)
# plot top 2 PCA dimensions
plotDimRed(spe, plot_type = "PCA")
# plot top 2 UMAP dimensions
plotDimRed(spe, plot_type = "UMAP")

#CLUSTERING
library(SpatialExperiment)
spe <- readRDS("spe_reduceddims.rds")
spe_full <- readRDS("spe_logcounts.rds")
#NON-SPATIAL CLUSTERING
library(scran)
# graph-based clustering
set.seed(123)
k <- 10
g <- buildSNNGraph(spe, k = k, use.dimred = "PCA")
g_walk <- igraph::cluster_walktrap(g)
clus <- g_walk$membership
table(clus)
# store cluster labels in column 'label' in colData
colLabels(spe) <- factor(clus)
library(ggspavis)
# plot clusters in spatial x-y coordinates
plotSpots(spe, annotate = "label", 
          pal = "libd_layer_colors")
# plot ground truth labels in spatial coordinates
plotSpots(spe, annotate = "ground_truth", 
          pal = "libd_layer_colors")
# plot clusters in PCA reduced dimensions
plotDimRed(spe, plot_type = "PCA", 
           annotate = "label", pal = "libd_layer_colors")
# plot clusters in UMAP reduced dimensions
plotDimRed(spe, plot_type = "UMAP", 
           annotate = "label", pal = "libd_layer_colors")
# save object(s)
saveRDS(spe, file = "spe_clustering.rds")
library(nnSVG)
# subsample spots for faster runtime in this example
# note: skip this step in full analysis
n <- 100
set.seed(123)
ix <- sample(seq_len(n), n)
spe_nnSVG <- spe_full[, ix]  ## note: using full object from logcounts step

# filter low-expressed and mitochondrial genes
# using stringent filtering for faster runtime in this example
# note: use default filtering in full analysis
spe_nnSVG <- filter_genes(
  spe_nnSVG, filter_genes_ncounts = 10, filter_genes_pcspots = 3
)
# re-calculate logcounts after filtering
spe_nnSVG <- logNormCounts(spe_nnSVG)
# run nnSVG
set.seed(123)
spe_nnSVG <- nnSVG(spe_nnSVG)
# select top SVGs
# note: using small subset in this example
# use larger set (e.g. top 1000 genes) in full analysis
n_top <- 50
ix_top <- order(rowData(spe_nnSVG)$rank)[1:n]
top_svgs <- rowData(spe_nnSVG)[ix_top, "gene_id"]
library(scater)
library(scran)
# dimensionality reduction
# compute PCA
# note: using small number of components in this example
# use larger number (e.g. ncomponents = 50) in full analysis
set.seed(123)
spe_nnSVG <- runPCA(spe_nnSVG, ncomponents = 10, subset_row = top_svgs)
# graph-based clustering
set.seed(123)
k <- 10
g <- buildSNNGraph(spe_nnSVG, k = k, use.dimred = "PCA")
g_walk <- igraph::cluster_walktrap(g)
clus <- g_walk$membership
table(clus)
# store cluster labels in column 'label' in colData
colLabels(spe_nnSVG) <- factor(clus)