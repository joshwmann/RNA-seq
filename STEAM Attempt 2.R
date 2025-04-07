##############################################################################
# 1) Install/Load Necessary Packages
##############################################################################
# install.packages("devtools")
# devtools::install_github("fanzhanglab/STEAM")
# install.packages("Seurat")
# install.packages("caret")

library(splatter)
library(scater)
library(scran)
library(igraph)
library(mclust)
library(mbkmeans)
library(DelayedMatrixStats)
library(STEAM)
library(Seurat)
library(caret)

##############################################################################
# 2) Simulate Data with Splatter
##############################################################################
set.seed(123)

# Increase cells to 2000 and create 10 distinct groups
params <- newSplatParams(batchCells = 2000,
                         group.prob = c(0.1, 0.15, 0.05, 0.07, 0.03,
                                        0.2, 0.06, 0.04, 0.16, 0.14))
sim <- splatSimulate(params, method = "groups", verbose = FALSE)

cat("Dimensions of simulated dataset:\n")
print(dim(sim))  # Genes x Cells

##############################################################################
# 3) Normalize, Find HVGs, Run PCA
##############################################################################
sim <- logNormCounts(sim)
dec <- modelGeneVar(sim)
top_hvgs <- getTopHVGs(dec, prop = 0.1)  # top 10% HVGs

# Run PCA on the HVGs
sim <- runPCA(sim, subset_row = top_hvgs)

##############################################################################
# 4) Clustering Methods
##############################################################################
set.seed(123)
graph_orig <- buildSNNGraph(sim, k = 10, use.dimred = "PCA")
clus_orig  <- igraph::cluster_walktrap(graph_orig)$membership

set.seed(123)
res_mbk  <- mbkmeans(sim[top_hvgs, ],
                     clusters = 10,
                     reduceMethod = NA,
                     whichAssay   = "logcounts")
clus_mbk <- res_mbk$Clusters

##############################################################################
# 5) Create Dummy (All-Zero) Spatial Coordinates
##############################################################################
n_cells <- ncol(sim)
dummy_spatial <- data.frame(
  x = rep(0, n_cells),
  y = rep(0, n_cells),
  row.names = colnames(sim)
)

##############################################################################
# 6) Prepare Gene Expression Matrix for STEAM
##############################################################################
expr_matrix <- assay(sim, "logcounts")

##############################################################################
# 7) Use STEAM to Evaluate Graph-Based Clustering Labels
##############################################################################
steam_graph <- LoadSTEAM(
  count_exp   = expr_matrix,
  spatial     = dummy_spatial,  # dummy coords for every cell
  labels      = clus_orig,
  Seurat.obj  = NULL
)

steam_graph <- RunSTEAM(
  steam_graph,
  train.ratio     = 0.8,
  n.size          = 5,   # STEAM will see every cell as in the same “zone”
  seed            = 123,
  cv.folds        = 3,   # reduce for demonstration
  cv.repeats      = 1,   # reduce for demonstration
  trainval.ratio  = 0.8,
  model           = "rf",
  n.tree          = 50,  # reduce for speed; typically ~100 or more
  kernel          = 'linear',
  train.folder.name = 'train.graph',
  allowParallel   = TRUE
)

cat("\n=== STEAM Performance (Graph-Based Labels) ===\n")
ViewMetrics(steam_graph)  # Check Kappa, F1, Accuracy, etc.

##############################################################################
# 8) Use STEAM to Evaluate mbkmeans Clustering Labels
##############################################################################
steam_mbk <- LoadSTEAM(
  count_exp   = expr_matrix,
  spatial     = dummy_spatial,   # same dummy coords
  labels      = clus_mbk,
  Seurat.obj  = NULL
)

steam_mbk <- RunSTEAM(
  steam_mbk,
  train.ratio     = 0.8,
  n.size          = 5,
  seed            = 123,
  cv.folds        = 3,
  cv.repeats      = 1,
  trainval.ratio  = 0.8,
  model           = "rf",
  n.tree          = 50,
  kernel          = 'linear',
  train.folder.name = 'train.mbk',
  allowParallel   = TRUE
)

cat("\n=== STEAM Performance (mbkmeans Labels) ===\n")
ViewMetrics(steam_mbk)
