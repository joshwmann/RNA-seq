##############################################################################
# 1) Simulate Data
##############################################################################
library(splatter)
library(scater)
library(scran)
library(igraph)
library(mclust)
library(mbkmeans)
library(DelayedMatrixStats)

set.seed(123)
params <- newSplatParams(batchCells = 200, group.prob = c(0.7, 0.3)) #200 cells, 70% one group, 30% another
sim <- splatSimulate(params, method = "groups", verbose = FALSE)

##############################################################################
# 2) Normalize, HVG selection, PCA
##############################################################################
sim <- logNormCounts(sim)
dec <- modelGeneVar(sim)
top_hvgs <- getTopHVGs(dec, prop = 0.1)
sim <- runPCA(sim, subset_row = top_hvgs)

##############################################################################
# 3A) Original Graph-Based Clustering
##############################################################################
graph_orig <- buildSNNGraph(sim, k = 10, use.dimred = "PCA")
clus_orig  <- igraph::cluster_walktrap(graph_orig)$membership

##############################################################################
# 3B) mbkmeans Clustering
##############################################################################
res_mbk  <- mbkmeans(sim[top_hvgs, ], clusters = 2, reduceMethod = NA, whichAssay = "logcounts")
clus_mbk <- res_mbk$Clusters

##############################################################################
# 4) ARI Computations
##############################################################################
true_labels <- colData(sim)$Group

ari_orig_vs_truth <- adjustedRandIndex(clus_orig,  true_labels)
ari_mbk_vs_truth  <- adjustedRandIndex(clus_mbk,   true_labels)
ari_orig_vs_mbk   <- adjustedRandIndex(clus_orig,  clus_mbk)

cat("ARI (Original vs. Truth):      ", ari_orig_vs_truth, "\n")
cat("ARI (mbkmeans vs. Truth):      ", ari_mbk_vs_truth,  "\n")
cat("ARI (Original vs. mbkmeans):   ", ari_orig_vs_mbk,   "\n")
