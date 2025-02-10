## ----eval=FALSE---------------------------------------------------------------
#  if (!requireNamespace("BiocManager"))
 #     install.packages("BiocManager")
  #BiocManager::install("mbkmeans")
  #BiocManager::install("TENxPBMCData")
  #BiocManager::install("DelayedMatrixStats")
## ----options, include=FALSE, message=FALSE, warning=FALSE---------------------
knitr::opts_chunk$set(cache=FALSE, error=FALSE, message=FALSE, warning=FALSE)

## ----packages, message=FALSE, warning=FALSE-----------------------------------
library(TENxPBMCData)
library(scater)
library(SingleCellExperiment)
library(mbkmeans)
library(DelayedMatrixStats)

## -----------------------------------------------------------------------------
tenx_pbmc4k <- TENxPBMCData(dataset = "pbmc4k")

set.seed(1034)
idx <- sample(seq_len(NCOL(tenx_pbmc4k)), 100)
sce <- tenx_pbmc4k[, idx]

#normalization
sce <- logNormCounts(sce)

vars <- rowVars(logcounts(sce))
names(vars) <- rownames(sce)
vars <- sort(vars, decreasing = TRUE)

sce1000 <- sce[names(vars)[1:1000],]

sce1000

## -----------------------------------------------------------------------------
res <- mbkmeans(sce1000, clusters = 5,
                reduceMethod = NA,
                whichAssay = "logcounts")

## -----------------------------------------------------------------------------
batchsize <- blocksize(sce1000)
batchsize

## -----------------------------------------------------------------------------
res_random <- mbkmeans(sce1000, clusters = 5, 
                       reduceMethod = NA,
                       whichAssay = "logcounts",
                       initializer = "random")
table(res$Clusters, res_random$Clusters)

## -----------------------------------------------------------------------------
sce1000 <- runPCA(sce1000, ncomponents=20)

ks <- seq(5, 15)
res <- lapply(ks, function(k) {
  mbkmeans(sce1000, clusters = k,
           reduceMethod = "PCA",
           calc_wcss = TRUE, num_init=10)
})

wcss <- sapply(res, function(x) sum(x$WCSS_per_cluster))
plot(ks, wcss, type = "b")

## -----------------------------------------------------------------------------
res_full <- mbkmeans(sce1000, clusters = 5,
                     reduceMethod = NA,
                     whichAssay = "logcounts",
                     initializer = "random",
                     batch_size = ncol(sce1000))
res_classic <- kmeans(t(logcounts(sce1000)), 
                      centers = 5)
table(res_full$Clusters, res_classic$cluster)

## -----------------------------------------------------------------------------
mat <- reducedDim(sce1000, "PCA")
dim(mat)

## -----------------------------------------------------------------------------
library(bluster)
clusterRows(mat, MbkmeansParam(centers=5))

## -----------------------------------------------------------------------------
clusterRows(mat, MbkmeansParam(centers=5, batch_size=10))

## -----------------------------------------------------------------------------
clusterRows(mat, MbkmeansParam(centers=5, batch_size=10), full = TRUE)

## -----------------------------------------------------------------------------
logcounts(sce1000)
clusterRows(t(logcounts(sce1000)), MbkmeansParam(centers=4))

