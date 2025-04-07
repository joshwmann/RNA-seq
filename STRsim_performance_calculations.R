##############################################################################
# Step 1: Run SRTsim to Simulate Spatial Transcriptomics Data (Reference-Based)
##############################################################################

##############################################################################
# 1) Load Required Packages
##############################################################################
library(STexampleData)     
library(SpatialExperiment) 
library(SRTsim)            
library(scater)            
library(scran)             
library(ggspavis)          

##############################################################################
# 2) Load Visium_humanDLPFC Data as Reference
##############################################################################
#Dataset is from the STexampleData package
spe <- Visium_humanDLPFC()

# Inspect the data
dim(spe)                  
assayNames(spe)           
head(rowData(spe))        
head(colData(spe))        
head(spatialCoords(spe))  

##############################################################################
# 3) Prepare Data for SRTsim
##############################################################################
# Extract gene expression and spatial metadata
example_count <- counts(spe)
example_count[1:5, 1:5]

# Extract x, y coordinates and ground truth labels
colnames(colData(spe))
example_loc <- cbind(spatialCoords(spe), label = as.numeric(as.factor(colData(spe)$ground_truth)))
colnames(example_loc) <- c("x", "y", "label")
head(example_loc)

head(example_loc)
# Create SRTsim object

#Ensuring Data Compatibility for SRTsim
all(colnames(example_count) == rownames(example_loc))  # Should return TRUE
#install.packages("Matrix")
library(Matrix)  # Ensure Matrix library is loaded
#Filter out genes with low expression before running
keep_genes <- rowSums(example_count > 10) > 20  # Keep genes detected in at least 5 spots
table(keep_genes)
example_count <- example_count[keep_genes, ]
example_count <- example_count[, colData(spe)$in_tissue == 1] #Keep spots over tissue
dim(example_count)  # Check number of remaining genes
# Convert example_count to a standard dense matrix
example_count <- as.matrix(example_count)
example_count[1:5, 1:5]
example_loc <- example_loc[colData(spe)$in_tissue == 1, ]
dim(example_count)
dim(example_loc)
head(example_loc)
all(colnames(example_count) == rownames(example_loc))
simSRT <- createSRT(count_in = example_count, loc_in = example_loc)
simSRT

##############################################################################
# 4) Reference-Based Simulation (Tissue & Domain Schemes)
##############################################################################
set.seed(123)  # to ensure reproducibility


# (A) Tissue-Based Simulation
simSRT_tissue <- srtsim_fit(simSRT, sim_schem = "tissue")  # estimate parameters
simSRT_tissue <- srtsim_count(simSRT_tissue)               # generate synthetic data

# (B) Domain-Based Simulation
#simSRT_domain <- srtsim_fit(simSRT, sim_schem = "domain")  # estimate parameters
#simSRT_domain <- srtsim_count(simSRT_domain)               # generate synthetic data


##############################################################################
# 5) Inspect and Visualize Synthetic Data
##############################################################################
# Check part of the simulated count matrix for tissue-based data
simCounts(simSRT_tissue)[1:5, 1:5]

# (Optional) Compare synthetic data to reference
simSRT_tissue <- compareSRT(simSRT_tissue)
visualize_metrics(simSRT_tissue)

# Visualize a specific gene’s expression in the simulated data
visualize_gene(simsrt = simSRT_tissue, plotgn = "ENSG00000117632", rev_y = TRUE)
#visualize_gene(simsrt = simSRT_domain, plotgn = "ENSG00000175130", rev_y = TRUE)


##############################################################################
# 6) Run STEAM (Tutorial Style) on SRTsim Output
##############################################################################

# 1) Load the tutorial libraries
library(Seurat)
library(caret)
library(STEAM)

# 2) Extract the simulated expression and metadata from simSRT_tissue
counts_sim <- simCounts(simSRT_tissue)            # get expression (sparse matrix)
counts_sim <- as.matrix(counts_sim)               # convert to base R matrix
counts_sim <- log2(counts_sim + 1)                # mimic the "scaled" data from tutorial

meta_sim <- simcolData(simSRT_tissue)             # get spot-level data (x, y, label)

# 3) Create the STEAM object exactly like the tutorial
STEAM.obj <- LoadSTEAM(
  count_exp  = counts_sim,                      # expression matrix
  spatial    = as.matrix(meta_sim[, c("x", "y")]),  # x,y coordinates
  labels     = meta_sim$label,                  # cluster labels
  Seurat.obj = NULL
)

str(STEAM.obj)

# 4) Run STEAM using the same parameters as in the tutorial
STEAM.obj <- RunSTEAM(
  STEAM.obj,
  train.ratio      = 0.8,      # same as tutorial
  n.size           = 5,
  seed             = 123,
  cv.folds         = 10,
  cv.repeats       = 3,
  trainval.ratio   = 0.8,
  model            = "rf",     # random forest, as in tutorial
  n.tree           = 100,
  kernel           = 'linear',
  train.folder.name = 'train.out',
  allowParallel    = TRUE
)

# 5) View performance metrics and do the same visualizations as tutorial
ViewMetrics(STEAM.obj)
STEAM:::feature_importance(STEAM.obj, top_n = 10, title = "Top Features by Importance")
STEAM::feature_expression(STEAM.obj, feature_name = "ENSG00000117632", title = "Expression Across Layers")
STEAM::plot_misclassified_cells(STEAM.obj)
