#install.packages("devtools")
library(devtools)


#devtools::install_github("fanzhanglab/STEAM")
#install.packages("Seurat")
#install.packages("caret")

# Load libraries

library(STEAM)
library(Seurat)
library(caret)

# Load data

# DLPFC.RDS is sample 151669 from the DLPFC dataset from Visium by 10X Genomics.
# This data is already preprocessed and includes a `SCTransform` scaled data matrix by `Seurat`, 
# spatial coordinates, and pre-identified cluster labels.

data(DLPFC)

# Create the STEAM object

# Create a STEAM object using the LoadSTEAM function.
STEAM.obj <- LoadSTEAM(
  count_exp = DLPFC$matrix,  # Scaled gene expression matrix
  spatial = DLPFC$coordinates,  # Spatial coordinates of cells
  labels = DLPFC$labels,  # Cluster labels
  Seurat.obj = NULL  # No additional Seurat object is passed
)

# Run STEAM on the STEAM object using specified parameters

STEAM.obj <- RunSTEAM(
  STEAM.obj,
  train.ratio = 0.8,  # Proportion of data for training
  n.size = 5,  # Neighborhood size for spatial averaging
  seed = 123,  # Random seed for reproducibility
  cv.folds = 10,  # Number of cross-validation folds
  cv.repeats = 3,  # Number of cross-validation repetitions
  trainval.ratio = 0.8,  # Ratio of training set split into train/validation
  model = "rf",  # Random Forest model (other options: "svm", "xgb", "multinom")
  n.tree = 100,  # Number of trees for Random Forest
  kernel = 'linear',  # Kernel type for SVM
  train.folder.name = 'train.out',  # Folder for storing training outputs
  allowParallel = TRUE  # Enable parallel computing for faster execution
)

# View the Performance Metrics from STEAM
ViewMetrics(STEAM.obj)  # Generates Kappa, F1-score, Accuracy, etc.

# View the Feature Importance

# Feature importance is calculated depending on the model used in STEAM.obj:
# - rf: Mean decrease in Gini from varImp() in `caret`
# - xgb: Gain-based importance from varImp() in `caret`
# - svm: Absolute magnitude of SVM coefficients
# - multinom: Sum of absolute coefficients across all classes

STEAM:::feature_importance(STEAM.obj, top_n = 10, title = "Top Features by Importance")

# View Gene Expression Distribution across Clusters

# Extracts and plots expression of a specified gene across different clusters/layers.
STEAM::feature_expression(STEAM.obj, feature_name = "TMSB10", title = "Expression Across Layers")

# View the Misclassifications

# Plots misclassified cells in spatial coordinates (black dots indicate errors).
STEAM::plot_misclassified_cells(STEAM.obj)
