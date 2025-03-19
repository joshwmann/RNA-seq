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

# Extract x, y coordinates and ground truth labels
colnames(colData(spe))
example_loc <- data.frame(
  x = colData(spe)$array_col,   # Spatial column coordinate
  y = colData(spe)$array_row,   # Spatial row coordinate
  label = colData(spe)$ground_truth  # Ground truth annotations
)
colnames(example_loc) <- c("x", "y", "label")

# Convert spatial coordinates to numeric
example_loc$x <- as.numeric(example_loc$x)
example_loc$y <- as.numeric(example_loc$y)
example_loc$label <- as.character(example_loc$label)

# Create SRTsim object

#Ensuring Data Compatibility for SRTsim
colnames(example_count)[1:5]  # Check first 5 column names of count matrix
rownames(example_loc)[1:5]    # Check first 5 row names of spatial metadata
rownames(example_loc) <- colData(spe)$barcode_id  # Set row names to barcode IDs
all(colnames(example_count) == rownames(example_loc))  # Should return TRUE
#install.packages("Matrix")
library(Matrix)  # Ensure Matrix library is loaded
# Convert example_count to a standard dense matrix
example_count <- as.matrix(example_count) 
#Filter out genes with low expression before running
keep_genes <- rowSums(example_count > 0) > 5  # Keep genes detected in at least 5 spots
example_count <- example_count[keep_genes, ]
dim(example_count)  # Check number of remaining genes
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
simSRT_domain <- srtsim_fit(simSRT, sim_schem = "domain")  # estimate parameters
simSRT_domain <- srtsim_count(simSRT_domain)               # generate synthetic data


##############################################################################
# 5) Inspect and Visualize Synthetic Data
##############################################################################
# Check part of the simulated count matrix for tissue-based data
simCounts(simSRT_tissue)[1:5, 1:5]

# (Optional) Compare synthetic data to reference
simSRT_tissue <- compareSRT(simSRT_tissue)
visualize_metrics(simSRT_tissue)

# Visualize a specific gene’s expression in the simulated data
visualize_gene(simsrt = simSRT_tissue, plotgn = "ENSG00000183036", rev_y = TRUE)
visualize_gene(simsrt = simSRT_domain, plotgn = "ENSG00000168314", rev_y = TRUE)
