#Installation
#BiocManager::install("scMultiSim")
library(scMultiSim)
library(tidyverse)


#Running Simulation
data(GRN_params_100)

results <- sim_true_counts(list(
  # required options
  GRN = GRN_params_100,
  tree = Phyla3(),
  num.cells = 500,
  # optional options
  num.cif = 20,
  discrete.cif = F,
  cif.sigma = 0.1
  # ... other options
))

#Shiny app
run_shiny()
res <- sim_true_counts(opt)

#Technical noise and batch effect
add_expr_noise(results)
divide_batches(results, nbatch = 2)

#Visualization
plot_tsne(results$counts, results$cell_meta$pop)

