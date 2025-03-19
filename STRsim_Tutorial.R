#install.packages('devtools')
devtools::install_github('xzhoulab/SRTsim')
library(SRTsim)
#Reference Based
str(exampleLIBD)
example_count   <- exampleLIBD$count
example_loc     <- exampleLIBD$info[,c("imagecol","imagerow","layer")]
colnames(example_loc) <- c("x","y","label")
simSRT  <- createSRT(count_in=example_count,loc_in =example_loc)
simSRT
## Set a seed for reproducible simulation
set.seed(1)
## Estimate model parameters for data generation
simSRT1 <- srtsim_fit(simSRT,sim_schem="tissue")
## Generate synthetic data with estimated parameters
simSRT1 <- srtsim_count(simSRT1)
## Set a seed for reproducible simulation
set.seed(1)
## Estimate model parameters for data generation
simSRT2 <- srtsim_fit(simSRT,sim_schem="domain")
## Generate synthetic data with estimated parameters
simSRT2 <- srtsim_count(simSRT2)
simCounts(simSRT1)[1:3,1:3]
## Compute metrics 
simSRT1   <- compareSRT(simSRT1)
## Visualize Metrics
visualize_metrics(simSRT1)
visualize_gene(simsrt=simSRT1,plotgn = "ENSG00000183036",rev_y=TRUE)
visualize_gene(simsrt=simSRT2,plotgn = "ENSG00000168314",rev_y=TRUE)

#Reference Free
library('SRTsim')
shinySRT1 <- SRTsim_shiny()
