## Spatial occupancy models 


## Description -----------------------------------------------------------------

#Based on the spatiotemporal occupancy vignette from spOccupancy:
#https://doserlab.com/files/spoccupancy-web/articles/spacetimemodelshtml

#This code was written to be called from the command line (e.g., on an HPC).
#This is highly recommended as the models will take a long time to run.
#Usage example: `Rscript code/02a_run_model.R <species> <n.batch> <output_dir>`
#To run interactively, set the arguments manually below.


## Set arguments and create output directory -----------------------------------

### If calling from command line:
# args <- commandArgs(trailingOnly = TRUE)
# species <- as.character(args[1])
# n.batch <- as.numeric(args[2])
# output_dir <- as.character(args[3]) 

### If running interactively:
species <- 'elephant'
n.batch <- 20000
output_dir <- 'outputs/model01'

#Create output directory if it doesn't exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  cat("\nCreated directory:", output_dir, "\n")
} else {
  cat("\nDirectory already exists:", output_dir, "\n")
}


## Hard-code model name --------------------------------------------------------

mod_name <- 'model01' #suggest not using underscores for ease of parsing later


## Print run info --------------------------------------------------------------

cat('\n##########################\n')
cat('Processing:', species, '\n')
cat("Total number of samples:", n.batch * 25, "\n") #hard coded with batch.length=25
cat('##########################\n')


## Load libraries --------------------------------------------------------------

library(tidyverse)
library(spOccupancy)


## Read in input objects -------------------------------------------------------

inputs <- readRDS('data/cleaned/sp_inputs.RDS')


## Specify formulas ------------------------------------------------------------

occ_formula <- ~fence_south*factor(year) +
                season_wet + season_dry +
                elev_scaledNS + season_wet:elev_scaledNS + season_dry:elev_scaledNS +
                slope_scaledNS + season_wet:slope_scaledNS + season_dry:slope_scaledNS +
                tree_vol_scaledNS + season_wet:tree_vol_scaledNS + season_dry:tree_vol_scaledNS +
                dist_boundary_log_scaledNS + season_wet:dist_boundary_log_scaledNS + season_dry:dist_boundary_log_scaledNS +
                dist_dambo_log_scaledNS + season_wet:dist_dambo_log_scaledNS + season_dry:dist_dambo_log_scaledNS +
                dist_rivers_log_scaledNS + season_wet:dist_rivers_log_scaledNS + season_dry:dist_rivers_log_scaledNS

det_formula <- ~fence_south + year_num + fence_south:year_num + 
                season_wet + season_dry +
                (1 | site)


## Model specifications --------------------------------------------------------

#spatial random effects
cov.model <- 'exponential'
n.neighbors <- 15

#temporal random effects
ar1 <- FALSE

#MCMC values (batch number is specified above in 'args')
n.burn <- 15000
n.thin <- 25
batch.length <- 25 #keep length=25 and increase n.batches above
n.chains <- 3 #all chains together in this version
n.report <- (n.batch) / 10 #report every 10%
n.threads <- 1 #found no performance gains with >1 using OpenMP


## Prep model ------------------------------------------------------------------

#get data for this species
inputs_sp <- inputs[[species]]

#compute pairwise distance between all sites
pw_dist <- dist(inputs_sp$coords)

#get 'z' states for initial values
z.inits <- apply(inputs_sp$y, c(1, 2), function(a) as.numeric(sum(a, na.rm = TRUE) > 0))

#no need when running all chains together, but still set seed for reproducibility
set.seed(7348)

#set inits
inits <- list(beta = 0,
              alpha = 0,
              z = z.inits,
              sigma.sq = 1,
              phi = (3 / mean(pw_dist)),
              phi = runif(1, 3 / max(pw_dist), 3 / min(pw_dist)),
              sigma.sq.t = 1.5,
              rho = 0.2)

#set priors
priors <- list(beta.normal = list(mean = 0, var = 2.72),
               alpha.normal = list(mean = 0, var = 2.72),
               # sigma.sq.t.ig = c(2, 0.5), #only if ar1=TRUE
               # rho.unif = c(-1, 1), #only if ar1=TRUE
               sigma.sq.ig = c(2, 1),
               sigma.sq.p.ig = list(2, 1),
               phi.unif = c(3 / max(pw_dist), 3 / min(pw_dist)))


## Run model -------------------------------------------------------------------

output <- stPGOcc(occ.formula = occ_formula,
                  det.formula = det_formula,
                  data = inputs_sp,
                  inits = inits,
                  priors = priors,
                  cov.model = cov.model,
                  n.neighbors = n.neighbors,
                  n.batch = n.batch,
                  batch.length = batch.length,
                  verbose = TRUE,
                  ar1 = ar1,
                  n.report = n.report,
                  n.burn = n.burn,
                  n.thin = n.thin,
                  n.chains = n.chains,
                  n.omp.threads = n.threads
                  )

cat('##########################\n')
cat('Finished sampling\n')


## Save ------------------------------------------------------------------------

output_path <- paste0(output_dir, '/', mod_name, '_', species, '.RDS')
saveRDS(output, output_path)

cat('Results saved to:', output_path, '\n')
cat('##########################\n')
