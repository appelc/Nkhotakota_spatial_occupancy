## Generate predictions at new sites from stPGOcc model


## Description -----------------------------------------------------------------

#This code was written to be called from the command line (e.g., on an HPC).
#This is highly recommended as predictions will take a long time.
#Usage example: `Rscript code/03a_prediction.R <species> <model_name> <res_meters>`
#To run interactively, set the arguments manually below.

#Note that in the manuscript I used a 100-m resolution for predictions.
#For the sake of file size, the 1000-m resolution predictions are shown here.


## Set arguments and create output directory -----------------------------------

# ### If calling from command line:
# args <- commandArgs(trailingOnly = TRUE)
# species <- as.character(args[1])
# model_name <- as.character(args[2]) 
# res_meters <- as.numeric(args[3]) 

### If running interactively:
species <- 'elephant'
model_name <- 'model01'
res_meters <- 1000 #resolution of covariates for prediction (30, 100, or 1000) in meters

#Create output directory 
output_dir <- file.path('predictions', model_name, species)
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  cat("\n##### Created directory:", output_dir, "\n")
} else {
  cat("\n##### Directory already exists:", output_dir, "\n")
}


## Print run info --------------------------------------------------------------

cat('\n##########################\n')
cat('Creating prediction plots for:', species, '      from:', model_name, '\n')
cat('##########################\n')


## Load libraries --------------------------------------------------------------

library(tidyverse)
library(spOccupancy)
library(data.table)
library(stars)
library(patchwork)


## Read in data ----------------------------------------------------------------

#model output
(mod_path <- paste0('outputs/', model_name, '/', model_name, '_', species, '.RDS'))
mod <- readRDS(mod_path)

if(!exists('mod')) {message('Could not find mod file') } else{
    message('Loaded mod file')
}

#read in covariate values for prediction
if (res_meters == 30) {cov_df <- readRDS('data/raw/covar_for_prediction_30m.RDS')}
if (res_meters == 100) {cov_df <- readRDS('data/raw/covar_for_prediction_100m.RDS')}
if (res_meters == 1000) {cov_df <- readRDS('data/raw/covar_for_prediction_1000m.RDS')}

if(!exists('cov_df')) { message('No covariate data found. res_meters must be one of 30, 100, 1000') } else {
  message('Found covariate data')
}

#coordinates:
if (res_meters == 30) {cov_coords <- fread('data/raw/covariates_resampled_30.csv')}
if (res_meters == 100) {cov_coords <- fread('data/raw/covariates_resampled_100m.csv')}
if (res_meters == 1000) {cov_coords <- fread('data/raw/covariates_resampled_1000m.csv')}

if(!exists('cov_coords')) { message('No coordinates found. res_meters must be one of 30, 100, 1000') } else {
  message('Found coordinates')
}

    
# Read in original input --------------------------------------------------------

inputs <- readRDS('data/cleaned/sp_inputs.RDS')
input_sp <- inputs[[species]]
  rm(inputs) #remove full inputs
  season_names <- input_sp$occ.covs$year_season[1,]
  year_names <- input_sp$occ.covs$year[1,]
  

# Specify prediction params ----------------------------------------------------
  
#set dimensions
(J.pred <- nrow(cov_df)) #number of sites to predict on
(n.seasons <- dim(input_sp$y)[2]) #number of primary periods
(p.occ <- ncol(mod$beta.samples)) #number of predictors (incl. intercept) -- depends on model structure
  # colnames(mod$beta.samples)

#extract coords
coords.0 <- as.matrix(cov_coords[, c('x', 'y')])
  dim(coords.0)
  
#indicate which primary time periods (season_years) we are predicting for (HARD-CODE FOR NOW)
#best to choose >1 (dim 1 causes problems further down)
predict.seasons <- c('2019_cool','2024_cool')
(t.cols <- which(season_names %in% predict.seasons))
  
#and subset covariate array for these time periods only (otherwise it will predict on all)
X.0_full <- cov_df
X.0 <- X.0_full[, t.cols, ]

#clean up before predicting
rm(input_sp)
rm(cov_df)
  
#number of sites to predict at a time 
#(most efficient to do in batches; see https://groups.google.com/g/spocc-spabund-users/c/uoRUji5vkUQ/m/Str0t4BKAgAJ)
chunk_size <- 500
vals <- split(1:J.pred, ceiling(seq_along(1:J.pred) / chunk_size)) #creates list with indices of each batch
length(vals) #number of batches


# Run prediction ---------------------------------------------------------------
 
#initialize object
psi.quants <- array(NA, dim = c(5, J.pred, length(predict.seasons))) #5: to save mean and 4 quantile values
# dim(psi.quants) #5 quantiles, n sites, n periods for prediction

#loop thru chunks, generating and saving predictions and then deleting 'pred' object
for (jj in 1:length(vals)) {
  print(paste("Currently on set ", jj, " out of ", length(vals), sep = ''))
  curr.indx <- vals[[jj]]
  out.pred <- predict(mod, X.0[curr.indx, , , drop = FALSE], coords.0[curr.indx, ],
                      t.cols = t.cols, n.omp.threads = 3, verbose = TRUE, 
                      type = 'occupancy', ignore.RE = TRUE) #ignore RE?
  quantiles.jj <- apply(out.pred$psi.0.samples, c(2, 3), 
                        quantile, c(0.025, 0.25, 0.5, 0.75, 0.975))
  rm(out.pred)
  psi.quants[, curr.indx, ] <- quantiles.jj
  gc()
}

dimnames(psi.quants) <- list('quantile' = c('lci_025','lci_25','median','uci_75','uci_975'),
                             'site' = 1:(dim(psi.quants)[2]),
                             'season' = predict.seasons)

saveRDS(psi.quants, 
     paste0(output_dir, '/predictions_', species, '_', res_meters, 'm.RDS'))


## Use script 03.... to make maps (reads in predictions files for all species together)

