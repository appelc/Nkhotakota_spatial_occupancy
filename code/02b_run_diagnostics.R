## Evaluate spOcc models using saved RDS objects (print summaries, trace plots, etc.)

#This code was written to be called from the command line (e.g., on an HPC).
#Usage example: `Rscript code/02b_run_diagnostics.R <species> <model_dir>`
#To run interactively, set the arguments manually below.


## Set arguments and create output directory ----------------------------------

### If calling from command line:
# args <- commandArgs(trailingOnly = TRUE)
# species <- as.character(args[1])
# model_dir <- as.character(args[2])

### If running interactively:
species <- 'elephant'
model_dir <- 'model01'

#Create output directory
new_dir <- paste('evaluation', model_dir, sep = '/')
if (!dir.exists(new_dir)) {
  dir.create(new_dir, recursive = TRUE)
  cat("\n##### Created directory:", new_dir, "\n")
} else {
  cat("\n##### Directory already exists:", new_dir, "\n")
}


## Load libraries --------------------------------------------------------------

library(tidyverse)
library(spOccupancy)
library(MCMCvis)
library(tools)


## Write function to perform evaluation ---------------------------------------

#helper function to extract statistics
extract_stat <- function(x, name) {
  if (is.null(dim(x))) {
    val <- unname(x[name]) #if vector
  } else {
    val <- x[, name] #if matrix
  }
  as.numeric(val)
}

#main function
evaluate_model <- function(rds_path, new_dir = '.') {

  #DELETE:
  rds_path <- species_file
  new_dir
  
  #for all chains:
  mod <- readRDS(rds_path)
  (model_name <- tools::file_path_sans_ext(basename(rds_path)))
  (species_name <- gsub(paste0(model_dir,'_'), '', model_name))
  
  print(rds_path)
  print(species_name)
  str(mod)

  #create output sub-directory
  output_dir <- file.path(new_dir, species_name)
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    cat("\nCreated directory:", output_dir, "\n")
  } else {
    cat("\nDirectory already exists:", output_dir, "\n")
  }

  #save txt file with model summary
  capture.output(summary(mod),
                 file = file.path(output_dir, paste0(species_name, '_model_summary.txt')))

  #compute summary values and store as dataframe:

  #beta (occupancy)
  sum_beta <- summary(mod$beta.samples)
  param_names_beta <- colnames(as.matrix(mod$beta.samples))

  beta_summary <- data.frame(
    parameter = param_names_beta,
    mean = extract_stat(sum_beta$statistics, "Mean"),
    sd = extract_stat(sum_beta$statistics, "SD"),
    q2.5 = extract_stat(sum_beta$quantiles, "2.5%"),
    # q10 = extract_stat(sum_beta$quantiles, "10%"),
    q50 = extract_stat(sum_beta$quantiles, "50%"),
    # q90 = extract_stat(sum_beta$quantiles, "90%"),
    q97.5 = extract_stat(sum_beta$quantiles, "97.5%"),
    Rhat = mod$rhat$beta,
    ESS = mod$ESS$beta,
    parameter_type = "occupancy"
  )
  
  #sigma.sq.psi (occupancy RE variance)
  if (!is.null(mod$sigma.sq.psi.samples)){
    sum_sigmaSqPsi <- summary(mod$sigma.sq.psi.samples)
    param_names_sigmaSqPsi <- colnames(as.matrix(mod$sigma.sq.psi.samples))
    
    sigmaSqPsi_summary <- data.frame(
      parameter = param_names_sigmaSqPsi,
      mean = extract_stat(sum_sigmaSqPsi$statistics, "Mean"),
      sd = extract_stat(sum_sigmaSqPsi$statistics, "SD"),
      q2.5 = extract_stat(sum_sigmaSqPsi$quantiles, "2.5%"),
      # q10 = extract_stat(sum_sigmaSqPsi$quantiles, "10%"),
      q50 = extract_stat(sum_sigmaSqPsi$quantiles, "50%"),
      # q90 = extract_stat(sum_sigmaSqPsi$quantiles, "90%"),
      q97.5 = extract_stat(sum_sigmaSqPsi$quantiles, "97.5%"),
      Rhat = mod$rhat$sigma.sq.psi,
      ESS = mod$ESS$sigma.sq.psi,
      parameter_type = "occupancy_RE_variance"
    )
  }
  
  #alpha (detection)
  sum_alpha <- summary(mod$alpha.samples)
  param_names_alpha <- colnames(as.matrix(mod$alpha.samples))

  alpha_summary <- data.frame(
    parameter = param_names_alpha,
    mean = extract_stat(sum_alpha$statistics, "Mean"),
    sd = extract_stat(sum_alpha$statistics, "SD"),
    q2.5 = extract_stat(sum_alpha$quantiles, "2.5%"),
    # q10 = extract_stat(sum_alpha$quantiles, "10%"),
    q50 = extract_stat(sum_alpha$quantiles, "50%"),
    # q90 = extract_stat(sum_alpha$quantiles, "90%"),
    q97.5 = extract_stat(sum_alpha$quantiles, "97.5%"),
    Rhat = mod$rhat$alpha,
    ESS = mod$ESS$alpha,
    parameter_type = "detection"
  )

  #sigma.sq.p (detection RE variance)
  if (!is.null(mod$sigma.sq.p.samples)){
    sum_sigmaSqDet <- summary(mod$sigma.sq.p.samples)
    param_names_sigmaSqDet <- colnames(as.matrix(mod$sigma.sq.p.samples))
    
    sigmaSqDet_summary <- data.frame(
      parameter = param_names_sigmaSqDet,
      mean = extract_stat(sum_sigmaSqDet$statistics, "Mean"),
      sd = extract_stat(sum_sigmaSqDet$statistics, "SD"),
      q2.5 = extract_stat(sum_sigmaSqDet$quantiles, "2.5%"),
      # q10 = extract_stat(sum_sigmaSqDet$quantiles, "10%"),
      q50 = extract_stat(sum_sigmaSqDet$quantiles, "50%"),
      # q90 = extract_stat(sum_sigmaSqDet$quantiles, "90%"),
      q97.5 = extract_stat(sum_sigmaSqDet$quantiles, "97.5%"),
      Rhat = mod$rhat$sigma.sq.p,
      ESS = mod$ESS$sigma.sq.p,
      parameter_type = "detection_RE_variance"
    )
  } else {
    message('Skipping sigma.sq.p')
  }
  
  #theta (spatial covariance and temporal covariance parameters)
  sum_theta <- summary(mod$theta.samples)
  param_names_theta <- colnames(as.matrix(mod$theta.samples))

  theta_summary <- data.frame(
    parameter = param_names_theta,
    mean = extract_stat(sum_theta$statistics, "Mean"),
    sd = extract_stat(sum_theta$statistics, "SD"),
    q2.5 = extract_stat(sum_theta$quantiles, "2.5%"),
    # q10 = extract_stat(sum_theta$quantiles, "10%"),
    q50 = extract_stat(sum_theta$quantiles, "50%"),
    # q90 = extract_stat(sum_theta$quantiles, "90%"),
    q97.5 = extract_stat(sum_theta$quantiles, "97.5%"),
    Rhat = mod$rhat$theta,
    ESS = mod$ESS$theta,
    parameter_type = "covariance"
  )

  #eta (RE for each year) -- only with AR(1)
  if (!is.null(mod$eta.samples)){
    sum_eta <- summary(mod$eta.samples)
    param_names_eta <- colnames(as.matrix(mod$eta.samples))
    
    eta_summary <- data.frame(
      parameter = param_names_eta,
      mean = extract_stat(sum_eta$statistics, "Mean"),
      sd = extract_stat(sum_eta$statistics, "SD"),
      q2.5 = extract_stat(sum_eta$quantiles, "2.5%"),
      # q10 = extract_stat(sum_eta$quantiles, "10%"),
      q50 = extract_stat(sum_eta$quantiles, "50%"),
      # q90 = extract_stat(sum_eta$quantiles, "90%"),
      q97.5 = extract_stat(sum_eta$quantiles, "97.5%"),
      Rhat = NA,
      ESS = NA,
      parameter_type = "AR1_RE"
    )
  } else {
    message('Skipping eta')
  }

  #beta star (occupancy RE)
  if (!is.null(mod$beta.star.samples)){
    sum_betaStar <- summary(mod$beta.star.samples)
    param_names_betaStar <- colnames(as.matrix(mod$beta.star.samples))

    cat('\n#########\n')
    print(param_names_betaStar)
    print(extract_stat(sum_betaStar$statistics, "Mean"))
    print(extract_stat(sum_betaStar$statistics, "SD"))
    print(extract_stat(sum_betaStar$quantiles, "2.5%"))
    print(extract_stat(sum_betaStar$quantiles, "50%"))
    print(extract_stat(sum_betaStar$quantiles, "97.5%"))
    print(mod$rhat$beta.star)
    print(mod$ESS$beta.star)

    betaStar_summary <- data.frame(
      parameter = param_names_betaStar,
      mean = extract_stat(sum_betaStar$statistics, "Mean"),
      sd = extract_stat(sum_betaStar$statistics, "SD"),
      q2.5 = extract_stat(sum_betaStar$quantiles, "2.5%"),
      # q10 = extract_stat(sum_betaStar$quantiles, "10%"),
      q50 = extract_stat(sum_betaStar$quantiles, "50%"),
      # q90 = extract_stat(sum_betaStar$quantiles, "90%"),
      q97.5 = extract_stat(sum_betaStar$quantiles, "97.5%"),
      # Rhat = mod$rhat$beta.star,
      # ESS = mod$ESS$beta.star,
      parameter_type = "occupancy_RE"
    )
  } else {
    message('Skipping beta.star')
  }
  
  #alpha star (detection RE)
  if (!is.null(mod$alpha.star.samples)){
    sum_alphaStar <- summary(mod$alpha.star.samples)
    param_names_alphaStar <- colnames(as.matrix(mod$alpha.star.samples))
    
    alphaStar_summary <- data.frame(
      parameter = param_names_alphaStar,
      mean = extract_stat(sum_alphaStar$statistics, "Mean"),
      sd = extract_stat(sum_alphaStar$statistics, "SD"),
      q2.5 = extract_stat(sum_alphaStar$quantiles, "2.5%"),
      # q10 = extract_stat(sum_alphaStar$quantiles, "10%"),
      q50 = extract_stat(sum_alphaStar$quantiles, "50%"),
      # q90 = extract_stat(sum_alphaStar$quantiles, "90%"),
      q97.5 = extract_stat(sum_alphaStar$quantiles, "97.5%"),
      # Rhat = mod$rhat$alpha.star,
      # ESS = mod$ESS$alpha.star,
      parameter_type = "detection_RE"
    )
  } else {
    message('Skipping alpha.star')
  }

  #bind and save
  param_summaries <- beta_summary %>% 
    bind_rows(
      if (exists('sigmaSqPsi_summary')) sigmaSqPsi_summary else NULL,
      if (exists('betaStar_summary')) betaStar_summary else NULL,
      if (exists('alpha_summary')) alpha_summary else NULL, 
      if (exists('alphaStar_summary')) alphaStar_summary else NULL,
      if (exists('sigmaSqDet_summary')) sigmaSqDet_summary else NULL,
      if (exists('theta_summary')) theta_summary else NULL, 
      if (exists('eta_summary')) eta_summary else NULL) %>%
    mutate(species = species_name)

  write.csv(param_summaries, paste0(output_dir, '/', species_name, '_parameter_summaries.csv'))


  #trace plots
  pdf(file.path(output_dir, paste0(species_name, '_traceplots_beta.pdf')))
  print(plot(mod, 'beta'))
  dev.off()

  pdf(file.path(output_dir, paste0(species_name, '_traceplots_alpha.pdf')))
  print(plot(mod, 'alpha'))
  dev.off()

  pdf(file.path(output_dir, paste0(species_name, '_traceplots_theta.pdf')))
  print(plot(mod, 'theta'))
  dev.off()


  #plot param estimates
  pdf(file.path(output_dir, paste0(species_name, '_parameter_plots.pdf')), width = 8, height = 10)
  beta_plot <- MCMCvis::MCMCplot(mod$beta.samples, ref_ovl = TRUE,
                    guide_axis = F, main = 'Occupancy')
  if(!is.null(mod$sigma.sq.psi.samples)) {
    sigmaPsi_plot <- MCMCvis::MCMCplot(mod$sigma.sq.psi.samples, ref_ovl = TRUE,
                                    guide_axis = F, main = 'Occupancy RE variance')
  }
  if(!is.null(mod$beta.star.samples)) {
    betaStar_plot <- MCMCvis::MCMCplot(mod$beta.star.samples, ref_ovl = TRUE,
                                    guide_axis = F, main = 'Occupancy RE')
  }
  alpha_plot <- MCMCvis::MCMCplot(mod$alpha.samples, ref_ovl = TRUE,
                    guide_axis = F, main = 'Detection')
  if(!is.null(mod$sigma.sq.p.samples)) {
    sigmaDet_plot <- MCMCvis::MCMCplot(mod$sigma.sq.p.samples, ref_ovl = TRUE,
                                    guide_axis = F, main = 'Detection RE')
  }
  theta_plot <- MCMCvis::MCMCplot(mod$theta.samples, ref_ovl = TRUE, 
                    guide_axis = F, main = 'Spatial covariance')
  if(!is.null(mod$eta.samples)) {
    eta_plot <- MCMCvis::MCMCplot(mod$eta.samples, ref_ovl = TRUE, 
                                  #labels = c('2019','2020','2021','2022','2023','2024'), 
                                  labels = c('2019_wet','2019_cool','2019_dry','2020_wet','2020_cool','2020_dry','2021_wet','2021_cool','2021_dry','2022_wet','2022_cool','2022_dry','2023_wet','2023_cool','2023_dry','2024_wet','2024_cool','2024_dry'),
                                  guide_lines = F, guide_axis = F, main = 'Temporal random effect')
  }
  dev.off()

  #posterior predictive checks (GoF) -- these have their own script now (04_ppc.R)

  #clean up
  rm(mod)
  gc()

}


## Perform evaluation ----------------------------------------------------------

#list all models
(model_path <- paste0('outputs/', model_dir))
(model_files <- list.files(model_path, full.names = TRUE, pattern = 'RDS'))

#extract file for this species
modname_part1 <- sapply(strsplit(model_dir, '\\_'), '[', 1)
file_name <- paste(modname_part1, species, sep = '_')
species_file <- model_files[grepl(file_name, model_files)]

print(model_files)
print(file_name)
print(species_file)

cat('\n###############\n')
cat('Evaluating:     ', species, '\n')

#run function
evaluate_model(species_file, new_dir)

cat('Done!\n')
cat('\n###############\n')


## Now proceed to 02c_run_ppc.R