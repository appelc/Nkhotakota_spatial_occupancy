## Posterior predictive checks for spOccupancy models using a sample of posterior values

#This code was written to be called from the command line (e.g., on an HPC).
#Usage example: `Rscript code/02c_run_ppc.R <species> <model_name> <n_samples>`
#To run interactively, set the arguments manually below.


## Set arguments and create output directory -----------------------------------

### If calling from command line:
# args <- commandArgs(trailingOnly = TRUE)
# species <- as.character(args[1])
# model_name <- as.character(args[2]) 
# n_samples <- as.numeric(args[3])

### If running interactively:
species <- 'elephant'
model_name <- 'model01'
n_samples <- 1000 #number of posterior samples to use (>10k?)

#Construct output directory path (these folders should already exist)
output_dir <- file.path('evaluation', model_name, species)
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  cat("\n##### Created directory:", output_dir, "\n")
} else {
  cat("\n##### Directory already exists:", output_dir, "\n")
}

## Print run info --------------------------------------------------------------

cat('\n##########################\n')
cat('Conducting posterior predictive checks for:', species, 'using', n_samples, 'samples', '\n')
cat(' \n')


## Load libraries --------------------------------------------------------------

library(tidyverse)
library(spOccupancy)
library(coda)


## Read in model output --------------------------------------------------------

file_name <- paste(model_name, species, sep = '_')
mod_path <- paste0('outputs/', model_name, '/', file_name, '.RDS')
mod <- readRDS(mod_path)

cat('Loaded output for', model_name, '\n')

#Read in input file for site/season names (must be an easier way)
inputs <- readRDS('data/cleaned/sp_inputs.RDS')
input_sp <- inputs[[species]]
site_names <- input_sp$occ.covs$site
season_names <- input_sp$occ.covs$year_season[1, ]
rm(inputs)


## Subset posterior samples ----------------------------------------------------

# mod$n.post #number of posterior samples
# dim(mod$psi.samples) #number of samples, number of sites, number of primary time periods

#how many to keep?
cat('Keeping', n_samples, 'posterior samples\n')

#subsample
set.seed(123)
keep_ids <- sample(1:dim(mod$psi.samples)[1], size = n_samples)

mod_reduced <- mod

#subset arrays
mod_reduced$z.samples    <- mod$z.samples[keep_ids, , , drop = FALSE]
mod_reduced$psi.samples  <- mod$psi.samples[keep_ids, , , drop = FALSE]

#subset MCMC matrices
mod_reduced$beta.samples  <- as.mcmc(as.matrix(mod$beta.samples)[keep_ids, ])
mod_reduced$alpha.samples <- as.mcmc(as.matrix(mod$alpha.samples)[keep_ids, ])
mod_reduced$theta.samples <- as.mcmc(as.matrix(mod$theta.samples)[keep_ids, ])
mod_reduced$w.samples     <- as.mcmc(as.matrix(mod$w.samples)[keep_ids, ])

if (!is.null(mod_reduced$eta.samples)){
   mod_reduced$eta.samples   <- as.mcmc(as.matrix(mod$eta.samples)[keep_ids, ])
} else {
  message('Skipping eta')
}

if (!is.null(mod_reduced$sigma.sq.psi.samples)){
   mod_reduced$sigma.sq.psi.samples   <- as.mcmc(as.matrix(mod$sigma.sq.psi.samples)[keep_ids, ])
} else {
  message('Skipping siqma.sq.psi')
}

if (!is.null(mod_reduced$sigma.sq.p.samples)){
   mod_reduced$sigma.sq.p.samples   <- as.mcmc(as.matrix(mod$sigma.sq.p.samples)[keep_ids, ])
} else {
  message('Skipping sigma.sq.p')
}

if (!is.null(mod_reduced$beta.star.samples)){
   mod_reduced$beta.star.samples   <- as.mcmc(as.matrix(mod$beta.star.samples)[keep_ids, ])
} else {
  message('Skipping beta.star')
}

if (!is.null(mod_reduced$alpha.star.samples)){
   mod_reduced$alpha.star.samples   <- as.mcmc(as.matrix(mod$alpha.star.samples)[keep_ids, ])
} else {
  message('Skipping alpha.star')
}

#update metadata
mod_reduced$n.post <- n_samples
mod_reduced$n.chains <- 1
mod_reduced$n.samples <- n_samples

#set class and n.posterior samples of new object
class(mod_reduced) <- 'stPGOcc'

#remove full model object to free up space
rm(mod)


## Posterior predictive checks -------------------------------------------------

cat('Starting PPC...\n')

#does the model adequately represent variation in occurrence/detection across space? (group = 1 for sites/rows)
ppc_sites_ft <- ppcOcc(mod_reduced, fit.stat = 'freeman-tukey', group = 1)
ppc_sites_cs <- ppcOcc(mod_reduced, fit.stat = 'chi-squared', group = 1)

  #calcuate bayesian p-values
  bpv_sites_ft <- mean(ppc_sites_ft$fit.y.rep > ppc_sites_ft$fit.y)
  print(bpv_sites_ft)
  bpv_sites_cs <- mean(ppc_sites_cs$fit.y.rep > ppc_sites_cs$fit.y)

  #print summaries
  capture.output(c(print(summary(ppc_sites_ft)), print(summary(ppc_sites_cs))),
                file = file.path(output_dir, paste0(species, '_ppc_sites.txt')))

  #plot diff between replicate and actual data
  sites_ft_df <- data.frame(ppc_sites_ft$fit.y.rep.group.quants[3, , ] - ppc_sites_ft$fit.y.group.quants[3, , ])
  sites_ft_df$site <- rep(1:nrow(sites_ft_df))
  sites_ft_df_long <- sites_ft_df %>% pivot_longer(cols = starts_with('X'), names_to = 'season', values_to = 'diff')
  sites_ft_df_long$site_name <- rep(site_names, length(unique(sites_ft_df_long$season)))
  sites_ft_df_long$season_name <- rep(season_names, length(unique(sites_ft_df_long$site)))
  labeled_sites <- sites_ft_df_long %>% filter(abs(diff) > 2)
  sites_ft_plot <- ggplot(sites_ft_df_long, aes(x = site_name, y = diff)) +
      geom_point(alpha = 0.3) + # aes(color = season_name)
      geom_text(data = labeled_sites, aes(label = paste(site_name, season_name)), size = 4) +
      theme_minimal() +
      labs(y = "Replicate - Observed Discrepancy", x = "Site",
           title = "PPC (groups = site, stat = freeman-tukey)")
  ggsave(paste0(output_dir, '/', species, '_ppc_sites_ft.pdf'), 
         sites_ft_plot, dpi = 300)

  sites_cs_df <- data.frame(ppc_sites_cs$fit.y.rep.group.quants[3, , ] - ppc_sites_cs$fit.y.group.quants[3, , ])
  sites_cs_df$site <- rep(1:nrow(sites_cs_df))
  sites_cs_df_long <- sites_cs_df %>% pivot_longer(cols = starts_with('X'), names_to = 'season', values_to = 'diff')
  sites_cs_df_long$site_name <- rep(site_names, length(unique(sites_cs_df_long$season)))
  sites_cs_df_long$season_name <- rep(season_names, length(unique(sites_cs_df_long$site)))
  labeled_sites <- sites_cs_df_long %>% filter(abs(diff) > 2)
  sites_cs_plot <- ggplot(sites_cs_df_long, aes(x = site_name, y = diff)) +
      geom_point(alpha = 0.3) + # aes(color = season_name)
      geom_text(data = labeled_sites, aes(label = paste(site_name, season_name)), size = 4) +
      theme_minimal() +
      labs(y = "Replicate - Observed Discrepancy", x = "Site",
           title = "PPC (groups = site, stat = chi-squared)")
ggsave(paste0(output_dir, '/', species, '_ppc_sites_cs.pdf'), 
       sites_cs_plot, dpi = 300)

#does the model adequately represent variation in occurrence/detection across time? (group = 2 for columns/replicates -- or years here...?)
ppc_time_ft <- ppcOcc(mod_reduced, fit.stat = 'freeman-tukey', group = 2)
ppc_time_cs <- ppcOcc(mod_reduced, fit.stat = 'chi-squared', group = 2)

  #calculate bayesian p-values
  bpv_time_ft <- mean(ppc_time_ft$fit.y.rep > ppc_time_ft$fit.y)
  bpv_time_cs <- mean(ppc_time_cs$fit.y.rep > ppc_time_cs$fit.y)

  #print summaries
  capture.output(c(print(summary(ppc_time_ft)), print(summary(ppc_time_cs))),
                file = file.path(output_dir, paste0(species, '_ppc_time.txt')))

  #plot diff between replicate and actual data
  time_ft_df <- data.frame(ppc_time_ft$fit.y.rep.group.quants[3, , ] - ppc_time_ft$fit.y.group.quants[3, , ])
  time_ft_df$season <- season_names
  time_ft_df_long <- time_ft_df %>% pivot_longer(cols = starts_with('X'), names_to = 'occasion', values_to = 'diff')
  time_ft_df_long$occasion <- as.numeric(gsub('X','',time_ft_df_long$occasion))
  time_ft_df_long$occasion_name <- paste('wk', time_ft_df_long$occasion, sep = '')
  labeled_occasions <- time_ft_df_long %>% filter(abs(diff) > 2)
  time_ft_plot <- ggplot(time_ft_df_long, aes(x = occasion, y = diff)) +
  geom_point(alpha = 0.3) + # aes(color = season_name)
  geom_text(data = labeled_occasions, aes(label = paste(season, occasion_name, sep = '_')), size = 4) +
  theme_minimal() +
  labs(y = "Replicate - Observed Discrepancy", x = "Occasion",
       title = "PPC (groups = time, stat = freeman-tukey)")
  ggsave(paste0(output_dir, '/', species, '_ppc_time_ft.pdf'), 
         time_ft_plot, dpi = 300)

  time_cs_df <- data.frame(ppc_time_cs$fit.y.rep.group.quants[3, , ] - ppc_time_cs$fit.y.group.quants[3, , ])
  time_cs_df$season <- season_names
  time_cs_df_long <- time_cs_df %>% pivot_longer(cols = starts_with('X'), names_to = 'occasion', values_to = 'diff')
  time_cs_df_long$occasion <- as.numeric(gsub('X','',time_cs_df_long$occasion))
  time_cs_df_long$occasion_name <- paste('wk', time_cs_df_long$occasion, sep = '')
  labeled_occasions <- time_cs_df_long %>% filter(abs(diff) > 2)
  time_cs_plot <- ggplot(time_cs_df_long, aes(x = occasion, y = diff)) +
  geom_point(alpha = 0.3) + # aes(color = season_name)
  geom_text(data = labeled_occasions, aes(label = paste(season, occasion_name, sep = '_')), size = 4) +
  theme_minimal() +
  labs(y = "Replicate - Observed Discrepancy", x = "Occasion",
       title = "PPC (groups = time, stat = chi-squared)")
  ggsave(paste0(output_dir, '/', species, '_ppc_time_cs.pdf'), 
         time_cs_plot, dpi = 300)

#save BPV
bayesian.p.vals <- data.frame('species' = species,
                              'model' = model_name,
                              'n_samples' = n_samples,
                              'group' = c('sites','sites','time','time'),
                              'fit_stat' = c('freeman-tukey','chi-squared','freeman-tukey','chi-squared'),
                              'bpv' = c(bpv_sites_ft, bpv_sites_cs, bpv_time_ft, bpv_time_cs))

write.csv(bayesian.p.vals, paste0(output_dir, '/', species, '_bayesian_pval.csv'))

cat('Done. Results saved to:', output_dir, '\n')
cat('##########################\n')

#clean up
gc()
