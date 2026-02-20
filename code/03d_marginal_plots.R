# 2. Create marginal plots for spOccupancy models

#Based on code from Kovalenko et al. https://github.com/vladkov88/nutcracker/blob/main/code/summary.R

library(tidyverse)
library(spOccupancy)
library(data.table)
library(patchwork)
library(cowplot) #to get legend


## Load data --------------------------------------------------------------------

#get covariate summary values (min/max, mean/sd scaled for N/S)
covar_data <- fread('H_sp_occ/covar_mean_sd_scaledNS.csv')

  #pull out ones I want
  covar_means <- list()
  sort(unique(covar_data$variable))
  covar_means[['elev']] <- covar_data[covar_data$variable == 'elev_mean_3x3',]
  covar_means[['slope']] <- covar_data[covar_data$variable == 'pct_slope_mean_3x3',]
  covar_means[['tree_vol']] <- covar_data[covar_data$variable == 'tree_vol_mean_3x3',]
  covar_means[['boundary']] <- covar_data[covar_data$variable == 'distance_to_perimeter',]
  covar_means[['boundary_log']] <- covar_data[covar_data$variable == 'dist_perimeter_log',]
  covar_means[['dambo']] <- covar_data[covar_data$variable == 'dist_dambo_mean_3x3',]
  covar_means[['dambo_log']] <- covar_data[covar_data$variable == 'dist_dambo_log',]
  covar_means[['rivers']] <- covar_data[covar_data$variable == 'dist_rivers_mean_3x3',]
  covar_means[['rivers_log']] <- covar_data[covar_data$variable == 'dist_rivers_log',]
  
#list model files
file_list <- list.files('H_sp_occ/outputs/model17MSversion/', pattern = 'RDS', full.names = TRUE)

#read one as an example to set dimensions of design matrix
mod_example <- readRDS(file_list[1])

#create an empty design matrix
(n_pred <- 500)                    #number of values to predict on
(n_covar <- dim(mod_example$X)[3]) #number of covariates for model11/model16/model17
(n_post <- mod_example$n.post *    #number of posterior samples
    mod_example$n.chains)

X.template <- matrix(0, n_pred, n_covar)    #design matrix
pred.template <- matrix(NA, n_post, n_pred) #predictions template

  dim(X.template)    #values * covariates
  dim(pred.template) #samples * values

#add parameters as column names
(colnames(X.template) <- colnames(mod_example$beta.samples))

#clean up
rm(mod_example)
rm(covar_data)
  
  
## Write function to predict ---------------------------------------------------

generate_preds <- function(covar_name, species_list) {
  
  #get mean/SD and create values to predict on
  covar_summary <- covar_means[[covar_name]] %>% filter(fence == 'North')
  occ_vals <- seq(from = covar_summary$min, to = covar_summary$max, length.out = n_pred)
  
  #apply log-transformation if necessary, and now get mean/SD of log values for next steps
  if (grepl('boundary|dambo|river', covar_name)){
    covar_summary <- covar_means[[paste0(covar_name, '_log')]] %>% filter(fence == 'North')
    occ_vals <- log(occ_vals)
    message('log transformation performed')
  }
  
  #initialize empty design matrix
  X.pred <- X.template
  
  #set Intercept
  X.pred[, '(Intercept)'] <- 1
  
  #find matching column name (just the main effect, NOT season interaction yet)
  covar_column <- colnames(X.pred)[grepl(covar_name, colnames(X.pred)) & !grepl('season', colnames(X.pred))]
  cat('adding to:', covar_column)
  
  #populate with scaled values for prediction
  X.pred[, covar_column] <- (occ_vals - covar_summary$mean) / covar_summary$sd
  
  #keep 'fence' at North, but if I want South:
  # X.pred[, 'fence_south] <- rep(1, 500)  
  
  #QC
  # head(X.pred) #default: north, 2019, cool, other covar at means
  
  #initialize object to store results
  predictions_df <- NULL
  
  #loop thru species
  for (pp in species_list){
    
    #get model file for this species
    (model_file <- file_list[grepl(pp, file_list)])
    mod <- readRDS(model_file)
    
    #loop thru seasons
    for (ss in c('cool','wet','dry')) {
      
      cat('\nPredicting', covar_name, 'for:', pp, ss)
      
      #populate corresponding season and interaction columns
      
      #cool: keep the same
      if(ss == 'cool') {X.pred.ss <- X.pred}
      
      #wet:
      colname_wet <- colnames(X.pred.ss)[grepl(covar_name, colnames(X.pred.ss)) & grepl('wet', colnames(X.pred.ss))]
      if(ss == 'wet') {X.pred.ss <- X.pred}
      if(ss == 'wet') {X.pred.ss[, 'season_wet'] <- rep(1, n_pred)}
      if(ss == 'wet') {X.pred.ss[, colname_wet] <- X.pred[, covar_column]}
      
      colname_dry <- colnames(X.pred.ss)[grepl(covar_name, colnames(X.pred.ss)) & grepl('dry', colnames(X.pred.ss))]
      if(ss == 'dry') {X.pred.ss <- X.pred}
      if(ss == 'dry') {X.pred.ss[, 'season_dry'] <- rep(1, n_pred)}
      if(ss == 'dry') {X.pred.ss[, colname_dry] <- X.pred[, covar_column]}
      
      #qc
      # head(X.pred.ss)
      
      #reset template each time
      psi.pred <- pred.template 
      
      #predict: multiply design matrix by each row of beta posterior samples
      for (i in 1:n_post) {
        psi.pred[i, ] <- plogis(X.pred.ss %*% as.matrix(mod$beta.samples[i, ])) 
      }  
      
      #convert to dataframe
      preds_df <- data.frame(mean = apply(psi.pred, 2, mean),
                             lci = apply(psi.pred, 2, quantile, 0.025),
                             uci = apply(psi.pred, 2, quantile, 0.975),
                             covar_values = occ_vals,
                             species = pp,
                             season = ss,
                             # covar = gsub('_scaledNS', '', covar_name),
                             covar = covar_name)
      
      #store
      predictions_df <- rbind(predictions_df, preds_df)
    }
  }
  
  return(predictions_df)
  
}


## Apply function for each covariate -------------------------------------------

#list all species [or run function separately for covariates and input only ones w/ sig. effects]
species <- c('buffalo','eland','elephant','impala','kudu','sable','warthog','waterbuck','zebra')
names(covar_means)

#predict (~15 min)
all_predictions <- list()
for (cc in c('elev','slope','tree_vol','boundary','dambo','rivers')){
  preds_cc <- generate_preds(covar_name = cc, species_list = species)
  all_predictions[[cc]] <- preds_cc
  write.csv(preds_cc, paste0('H_sp_occ/figures/marginal_plots/marginal_data/model17ms_', cc, '.csv'))
}

str(all_predictions)
saveRDS(all_predictions, 'H_sp_occ/figures/marginal_plots/marginal_data/model17ms_all_predictions.RDS')

  
## Plot ------------------------------------------------------------------------

#set colors
color_palette <- c(dry = "#E69F00", wet = "#56B4E9", cool = "#009E73") 
  
#read input file to get actual values for rugs (get raw, not scaled, values)
inputs <- readRDS('H_sp_occ/sp_inputs.RDS')
rug_values <- list(elev = inputs$aardvark$occ.covs$elev,
                   slope_rug = inputs$aardvark$occ.covs$slope,
                   treeVol_rug = inputs$aardvark$occ.covs$tree_vol,
                   distBoundary_rug = inputs$aardvark$occ.covs$dist_boundary,
                   distDambos_rug = inputs$aardvark$occ.covs$dist_dambo,
                   distRivers_rug = inputs$aardvark$occ.covs$dist_rivers)
  
#get effects table so I can find ones with >0.9 probability
all_effects <- fread('H_sp_occ/results/model17_MSversion/covariate_effects_all_species.csv')


### ELEVATION (eland wet, sable wet, waterbuck wet) ----------------------------
  
  #read in OR select
  # preds_elevation <- fread('H_sp_occ/figures/marginal_plots/marginal_data/model17')
  preds_elevation <- all_predictions[['elev']]  

  #select seasons I actually want to plot for each species
  all_effects %>% filter(variable == 'elev_scaledNS',
                         prob_positive >= 0.9 | prob_negative >= 0.9) %>% pull(season, species)

  elev_df <- preds_elevation %>% 
                filter(
                  (species == 'eland' & season == 'wet') |
                    (species == 'sable' & season == 'wet') |
                    (species == 'zebra' & season == 'cool')
                  )
  elev_df %>% select(species, season, covar) %>% ftable()
  
  #get rug values
  elev_rug <- data.frame('values' = rug_values$elev)
  
  #plot
  plot_elev <- ggplot(elev_df, aes(x = covar_values, y = mean, color = season, fill = season)) +
    # geom_ribbon(aes(ymin = lci, ymax = uci), alpha = 0.4, color = NA) +
    geom_ribbon(aes(ymin = lci, ymax = uci), fill = NA, alpha = 0.4, linetype = 'dashed') +
    geom_line(linewidth = 1) +
    scale_color_manual(values = color_palette, drop = FALSE) +
    scale_fill_manual(values = color_palette, drop = FALSE) +
    scale_y_continuous(limits = c(0, 1)) +
    facet_wrap(~species, scales = 'free', nrow = 1) +
    geom_rug(data = elev_rug, mapping = aes(x = values), inherit.aes = FALSE) +
    labs(x = 'Elevation (m)', #unique(plot_df$covar), 
         y = 'ψ ± 95% CI', 
         # title = unique(plot_df$species),
         # subtitle = unique(plot_df$season)
         ) + 
    theme_bw(base_size = 12) + theme(strip.background = element_rect(fill = 'transparent'))
  plot_elev

  # ggsave(paste0('H_sp_occ/figures/marginal_plots/', 'model17_', 'elevation', '.tif'),
  #        plot_elev, width = 6, height = 3) #3 cols, 1 row
      #width=2 per panel (column)
      #height=3 per panel (row)


### SLOPE (eland wet, elephant wet/dry, kudu cool, warthog all, waterbuck all) -------------------------------
  
  #read in OR select
  # preds_elevation <- fread('H_sp_occ/figures/marginal_plots/marginal_data/model17')
  preds_slope <- all_predictions[['slope']]  
  
  #select seasons I actually want to plot for each species
  all_effects %>% filter(variable == 'slope_scaledNS',
                         prob_positive >= 0.9 | prob_negative >= 0.9) %>% pull(season, species)
  
  #select seasons I actually want to plot for each species
  slope_df <- preds_slope %>% 
    filter(
      (species == 'buffalo' & season == 'wet') |
        (species == 'eland' & season == 'wet') |
        (species == 'elephant' & season == 'wet') |
        (species == 'elephant' & season == 'dry') |
        (species == 'kudu' & season == 'cool') |
        (species == 'warthog' & season == 'wet') |
        (species == 'warthog' & season == 'dry') |
        (species == 'waterbuck' & season == 'cool') |
        (species == 'waterbuck' & season == 'wet')
  )
  slope_df %>% select(species, season, covar) %>% ftable()
  
  #get rug values
  slope_rug <- data.frame('values' = rug_values$slope_rug)
  
  #plot
  plot_slope <- ggplot(slope_df, aes(x = covar_values, y = mean, color = season, fill = season)) +
    # geom_ribbon(aes(ymin = lci, ymax = uci), color = NA, alpha = 0.2) +
    geom_ribbon(aes(ymin = lci, ymax = uci), fill = NA, alpha = 0.4, linetype = 'dashed') +
    geom_line(linewidth = 1) +
    scale_color_manual(values = color_palette, drop = FALSE) +
    scale_fill_manual(values = color_palette, drop = FALSE) +
    scale_y_continuous(limits = c(0, 1)) +
    facet_wrap(~species, scales = 'free', nrow = 1) +
    geom_rug(data = slope_rug, mapping = aes(x = values), inherit.aes = FALSE) +
    labs(x = 'Slope (degrees)', #unique(plot_df$covar), 
         y = 'ψ ± 95% CI', 
         # title = unique(plot_df$species),
         # subtitle = unique(plot_df$season)
    ) + 
    theme_bw(base_size = 12) + theme(strip.background = element_rect(fill = 'transparent'))
  plot_slope
  
  # ggsave(paste0('H_sp_occ/figures/marginal_plots/', 'model16_', 'slope', '.tif'),
  #        plot_slope, width = 10, height = 6)
  #   #5 cols, 3 rows
  
  
### TREEVOL (elephant dry, kudu wet) -------------------------------
  
  #read in OR select
  # preds_elevation <- fread('H_sp_occ/figures/marginal_plots/marginal_data/model17')
  preds_treeVol <- all_predictions[['tree_vol']]  
  
  #select seasons I actually want to plot for each species
  all_effects %>% filter(variable == 'tree_vol_scaledNS',
                         prob_positive >= 0.9 | prob_negative >= 0.9) %>% pull(season, species)
  
  #select seasons I actually want to plot for each species
  treeVol_df <- preds_treeVol %>% 
    filter(
      (species == 'elephant' & season == 'dry') |
        (species == 'kudu' & season == 'wet') |
        (species == 'warthog' & season == 'cool') |
        # (species == 'zebra' & season == 'wet')
        (species == 'eland' & season == 'dry')
      
    )
  treeVol_df %>% select(species, season, covar) %>% ftable()
  
  #get rug values
  treeVol_rug <- data.frame('values' = rug_values$treeVol_rug)
  
  #plot
  plot_treeVol <- ggplot(treeVol_df, aes(x = covar_values, y = mean, color = season, fill = season)) +
    # geom_ribbon(aes(ymin = lci, ymax = uci), color = NA, alpha = 0.2) +
    geom_ribbon(aes(ymin = lci, ymax = uci), fill = NA, alpha = 0.4, linetype = 'dashed') +
    geom_line(linewidth = 1) +
    scale_color_manual(values = color_palette, drop = FALSE) +
    scale_fill_manual(values = color_palette, drop = FALSE) +
    scale_y_continuous(limits = c(0, 1)) +
    facet_wrap(~species, scales = 'free', nrow = 1) +
    geom_rug(data = treeVol_rug, mapping = aes(x = values), inherit.aes = FALSE) +
    labs(x = 'Tree volume (trees/ha^3)', #unique(plot_df$covar), 
         y = 'ψ ± 95% CI', 
         # title = unique(plot_df$species),
         # subtitle = unique(plot_df$season)
    ) + 
    theme_bw(base_size = 12) + theme(strip.background = element_rect(fill = 'transparent'))
  plot_treeVol
  
  # ggsave(paste0('H_sp_occ/figures/marginal_plots/', 'model16_', 'treeVol', '.tif'),
  #        plot_treeVol, width = 6, height = 3)
  #     #2 cols, 2 rows
  
  
### BOUNDARY (eland cool/dry, elephant dry, sable cool/wet, warthog all) -------------------------------
  
  #read in OR select
  # preds_elevation <- fread('H_sp_occ/figures/marginal_plots/marginal_data/model17')
  preds_distBoundary <- all_predictions[['boundary']]  
  
  #select seasons I actually want to plot for each species
  all_effects %>% filter(variable == 'dist_boundary_log_scaledNS',
                         prob_positive >= 0.9 | prob_negative >= 0.9) %>% pull(season, species)
  
  #select seasons I actually want to plot for each species
  boundary_df <- preds_distBoundary %>% 
    filter(
      (species == 'eland' & season == 'cool') |
        (species == 'eland' & season == 'dry') |
        # (species == 'elephant' & season == 'dry') |
        (species == 'sable' & season == 'cool') |
        (species == 'sable' & season == 'wet') |
        (species == 'warthog' & season == 'cool') |
        (species == 'warthog' & season == 'dry') |
        (species == 'zebra' & season == 'cool') |
        (species == 'zebra' & season == 'wet') |
        (species == 'zebra' & season == 'dry')
    )
  boundary_df %>% select(species, season, covar) %>% ftable()
  
  #get rug values
  boundary_rug <- data.frame('values' = rug_values$distBoundary_rug)
  
  #back-transform log values
  boundary_df <- boundary_df %>% mutate(covar_values_bt = exp(covar_values))
  
  #plot
  plot_distBoundary <- ggplot(boundary_df, aes(x = covar_values_bt/1000, y = mean, color = season, fill = season)) +
    # geom_ribbon(aes(ymin = lci, ymax = uci), color = NA, alpha = 0.2) +
    geom_ribbon(aes(ymin = lci, ymax = uci), fill = NA, alpha = 0.4, linetype = 'dashed') +
    geom_line(linewidth = 1) +
    scale_color_manual(values = color_palette, drop = FALSE) +
    scale_fill_manual(values = color_palette, drop = FALSE) +
    scale_y_continuous(limits = c(0, 1)) +
    facet_wrap(~species, scales = 'free', nrow = 1) +
    geom_rug(data = boundary_rug, mapping = aes(x = values/1000), inherit.aes = FALSE) +
    labs(x = 'Distance to boundary (km)', #unique(plot_df$covar), 
         y = 'ψ ± 95% CI', 
         # title = unique(plot_df$species),
         # subtitle = unique(plot_df$season)
    ) + 
    theme_bw(base_size = 12) + theme(strip.background = element_rect(fill = 'transparent'))
  plot_distBoundary
  
  # ggsave(paste0('H_sp_occ/figures/marginal_plots/', 'model16_', 'distBoundary', '.tif'),
  #        plot_distBoundary, width = 8, height = 6)
      #4 cols, 3 rows
  

### DAMBOS (eland dry, elephant dry, waterbuck cool/dry) -------------------------------
  
  #read in OR select
  # preds_elevation <- fread('H_sp_occ/figures/marginal_plots/marginal_data/model17')
  preds_distDambos <- all_predictions[['dambo']]
  
  #select seasons I actually want to plot for each species
  all_effects %>% filter(variable == 'dist_dambo_log_scaledNS',
                         prob_positive >= 0.9 | prob_negative >= 0.9) %>% pull(season, species)
  
  #select seasons I actually want to plot for each species
  dambos_df <- preds_distDambos %>% 
    filter(
      (species == 'buffalo' & season == 'cool') |
        (species == 'buffalo' & season == 'wet') |
        (species == 'buffalo' & season == 'dry') |
        (species == 'eland' & season == 'dry') |
        (species == 'elephant' & season == 'cool') |
        (species == 'elephant' & season == 'dry') |
        (species == 'warthog' & season == 'cool') |
        (species == 'waterbuck' & season == 'cool') |
        (species == 'zebra' & season == 'cool')
    )
  dambos_df %>% select(species, season, covar) %>% ftable()
  
  #get rug values
  dambos_rug <- data.frame('values' = rug_values$distDambos_rug)
  
  #back-transform log values
  dambos_df <- dambos_df %>% mutate(covar_values_bt = exp(covar_values))
  
  #plot
  plot_distDambos <- ggplot(dambos_df, aes(x = covar_values_bt/1000, y = mean, color = season, fill = season)) +
    # geom_ribbon(aes(ymin = lci, ymax = uci), color = NA, alpha = 0.2) +
    geom_ribbon(aes(ymin = lci, ymax = uci), fill = NA, alpha = 0.4, linetype = 'dashed') +
    geom_line(linewidth = 1) +
    scale_color_manual(values = color_palette, drop = FALSE) +
    scale_fill_manual(values = color_palette, drop = FALSE) +
    scale_y_continuous(limits = c(0, 1)) +
    facet_wrap(~species, scales = 'free', nrow = 1) +
    geom_rug(data = dambos_rug, mapping = aes(x = values/1000), inherit.aes = FALSE) +
    labs(x = 'Distance to dambos (km)', #unique(plot_df$covar), 
         y = 'ψ ± 95% CI', 
         # title = unique(plot_df$species),
         # subtitle = unique(plot_df$season)
    ) + 
    theme_bw(base_size = 12) + theme(strip.background = element_rect(fill = 'transparent'))
  plot_distDambos
  
  # ggsave(paste0('H_sp_occ/figures/marginal_plots/', 'model16_', 'distDambos', '.tif'),
  #        plot_distDambos, width = 10, height = 3)
  #3 cols, 3 rows
  
  
### RIVERS (eland wet/dry, elephant dry) -------------------------------
  
  #read in OR select
  # preds_elevation <- fread('H_sp_occ/figures/marginal_plots/marginal_data/model17')
  preds_distRivers <- all_predictions[['rivers']]
  
  #select seasons I actually want to plot for each species
  all_effects %>% filter(variable == 'dist_rivers_log_scaledNS',
                         prob_positive >= 0.9 | prob_negative >= 0.9) %>% pull(season, species)
  
  #select seasons I actually want to plot for each species
  rivers_df <- preds_distRivers %>% 
    filter(
      (species == 'elephant' & season == 'dry') |
        (species == 'waterbuck' & season == 'dry')
    )
  rivers_df %>% select(species, season, covar) %>% ftable()
  
  #get rug values
  rivers_rug <- data.frame('values' = rug_values$distRivers_rug)
  
  #back-transform log values
  rivers_df <- rivers_df %>% mutate(covar_values_bt = exp(covar_values))
  
  #plot
  plot_distRivers <- ggplot(rivers_df, aes(x = covar_values_bt/1000, y = mean, color = season, fill = season)) +
    # geom_ribbon(aes(ymin = lci, ymax = uci), color = NA, alpha = 0.2) +
    geom_ribbon(aes(ymin = lci, ymax = uci), fill = NA, alpha = 0.4, linetype = 'dashed') +
    geom_line(linewidth = 1) +
    scale_color_manual(values = color_palette, drop = FALSE) +
    scale_fill_manual(values = color_palette, drop = FALSE) +
    scale_y_continuous(limits = c(0, 1)) +
    facet_wrap(~species, scales = 'free', nrow = 1) +
    geom_rug(data = rivers_rug, mapping = aes(x = values/1000), inherit.aes = FALSE) +
    labs(x = 'Distance to rivers (km)', #unique(plot_df$covar), 
         y = 'ψ ± 95% CI', 
         # title = unique(plot_df$species),
         # subtitle = unique(plot_df$season)
    ) + 
    theme_bw(base_size = 12) + theme(strip.background = element_rect(fill = 'transparent'))
  plot_distRivers
  
  # ggsave(paste0('H_sp_occ/figures/marginal_plots/', 'model16_', 'distRivers', '.tif'),
  #        plot_distRivers, width = 6, height = 3)
  # #2 cols, 2 rows
  
  
## COMBINE PLOTS ---------------------------------------------------------------
  
row1 <- (plot_treeVol + theme(legend.position = 'none', 
                                   strip.background = element_rect(color = NA, linewidth = 2), #color = 'black'
                                   strip.text = element_text(face = 'bold', size = 12),
                                   plot.background = element_rect(fill = NA, color = "black", linewidth = 0.5)) +
           ggtitle('A. Tree Volume')) +
  # plot_spacer() +
  plot_layout(widths = c(4,2))

row2 <- (plot_distBoundary + theme(legend.position = 'none', 
                           strip.background = element_rect(color = NA, linewidth = 2), #darkorange
                           strip.text = element_text(face = 'bold', size = 12),
                           plot.background = element_rect(fill = "grey95", color = "black", linewidth = 0.5)) +
           ggtitle('B. Distance to Boundary')) +
  plot_layout(widths = c(4, 2))
row2

row3 <- (plot_elev + theme(legend.position = 'none', 
                            strip.background = element_rect(color = NA, linewidth = 2), #purple2
                            strip.text = element_text(face = 'bold', size = 12),
                            plot.background = element_rect(fill = NA, color = "black", linewidth = 0.5)) +
           ggtitle('C. Elevation')) +
        (plot_distRivers + theme(legend.position = 'none', 
                           strip.background = element_rect(color = NA, linewidth = 2), #purple2
                           strip.text = element_text(face = 'bold', size = 12),
                           plot.background = element_rect(fill = "grey95", color = "black", linewidth = 0.5)) +
           ggtitle('D. Distance to Rivers')) +
  plot_layout(widths = c(3,2,1), 
              axis_titles = 'collect')

row4 <- (plot_distDambos + theme(legend.position = 'none', 
                                 strip.background = element_rect(color = NA, linewidth = 2), #'lightblue'
                                 strip.text = element_text(face = 'bold', size = 12),
                                 plot.background = element_rect(fill = "grey95", color = "black", linewidth = 0.5)) +
           ggtitle('E. Distance to Dambos')) + 
  plot_layout(widths = c(6))   

row5 <- (plot_slope + theme(legend.position = 'none', 
                                 strip.background = element_rect(color = NA, linewidth = 2), #'lightblue'
                                 strip.text = element_text(face = 'bold', size = 12),
                                 plot.background = element_rect(fill = NA, color = "black", linewidth = 0.5)) +
           ggtitle('F. Slope')) + 
  plot_layout(widths = c(6)) 

# #modify legend using an example plot  
# legend <- get_legend(plot_distBoundary)
plot_distBoundary2 <- plot_distBoundary +
  theme(legend.text = element_text(size = 16),
        legend.direction = "horizontal",
        legend.box.spacing = unit(0.5, "cm"),
        # legend.position = 'bottom',
        legend.title.position = 'top',
        legend.justification = "center",
        legend.title = element_text(size = 18, hjust = 0.5)
  )
legend <- get_legend(plot_distBoundary2)

combined_marginals <- ((row1 + legend) / (row2 + plot_spacer()) / (row3 + plot_spacer()) / row4 / row5)
combined_marginals

#Save  
ggsave('H_sp_occ/figures/marginal_plots/model17ms_allMarginals.tif',
       combined_marginals,
       width = 10,
       height = 10)

# Another way:
# layout2 <- "
# AAABBB
# CCCCCC
# EEEEEE
# FFFFFF
# GGGGGG
# "
# 
# row1b <- (plot_elev + theme(legend.position = 'none', 
#                             strip.background = element_rect(color = NA, linewidth = 2), #purple2
#                             strip.text = element_text(face = 'bold', size = 12),
#                             plot.background = element_rect(fill = NA, color = "black", linewidth = 0.5)) +
#             ggtitle('A. Elevation')) +
#   plot_layout(widths = c(3,3), 
#               axis_titles = 'collect')
# 
# row2b <- (plot_distBoundary + theme(legend.position = 'none', 
#                                    strip.background = element_rect(color = NA, linewidth = 2), #darkorange
#                                    strip.text = element_text(face = 'bold', size = 12),
#                                    plot.background = element_rect(fill = "grey95", color = "black", linewidth = 0.5)) +
#            ggtitle('B. Distance to Boundary')) +
#   plot_layout(widths = c(4, 2))
# 
# row3b <- (plot_treeVol + theme(legend.position = 'none', 
#                                strip.background = element_rect(color = NA, linewidth = 2), #color = 'black'
#                                strip.text = element_text(face = 'bold', size = 12),
#                                plot.background = element_rect(fill = NA, color = "black", linewidth = 0.5)) +
#             ggtitle('C. Tree Volume')) +
#   (plot_distRivers + theme(legend.position = 'none', 
#                            strip.background = element_rect(color = NA, linewidth = 2), #purple2
#                            strip.text = element_text(face = 'bold', size = 12),
#                            plot.background = element_rect(fill = "grey95", color = "black", linewidth = 0.5)) +
#      ggtitle('D. Distance to Rivers')) +
#   plot_layout(widths = c(4,2), 
#               axis_titles = 'collect')
# 
# #modify legend using an example plot  
# plot_distBoundary2 <- plot_distBoundary + 
#   theme(legend.text = element_text(size = 16), 
#         legend.direction = "horizontal",
#         legend.box.spacing = unit(0.5, "cm"),
#         # legend.position = 'bottom',
#         legend.title.position = 'top',
#         legend.justification = "center",    
#         legend.title = element_text(size = 18, hjust = 0.5)
#   )
# legend2 <- get_legend(plot_distBoundary2)
# 
# #combine
# combined3 <- (row1b + legend2) + row2 + row3b + row4 + row5 + plot_layout(design = layout2)
# combined3
# 
# ggsave('H_sp_occ/figures/marginal_plots/model17ms_allMarginals_v2.tif',
#        combined3,
#        width = 10,
#        height = 10)


## PLOT ANOTHER WAY ------------------------------------------------------------
  
files <- list.files('H_sp_occ/figures/marginal_plots/marginal_data/', pattern = 'csv', full.names = TRUE)
file_list <- lapply(files, fread)

preds_df <- rbindlist(file_list)  

preds_df %>% select(species, season, covar) %>% ftable()

#remove any if I want to, with 'select'

#get rug values -- though this isn't working yet
df_rug <- imap_dfr(rug_values, ~ tibble(name = .y, values = .x))
  
#plot version 1
plot1 <- ggplot(preds_df, aes(x = covar_values, y = mean, color = species, linetype = season)) +
  # geom_ribbon(aes(ymin = lci, ymax = uci), fill = NA, linetype = 'dotted', alpha = 0.4) +
  geom_line(linewidth = 1) +    
  # scale_color_manual(values = color_palette, drop = FALSE) +
  # scale_fill_manual(values = color_palette, drop = FALSE) +
  scale_y_continuous(limits = c(0, 1)) +
  facet_wrap(~covar ~season, scales = 'free') +
  # geom_rug(data = df_rug, mapping = aes(x = values), inherit.aes = FALSE, color = 'gray40') +
  labs(
       y = 'Occupancy Probability (ψ ± 95% CI)', 
       # title = unique(plot_df$species),
       # subtitle = unique(plot_df$season)
  ) + 
  theme_bw(base_size = 12) + theme(strip.background = element_rect(fill = 'transparent'))
plot1





