## 1. Assess effects of environmental covariates from spOccupancy models

library(ggplot2)
library(data.table)
library(dplyr)
library(tidyr)


## Write function to calculate derived params ----------------------------------
calculate_params <- function(beta_samples) {
  
  # Define the variables and seasons
  continuous_vars <- c("elev_scaledNS", "slope_scaledNS", "tree_vol_scaledNS", 
                       "dist_boundary_log_scaledNS", "dist_dambo_log_scaledNS", "dist_rivers_log_scaledNS")
  seasons <- c("cool", "wet", "dry")
  
  # Initialize results data frame
  derived_params <- data.frame()
  
  for(var in continuous_vars) {
    for(season in seasons) {
      
      # Calculate the derived parameter for this variable-season combination
      if(season == "cool") {
        # For cool season, just the main effect
        effect_samples <- beta_samples[, var]
        param_name <- paste0(var, "/", season)
        
      } else if(season == "wet") {
        # For wet season, main effect + interaction
        main_effect <- beta_samples[, var]
        interaction_col <- paste0("season_wet:", var)
        
        if(interaction_col %in% colnames(beta_samples)) {
          interaction_effect <- beta_samples[, interaction_col]
          effect_samples <- main_effect + interaction_effect
        } else {
          effect_samples <- main_effect  # If no interaction column
        }
        param_name <- paste0(var, "/", season)
        
      } else if(season == "dry") {
        # For dry season, main effect + interaction
        main_effect <- beta_samples[, var]
        interaction_col <- paste0("season_dry:", var)
        
        if(interaction_col %in% colnames(beta_samples)) {
          interaction_effect <- beta_samples[, interaction_col]
          effect_samples <- main_effect + interaction_effect
        } else {
          effect_samples <- main_effect  # If no interaction column
        }
        param_name <- paste0(var, "/", season)
      }
      
      # Calculate summary statistics
      quantiles <- quantile(effect_samples, probs = c(0.025, 0.5, 0.975))
      prob_positive <- mean(effect_samples > 0)
      prob_negative <- mean(effect_samples < 0)
      
      # Store results
      derived_params <- rbind(derived_params, data.frame(
        parameter = param_name,
        variable = var,
        season = season,
        q025 = quantiles[1],
        median = quantiles[2],
        q975 = quantiles[3],
        mean = mean(effect_samples),
        sd = sd(effect_samples),
        prob_positive = prob_positive,
        prob_negative = prob_negative,
        # significant_positive = prob_positive > 0.95,
        # significant_negative = prob_negative > 0.95
        significant_positive = quantiles[1] > 0,  # 95% CI entirely above 0
        significant_negative = quantiles[3] < 0   # 95% CI entirely below 0
      ))
    }
  }
  
  # Clean up row names
  rownames(derived_params) <- NULL
  
  return(derived_params)
}


## Write function to print results file ----------------------------------------
print_params <- function(derived_params, species_name) {
  
  cat("=== DERIVED PARAMETER SUMMARY ===\n\n")
  cat("===", species_name, "===\n\n")
  cat("Effects are on the logit scale (occupancy model linear predictor)\n\n")
  
  # Group by variable for cleaner output
  vars <- unique(derived_params$variable)
  var_labels <- c(
    "elev_scaledNS" = "Elevation",
    "slope_scaledNS" = "Slope", 
    "tree_vol_scaledNS" = "Tree Volume",
    "dist_boundary_log_scaledNS" = "Distance to Boundary",
    "dist_dambo_log_scaledNS" = "Distance to Dambo",
    "dist_rivers_log_scaledNS" = "Distance to Rivers"
  )
  
  for(var in vars) {
    var_data <- derived_params[derived_params$variable == var, ]
    
    cat(sprintf("--- %s ---\n", var_labels[var]))
    
    for(i in 1:nrow(var_data)) {
      row <- var_data[i, ]
      significance <- ""
      if(row$significant_positive) significance <- "[95%CI ≠ 0]"
      if(row$significant_negative) significance <- "[95%CI ≠ 0]"
      
      f_score <- ""
      if(row$median >= 0) f_score <- row$prob_positive
      if(row$median < 0) f_score <- row$prob_negative
      
      f_dir <- ""
      if(row$median >= 0) f_dir <- "Pr effect is POS +"
      if(row$median < 0) f_dir <- "Pr effect is NEG -"
      
      f_level <- ""
      if((abs(row$prob_positive) >= 0.7) | (abs(row$prob_negative) >= 0.7)) f_level <- "*WEAK"
      if((abs(row$prob_positive) >= 0.8) | (abs(row$prob_negative) >= 0.8)) f_level <- "**MODERATE"
      if((abs(row$prob_positive) >= 0.9) | (abs(row$prob_negative) >= 0.9)) f_level <- "***STRONG"
      
      # cat(sprintf("  %s: %.3f [%.3f, %.3f] (P(>0)=%.2f)%s\n", 
      #             toupper(row$season), 
      #             row$median, row$q025, row$q975, 
      #             row$prob_positive, significance))
      
      cat(sprintf("  %s: %.2f %s %s %.3f [%.2f, %.2f] %s %s %s\n",
                  toupper(row$season), 
                  f_score, f_dir, '  |  ',
                  row$median, row$q025, row$q975, 
                  '  |  ', f_level, significance))
      
    }
    cat("\n")
  }
}


## Write function to create results table (S2 TABLE 1) -------------------------
create_table <- function(combined_results) {
  
  # Create nice variable labels
  var_labels <- c(
    "elev_scaledNS" = "Elevation",
    "slope_scaledNS" = "Slope", 
    "tree_vol_scaledNS" = "Tree Volume",
    "dist_boundary_log_scaledNS" = "Distance to Boundary",
    "dist_dambo_log_scaledNS" = "Distance to Dambo",
    "dist_rivers_log_scaledNS" = "Distance to Rivers"
  )
  
  # Create significance indicators based on probability thresholds
  combined_results$significance_stars <- case_when(
    # For positive effects
    combined_results$prob_positive >= 0.9 ~ " **",
    combined_results$prob_positive >= 0.8 ~ " *",
    # For negative effects  
    combined_results$prob_negative >= 0.9 ~ " **",
    combined_results$prob_negative >= 0.8 ~ " *",
    # No significance
    TRUE ~ ""
  )
  
  # Create formatted effect sizes with CIs and significance
  combined_results$effect_with_sig <- sprintf("%.2f\n[%.2f, %.2f]%s", 
                                              combined_results$median, 
                                              combined_results$q025, 
                                              combined_results$q975,
                                              combined_results$significance_stars)
  
  # Create species-season identifier
  combined_results$species_season <- paste(combined_results$species, toupper(combined_results$season))
  
  # Replace variable names with labels
  combined_results$variable_label <- var_labels[combined_results$variable]
  
  # Reshape to wide format
  summary_table <- combined_results %>%
    select(species, season, species_season, variable_label, effect_with_sig) %>%
    pivot_wider(names_from = variable_label, 
                values_from = effect_with_sig) %>%
    arrange(species, match(season, c("cool", "wet", "dry"))) %>%
    select(-season)
  
  return(summary_table)
}


## Calculate derived params ----------------------------------------------------

files <- list.files('H_sp_occ/outputs/model17MSversion/', pattern = 'RDS', full.names = TRUE)

#keep species for which this model converged
# files11_keep <- files11[!grepl('buffalo|impala|zebra',files11)]
files_keep <- files

#object to store results
derived_params <- list()
seasonal_params <- NULL

#loop thru model files (species)
for (ff in files_keep){
  
  #extract species name from file name
  species_name <- gsub('H_sp_occ/outputs/model17MSversion//model17_', '', ff)
  species_name <- gsub('_allChains.RDS', '', species_name)
  
  cat('## Predicting for', species_name, '\n')
  
  mod <- readRDS(ff)
  
  #get beta samples and calculate derived parameters
  beta_samples_ff <- mod$beta.samples
  params_ff <- calculate_params(beta_samples_ff)
  
  params_ff$species <- species_name
  derived_params[[species_name]] <- params_ff
  
  #generate summary and save to file
  # capture.output(print_params(params_ff, species_name),
  #                file = file.path('H_sp_occ/results/model17/', paste0(species_name, '_covariate_effects.txt')))
  
  #derive seasonal effects (reference level is Cool)
  dry_effect <- beta_samples_ff[,'season_dry']
  wet_effect <- beta_samples_ff[,'season_wet']
  
  dry_df <- c(quantile(dry_effect, c(0.025, 0.5, 0.975)), 
              'prob_neg' = mean(dry_effect < 0),
              'prob_pos' = mean(dry_effect > 0))
  wet_df <- c(quantile(wet_effect, c(0.025, 0.5, 0.975)), 
              'prob_neg' = mean(wet_effect < 0),
              'prob_pos' = mean(wet_effect > 0))
  
  seasonal_df <- bind_rows(dry = dry_df, wet = wet_df, .id = 'season') %>% 
    mutate(across(where(is.numeric), ~ round(.x, 4))) %>%
    mutate(species = species_name,
           f_score = ifelse(`50%`>0, prob_pos, prob_neg))
  
  seasonal_params <- rbind(seasonal_params, seasonal_df)
  
}

all_derived <- rbindlist(derived_params)
write.csv(all_derived, 'H_sp_occ/results/model17_MSversion/covariate_effects_all_species.csv')

#create summary table for Appendix S3 Table 1
summary_table <- create_table(all_derived)
write.csv(summary_table, 'H_sp_occ/results/model17_MSversion/formatted_covariate_table.csv')


#review seasonal effects for results
seasonal_params %>% filter(f_score >= 0.9)


## Parameter plot for FIG 5 ----------------------------------------------------

color_palette <- c(dry = "#E69F00", wet = "#56B4E9", cool = "#009E73") 

#read in if necessary
all_derived <- fread('H_sp_occ/results/model17_MSversion/covariate_effects_all_species.csv')
effects <- all_derived
effects$variable <- gsub('_scaledNS', '', effects$variable)

#reorder variables for facets and change some names
effects$variable <- factor(effects$variable,
                           levels = c('elev', 'slope', 'tree_vol',
                                      'dist_boundary_log', 'dist_dambo_log', 'dist_rivers_log'),
                           labels = c('Elevation', 'Slope', 'Tree Volume',
                                      'Distance to Boundary', 'Distance to Dambo', 'Distance to Rivers'))

color_palette <- 
  c('cool' = '#009E73', 'wet' = '#56B4E9', 'dry' = '#E69F00')

#note f-score >0.9
effects <- effects %>%
  mutate(high_prob = ifelse(prob_negative >= 0.9 | prob_positive >= 0.9, TRUE, FALSE))

#plot
effects_plot <- ggplot(effects, aes(species, y = median, color = season)) +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  geom_pointrange(aes(ymin = q025, ymax = q975, 
                      fill = ifelse(high_prob, 'white', as.character(season)),
                      shape = high_prob,
                      linetype = high_prob), 
                  size = 0.8,
                  position = position_dodge(width = 0.6)) +
  # scale_shape_manual(values = c(21, 24)) +
  scale_shape_manual(values = c('TRUE' = 24, 
                                'FALSE' = 21), 
                     labels = c('TRUE' = 'Probability ≥0.9', 
                                'FALSE' = 'Probability <0.9'),
                     drop = FALSE) +
  scale_color_manual(values = color_palette, drop = FALSE) +
  scale_fill_manual(values = c('white' = 'white', 'cool' = '#009E73', 'wet' = '#56B4E9', 'dry' = '#E69F00'), 
                    drop = FALSE) +
  # scale_fill_manual(values = color_palette, drop = FALSE) +
  facet_wrap(~variable, scales = 'free', nrow = 3) +
  ylim(c(-6,7)) +
  labs(y = 'Effect size (coefficient median ± 95% CI)',
       color = 'Season',
       fill = 'Probability ≥0.9',
       shape = 'Probability ≥0.9',
       linetype = 'Probability ≥0.9') +
  guides(linetype = 'none',
         fill = 'none') +
  theme_bw(base_size = 14) + theme(axis.text.x = element_text(angle = 45, hjust = 1),
                                   axis.title.x = element_blank(),
                                   legend.position = 'top',
                                   legend.title = element_blank(),
                                   legend.text = element_text(size = 14),
                                   strip.background = element_rect(fill = NA))

effects_plot

ggsave('H_sp_occ/figures/model17ms_effect_sizes_all2.tif',
       effects_plot, width = 8, height = 8)

## Report some
effects %>%
  filter(high_prob == TRUE) %>%
  select(species, variable, season, median, q025, q975, prob_positive, prob_negative) %>%
  arrange(variable, species)

## Table for Appendix
#create table with rows for species and season, columns for variable, and values are prob_positive or prob_negative with + - indicators
effects_table <- effects %>%
  mutate(effect_direction = ifelse(median >= 0, paste0(round(prob_positive, 2), ' +'), 
                                   paste0(round(prob_negative, 2), ' -'))) %>%
  select(species, season, variable, effect_direction) %>%
  pivot_wider(names_from = variable, values_from = effect_direction) %>%
  arrange(species, match(season, c('cool', 'wet', 'dry')))

effects_tableplo
write.csv(effects_table, 'H_sp_occ/results/model17_MSversion/effect_direction_table.csv', row.names = FALSE)


## Now go to marginal plots script.


## OPTIONAL: calculate seasonal effects. how to report this? --------------------------------------------------





cat("=== DERIVED PARAMETER SUMMARY ===\n\n")
cat("===", species_name, "===\n\n")
cat("Effects are on the logit scale (occupancy model linear predictor)\n\n")

# Group by variable for cleaner output
vars <- unique(derived_params$variable)
var_labels <- c(
  "elev_scaledNS" = "Elevation",
  "slope_scaledNS" = "Slope", 
  "tree_vol_scaledNS" = "Tree Volume",
  "dist_boundary_log_scaledNS" = "Distance to Boundary",
  "dist_dambo_log_scaledNS" = "Distance to Dambo",
  "dist_rivers_log_scaledNS" = "Distance to Rivers"
)

for(var in vars) {
  var_data <- derived_params[derived_params$variable == var, ]
  
  cat(sprintf("--- %s ---\n", var_labels[var]))
  
  for(i in 1:nrow(var_data)) {
    row <- var_data[i, ]
    significance <- ""
    if(row$significant_positive) significance <- "[95%CI ≠ 0]"
    if(row$significant_negative) significance <- "[95%CI ≠ 0]"
    
    f_score <- ""
    if(row$median >= 0) f_score <- row$prob_positive
    if(row$median < 0) f_score <- row$prob_negative
    
    f_dir <- ""
    if(row$median >= 0) f_dir <- "Pr effect is POS +"
    if(row$median < 0) f_dir <- "Pr effect is NEG -"
    
    f_level <- ""
    if((abs(row$prob_positive) >= 0.7) | (abs(row$prob_negative) >= 0.7)) f_level <- "*WEAK"
    if((abs(row$prob_positive) >= 0.8) | (abs(row$prob_negative) >= 0.8)) f_level <- "**MODERATE"
    if((abs(row$prob_positive) >= 0.9) | (abs(row$prob_negative) >= 0.9)) f_level <- "***STRONG"
    
    # cat(sprintf("  %s: %.3f [%.3f, %.3f] (P(>0)=%.2f)%s\n", 
    #             toupper(row$season), 
    #             row$median, row$q025, row$q975, 
    #             row$prob_positive, significance))
    
    cat(sprintf("  %s: %.2f %s %s %.3f [%.2f, %.2f] %s %s %s\n",
                toupper(row$season), 
                f_score, f_dir, '  |  ',
                row$median, row$q025, row$q975, 
                '  |  ', f_level, significance))
    
  }
  cat("\n")
}


#season effects (reference level is Cool)
dry_effect <- mod$beta.samples[,'season_dry']
wet_effect <- mod$beta.samples[,'season_wet']

#dry season effect:
quantile(dry_effect, c(0.025, 0.5, 0.975))
mean(dry_effect < 0)

#wet season effect:
quantile(wet_effect, c(0.025, 0.5, 0.975))
mean(wet_effect > 0)





## Only show ones with ≥0.9 f-score:

effects_sig <- effects %>% filter(prob_positive >= 0.9 | prob_negative >= 0.9)

effects_sig_plot <- ggplot(effects_sig, aes(species, y = median, color = season)) +
  geom_pointrange(aes(ymin = q025, ymax = q975), position = position_dodge(width = 0.5)) +
  # geom_point(position = position_dodge(width = 0.5)) +
  scale_color_manual(values = color_palette, drop = FALSE) +
  # scale_fill_manual(values = color_palette, drop = FALSE) +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  facet_wrap(~variable, scales = 'free') +
  theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1))

effects_sig_plot
# ggsave('H_sp_occ/figures/model11_effect_sizes_sigOnly.tif',
#        effects_sig_plot, width = 8, height = 6)
ggsave('H_sp_occ/figures/model16_effect_sizes_sigOnly.tif',
       effects_sig_plot, width = 8, height = 6)


