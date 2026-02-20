## 3. perform linear regression on yearly 'psi' estimates

library(data.table)
library(tidyverse)
library(broom)


# Create data for predicting ---------------------------------------------------

# Define factor levels
fence_levels <- c("north", "south")  # assuming north is reference level
season_levels <- c("cool", "wet", "dry")  # assuming cool is reference level
year_levels <- 2019:2024

# Create all combinations of fence, season, year
pred_data <- expand.grid(
  fence = fence_levels,
  season = season_levels,
  year = year_levels,
  # Set all continuous variables to their means (0 after scaling)
  elev_scaledNS = 0,
  slope_scaledNS = 0,
  tree_vol_scaledNS = 0,
  dist_boundary_log_scaledNS = 0,
  dist_dambo_log_scaledNS = 0,
  dist_rivers_log_scaledNS = 0
)

# Format factor variables to match
pred_data$fence_south <- ifelse(pred_data$fence == "south", 1, 0)
pred_data$season_wet <- ifelse(pred_data$season == "wet", 1, 0)
pred_data$season_dry <- ifelse(pred_data$season == "dry", 1, 0)
pred_data$year_factor <- factor(pred_data$year)

head(pred_data)

#save for temp 
write.csv(pred_data, 'H_sp_occ/outputs/prediction_data.csv', row.names = FALSE)


# Function to create design matrix for predictions -----------------------------
create_design_matrix <- function(data) {
  
  #for testing
  # data <- pred_data
  
  # number of covar combinations I'm predicting on
  n <- nrow(data)
  
  # Initialize design matrix with intercept, and set column names
  X <- matrix(0, nrow = n, ncol = 32)  # 32 columns based on the beta matrix
  colnames(X) <- c(
    "(Intercept)", "fence_south", "factor(year)2020", "factor(year)2021",
    "factor(year)2022", "factor(year)2023", "factor(year)2024", "season_wet",
    "season_dry", "elev_scaledNS", "slope_scaledNS", "tree_vol_scaledNS",
    "dist_boundary_log_scaledNS", "dist_dambo_log_scaledNS", "dist_rivers_log_scaledNS",
    "fence_south:factor(year)2020", "fence_south:factor(year)2021", 
    "fence_south:factor(year)2022", "fence_south:factor(year)2023", 
    "fence_south:factor(year)2024", "season_wet:elev_scaledNS", 
    "season_dry:elev_scaledNS", "season_wet:slope_scaledNS", 
    "season_dry:slope_scaledNS", "season_wet:tree_vol_scaledNS", 
    "season_dry:tree_vol_scaledNS", "season_wet:dist_boundary_log_scaledNS", 
    "season_dry:dist_boundary_log_scaledNS", "season_wet:dist_dambo_log_scaledNS", 
    "season_dry:dist_dambo_log_scaledNS", "season_wet:dist_rivers_log_scaledNS", 
    "season_dry:dist_rivers_log_scaledNS"
  )
  
  # Fill in the design matrix -- iterate through each combination of covar values
  for(i in 1:n) {
    
    # Intercept
    X[i, "(Intercept)"] <- 1
    
    # Fence effect
    X[i, "fence_south"] <- data$fence_south[i]
    
    # Year effects (2019 is reference)
    if(data$year[i] == 2020) X[i, "factor(year)2020"] <- 1
    if(data$year[i] == 2021) X[i, "factor(year)2021"] <- 1
    if(data$year[i] == 2022) X[i, "factor(year)2022"] <- 1
    if(data$year[i] == 2023) X[i, "factor(year)2023"] <- 1
    if(data$year[i] == 2024) X[i, "factor(year)2024"] <- 1
    
    # Season effects
    X[i, "season_wet"] <- data$season_wet[i]
    X[i, "season_dry"] <- data$season_dry[i]
    
    # Continuous variables (all set to 0)
    X[i, "elev_scaledNS"] <- data$elev_scaledNS[i]
    X[i, "slope_scaledNS"] <- data$slope_scaledNS[i]
    X[i, "tree_vol_scaledNS"] <- data$tree_vol_scaledNS[i]
    X[i, "dist_boundary_log_scaledNS"] <- data$dist_boundary_log_scaledNS[i]
    X[i, "dist_dambo_log_scaledNS"] <- data$dist_dambo_log_scaledNS[i]
    X[i, "dist_rivers_log_scaledNS"] <- data$dist_rivers_log_scaledNS[i]
    
    # Interaction terms - fence:year
    if(data$fence_south[i] == 1 && data$year[i] == 2020) X[i, "fence_south:factor(year)2020"] <- 1
    if(data$fence_south[i] == 1 && data$year[i] == 2021) X[i, "fence_south:factor(year)2021"] <- 1
    if(data$fence_south[i] == 1 && data$year[i] == 2022) X[i, "fence_south:factor(year)2022"] <- 1
    if(data$fence_south[i] == 1 && data$year[i] == 2023) X[i, "fence_south:factor(year)2023"] <- 1
    if(data$fence_south[i] == 1 && data$year[i] == 2024) X[i, "fence_south:factor(year)2024"] <- 1
    
    # Interaction terms - season:continuous (all continuous vars are 0, so these will be 0)
    X[i, "season_wet:elev_scaledNS"] <- data$season_wet[i] * data$elev_scaledNS[i]
    X[i, "season_dry:elev_scaledNS"] <- data$season_dry[i] * data$elev_scaledNS[i]
    X[i, "season_wet:slope_scaledNS"] <- data$season_wet[i] * data$slope_scaledNS[i]
    X[i, "season_dry:slope_scaledNS"] <- data$season_dry[i] * data$slope_scaledNS[i]
    X[i, "season_wet:tree_vol_scaledNS"] <- data$season_wet[i] * data$tree_vol_scaledNS[i]
    X[i, "season_dry:tree_vol_scaledNS"] <- data$season_dry[i] * data$tree_vol_scaledNS[i]
    X[i, "season_wet:dist_boundary_log_scaledNS"] <- data$season_wet[i] * data$dist_boundary_log_scaledNS[i]
    X[i, "season_dry:dist_boundary_log_scaledNS"] <- data$season_dry[i] * data$dist_boundary_log_scaledNS[i]
    X[i, "season_wet:dist_dambo_log_scaledNS"] <- data$season_wet[i] * data$dist_dambo_log_scaledNS[i]
    X[i, "season_dry:dist_dambo_log_scaledNS"] <- data$season_dry[i] * data$dist_dambo_log_scaledNS[i]
    X[i, "season_wet:dist_rivers_log_scaledNS"] <- data$season_wet[i] * data$dist_rivers_log_scaledNS[i]
    X[i, "season_dry:dist_rivers_log_scaledNS"] <- data$season_dry[i] * data$dist_rivers_log_scaledNS[i]
  }
  
  return(X)
}


# Create design matrix ---------------------------------------------------------

X_pred <- create_design_matrix(pred_data)
str(X_pred)
colnames(X_pred)
head(X_pred)

#save temporarily for practice
write.csv(X_pred, 'H_sp_occ/outputs/prediction_design_matrix.csv', row.names = FALSE)


# Function to make predictions given beta estimates ----------------------------
make_predictions <- function(beta_samples, X_design, pred_data) {
  # beta_samples should be a matrix where each row is a posterior sample
  # and columns correspond to the parameters

  ##for testing, first read in files below
  beta_samples <- beta_samples
  X_design <- X_pred
  pred_data <- pred_data
  ##
    
  # Calculate linear predictor for each posterior sample
  linear_pred <- X_design %*% t(beta_samples)

  # Transform to probability scale using logistic function
  psi_pred <- plogis(linear_pred)
  
  # Calculate summary statistics
  psi_mean <- rowMeans(psi_pred)
  psi_lower <- apply(psi_pred, 1, quantile, probs = 0.025)
  psi_upper <- apply(psi_pred, 1, quantile, probs = 0.975)
  
  # Add predictions to data frame
  results <- cbind(pred_data[, c("fence", "season", "year")], 
                   psi_mean = psi_mean,
                   psi_lower = psi_lower,
                   psi_upper = psi_upper)
  
  return(results)
}


# PREDICT (mod11) -------------------------------------------------------------

files11 <- list.files('H_sp_occ/outputs/model11/', pattern = 'RDS', full.names = TRUE)

pred_list <- list()
for (ff in files11){
  
  species_name <- gsub('H_sp_occ/outputs/model11//model11_', '', ff)
  species_name <- gsub('_allChains.RDS', '', species_name)
  
  cat('## Predicting for', species_name, '\n')
  
  mod <- readRDS(ff)
  
  beta_samples <- mod$beta.samples
  predictions <- make_predictions(beta_samples, X_pred, pred_data)

  predictions$species <- species_name
  pred_list[[species_name]] <- predictions
    
}

all_preds <- rbindlist(pred_list)
# write.csv(all_preds, 'H_sp_occ/results/model11/mod11_trend_predictions.csv')


# PREDICT (mod13) -------------------------------------------------------------

#this model is like the previous but with no envir variables (elev, slope, etc.)

files13 <- list.files('H_sp_occ/outputs/model13/', pattern = 'RDS', full.names = TRUE)

#modify prediction data
pred_data_reduced <- pred_data %>% select(fence, season, year, fence_south, season_wet, season_dry, year_factor)

#modify design matrix
variables_to_keep <- c('(Intercept)', 'fence_south', 'factor(year)2020', 'factor(year)2021',
                       'factor(year)2022', 'factor(year)2023', 'factor(year)2024', 'season_wet',
                       'season_dry', 'fence_south:factor(year)2020', 'fence_south:factor(year)2021', 
                       'fence_south:factor(year)2022', 'fence_south:factor(year)2023', 
                       'fence_south:factor(year)2024')

X_pred_reduced <- X_pred[, variables_to_keep]

#predict
pred_reduced_list <- list()

for (rr in files13){
  
  species_name <- gsub('H_sp_occ/outputs/model13//model13_', '', rr)
  species_name <- gsub('_allChains.RDS', '', species_name)
  
  cat('## Predicting for', species_name, '\n')
  
  mod <- readRDS(rr)
  
  beta_samples <- mod$beta.samples
  predictions <- make_predictions(beta_samples, X_pred_reduced, pred_data_reduced)
  
  predictions$species <- species_name
  pred_reduced_list[[species_name]] <- predictions
  
}

all_preds_reduced <- rbindlist(pred_reduced_list)
write.csv(all_preds_reduced, 'H_sp_occ/results/model13/mod13_trend_predictions.csv')


# PREDICT (mod17) -------------------------------------------------------------

files17 <- list.files('H_sp_occ/outputs/model17/', pattern = 'RDS', full.names = TRUE)

pred_list <- list()
for (ff in files17){
  
  species_name <- gsub('H_sp_occ/outputs/model17//model17_', '', ff)
  species_name <- gsub('_allChains.RDS', '', species_name)
  
  cat('## Predicting for', species_name, '\n')
  
  mod <- readRDS(ff)
  
  beta_samples <- mod$beta.samples
  predictions <- make_predictions(beta_samples, X_pred, pred_data)
  
  predictions$species <- species_name
  pred_list[[species_name]] <- predictions
  
}

all_preds <- rbindlist(pred_list)
write.csv(all_preds, 'H_sp_occ/results/model17/mod17_trend_predictions.csv')


## Function to perform regression analysis -------------------------------------
analyze_yearly_trends <- function(predictions) {
  
  # Split by fence and season combinations for separate trend analyses
  trend_results <- predictions %>%
    group_by(fence, season) %>%
    do({
      model <- lm(psi_mean ~ year, data = .)
      model_summary <- summary(model)
      
      data.frame(
        slope = coef(model)[2],
        intercept = coef(model)[1],
        r_squared = model_summary$r.squared,
        p_value = model_summary$coefficients[2, 4],
        slope_se = model_summary$coefficients[2, 2]
      )
    }) %>%
    mutate(
      # Calculate trend over study period
      change_total = slope * (max(predictions$year) - min(predictions$year)),
      change_percent = (change_total / intercept) * 100,
      # Add significance and interpretation
      significant = p_value < 0.05,
      trend_direction = ifelse(slope > 0, "Increasing", "Decreasing"),
      trend_strength = case_when(
        abs(slope) < 0.01 ~ "Weak",
        abs(slope) < 0.02 ~ "Moderate", 
        TRUE ~ "Strong"
      )
    )
  
  return(trend_results)
}


# PERFORM REGRESSION ------------------------------------------------------------------

#read in if necessary
# all_preds <- fread('H_sp_occ/results/model11/mod11_trend_predictions.csv')
# all_preds_reduced <- fread('H_sp_occ/results/model13/mod13_trend_predictions.csv')
# all_preds <- fread('H_sp_occ/results/model17/mod17_trend_predictions.csv')

#perform regression on mod11 data 
trend_list <- list()
for (ss in unique(all_preds$species)) {

  preds_ss <- all_preds %>% filter(species == ss)
  trend_ss <- analyze_yearly_trends(preds_ss)
  trend_ss$species <- ss
  trend_list[[ss]] <- trend_ss
  
}

all_trends <- rbindlist(trend_list)
# write.csv(all_trends, 'H_sp_occ/results/model11/mod11_regression.csv')
write.csv(all_trends, 'H_sp_occ/results/model17/mod17_regression.csv')

## TO REPORT IN APPENDIX:
all_trends <- fread('H_sp_occ/results/model17/mod17_regression.csv')

#COOL
all_trends %>% select(fence, season, species, slope, r_squared, p_value) %>%
  filter(fence == 'south', season == 'cool') %>%
  pivot_wider(names_from = fence, values_from = c(slope, r_squared, p_value)) %>%
  mutate(across(where(is.numeric), ~ round(.x, 2)))

all_trends %>% select(fence, season, species, slope, r_squared, p_value) %>%
  filter(fence == 'north', season == 'cool') %>%
  pivot_wider(names_from = fence, values_from = c(slope, r_squared, p_value)) %>%
  mutate(across(where(is.numeric), ~ round(.x, 2)))

#WET
all_trends %>% select(fence, season, species, slope, r_squared, p_value) %>%
  filter(fence == 'south', season == 'wet') %>%
  pivot_wider(names_from = fence, values_from = c(slope, r_squared, p_value)) %>%
  mutate(across(where(is.numeric), ~ round(.x, 2)))

all_trends %>% select(fence, season, species, slope, r_squared, p_value) %>%
  filter(fence == 'north', season == 'wet') %>%
  pivot_wider(names_from = fence, values_from = c(slope, r_squared, p_value)) %>%
  mutate(across(where(is.numeric), ~ round(.x, 2)))

#DRY
all_trends %>% select(fence, season, species, slope, r_squared, p_value) %>%
  filter(fence == 'south', season == 'dry') %>%
  pivot_wider(names_from = fence, values_from = c(slope, r_squared, p_value)) %>%
  mutate(across(where(is.numeric), ~ round(.x, 3)))

all_trends %>% select(fence, season, species, slope, r_squared, p_value) %>%
  filter(fence == 'north', season == 'dry') %>%
  pivot_wider(names_from = fence, values_from = c(slope, r_squared, p_value)) %>%
  mutate(across(where(is.numeric), ~ round(.x, 2)))



# Perform regression on model13 data
# trend_list_reduced <- list()
# for (tt in unique(all_preds_reduced$species)) {
#   
#   preds_tt <- all_preds_reduced %>% filter(species == tt)
#   trend_tt <- analyze_yearly_trends(preds_tt)
#   trend_tt$species <- tt
#   trend_list_reduced[[tt]] <- trend_tt
#   
# }
# 
# all_trends_reduced <- rbindlist(trend_list_reduced)
# write.csv(all_trends_reduced, 'H_sp_occ/results/model13/mod13_regression.csv')


# PERFORM POOLED REGRESSION ----------------------------------------------------

#read in if necessary
# all_preds <- fread('H_sp_occ/results/model11/mod11_trend_predictions.csv')
# all_preds_reduced <- fread('H_sp_occ/results/model13/mod13_trend_predictions.csv')

#inspect
nrow(all_preds)
  9*6*3*2 #9 species, 6 years, 3 seasons, 2 fence regions

#pool seasons (full model)
trends_pooled <- all_preds %>%
    group_by(fence, species) %>%
    do({
      model <- lm(psi_mean ~ year, data = .)
      model_summary <- summary(model)
      
      data.frame(
        slope = coef(model)[2],
        intercept = coef(model)[1],
        r_squared = model_summary$r.squared,
        p_value = model_summary$coefficients[2, 4],
        slope_se = model_summary$coefficients[2, 2]
      )
    })
trends_pooled
# write.csv(trends_pooled, 'H_sp_occ/results/model11/mod11_regression_pooledSeasons.csv')
write.csv(trends_pooled, 'H_sp_occ/results/model17/mod17_regression_pooledSeasons.csv')

#pool seasons (reduced model)
# trends_pooled_reduced <- all_preds_reduced %>%
#   group_by(fence, species) %>%
#   do({
#     model <- lm(psi_mean ~ year, data = .)
#     model_summary <- summary(model)
#     
#     data.frame(
#       slope = coef(model)[2],
#       intercept = coef(model)[1],
#       r_squared = model_summary$r.squared,
#       p_value = model_summary$coefficients[2, 4],
#       slope_se = model_summary$coefficients[2, 2]
#     )
#   })
# trends_pooled_reduced
# write.csv(trends_pooled_reduced, 'H_sp_occ/results/model13/mod13_regression_pooledSeasons.csv')

# trends_pooled_reduced %>% filter(species %in% c('buffalo','impala','zebra')) %>%
#                   select(fence, species, slope, r_squared, p_value)


##TO REPORT IN TABLE 2:
trends_pooled <- fread('H_sp_occ/results/model17/mod17_regression_pooledSeasons.csv')

#SOUTH
trends_pooled %>% select(fence, species, slope, r_squared, p_value) %>%
  filter(fence == 'south') %>%
  pivot_wider(names_from = fence, values_from = c(slope, r_squared, p_value)) %>%
  mutate(across(where(is.numeric), ~ round(.x, 2)))

#NORTH
trends_pooled %>% select(fence, species, slope, r_squared, p_value) %>%
  filter(fence == 'north') %>%
  pivot_wider(names_from = fence, values_from = c(slope, r_squared, p_value)) %>%
  mutate(across(where(is.numeric), ~ round(.x, 2)))


# PLOT -------------------------------------------------------------------------

#read in if necessary
# all_preds <- fread('H_sp_occ/results/model11/mod11_trend_predictions.csv')
# all_preds_reduced <- fread('H_sp_occ/results/model13/mod13_trend_predictions.csv')
all_preds <- fread('H_sp_occ/results/model17/mod17_trend_predictions.csv')

#combine
# all_preds$model <- 'model11'
# all_preds_reduced$model <- 'model13'
# all_preds_combined <- bind_rows(all_preds, all_preds_reduced)
all_preds$model <- 'model17'

#keep the correct model for each species
# table(all_preds_combined$model, all_preds_combined$species)
# mod11_sp <- c('eland','elephant','kudu','sable','warthog','waterbuck')
# mod13_sp <- c('buffalo','impala','zebra')
# all_preds_filtered <- all_preds_combined %>% filter((model == 'model11' & species %in% mod11_sp) |
#                                                     (model == 'model13' & species %in% mod13_sp))
# table(all_preds_filtered$model, all_preds_filtered$species)

all_preds_filtered <- all_preds
table(all_preds_filtered$model, all_preds_filtered$species)

#make dataframe to plot vertical line at 2022 for all sp. except elephants
vline_data <- data.frame(x = 2022,
                         species = c('buffalo','eland','impala','kudu','sable','warthog','waterbuck','zebra'))


## Note these are just for one season (see below for pooled)

## POINT/LINE PLOTS (POSTERIOR QUANTILES)
# panel_preds <- ggplot(all_preds_filtered %>% filter(season == 'cool'), 
#        aes(x = year, y = psi_mean, color = fence)) +
#   geom_vline(data = vline_data, aes(xintercept = x), color = 'gray40', linetype = 'dotdash') + #mark 2022 for secondary translocations
#   geom_line(aes(linetype = fence), position = position_dodge(width = 0.5)) +
#   geom_pointrange(aes(ymin = psi_lower, ymax = psi_upper, shape = fence),  #, linetype = fence
#                   alpha = 0.5,
#                   position = position_dodge(width = 0.6)) +
#   geom_point(aes(shape = fence), size = 3,
#              position = position_dodge(width = 0.6)) +
#   scale_color_manual(values = c("south" = "#008080", "north" = "#e69f00")) +
#   # scale_color_manual(values = c("south" = "#1b9e77", "north" = "#d95f02")) +
#   scale_x_continuous(breaks = c(2019, 2020, 2021, 2022, 2023, 2024)) +
#   facet_wrap(~species, nrow = 3) +
#   labs(x = NULL, y = "Occupancy Probability (ψ ± 95% CI)") +
#   theme_bw(base_size = 14) + theme(axis.text.x = element_text(angle = 45, hjust = 1),
#                                    # legend.position = 'none'
#                                    legend.title = element_blank())
# panel_preds

# ggsave('H_sp_occ/figures/psi_trends_panel_dry.png', panel_preds, width = 7, height = 5)
# ggsave('H_sp_occ/figures/psi_trends_panel_cool_mod17.png', panel_preds, width = 7, height = 5)


## REGRESSION PLOTS
panel_lm <- ggplot(all_preds_filtered %>% filter(season == 'cool'), 
                   aes(x = year, y = psi_mean, color = fence, shape = fence)) + #linetype = fence, 
  geom_vline(data = vline_data, aes(xintercept = x), color = 'gray40', 
             linetype = 'dotdash', linewidth = 1) + #mark 2022 for secondary translocations
  geom_pointrange(aes(ymin = psi_lower, ymax = psi_upper, shape = fence),  #, linetype = fence
                  alpha = 0.5,
                  position = position_dodge(width = 0.6),
                  show.legend = FALSE) +
  geom_point(size = 3, position = position_dodge(width = 0.6)) +
  geom_smooth(method = "lm", se = TRUE, level = 0.95, alpha = 0.2,
              aes(fill = fence, linetype = fence)) +
  scale_x_continuous(breaks = 2019:2024) +
  scale_color_manual(values = c(south = "#008080", north = "#e69f00")) +
  scale_fill_manual(values = c(south = "#008080", north = "#e69f00")) +
  facet_wrap(~species, nrow = 3) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = NULL, y = "Occupancy Probability (ψ ± 95% CI)") +
  theme_bw(base_size = 14) + theme(axis.text.x = element_text(angle = 45, hjust = 1),
                          # legend.position = 'none'
                          legend.title = element_blank(),
                     strip.background = element_rect(fill = NA))
panel_lm

# ggsave('H_sp_occ/figures/psi_lm_panel_cool.png', panel_lm, width = 7, height = 5)
ggsave('H_sp_occ/figures/psi_lm_panel_cool_mod17.png', panel_lm, width = 7, height = 5)


## POOLED REGRESSION PLOTS -- IS THIS RIGHT? I just didn't filter by season.
panel_lm_pooled <- ggplot(all_preds_filtered, 
                   aes(x = year, y = psi_mean, color = fence, shape = fence)) + #linetype = fence, 
  geom_vline(data = vline_data, aes(xintercept = x), color = 'gray40', 
             linetype = 'dotdash', linewidth = 1) + #mark 2022 for secondary translocations
  geom_pointrange(aes(ymin = psi_lower, ymax = psi_upper, shape = fence),  #, linetype = fence
                  alpha = 0.5,
                  position = position_dodge(width = 0.6),
                  show.legend = FALSE) +
  geom_point(size = 3, position = position_dodge(width = 0.6)) +
  geom_smooth(method = "lm", se = TRUE, level = 0.95, alpha = 0.2,
              aes(fill = fence, linetype = fence)) +
  scale_x_continuous(breaks = 2019:2024) +
  scale_color_manual(values = c(south = "#008080", north = "#e69f00")) +
  scale_fill_manual(values = c(south = "#008080", north = "#e69f00")) +
  facet_wrap(~species, nrow = 3) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = NULL, y = "Occupancy Probability (ψ ± 95% CI)") +
  theme_bw(base_size = 14) + theme(axis.text.x = element_text(angle = 45, hjust = 1),
                                   # legend.position = 'none'
                                   legend.title = element_blank(),
                                   strip.background = element_rect(fill = NA))
panel_lm_pooled

# ggsave('H_sp_occ/figures/psi_lm_panel_pooled.png', panel_lm_pooled, width = 7, height = 5)
ggsave('H_sp_occ/figures/psi_lm_panel_pooled_mod17.png', panel_lm_pooled, width = 7, height = 5)



## TRY TO ADD TABLE TO PLOT? NAH -----------------------------------------------

library(gt)

trends_pooled <- fread('H_sp_occ/results/model11/mod11_regression_pooledSeasons.csv')
trends_pooled_reduced <- fread('H_sp_occ/results/model13/mod13_regression_pooledSeasons.csv')

#simplify
mod11_simple <- trends_pooled %>% filter(!species %in% c('buffalo','impala','zebra')) %>%
  select(fence, species, slope, r_squared, p_value)

mod13_simple <- trends_pooled_reduced %>% filter(species %in% c('buffalo','impala','zebra')) %>%
  select(fence, species, slope, r_squared, p_value)

#combine and reshape
trends_df <- bind_rows(mod11_simple, mod13_simple)
trends_df <- trends_df %>% pivot_wider(names_from = fence,
                                       values_from = c(slope, r_squared, p_value)) %>%
                           select(species, slope_south, r_squared_south, p_value_south,
                                  slope_north, r_squared_north, p_value_north,)

#clean up?
colnames(trends_df) <- gsub('r_squared','R²',colnames(trends_df))
colnames(trends_df) <- gsub('p_value','p',colnames(trends_df))
colnames(trends_df) <- gsub('_','\n',colnames(trends_df))

#convert to gt object
trends_table <- trends_df %>%
  gt(rowname_col = "species") %>%
  fmt_number(columns = everything(), decimals = 2)

trends_table

panel_w_table <- (panel_lm_pooled + theme(legend.position = 'top')) + 
  wrap_table(trends_table, panel = "full", space = "fixed")
ggsave('H_sp_occ/figures/psi_lm_panel_pooled_wTable.png', panel_w_table, width = 12, height = 5)


#cool!

#testing other ways:



#testing to clean up labels:
#testing
trends_df2 <- trends_df %>%
  rename_with(~sub("_north$", "", .x), ends_with("north")) %>%
  rename_with(~sub("_south$", "", .x), ends_with("south"))

trends_table2 <- trends_df %>%
  gt(rowname_col = "species") %>%
  # Add spanner headers for Region
  tab_spanner(
    label = "north",
    columns = c(slope, r_squared, p_value)
  ) %>%
  tab_spanner(
    label = "south",
    columns = c(Slope1, R21, `p-value1`)
  ) %>%
  # clean up labels
  cols_label(
    slope = "slope",
    r_squared = "R²",
    p_value = "p",
    slope1 = "slope",
    r_squared1 = "R²",
    p_value1 = "p"
  ) %>%
  fmt_number(columns = everything(), decimals = 3)

