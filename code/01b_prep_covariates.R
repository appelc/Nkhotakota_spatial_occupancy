## Format detection data for spOccupancy models (specifically multi-season)

library(data.table)
library(tidyverse)
library(reshape2)


## Read in data ----------------------------------------------------------------

#detection arrays
dh_arrays <- readRDS('data/cleaned/detection_arrays_yearly.RDS')

#site-level envir covariates
env_covar <- fread('data/raw/covariates.csv')


## Format structural covariates (OCCUPANCY) ------------------------------------

dim(dh_arrays[[1]]) #sites, year-seasons, weeks (max)

#site effect (vector w/ length n_sites)
(site_covar <- rownames(dh_arrays$elephant[,,1]))
  head(site_covar,10)
  (n_sites <- length(site_covar))

#year-season effect (n_sites*n_year_seasons matrix)
(year_seasons <- colnames(dh_arrays$elephant[,,1]))
year_season_covar <- matrix(rep(year_seasons, times = n_sites), nrow = n_sites, byrow = TRUE)
  head(year_season_covar); dim(year_season_covar)
  (n_year_seasons <- ncol(year_season_covar)) #remember 'year_seasons' are primary periods

#season (n_sites*n_year_seasons matrix)
(seasons <- sapply(strsplit(year_seasons, '_'), '[', 2))
season_covar <- matrix(rep(seasons, times = n_sites), nrow = n_sites, byrow = TRUE)
  head(season_covar); dim(season_covar)
  (n_seasons <- length(unique(seasons)))
  
#year effect (n_sites*n_year_seasons matrix)
(years <- sapply(strsplit(year_seasons, '_'), '[', 1))
year_covar <- matrix(rep(years, times = n_sites), nrow = n_sites, byrow = TRUE)
  head(year_covar); dim(year_covar)
  (n_years <- length(unique(years)))


## Format environmental covariates (OCCUPANCY) ---------------------------------

#these need to be vectors w/ length n_sites

names(env_covar)

#which to keep?
(covar_names <- c('site','elev_mean_3x3','pct_slope_mean_3x3','tree_vol_mean_3x3',
                  'dist_dambo_mean_3x3','dist_rivers_mean_3x3','distance_to_perimeter','fence'))
head(env_covar[,..covar_names]) 
#these are NOT standardized (we'll do it in the model specification via scale() -- or below for N/S separately)

#keep only the select covariates
env_covar <- env_covar %>% select(all_of(covar_names))    
length(unique(env_covar$site)) #all sites with covariates here

#keep only sites that have data in det histories
site_covar
env_covar <- env_covar %>% filter(site %in% site_covar)
length(unique(env_covar$Site)) #should match n_sites now

#use North/South for 'fence' instead of 0/1
env_covar$fence <- ifelse(env_covar$fence == 0, 'North', 'South')

#sort by site
env_covar <- env_covar %>% arrange(site)


## Log-transform the distance covariates ---------------------------------------

head(env_covar)

env_covar <- env_covar %>% mutate(dist_perimeter_log = log(distance_to_perimeter),
                                  dist_dambo_log = log(dist_dambo_mean_3x3),
                                  dist_rivers_log = log(dist_rivers_mean_3x3))

hist(env_covar$distance_to_perimeter)
hist(env_covar$dist_perimeter_log)

hist(env_covar$dist_dambo_mean_3x3)
hist(env_covar$dist_dambo_log)

hist(env_covar$dist_rivers_mean_3x3)
hist(env_covar$dist_rivers_log)


## Scale covariates ------------------------------------------------------------

head(env_covar)

#use the mean/SD from each respective region (north/south)
env_covar_scaled <- env_covar %>% group_by(fence) %>%
                      mutate(elev_mean_scaled = (elev_mean_3x3 - mean(elev_mean_3x3, na.rm = TRUE)) / sd(elev_mean_3x3, na.rm = TRUE),
                             pct_slope_mean_scaled = (pct_slope_mean_3x3 - mean(pct_slope_mean_3x3)) / sd(pct_slope_mean_3x3),
                             tree_vol_mean_scaled = (tree_vol_mean_3x3 - mean(tree_vol_mean_3x3)) / sd(tree_vol_mean_3x3),
                             dist_dambo_mean_scaled = (dist_dambo_mean_3x3 - mean(dist_dambo_mean_3x3)) / sd(dist_dambo_mean_3x3),
                             dist_rivers_mean_scaled = (dist_rivers_mean_3x3 - mean(dist_rivers_mean_3x3)) / sd(dist_rivers_mean_3x3),
                             dist_perimeter_scaled = (distance_to_perimeter - mean(distance_to_perimeter)) / sd(distance_to_perimeter),
                             dist_dambo_log_scaled = (dist_dambo_log - mean(dist_dambo_log)) / sd(dist_dambo_log),
                             dist_rivers_log_scaled = (dist_rivers_log - mean(dist_rivers_log)) / sd(dist_rivers_log),
                             dist_perimeter_log_scaled = (dist_perimeter_log - mean(dist_perimeter_log)) / sd(dist_perimeter_log)
                             ) %>% 
                      ungroup()

env_covar_scaled

#QC (elevation as example)
env_covar %>% summarise(mean = mean(elev_mean_3x3),
                        SD = sd(elev_mean_3x3))
env_covar %>% group_by(fence) %>% summarise(mean = mean(elev_mean_3x3),
                                            SD = sd(elev_mean_3x3))

#check a north value
head(env_covar_scaled)
(1055 - 829) / 186 #using mean/SD from all
(1055 - 898) / 174 #using mean/SD from N

#check a south value
tail(env_covar_scaled)
(589 - 829) / 186 #using mean/SD from all
(589 - 736) / 161 #using mean/SD from N

#they're close to mean=0 and sd=1 but not exactly
env_covar_scaled %>% group_by(fence) %>% summarise(mean = mean(elev_mean_scaled),
                                                   SD = sd(elev_mean_scaled))

#view log/scaled
hist(env_covar_scaled$dist_dambo_mean_3x3)
hist(env_covar_scaled$dist_dambo_mean_scaled)
hist(env_covar_scaled$dist_dambo_log)
hist(env_covar_scaled$dist_dambo_log_scaled)

#combine and save
occ_covar_scaled <- list('site' = site_covar, 
                         'year_season' = year_season_covar,
                         'season' = season_covar, 
                         'year' = year_covar,
                         'elev' = env_covar_scaled$elev_mean_3x3,
                         'elev_scaledNS' = env_covar_scaled$elev_mean_scaled,
                         'slope' = env_covar_scaled$pct_slope_mean_3x3,
                         'slope_scaledNS' = env_covar_scaled$pct_slope_mean_scaled,
                         'tree_vol' = env_covar_scaled$tree_vol_mean_3x3,
                         'tree_vol_scaledNS' = env_covar_scaled$tree_vol_mean_scaled,
                         'dist_dambo' = env_covar_scaled$dist_dambo_mean_3x3,
                         'dist_dambo_scaledNS' = env_covar_scaled$dist_dambo_mean_scaled,
                         'dist_dambo_log' = env_covar_scaled$dist_dambo_log,
                         'dist_dambo_log_scaledNS' = env_covar_scaled$dist_dambo_log_scaled,
                         'dist_rivers' = env_covar_scaled$dist_rivers_mean_3x3,
                         'dist_rivers_scaledNS' = env_covar_scaled$dist_rivers_mean_scaled,
                         'dist_rivers_log' = env_covar_scaled$dist_rivers_log,
                         'dist_rivers_log_scaledNS' = env_covar_scaled$dist_rivers_log_scaled,
                         'dist_boundary' = env_covar_scaled$distance_to_perimeter,
                         'dist_boundary_scaledNS' = env_covar_scaled$dist_perimeter_scaled,
                         'dist_boundary_log' = env_covar_scaled$dist_perimeter_log,
                         'dist_boundary_log_scaledNS' = env_covar_scaled$dist_perimeter_log_scaled,
                         'fence' = env_covar_scaled$fence
)
str(occ_covar_scaled) #make sure they are all either vectors (length n_sites) or matrices (n_sites * n_year_seasons)

#save
saveRDS(occ_covar_scaled, 'data/cleaned/occ_covar.RDS')

#also save mean/SD and min/max (we'll need them for prediction)
covar_means_grouped <- env_covar %>% 
        pivot_longer(cols = c(elev_mean_3x3, pct_slope_mean_3x3, tree_vol_mean_3x3,
                              dist_dambo_mean_3x3, dist_rivers_mean_3x3, distance_to_perimeter,
                              dist_dambo_log, dist_rivers_log, dist_perimeter_log),
                     names_to = 'variable', values_to = 'value') %>%
        group_by(fence, variable) %>%
        summarise(mean = mean(value, na.rm = TRUE),
                  sd = sd(value, na.rm = TRUE),
                  min = min(value, na.rm = TRUE),
                  max = max(value, na.rm = TRUE),
                  .groups = 'drop')

write.csv(covar_means_grouped, 'data/cleaned/covar_mean_sd_scaledNS.csv')


## Format coordinates ----------------------------------------------------------

plot_coords <- fread('data/raw/coordinates.csv')
head(plot_coords)
length(unique(plot_coords$site)) #has all 210

#keep only sites that have data in det histories
site_covar
plot_coords <- plot_coords %>% filter(site %in% site_covar)
length(unique(plot_coords$site)) #now matches n_sites

#make sure it's sorted by site
plot_coords <- plot_coords %>% arrange(site)

#save
saveRDS(plot_coords, 'data/cleaned/coords.RDS')


## Now go to script '01c_prep_inputs.R'
