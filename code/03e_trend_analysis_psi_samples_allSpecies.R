# Deriving yearly occupancy using psi.samples directly

library(data.table)
library(tidyverse)
library(spOccupancy)
library(patchwork)


## Read in input data ----------------------------------------------------------
sp_inputs <- readRDS('H_sp_occ/sp_inputs.RDS')
sp_inputs_example <- sp_inputs$elephant #only need one example to get covar names
rm(sp_inputs)

#get site names, time period names, and north/south from input file
site_names <- sp_inputs_example$occ.covs$site
yearSeason_names <- sp_inputs_example$occ.covs$year_season[1,]
fence_names <- data.frame(site = sp_inputs_example$occ.covs$site,
                          fence = sp_inputs_example$occ.covs$fence)

north_sites <- fence_names$site[fence_names$fence == "North"]
south_sites <- fence_names$site[fence_names$fence == "South"]

cool_names <- yearSeason_names[grepl("_cool$", yearSeason_names)]
dry_names <- yearSeason_names[grepl("_dry$", yearSeason_names)]
wet_names <- yearSeason_names[grepl("_wet$", yearSeason_names)]


## Read in model outputs -------------------------------------------------------
output_files <- list.files('H_sp_occ/outputs/model17MSversion/', pattern = 'RDS', full.names = TRUE)

#create object to store formatted psi samples
psi_matrices <- list()

#read in each RDS, extract, and format psi samples
for (ff in output_files){
  
  #extract species name
  sp_name <- strsplit(basename(ff), '_')[[1]][2]
  cat('Reading data for:', sp_name, '\n')
  
  #read output file and extract psi samples
  output <- readRDS(ff)
  psi_samples <- output$psi.samples
  
  dim(output$psi.samples)
  
  #set dimension names
  dimnames(psi_samples) <- list(
    iter = seq_len(dim(psi_samples)[1]),
    site = site_names,
    time = yearSeason_names)
  
  #subset north/south
  psi_samples_north <- psi_samples[, north_sites, ]
  psi_samples_south <- psi_samples[, south_sites, ]
  
  #further subset by season
  psi_samples_north_cool <- psi_samples_north[, , cool_names]
  psi_samples_north_dry <- psi_samples_north[, , dry_names]
  psi_samples_north_wet <- psi_samples_north[, , wet_names]
  psi_samples_south_cool <- psi_samples_south[, , cool_names]
  psi_samples_south_dry <- psi_samples_south[, , dry_names]
  psi_samples_south_wet <- psi_samples_south[, , wet_names]

  #for each array, for each iteration, take the mean across sites for each year_season
  yearly_matrices <- list(
    north_cool = apply(psi_samples_north_cool, c(1, 3), mean),
    north_dry = apply(psi_samples_north_dry, c(1, 3), mean),
    north_wet = apply(psi_samples_north_wet, c(1, 3), mean),
    south_cool = apply(psi_samples_south_cool, c(1, 3), mean),
    south_dry = apply(psi_samples_south_dry, c(1, 3), mean),
    south_wet = apply(psi_samples_south_wet, c(1, 3), mean)
  )
  
  #store
  psi_matrices[[sp_name]] <- yearly_matrices
  
  #clear unnecessary objects before moving to next species
  rm(output, psi_samples, 
     psi_samples_north, psi_samples_south,
     psi_samples_north_cool, psi_samples_north_dry, psi_samples_north_wet,
     psi_samples_south_cool, psi_samples_south_dry, psi_samples_south_wet,
     yearly_matrices)
  
}

names(psi_matrices)
str(psi_matrices$buffalo) #each species has 6 matrices (fence*seasons)
dim(psi_matrices$buffalo$north_cool) #each matrix has dimensions iterations*years


## Run post-hoc trend analysis -------------------------------------------------

#create covariate dataframe for regression
(covs <- data.frame(year = 0:ncol(psi_matrices$elephant[[1]])))

#set priors and initial values (I don't have RE so no need for sigma.sq)
inits_list <- list(beta = 0, tau.sq = 1)
priors_list <- list(beta.normal = list(mean = 0, var = 10000),
               tau.sq.ig = c(0.001, 0.001))

#initialize objects to store results
psi_long_df <- data.frame()
posthoc_results_df <- data.frame()
posthoc_results_list <- list()

#for each species,
for (ss in names(psi_matrices)){
  yearly_matrices <- psi_matrices[[ss]]
  cat('Processing:', ss, '\n')
  
  #for each of the 6 fence_season combinations, run postHocLM and store results:
  for(aa in names(yearly_matrices)){
    data <- yearly_matrices[[aa]]
    
    #summarize and store quantiles across iterations for each year_season (but we'll do regression on full samples below)
    psi_quants <- apply(data, 2, quantile, c(0.025, 0.5, 0.975))
    psi_long_ss <- as.data.frame(t(psi_quants)) %>%
      mutate(species = ss,
             year_season = rownames(.),
             fence_season = aa) %>%
      separate(year_season, into = c("year", "season"), sep = "_", remove = FALSE) %>%
      separate(fence_season, into = c("fence", "season2"), sep = "_", remove = FALSE)
    psi_long_df <- rbind(psi_long_df,
                         psi_long_ss)
    
    #create input list for regression
    input <- list(y = data,
                  covs = covs)
    
    #run post-hoc regression on psi samples
    posthoc_result <- postHocLM(
      formula = ~ year,
      data = input,
      inits = inits_list,
      priors = priors_list,
      n.chains = 3,
      # n.report = 10000,
      verbose = FALSE
    )
    print(summary(posthoc_result))
    
    #summarize beta samples (slope and intercept) from regression output
    beta_samples <- as.data.frame(posthoc_result$beta.samples)
    trend_aa <- ifelse(mean(beta_samples[, 2]) > 0, "increasing", "decreasing")
    f_score_aa <- ifelse(trend_aa == "increasing",
                         mean(beta_samples[, 2] > 0),
                         mean(beta_samples[, 2] < 0))
    
    #store results
    fence_aa <- strsplit(aa, "_")[[1]][1]
    season_aa <- strsplit(aa, "_")[[1]][2]
    
    posthoc_results_df <- rbind(posthoc_results_df,
                                data.frame(
                                  species = ss,
                                  fence = fence_aa,
                                  season = season_aa,
                                  trend = trend_aa,
                                  f_score = f_score_aa,
                                  intercept_mean = mean(posthoc_result$beta.samples[,1]),
                                  intercept_median = quantile(posthoc_result$beta.samples[,1], 0.5),
                                  intercept_lower_CI = quantile(posthoc_result$beta.samples[,1], 0.025),
                                  intercept_upper_CI = quantile(posthoc_result$beta.samples[,1], 0.975),
                                  slope_mean = mean(posthoc_result$beta.samples[,2]),
                                  slope_median = median(posthoc_result$beta.samples[,2]),
                                  slope_lower_CI = quantile(posthoc_result$beta.samples[,2], 0.025),
                                  slope_upper_CI = quantile(posthoc_result$beta.samples[,2], 0.975),
                                  r2 = mean(posthoc_result$bayes.R2)
                                ))
    
    posthoc_results_list[[ss]][[aa]] <- posthoc_result
  }
  
}

#QC
summary(posthoc_results_list$elephant$north_cool)
quantile(posthoc_results_list$elephant$north_cool$beta.samples[,2], probs = c(0.025, 0.5, 0.975))
posthoc_results_df %>% filter(species == 'elephant', season == 'cool', fence == 'north')

dim(posthoc_results_list$elephant$north_cool$y.hat.samples) #(28200 samples * 3 chains) x 6 year_seasons

head(posthoc_results_df) #regression slope, intercept, probability of increase, R2, etc.
head(psi_long_df) #predicted occupancy quantiles by year and season


## Plot results ----------------------------------------------------------------

#Option 1: plot slope values colored by F-scores
plot1 <- ggplot(posthoc_results_df, aes(x = season, y = slope_mean, 
                          ymin = slope_lower_CI, ymax = slope_upper_CI,
                          color = f_score, shape = fence)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_pointrange(size = 0.8, position = position_dodge(width = 0.5)) +
  facet_wrap(~species) +
  scale_color_gradient2(
    low = "red", mid = "gray80", high = "blue", midpoint = 0.5,
    name = "P(increasing)") +
  labs(
    x = "Season",
    y = "Occupancy trend (slope per year)",
    title = "Temporal occupancy trends") +
  theme_bw(base_size = 14) + theme(legend.position = "right",
                                   panel.grid.minor = element_blank())
plot1

#ordered by trend
slope_order <- posthoc_results_df %>%
  filter(season == "cool", fence == "south") %>%
  arrange(slope_mean) %>%
  pull(species)
posthoc_results_df$species <- factor(posthoc_results_df$species,
                                     levels = slope_order)
posthoc_results_df$fence <- factor(posthoc_results_df$fence,
                                   levels = c('south','north'))

plot1b_cool <- ggplot(posthoc_results_df[posthoc_results_df$season == 'cool',], 
                 aes(x = species, y = slope_mean, 
                     ymin = slope_lower_CI, ymax = slope_upper_CI,
                     color = f_score, shape = fence)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_pointrange(size = 0.8, position = position_dodge(width = 0.5)) +
  facet_wrap(~fence) +
  scale_color_gradient2(
    low = "#f0f0f0", high = "#1b7837",
    midpoint = 0.45,
    name = "Probability") +
  labs(
    y = "Occupancy trend (slope per year)",
    title = "Temporal occupancy trends") +
  guides(shape = 'none') +
  coord_flip() +
  theme_bw(base_size = 14) + theme(legend.position = "right",
                                   axis.title.y = element_blank(),
                                   strip.background = element_rect(fill = NA),
                                   panel.grid.minor = element_blank())
plot1b_cool

#icons a different way
plot1c_cool <- ggplot(posthoc_results_df[posthoc_results_df$season == 'cool',], 
                      aes(x = species, y = slope_mean, 
                          ymin = slope_lower_CI, ymax = slope_upper_CI,
                          shape = fence)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_pointrange(aes(color = slope_mean > 0, alpha = f_score, size = f_score),  
                  # size = 1,
                  position = position_dodge(width = 0.5)) +
  facet_wrap(~fence) +
  # scale_shape_manual(values = c("north" = 24, "south" = 25)) + 
  # scale_fill_manual(
  #   values = c("decreasing" = 'white'),
  #   labels = c("Decreasing"),
  #   name = "Trend direction"
  # ) +
  scale_alpha_continuous(range = c(0.3, 1),
                         guide = 'none') +
  scale_color_manual(
    values = c("FALSE" = "#8B4513", "TRUE" = "#1f78b4"),
    labels = c("Decreasing", "Increasing")
  ) +
  scale_size_continuous(
    range = c(0.5, 2),
  ) +
  labs(
    y = "Occupancy trend (slope per year)",
    # title = "Temporal occupancy trends",
    color = "Trend direction",
    alpha = 'Probability',
    size = 'Probability') +
  guides(shape = 'none') +
  coord_flip() +
  theme_bw(base_size = 14) + theme(legend.position = "right",
                                   legend.text = element_text(size = 12),
                                   axis.title.y = element_blank(),
                                   strip.background = element_rect(fill = NA),
                                   panel.grid.minor = element_blank())
plot1c_cool

ggsave('H_sp_occ/figures/111025/psi_trend_mod17ms_cool.png', plot1c_cool, width = 7, height = 5)


#Option 2: F-scores colored by slope values (like Doser et al. 2022)
species_order_cool <- posthoc_results_df %>%
  filter(season == "cool", fence == "south") %>%
  arrange(f_score) %>%
  pull(species)

posthoc_results_df$species <- factor(posthoc_results_df$species,
                                     levels = species_order_cool)

plot2 <- ggplot(posthoc_results_df, aes(x = species, y = f_score, color = slope_mean, shape = fence)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_point(size = 3) +
  facet_wrap(~season) +
  scale_color_gradient2(low = "red", mid = "gray80", high = "blue",
                        name = "Mean slope") +
  labs(x = "Season",
       y = "Probability of increasing trend",
       title = "Temporal occupancy trends") +
  coord_flip(ylim = c(0, 1)) +
  theme_bw(base_size = 14) +
  theme(
    legend.position = "right",
    panel.grid.minor = element_blank()
  )
plot2

ggsave('H_sp_occ/figures/111025/psi_Fscore_mod17ms.png', plot2, width = 7, height = 5)


#Option 3: plot the occupancy estimates directly with their uncertainty 
# with geom_smooth here (but it's NOT the same as the results from postHocLM)
vline_data <- data.frame(x = 2022,
                         species = c('buffalo','eland','impala','kudu','sable','warthog','waterbuck','zebra'))

plot3_wet <- ggplot(psi_long_df[psi_long_df$season == 'wet',], 
                    aes(x = as.numeric(year), y = `50%`, color = fence, fill = fence, shape = fence)) +
  geom_vline(data = vline_data, aes(xintercept = x), color = 'gray40', 
             linetype = 'dotdash', linewidth = 1) + #mark 2022 for secondary translocations
  geom_pointrange(aes(ymin = `2.5%`, ymax = `97.5%`),
                  position = position_dodge(width = 0.6),
                  show.legend = FALSE) +
  geom_point(size = 3, position = position_dodge(width = 0.6)) +
  # geom_ribbon(aes(ymin = `2.5%`, ymax = `97.5%`), position = position_dodge(width = 0.6), alpha = 0.2, color = NA) +
  geom_smooth(aes(y = `50%`, linetype = fence), method = "lm", se = TRUE, level = 0.95, 
              position = position_dodge(width = 0.6), alpha = 0.2) +
  scale_x_continuous(breaks = 2019:2024) +
  facet_wrap(~species) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(
    x = "Year",
    y = "Occupancy probability (ψ ± 95% CI)",
    color = "Fence status",
    fill = 'Fence status',
    shape = 'Fence status',
    linetype = 'Fence status'
    # title = 'Elephant'
  ) +
  scale_color_manual(values = c(south = "#008080", north = "#e69f00")) +
  scale_fill_manual(values = c(south = "#008080", north = "#e69f00")) +
  theme_bw(base_size = 14) + theme(axis.text.x = element_text(angle = 45, hjust = 1),
                                   axis.title.x = element_blank(),
                                   legend.position = 'right',
                                   # legend.position.inside = c(0.9, 0.1),
                                   legend.title = element_blank(),
                                   legend.text = element_text(size = 14),
                                   legend.background = element_blank(),
                                   strip.background = element_rect(fill = NA))
plot3_wet

ggsave('H_sp_occ/figures/111025/psi_lm_geomSmooth_mod17ms_wet.png', plot3_wet, width = 7, height = 5)

#Option 4: plot the occupancy estimates directly with their uncertainty 
# but just the uncertainty in point estimates
plot4_cool <- ggplot(psi_long_df[psi_long_df$season == 'cool',], 
                    aes(x = as.numeric(year), y = `50%`, color = fence, fill = fence, shape = fence)) +
  geom_vline(data = vline_data, aes(xintercept = x), color = 'gray40', 
             linetype = 'dotdash', linewidth = 1) + #mark 2022 for secondary translocations
  geom_line(linewidth = 1, aes(linetype = fence), position = position_dodge(width = 0.6)) +
  geom_pointrange(aes(ymin = `2.5%`, ymax = `97.5%`),
                  position = position_dodge(width = 0.6),
                  show.legend = FALSE) +
  geom_point(size = 3, position = position_dodge(width = 0.6)) +
  # geom_ribbon(aes(ymin = `2.5%`, ymax = `97.5%`), position = position_dodge(width = 0.6), alpha = 0.2, color = NA) +
  scale_x_continuous(breaks = 2019:2024) +
  facet_wrap(~species) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(
    x = "Year",
    y = "Occupancy probability (ψ ± 95% CI)",
    color = "Fence status",
    fill = 'Fence status',
    shape = 'Fence status',
    linetype = 'Fence status'
    # title = 'Elephant'
  ) +
  scale_color_manual(values = c(south = "#008080", north = "#e69f00")) +
  scale_fill_manual(values = c(south = "#008080", north = "#e69f00")) +
  theme_bw(base_size = 14) + theme(axis.text.x = element_text(angle = 45, hjust = 1),
                                   axis.title.x = element_blank(),
                                   legend.position = 'right',
                                   # legend.position.inside = c(0.9, 0.1),
                                   legend.title = element_blank(),
                                   legend.text = element_text(size = 14),
                                   legend.background = element_blank(),
                                   strip.background = element_rect(fill = NA))
plot4_cool

ggsave('H_sp_occ/figures/111025/psi_est_CI_mod17ms_cool.png', plot4_cool, width = 7, height = 5) #with geompointrange but no ribbon
# ggsave('H_sp_occ/figures/111025/psi_est_ribbon_mod17.png', plot4_wet, width = 7, height = 5) #with ribbon but no geompointrange

## FOR PRESENTATION: INDIVIDUAL SPECIES
plot4_cool_sp <- ggplot(psi_long_df[psi_long_df$season == 'cool' & psi_long_df$species == 'waterbuck',], 
                     aes(x = as.numeric(year), y = `50%`, color = fence, fill = fence, shape = fence)) +
  geom_vline(data = vline_data, aes(xintercept = x), color = 'gray40', 
             linetype = 'dotdash', linewidth = 1) + #mark 2022 for secondary translocations
  geom_line(linewidth = 1, aes(linetype = fence), position = position_dodge(width = 0.6)) +
  geom_pointrange(aes(ymin = `2.5%`, ymax = `97.5%`),
                  position = position_dodge(width = 0.6),
                  show.legend = FALSE) +
  geom_point(size = 3, position = position_dodge(width = 0.6)) +
  # geom_ribbon(aes(ymin = `2.5%`, ymax = `97.5%`), position = position_dodge(width = 0.6), alpha = 0.2, color = NA) +
  scale_x_continuous(breaks = 2019:2024) +
  # facet_wrap(~species) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(
    x = "Year",
    # y = "Occupancy probability (ψ ± 95% CI)",
    y = NULL,
    color = "Fence status",
    fill = 'Fence status',
    shape = 'Fence status',
    linetype = 'Fence status',
    title = 'Waterbuck'
  ) +
  scale_color_manual(values = c(south = "#008080", north = "#e69f00")) +
  scale_fill_manual(values = c(south = "#008080", north = "#e69f00")) +
  theme_bw(base_size = 14) + theme(axis.text.x = element_text(angle = 45, hjust = 1),
                                   axis.title.x = element_blank(),
                                   legend.position = 'none',
                                   # legend.position.inside = c(0.9, 0.1),
                                   legend.title = element_blank(),
                                   legend.text = element_text(size = 14),
                                   legend.background = element_blank(),
                                   strip.background = element_rect(fill = NA))
plot4_cool_sp
ggsave('H_sp_occ/figures/111025/psi_est_postHocLM_mod17ms_WATERBUCK.png', plot4_cool_sp, width = 2, height = 2) #with ribbon but no geompointrange




#Option 5:, make some predictions from the LM results and plot them
year_range <- seq(min(as.numeric(psi_long_df$year)), max(as.numeric(psi_long_df$year)), by = 1)

trend_lines <- posthoc_results_df %>%
  tidyr::expand_grid(year_numeric = 0:(length(year_range)-1)) %>%
  mutate(pred = intercept_mean + slope_mean * year_numeric,
         pred_lower = intercept_lower_CI + slope_lower_CI * year_numeric,
         pred_upper = intercept_upper_CI + slope_upper_CI * year_numeric) %>%
  mutate(year = year_numeric + min(as.numeric(psi_long_df$year)))

plot5_wet <- ggplot(psi_long_df[psi_long_df$season == 'wet',], 
                    aes(x = as.numeric(year), y = `50%`, color = fence, fill = fence, shape = fence)) +
  geom_vline(data = vline_data, aes(xintercept = x), color = 'gray40', 
             linetype = 'dotdash', linewidth = 1) + #mark 2022 for secondary translocations
  # geom_line(size = 1, aes(linetype = fence)) +
  geom_point(size = 2) +
  # geom_smooth(se = FALSE, color = "gray60", linetype = "dashed") +
  geom_line(
    data = trend_lines[trend_lines$season == 'wet',],
    aes(y = pred, color = fence, group = interaction(fence, season)),
    linewidth = 1.2
  ) +
  geom_ribbon(
    data = trend_lines[trend_lines$season == 'wet',],
    aes(ymin = pred_lower, y = pred, ymax = pred_upper, fill = fence, group = interaction(fence, season)),
    alpha = 0.2,
    color = NA
  ) +
  # geom_smooth(aes(y = `50%`), method = "lm", se = TRUE, level = 0.95, alpha = 0.2) +
  facet_wrap(~species) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(
    x = "Year",
    y = "Predicted occupancy probability (ψ ± 95% CI)",
    color = "Fence status",
    fill = 'Fence status',
    linetype = 'Fence status',
    shape = 'Fence status'
  ) +
  scale_color_manual(values = c(south = "#008080", north = "#e69f00")) +
  scale_fill_manual(values = c(south = "#008080", north = "#e69f00")) +
  theme_bw(base_size = 14) + theme(axis.text.x = element_text(angle = 45, hjust = 1),
                                   axis.title.x = element_blank(),
                                   legend.position = 'right',
                                   # legend.position.inside = c(0.9, 0.1),
                                   legend.title = element_blank(),
                                   legend.text = element_text(size = 14),
                                   legend.background = element_blank(),
                                   strip.background = element_rect(fill = NA))
plot5_wet

ggsave('H_sp_occ/figures/111025/psi_est_postHocLM_mod17ms.png', plot5_wet, width = 7, height = 5) #with ribbon but no geompointrange

## Combine estimated psi + slope values
plot6 <- (plot4_cool + theme(legend.position = 'top')) + 
  (plot1c_cool + guides(size = 'none') + labs(y = 'Annual trend') +
     theme(legend.position = 'top', 
           legend.title = element_blank())) +
  plot_annotation(tag_levels = 'A') +
  plot_layout(widths = c(2, 1))
plot6

ggsave('H_sp_occ/figures/111025/mod17ms_psi_est_trend_combined_cool.png', plot6, width = 9, height = 6)


## TO REPORT FOR TABLE 2 (cool season) AND APPENDIX S3 (wet/dry) ---------------
#cool
posthoc_results_df %>% select(species, fence, season, trend, f_score, slope_mean, slope_lower_CI, slope_upper_CI, r2) %>%
  filter(season == 'cool', fence == 'south') %>%
  select(season, species, fence, slope_mean, slope_lower_CI, slope_upper_CI, r2, f_score, trend) %>%
  mutate(across(where(is.numeric), ~ round(., 2)))

posthoc_results_df %>% select(species, fence, season, trend, f_score, slope_mean, slope_lower_CI, slope_upper_CI, r2) %>%
  filter(season == 'cool', fence == 'north') %>%
  select(season, species, fence, slope_mean, slope_lower_CI, slope_upper_CI, r2, f_score, trend) %>%
  mutate(across(where(is.numeric), ~ round(., 2)))

#dry
posthoc_results_df %>% select(species, fence, season, trend, f_score, slope_mean, slope_lower_CI, slope_upper_CI, r2) %>%
  filter(season == 'dry', fence == 'south') %>%
  select(season, species, fence, slope_mean, slope_lower_CI, slope_upper_CI, r2, f_score, trend) %>%
  mutate(across(where(is.numeric), ~ round(., 2)))

posthoc_results_df %>% select(species, fence, season, trend, f_score, slope_mean, slope_lower_CI, slope_upper_CI, r2) %>%
  filter(season == 'dry', fence == 'north') %>%
  select(season, species, fence, slope_mean, slope_lower_CI, slope_upper_CI, r2, f_score, trend) %>%
  mutate(across(where(is.numeric), ~ round(., 2)))

#wet
posthoc_results_df %>% select(species, fence, season, trend, f_score, slope_mean, slope_lower_CI, slope_upper_CI, r2) %>%
  filter(season == 'wet', fence == 'south') %>%
  select(season, species, fence, slope_mean, slope_lower_CI, slope_upper_CI, r2, f_score, trend) %>%
  mutate(across(where(is.numeric), ~ round(., 2)))

posthoc_results_df %>% select(species, fence, season, trend, f_score, slope_mean, slope_lower_CI, slope_upper_CI, r2) %>%
  filter(season == 'wet', fence == 'north') %>%
  select(season, species, fence, slope_mean, slope_lower_CI, slope_upper_CI, r2, f_score, trend) %>%
  mutate(across(where(is.numeric), ~ round(., 2)))


## Dry and wet season plots for Appendix S3 ------------------------------------

#TRENDS

#Dry
species_order_dry <- posthoc_results_df %>%
  filter(season == "dry", fence == "south") %>%
  arrange(slope_mean) %>%
  pull(species)

posthoc_results_df$species <- factor(posthoc_results_df$species,
                                     levels = species_order_dry)

plot1c_dry <- ggplot(posthoc_results_df[posthoc_results_df$season == 'dry',], 
                      aes(x = species, y = slope_mean, 
                          ymin = slope_lower_CI, ymax = slope_upper_CI,
                          shape = fence)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_pointrange(aes(color = slope_mean > 0, alpha = f_score, size = f_score),  
                  position = position_dodge(width = 0.5)) +
  facet_wrap(~fence) +
  scale_alpha_continuous(range = c(0.3, 1),
                         guide = 'none') +
  scale_color_manual(
    values = c("FALSE" = "#8B4513", "TRUE" = "#1f78b4"),
    labels = c("Decreasing", "Increasing")
  ) +
  scale_size_continuous(
    range = c(0.5, 2),
  ) +
  labs(
    y = "Occupancy trend (slope per year)",
    color = "Trend direction",
    alpha = 'Probability',
    size = 'Probability') +
  guides(shape = 'none') +
  coord_flip() +
  theme_bw(base_size = 14) + theme(legend.position = "right",
                                   legend.text = element_text(size = 12),
                                   axis.title.y = element_blank(),
                                   strip.background = element_rect(fill = NA),
                                   panel.grid.minor = element_blank())

#Wet
species_order_wet <- posthoc_results_df %>%
  filter(season == "wet", fence == "south") %>%
  arrange(slope_mean) %>%
  pull(species)

posthoc_results_df$species <- factor(posthoc_results_df$species,
                                     levels = species_order_wet)

plot1c_wet <- ggplot(posthoc_results_df[posthoc_results_df$season == 'wet',], 
                     aes(x = species, y = slope_mean, 
                         ymin = slope_lower_CI, ymax = slope_upper_CI,
                         shape = fence)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_pointrange(aes(color = slope_mean > 0, alpha = f_score, size = f_score),  
                  position = position_dodge(width = 0.5)) +
  facet_wrap(~fence) +
  scale_alpha_continuous(range = c(0.3, 1),
                         guide = 'none') +
  scale_color_manual(
    values = c("FALSE" = "#8B4513", "TRUE" = "#1f78b4"),
    labels = c("Decreasing", "Increasing")
  ) +
  scale_size_continuous(
    range = c(0.5, 2),
  ) +
  labs(
    y = "Occupancy trend (slope per year)",
    color = "Trend direction",
    alpha = 'Probability',
    size = 'Probability') +
  guides(shape = 'none') +
  coord_flip() +
  theme_bw(base_size = 14) + theme(legend.position = "right",
                                   legend.text = element_text(size = 12),
                                   axis.title.y = element_blank(),
                                   strip.background = element_rect(fill = NA),
                                   panel.grid.minor = element_blank())

#yearly psi
plot4_dry <- ggplot(psi_long_df[psi_long_df$season == 'dry',], 
                     aes(x = as.numeric(year), y = `50%`, color = fence, fill = fence, shape = fence)) +
  geom_vline(data = vline_data, aes(xintercept = x), color = 'gray40', 
             linetype = 'dotdash', linewidth = 1) + #mark 2022 for secondary translocations
  geom_line(linewidth = 1, aes(linetype = fence), position = position_dodge(width = 0.6)) +
  geom_pointrange(aes(ymin = `2.5%`, ymax = `97.5%`),
                  position = position_dodge(width = 0.6),
                  show.legend = FALSE) +
  geom_point(size = 3, position = position_dodge(width = 0.6)) +
  scale_x_continuous(breaks = 2019:2024) +
  facet_wrap(~species) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(
    x = "Year",
    y = "Occupancy probability (ψ ± 95% CI)",
    color = "Fence status",
    fill = 'Fence status',
    shape = 'Fence status',
    linetype = 'Fence status'
  ) +
  scale_color_manual(values = c(south = "#008080", north = "#e69f00")) +
  scale_fill_manual(values = c(south = "#008080", north = "#e69f00")) +
  theme_bw(base_size = 14) + theme(axis.text.x = element_text(angle = 45, hjust = 1),
                                   axis.title.x = element_blank(),
                                   legend.position = 'right',
                                   # legend.position.inside = c(0.9, 0.1),
                                   legend.title = element_blank(),
                                   legend.text = element_text(size = 14),
                                   legend.background = element_blank(),
                                   strip.background = element_rect(fill = NA))

plot4_wet <- ggplot(psi_long_df[psi_long_df$season == 'wet',], 
                     aes(x = as.numeric(year), y = `50%`, color = fence, fill = fence, shape = fence)) +
  geom_vline(data = vline_data, aes(xintercept = x), color = 'gray40', 
             linetype = 'dotdash', linewidth = 1) + #mark 2022 for secondary translocations
  geom_line(linewidth = 1, aes(linetype = fence), position = position_dodge(width = 0.6)) +
  geom_pointrange(aes(ymin = `2.5%`, ymax = `97.5%`),
                  position = position_dodge(width = 0.6),
                  show.legend = FALSE) +
  geom_point(size = 3, position = position_dodge(width = 0.6)) +
  scale_x_continuous(breaks = 2019:2024) +
  facet_wrap(~species) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(
    x = "Year",
    y = "Occupancy probability (ψ ± 95% CI)",
    color = "Fence status",
    fill = 'Fence status',
    shape = 'Fence status',
    linetype = 'Fence status'
  ) +
  scale_color_manual(values = c(south = "#008080", north = "#e69f00")) +
  scale_fill_manual(values = c(south = "#008080", north = "#e69f00")) +
  theme_bw(base_size = 14) + theme(axis.text.x = element_text(angle = 45, hjust = 1),
                                   axis.title.x = element_blank(),
                                   legend.position = 'right',
                                   # legend.position.inside = c(0.9, 0.1),
                                   legend.title = element_blank(),
                                   legend.text = element_text(size = 14),
                                   legend.background = element_blank(),
                                   strip.background = element_rect(fill = NA))


#combine
plot6_dry <- (plot4_dry + theme(legend.position = 'top')) + 
  (plot1c_dry + guides(size = 'none') + labs(y = 'Annual trend') +
     theme(legend.position = 'top', 
           legend.title = element_blank())) +
  plot_annotation(tag_levels = 'A', title = 'Dry season') +
  plot_layout(widths = c(2, 1))
plot6_dry
ggsave('H_sp_occ/figures/111025/mod17ms_psi_est_trend_combined_dry.png', plot6_dry, width = 9, height = 6)

plot6_wet <- (plot4_wet + theme(legend.position = 'top')) + 
  (plot1c_wet + guides(size = 'none') + labs(y = 'Annual trend') +
     theme(legend.position = 'top', 
           legend.title = element_blank())) +
  plot_annotation(tag_levels = 'A', title = 'Wet season') +
  plot_layout(widths = c(2, 1))
plot6_wet
ggsave('H_sp_occ/figures/111025/mod17ms_psi_est_trend_combined_wet.png', plot6_wet, width = 9, height = 6)
