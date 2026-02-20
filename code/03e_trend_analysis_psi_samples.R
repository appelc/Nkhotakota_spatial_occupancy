# Deriving yearly occupancy using psi.samples directly

library(data.table)
library(tidyverse)
library(spOccupancy)


## Read in model outputs -------------------------------------------------------
# output_files <- list.files('H_sp_occ/outputs/model17/', pattern = 'RDS', full.names = TRUE)
output <- readRDS('H_sp_occ/outputs/model17/model17_elephant_allChains.RDS') #example for now
str(output)


## Read in input data ----------------------------------------------------------
sp_inputs <- readRDS('H_sp_occ/sp_inputs.RDS')

input <- sp_inputs$elephant #example for now
str(input)


## Extract and format psi samples ----------------------------------------------

#get site names, time period names, and north/south from input file
site_names <- input$occ.covs$site
yearSeason_names <- input$occ.covs$year_season[1,]
fence_names <- data.frame(site = input$occ.covs$site,
                          fence = input$occ.covs$fence)

#extract psi samples and set dimnames
psi_samples <- output$psi.samples
dim(psi_samples)  

dimnames(psi_samples) <- list(
  iter = seq_len(dim(psi_samples)[1]),
  site = site_names,
  time = yearSeason_names
)  

dimnames(psi_samples)
head(psi_samples)

#subset north/south
north_sites <- fence_names$site[fence_names$fence == "North"]
south_sites <- fence_names$site[fence_names$fence == "South"]

psi_samples_north <- psi_samples[, north_sites, ]
psi_samples_south <- psi_samples[, south_sites, ]

dim(psi_samples_north) #28200 samples, 115 sites x  17 year_seasons
dim(psi_samples_south) #28200 samples, 84 sites x  17 year_seasons


## -----------------------------------------------------------------------------

## THREE APPROACHES TO TREND ANALYSIS:
### 1) psi ~ year_season (separately for North/South) -- e.g., change from each year_season to the next
#### -first split psi samples by north/south
#### -then summarize by year_season (mean and quantiles across all sites for each year_season)

### 2) psi ~ year (separately for North/South and Wet/Cool/Dry) -- e.g., change from 2019_wet to 2020_wet to 2021_wet, etc.
#### -first split psi samples by north/south
#### -then split samples by year
#### -then summarize by year and season (mean and quantiles across all sites for each year and season)
#### -then plot trends separately for Wet/Cool/Dry

### 3) psi ~ year (separately for North/South, averaging across seasons within year) -- e.g., change from 2019 to 2020 to 2021, etc.
#### -first split psi samples by north/south
#### -then summarize by year (mean and quantiles across all sites for each year, averaging across seasons within year) **this part could be tricky**
#### -then plot trends


## 1) Summarize by year_season -------------------------------------------------

#for each mcmc sample, calculate the mean and quantiles across all sites for each year_season
psi_mean_north <- apply(psi_samples_north, c(1, 3), mean)
psi_mean_south <- apply(psi_samples_south, c(1, 3), mean)

dim(psi_mean_north); dim(psi_mean_south) #samples * year_seasons

psi_quants_north <- apply(psi_mean_north, 2, quantile, c(0.025, 0.5, 0.975))
psi_quants_south <- apply(psi_mean_south, 2, quantile, c(0.025, 0.5, 0.975))

dim(psi_quants_north); dim(psi_quants_south) #3 quantiles * year_seasons

psi_quants_north

#pivot long
psi_long_north <- as.data.frame(t(psi_quants_north)) %>%
  mutate(year_season = rownames(.))
psi_long_south <- as.data.frame(t(psi_quants_south)) %>%
  mutate(year_season = rownames(.))

#perform linear regression and store outputs to report

#with lm()

#try with postHocLM()



#combine for plotting
psi_long_combined <- bind_rows(
  psi_long_north %>% mutate(fence = "North"),
  psi_long_south %>% mutate(fence = "South")
)

#set year_season as factor with correct order
psi_long_combined$year_season <- factor(psi_long_combined$year_season, levels = yearSeason_names)

#make dataframe to plot vertical line at 2022 for all sp. except elephants
# vline_data <- data.frame(x = '2022_cool',
# species = c('buffalo','eland','impala','kudu','sable','warthog','waterbuck','zebra'))

#plot
plot1 <- ggplot(psi_long_combined, aes(x = year_season, y = `50%`, color = fence, shape = fence)) +
  # geom_vline(data = vline_data, aes(xintercept = x), color = 'gray40',
             # linetype = 'dotdash', linewidth = 1) + #mark 2022 for secondary translocations
  geom_vline(aes(xintercept = '2022_cool'), color = 'gray40',
             linetype = 'dotdash', linewidth = 1) + #mark 2022 for secondary translocations
  geom_pointrange(aes(ymin = `2.5%`, ymax = `97.5%`),  #, linetype = fence
                  alpha = 0.5,
                  position = position_dodge(width = 0.6),
                  show.legend = FALSE) +
  geom_point(size = 3, position = position_dodge(width = 0.6)) +
  geom_smooth(method = "lm", se = TRUE, level = 0.95, alpha = 0.2,
              aes(x = as.numeric(year_season), y = `50%`, fill = fence, linetype = fence)) +
  # scale_x_continuous(breaks = 2019:2024) +
  scale_color_manual(values = c(South = "#008080", North = "#e69f00")) +
  scale_fill_manual(values = c(South = "#008080", North = "#e69f00")) +
  # facet_wrap(~species, nrow = 3) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = NULL, y = "Occupancy Probability (ψ ± 95% CI)",
       title = 'Elephant') +
  theme_bw(base_size = 14) + theme(axis.text.x = element_text(angle = 45, hjust = 1),
                                   legend.position = 'inside',
                                   legend.position.inside = c(0.8, 0.1),
                                   legend.title = element_blank(),
                                   strip.background = element_rect(fill = NA))
plot1


## 2) Summarize by year (separate seasons) -------------------------------------

#get year_season names
cool_names <- yearSeason_names[grepl("_cool$", yearSeason_names)]
dry_names <- yearSeason_names[grepl("_dry$", yearSeason_names)]
wet_names <- yearSeason_names[grepl("_wet$", yearSeason_names)]

#subset by season AND fence
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
dim(yearly_matrices$north_cool) #each will be iter x years

#create covariate dataframe for regression
(covs <- data.frame(year = 0:ncol(yearly_matrices[[1]])))

#set priors and initial values (I don't have RE so no need for sigma.sq)
inits_list <- list(beta = 0, tau.sq = 1)
priors_list <- list(beta.normal = list(mean = 0, var = 2.7),
               tau.sq.ig = c(0.001, 0.001))

#initialize objects to store results
psi_long_df <- data.frame()
posthoc_results_list <- list()
posthoc_results_df <- data.frame()

#for each of the 6 fence_season combinations, run postHocLM and store results:
for(aa in names(yearly_matrices)){
  
  #extract the matrix
  data <- yearly_matrices[[aa]]
  
  #summarize and store quantiles across iterations for each year_season (but we'll do regression on full samples below)
  psi_quants <- apply(data, 2, quantile, c(0.025, 0.5, 0.975))
  psi_long <- as.data.frame(t(psi_quants)) %>%
    mutate(year_season = rownames(.),
           fence_season = aa) %>%
    separate(year_season, into = c("year", "season"), sep = "_", remove = FALSE) %>%
    separate(fence_season, into = c("fence", "season2"), sep = "_", remove = FALSE)
  psi_long_df <- rbind(psi_long_df,
                       psi_long)
  
  #create input list for regression (w/ matrix and covariates)
  input <- list(y = data,
                covs = covs)
  
  #run post-hoc regression on psi samples
  posthoc_result <- postHocLM(
    formula = ~ year,
    data = input,
    inits = inits_list,
    priors = priors_list,
    n.chains = 3,
    n.report = 5000,
    verbose = TRUE
  )
  # print(summary(posthoc_result))
  
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
                                fence = fence_aa,
                                season = season_aa,
                                trend = trend_aa,
                                f_score = f_score_aa,
                                intercept_mean = mean(posthoc_result$beta.samples[,1]),
                                slope_mean = mean(posthoc_result$beta.samples[,2]),
                                slope_median = median(posthoc_result$beta.samples[,2]),
                                slope_lower_CI = quantile(posthoc_result$beta.samples[,2], 0.025),
                                slope_upper_CI = quantile(posthoc_result$beta.samples[,2], 0.975),
                                r2 = mean(posthoc_result$bayes.R2)
                              ))
  
  posthoc_results_list[[aa]] <- posthoc_result
}

#QC
# summary(posthoc_results_list$north_cool)
# quantile(posthoc_results_list$north_cool$beta.samples[,2], probs = c(0.025, 0.5, 0.975))
# posthoc_results_df
# 
# dim(posthoc_results_list$north_cool$y.hat.samples) #(28200 samples * 3 chains) x 6 year_seasons

head(posthoc_results_df) #regression slope, intercept, probability of increase, R2, etc.
head(psi_long_df) #predicted occupancy quantiles by year and season

## Several ways to plot:

#plot slope values colored by F-scores
plot2 <- ggplot(posthoc_results_df, aes(x = season, y = slope_mean, 
                          ymin = slope_lower_CI, ymax = slope_upper_CI,
                          color = f_score)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_pointrange(size = 0.8, position = position_dodge(width = 0.5)) +
  facet_wrap(~ fence) +
  scale_color_gradient2(
    low = "red", mid = "gray80", high = "blue", midpoint = 0.5,
    name = "P(increasing)") +
  labs(
    x = "Season",
    y = "Occupancy trend (slope per year)",
    title = "Temporal occupancy trends") +
  theme_bw(base_size = 14) + theme(legend.position = "right",
                                   panel.grid.minor = element_blank())
plot2

#or F-scores colored by slope values (like Doser et al. 2022_
plot3 <- ggplot(posthoc_results_df, aes(x = season, y = f_score, color = slope_mean, shape = fence)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_point(size = 3) +
  facet_wrap(~ fence) +
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
plot3

#or plot the occupancy estimates directly with their uncertainty -- with or without geom_smooth here (that's NOT the same as the results from postHocLM)
plot4 <- ggplot(psi_long_df, aes(x = as.numeric(year), y = `50%`, color = fence, fill = fence)) +
  # geom_line(size = 1, aes(linetype = fence)) +
  # geom_point(size = 2) +
  geom_pointrange(aes(ymin = `2.5%`, ymax = `97.5%`),
                  alpha = 0.5,
                  position = position_dodge(width = 0.6),
                  show.legend = FALSE) +
  # geom_ribbon(aes(ymin = `2.5%`, ymax = `97.5%`), position = position_dodge(width = 0.6), alpha = 0.2, color = NA) +
  geom_smooth(aes(y = `50%`), method = "lm", se = TRUE, level = 0.95, alpha = 0.2) +
  facet_wrap(~ season) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(
    x = "Year",
    y = "Predicted occupancy probability (ψ ± 95% CI)",
    color = "Fence status",
    fill = 'Fence status',
    linetype = 'Fence status',
    title = 'Elephant'
  ) +
  scale_color_manual(values = c(south = "#008080", north = "#e69f00")) +
  scale_fill_manual(values = c(south = "#008080", north = "#e69f00")) +
  theme_bw(base_size = 14) + theme(axis.text.x = element_text(angle = 45, hjust = 1),
                                   axis.title.x = element_blank(),
                                   legend.position = 'inside',
                                   legend.position.inside = c(0.9, 0.1),
                                   legend.title = element_blank(),
                                   legend.text = element_text(size = 14),
                                   legend.background = element_blank(),
                                   strip.background = element_rect(fill = NA))
plot4

#OR, make some predictions from the LM results and plot them
year_range <- seq(min(as.numeric(psi_long_df$year)), max(as.numeric(psi_long_df$year)), by = 1)

trend_lines <- posthoc_results_df %>%
  tidyr::expand_grid(year_numeric = 0:(length(year_range)-1)) %>%
  mutate(pred = intercept_mean + slope_mean * year_numeric,
         pred_lower = intercept_mean + slope_lower_CI * year_numeric,
         pred_upper = intercept_mean + slope_upper_CI * year_numeric) %>%
  mutate(year = year_numeric + min(as.numeric(psi_long_df$year)))

plot5 <- ggplot(psi_long_df, aes(x = as.numeric(year), y = `50%`, color = fence, fill = fence)) +
  # geom_line(size = 1, aes(linetype = fence)) +
  geom_point(size = 2) +
  # geom_smooth(se = FALSE, color = "gray60", linetype = "dashed") +
  geom_line(
    data = trend_lines,
    aes(y = pred, color = fence, group = interaction(fence, season)),
    linewidth = 1.2
  ) +
  geom_ribbon(
    data = trend_lines,
    aes(ymin = pred_lower, y = pred, ymax = pred_upper, fill = fence, group = interaction(fence, season)),
    alpha = 0.2,
    color = NA
  ) +
  # geom_smooth(aes(y = `50%`), method = "lm", se = TRUE, level = 0.95, alpha = 0.2) +
  facet_wrap(~ season) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(
    x = "Year",
    y = "Predicted occupancy probability (ψ ± 95% CI)",
    color = "Fence status",
    fill = 'Fence status',
    linetype = 'Fence status'
  ) +
  scale_color_manual(values = c(south = "#008080", north = "#e69f00")) +
  scale_fill_manual(values = c(south = "#008080", north = "#e69f00")) +
  theme_bw(base_size = 14) + theme(axis.text.x = element_text(angle = 45, hjust = 1),
                                   axis.title.x = element_blank(),
                                   legend.position = 'inside',
                                   legend.position.inside = c(0.9, 0.1),
                                   legend.title = element_blank(),
                                   legend.text = element_text(size = 14),
                                   legend.background = element_blank(),
                                   strip.background = element_rect(fill = NA))
plot5




## 3) Summarize by year (averaging across seasons) --------------------------------

#need an extra level of summarizing to get mean across seasons within year





  
  
  
  