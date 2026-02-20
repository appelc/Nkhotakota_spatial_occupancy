## Create maps using predictions generated at new sites from stPGOcc model

library(tidyverse)
library(spOccupancy)
library(data.table)
library(stars)
library(patchwork)
library(rphylopic)
library(cowplot)

#set resolution for predictions (one of 30, 100, 1000) -- note that the manuscript used 100 m
res_meters <- '1000' 


## Prep ------------------------------------------------------------------------

#read in coordinates
if (res_meters == 30) {cov_coords <- fread('data/raw/covariates_resampled_30.csv')}
if (res_meters == 100) {cov_coords <- fread('data/raw/covariates_resampled_100m.csv')}
if (res_meters == 1000) {cov_coords <- fread('data/raw/covariates_resampled_1000m.csv')}

if(!exists('cov_coords')) { message('No coordinates found. res_meters must be one of 30, 100, 1000') } else {
  message('Found coordinates')
}

#extract coords
coords.0 <- as.matrix(cov_coords[, c('x', 'y')])

#indicate which primary time periods (season_years) we have predictions for
predict.seasons <- c('2019_cool','2024_cool')


## Add phylopic IDs ------------------------------------------------------------

phylopic_ids <- data.frame('species' = c('buffalo',
                                         'eland',
                                         'elephant',
                                         'impala',
                                         'kudu',
                                         'sable',
                                         'warthog',
                                         'waterbuck',
                                         'zebra'),
                           'id' = c('65c4a9b3-dcde-4f0f-9a1f-8d71e74be9ec',
                                    '7db171af-ac7a-4859-9c2b-66488a5a5c95',
                                    '62398ac0-f0c3-48f8-8455-53512a05fbc4',
                                    'e07d1491-1d85-4c47-9f7d-075ea57bf0c5',
                                    'a3c50373-4fc4-44b2-9c09-fb9959f5893e',
                                    'a6c8cc22-c035-48e1-bd2a-98a29446b392',
                                    'd594a0c4-2708-4cde-ba41-08fde8b1184f',
                                    'f93103f1-e2a0-4c73-b274-c7b51afe4db0',
                                    '81caf94e-5cbe-4e5e-8101-545abea2bfc0'),
                           'size_scaling' = c(0.3,
                                              0.4,
                                              0.4,
                                              0.3,
                                              0.5,
                                              0.4,
                                              0.3,
                                              0.5,
                                              0.4))


## Make maps -------------------------------------------------------------------
model_name <- 'model01'
season <- '2024_cool' #choose one here

#read in predicted values for all species
(files <- list.files(paste('predictions/', model_name, sep = ''), pattern = 'RDS', full.names = TRUE, recursive = TRUE))

#create object to store maps
plots_list <- list()

#loop thru each species:
for (ff in files){
  
  #species name
  (species_name <- sapply(strsplit(ff, '/'), function(x) x[3]))
  
  #read in predicted values
  psi.quants <- readRDS(ff)
  
  #create dataframe for plotting
  plot.data <- NULL
  stars.data <- list()
  for (ss in predict.seasons){
    psi.quants.season <- psi.quants[,,ss]
    plot.data.season <- data.frame(
      x = coords.0[,1],
      y = coords.0[,2],
      mean.psi = psi.quants.season[3,], #0.50 quantile
      lci.psi = psi.quants.season[1,],  #0.025 quantile
      uci.psi = psi.quants.season[5,],  #0.975 quantile
      season = ss #ideally match to actual name
    )
    plot.data <- rbind(plot.data, plot.data.season)
    stars.data[[ss]] <- st_as_stars(plot.data.season, dims = c('x','y'))
  }
  
  #get just the season I want and convert to stars data
  plot.data.2024_cool <- plot.data %>% filter(season == season)
  stars.2024_cool <- st_as_stars(plot.data.2024_cool, dims = c('x','y'))
  
  #plot
  plot.2024_cool <- ggplot() +
    geom_stars(data = stars.2024_cool, aes(x = x, y = y, fill = mean.psi)) +
    scale_fill_viridis_c(name = 'median psi', na.value = 'transparent', direction = -1,
                         option = 'D') + #D or E
    # add_phylopic(uuid = phylopic_ids[phylopic_ids$species == species_name,]$id,
    #              x = 594000, y = 8570000, 
    #              height = 15000) + #10000
    ggtitle(species_name) +
    
    guides(fill = guide_colorbar(
      title.position = "left",   # move title above
      title.hjust = 0.5,        # center title
      barwidth = 10,            # width of the bar
      barheight = 0.5,          # height of the bar
      direction = "horizontal"  # make it horizontal
    )) +
    
    theme_void() + theme(plot.title = element_text(size = 14, face = 'bold', hjust = 0.5),
                         legend.position = 'bottom')
  plot.2024_cool
  
  #store
  plots_list[[species_name]] <- plot.2024_cool
  
}


## Plot together ---------------------------------------------------------------

#adding phylopic icons manually for control over placement

combined_2024cool <- (((plots_list$buffalo + add_phylopic(uuid = '65c4a9b3-dcde-4f0f-9a1f-8d71e74be9ec',
                                                        x = 594000, y = 8610000, height = 11000)) +
                         (plots_list$eland + add_phylopic(uuid = '7db171af-ac7a-4859-9c2b-66488a5a5c95',
                                                        x = 594000, y = 8610000, height = 13000)) +
                         (plots_list$elephant + add_phylopic(uuid = '62398ac0-f0c3-48f8-8455-53512a05fbc4',
                                                             x = 594000, y = 8608000, height = 11000)) +
                         plot_layout(nrow = 1)) /
                        
                      ((plots_list$impala + add_phylopic(uuid = 'e07d1491-1d85-4c47-9f7d-075ea57bf0c5',
                                                        x = 594000, y = 8610000, height = 10000)) +
                         (plots_list$kudu + add_phylopic(uuid = 'a3c50373-4fc4-44b2-9c09-fb9959f5893e',
                                     x = 594000, y = 8610000, height = 15000)) +
                         (plots_list$sable + add_phylopic(uuid = 'a6c8cc22-c035-48e1-bd2a-98a29446b392',
                                      x = 594000, y = 8610000, height = 15000)) +
                         plot_layout(nrow = 1)) /
                        
                      ((plots_list$warthog + add_phylopic(uuid = 'd594a0c4-2708-4cde-ba41-08fde8b1184f',
                                          x = 594000, y = 8610000, height = 10000)) +
                         (plots_list$waterbuck + add_phylopic(uuid = 'f93103f1-e2a0-4c73-b274-c7b51afe4db0',
                                            x = 594000, y = 8610000, height = 15000)) +
                         (plots_list$zebra + add_phylopic(uuid = '81caf94e-5cbe-4e5e-8101-545abea2bfc0',
                                                              x = 594000, y = 8610000, height = 10000)) +
                         plot_layout(nrow = 1))) +
  
                    plot_layout(guides = 'collect') & theme(legend.position = 'bottom')

combined_2024cool

## FIGURE 4
ggsave(paste('predictions/map_panel_', model_name, '_', season, '_', res_meters, 'm.tif', sep = ''),
       combined_2024cool, width = 4, height = 6)

