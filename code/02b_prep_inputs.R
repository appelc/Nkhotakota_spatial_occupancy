## Prep inputs for spOccupancy

library(tidyverse)
library(spOccupancy)


## Read in detection data and covariates ---------------------------------------

dh_arrays <- readRDS('A_data/cleaned/detection_arrays_yearly_070225.RDS')
  str(dh_arrays)

occ_covar <- readRDS('H_sp_occ/occ_covar.RDS') #added scaled (N/S) versions 082525; added log and scaled-log versions 092025
  str(occ_covar)

coords <- readRDS('H_sp_occ/coords.RDS')
  str(coords)
  
  
## Remove 2018 dry season & 2024 dry season ------------------------------------

#remove from detections
dimnames(dh_arrays[[1]])[[2]]
head(dh_arrays$aardvark) #view example

dh_arrays <- lapply(dh_arrays, function(x) {
  x[, !dimnames(x)[[2]] %in% c('2018_dry','2024_dry'), , drop = FALSE]
})

dimnames(dh_arrays[[1]])[[2]]
head(dh_arrays$aardvark,10) #check example
str(dh_arrays) #should now be 17 columns

#remove from covariates
str(occ_covar)
occ_covar$year_season <- occ_covar$year_season[,-c(1,19)]
occ_covar$season <- occ_covar$season[,-c(1,19)]
occ_covar$year <- occ_covar$year[,-c(1,19)]
str(occ_covar) #should now be 17 columns for the time-varying ones


## Summarize -------------------------------------------------------------------

str(dh_arrays[[1]])
head(dh_arrays[[1]],5) #rows=sites, cols=season_years, slices=weeks

(nsites <- dim(dh_arrays[[1]])[1])
(nseasons <- 3)
(nyear_seasons <- dim(dh_arrays[[1]])[2])
(nweeks <- dim(dh_arrays[[1]])[3])


## Find sites with no data -----------------------------------------------------

#pick one species (NAs should be the same for all)
sp <- 'elephant'
y <- dh_arrays[[sp]]

#check which sites have all NAs across all years and replicates
all_na_sites <- apply(y, 1, function(site) all(is.na(site)))
cat("sites with all missing data:", sum(all_na_sites), "\n") 
which(all_na_sites)

#identify sites to keep (those that have at least one non-NA observation)
sites_to_keep <- apply(y, 1, function(site) !all(is.na(site)))
cat("Keeping", sum(sites_to_keep), "out of", length(sites_to_keep), "sites\n") 

## THERE IS ONE WITH NO DATA NOW THAT WE'VE REMOVED 2024_dry (E13-C)


## Format OCC covariates -------------------------------------------------------

str(occ_covar)

#create dummy variables for 'season' and 'fence'
occ_covar$season_wet  <- ifelse(occ_covar$season == 'wet', 1, 0)
occ_covar$season_cool <- ifelse(occ_covar$season == 'cool', 1, 0)
occ_covar$season_dry  <- ifelse(occ_covar$season == 'dry', 1, 0)

occ_covar$fence_south <- ifelse(occ_covar$fence == 'South', 1, 0)
occ_covar$fence_north <- ifelse(occ_covar$fence == 'North', 1, 0)

str(occ_covar)

#now subset to remove fully NA site-year_season rows
occ_covar_clean <- lapply(occ_covar, function(cov) {
    if(is.matrix(cov)) {
      cov[sites_to_keep, ]  #for site * year matrices
    } else if(is.vector(cov)) {
      cov[sites_to_keep]    #for site-only vectors
    }
  })
str(occ_covar_clean) #should be length 199 now

#create numeric versions of SITE and HEX (*082125: RE-CODING TO START AT 0, NOT 1*)
occ_covar_clean$site_num        <- as.integer(0:(length(occ_covar_clean$site)-1))
occ_covar_clean$hex_num         <- as.integer(0:(length(occ_covar_clean$hex)-1))

#create numeric versions of YEAR and YEAR_SEASON 
#(need to make sure order is correct first)
year_season_levels <- as.factor(unique(occ_covar_clean$year_season))
occ_covar_clean$year_season_num <- matrix(match(occ_covar_clean$year_season, year_season_levels), 
                                          nrow = nrow(occ_covar_clean$year_season)) - 1 #subtract 1 so it starts at 0

year_levels <- levels(as.factor(unique(occ_covar_clean$year)))
occ_covar_clean$year_num <- matrix(match(occ_covar_clean$year, year_levels), 
                                   nrow = nrow(occ_covar_clean$year)) - 1

str(occ_covar_clean) #should be length 199 now


## Format DET covariates -------------------------------------------------------

det_covar_clean = list('season' = occ_covar_clean$season,
                       'season_wet' = occ_covar_clean$season_wet,
                       'season_cool' = occ_covar_clean$season_cool,
                       'season_dry' = occ_covar_clean$season_dry,
                       'site' = as.numeric(occ_covar_clean$site_num),
                       'fence' = occ_covar_clean$fence,
                       'fence_south' = occ_covar_clean$fence_south,
                       'fence_north' = occ_covar_clean$fence_north,
                       'year_season_num' = occ_covar_clean$year_season_num,
                       'year_season' = occ_covar_clean$year_season,
                       'year_num' = occ_covar_clean$year_num,
                       'year' = occ_covar_clean$year) #added 092125

str(det_covar_clean) #should be length 199 already

#create dummy variables for 'year' (so we can use as factor)
det_covar_clean$year_2019 <- ifelse(det_covar_clean$year == '2019', 1, 0)
det_covar_clean$year_2020 <- ifelse(det_covar_clean$year == '2020', 1, 0)
det_covar_clean$year_2021 <- ifelse(det_covar_clean$year == '2021', 1, 0)
det_covar_clean$year_2022 <- ifelse(det_covar_clean$year == '2022', 1, 0)
det_covar_clean$year_2023 <- ifelse(det_covar_clean$year == '2023', 1, 0)
det_covar_clean$year_2024 <- ifelse(det_covar_clean$year == '2024', 1, 0)

#QC
head(det_covar_clean$year)
head(det_covar_clean$year_2019)
head(det_covar_clean$year_2020)
head(det_covar_clean$year_2021)
head(det_covar_clean$year_2022)
head(det_covar_clean$year_2023)
head(det_covar_clean$year_2024)

#create dummy variables for 'year_season' (so we can use as factor)
det_covar_clean$yearsn_2019wet <- ifelse(det_covar_clean$year_season == '2019_wet', 1, 0)
det_covar_clean$yearsn_2019cool <- ifelse(det_covar_clean$year_season == '2019_cool', 1, 0)
det_covar_clean$yearsn_2019dry <- ifelse(det_covar_clean$year_season == '2019_dry', 1, 0)
det_covar_clean$yearsn_2020wet <- ifelse(det_covar_clean$year_season == '2020_wet', 1, 0)
det_covar_clean$yearsn_2020cool <- ifelse(det_covar_clean$year_season == '2020_cool', 1, 0)
det_covar_clean$yearsn_2020dry <- ifelse(det_covar_clean$year_season == '2020_dry', 1, 0)
det_covar_clean$yearsn_2021wet <- ifelse(det_covar_clean$year_season == '2021_wet', 1, 0)
det_covar_clean$yearsn_2021cool <- ifelse(det_covar_clean$year_season == '2021_cool', 1, 0)
det_covar_clean$yearsn_2021dry <- ifelse(det_covar_clean$year_season == '2021_dry', 1, 0)
det_covar_clean$yearsn_2022wet <- ifelse(det_covar_clean$year_season == '2022_wet', 1, 0)
det_covar_clean$yearsn_2022cool <- ifelse(det_covar_clean$year_season == '2022_cool', 1, 0)
det_covar_clean$yearsn_2022dry <- ifelse(det_covar_clean$year_season == '2022_dry', 1, 0)
det_covar_clean$yearsn_2023wet <- ifelse(det_covar_clean$year_season == '2023_wet', 1, 0)
det_covar_clean$yearsn_2023cool <- ifelse(det_covar_clean$year_season == '2023_cool', 1, 0)
det_covar_clean$yearsn_2023dry <- ifelse(det_covar_clean$year_season == '2023_dry', 1, 0)
det_covar_clean$yearsn_2024wet <- ifelse(det_covar_clean$year_season == '2024_wet', 1, 0)
det_covar_clean$yearsn_2024cool <- ifelse(det_covar_clean$year_season == '2024_cool', 1, 0)

#QC
head(det_covar_clean$year_season)
head(det_covar_clean$yearsn_2019wet)
head(det_covar_clean$yearsn_2019cool)
head(det_covar_clean$yearsn_2019dry)
head(det_covar_clean$yearsn_2024wet)
head(det_covar_clean$yearsn_2024cool) #etc.

str(det_covar_clean)


## Format coordinates ----------------------------------------------------------

#site-level
head(coords)
keep <- names(sites_to_keep[sites_to_keep == TRUE]) #sites with data
coords_clean <- coords %>% filter(Site %in% keep) 

str(coords_clean) #199

#save modified
# saveRDS(coords_clean, 'H_sp_occ/coords_cleaned.RDS')

coords_clean <- coords_clean %>% select(UTMXcorrec, UTMYcorrec) %>%
                  rename(X = UTMXcorrec, Y = UTMYcorrec) %>% as.matrix()


# #hex-level
# coords_tmp <- coords %>% separate_wider_delim(Site, delim = '-', names = c('hex','cam'), cols_remove = FALSE)
# coords_hex_clean <- coords_hex %>% 
#                         rename(hex = GRID_ID_PADDED) %>% 
#                         arrange(hex) %>%
#                         filter(hex %in% coords_tmp$hex) %>%
#                         select(X, Y, hex) 
#   
# #also need an index matching each site with the coords to use (hex centroid)  
# coords_index <- match(coords_tmp$hex, coords_hex_clean$hex)
# str(coords_index)  



## Create input objects --------------------------------------------------------

inputs <- list()

#loop thru species
for (ss in names(dh_arrays)){
  
  #get data for this species
  y_ss <- dh_arrays[[ss]]
  
  #keep sites that are not all NAs
  y_ss <- y_ss[sites_to_keep, , ]
  
  #format input object
  inputs[[ss]] <- list(y = y_ss,
                       occ.covs = occ_covar_clean,
                       det.covs = det_covar_clean,
                       coords = coords_clean)
}

names(inputs)
str(inputs[[1]])


## Create input objects (hex-level coords) -------------------------------------

# inputs_hex <- list()
# 
# #loop thru species
# for (ss in names(dh_arrays)){
#   
#   #get data for this species
#   y_ss <- dh_arrays[[ss]]
#   
#   #format input object
#   inputs_hex[[ss]] <- list(y = y_ss,
#                        occ.covs = occ_covar,
#                        det.covs = list('season' = occ_covar$season,
#                                        'season_wet' = occ_covar$season_wet,
#                                        'season_cool' = occ_covar$season_cool,
#                                        'season_dry' = occ_covar$season_dry,
#                                        'site' = occ_covar$site_num,
#                                        'fence' = occ_covar$fence,
#                                        'fence_south' = occ_covar$fence_south,
#                                        'fence_north' = occ_covar$fence_north,
#                                        'year_season_num' = occ_covar$year_season_num,
#                                        'year_num' = occ_covar$year_num),
#                        coords = coords_hex_clean,
#                        grid.index = coords_index)
# }
# 
# names(inputs_hex)
# str(inputs_hex[[1]])


## EXPLORE ---------------------------------------------------------------------

#detections
head(inputs[[1]]$y,10)
dim(inputs[[1]]$y) #200 sites, 18 year-seasons, 22 weeks

#occupancy covariates
names(inputs[[1]]$occ.covs)
str(inputs[[1]]$occ.covs)

#years/seasons need to be structured as time-varying
head(inputs[[1]]$occ.covs$year)
head(inputs[[1]]$occ.covs$year_num) #0-5
head(inputs[[1]]$occ.covs$year_season)
head(inputs[[1]]$occ.covs$year_season_num) #0-17
head(inputs[[1]]$occ.covs$season)

#all others can just be site-varying 
head(inputs[[1]]$occ.covs$site); length(inputs[[1]]$occ.covs$site)
head(inputs[[1]]$occ.covs$elev)
head(inputs[[1]]$occ.covs$elev_scaledNS)
head(inputs[[1]]$occ.covs$fence)
#etc

#detection covariates
names(inputs[[1]]$det.covs)
str(inputs[[1]]$det.covs)
head(inputs[[1]]$det.covs$season)
head(inputs[[1]]$det.covs$site, 50)
head(inputs[[1]]$det.covs$year)

#coordinates
head(inputs[[1]]$coords)
dim(inputs[[1]]$coords)
# head(inputs_hex[[1]]$coords)

#ensure no NAs in occurrence covar
sum(is.na(inputs[[1]]$occ.covs$years))
sum(is.na(inputs[[1]]$occ.covs$season))
sum(is.na(inputs[[1]]$occ.covs$site))
sum(is.na(inputs[[1]]$occ.covs$elev))
sum(is.na(inputs[[1]]$occ.covs$tree_vol))
sum(is.na(inputs[[1]]$occ.covs$dist_dambo))
sum(is.na(inputs[[1]]$occ.covs$dist_boundary_log_scaledNS))
sum(is.na(inputs[[1]]$occ.covs$fence))
#etc.

#save
saveRDS(inputs, 'H_sp_occ/sp_inputs.RDS') 
  #saved 090525 without 2024_dry
  #saved 092025 with log versions of 'dist' covar
  #saved 092125 with 'year' dummy variables in det.covs to use as factor


# saveRDS(inputs_hex, 'H_sp_occ/sp_inputs_hexCoords.RDS')


