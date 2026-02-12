## Format detection data for spatiotemporal occupancy models

library(data.table)
library(tidyverse)
library(fuzzyjoin)
library(purrr)
library(reshape2)


## Read in long-format detection data ------------------------------------------

#This file should have one row per site-date-species and contain ALL combinations.
#For example, it should include all dates between the earliest and latest date of deployment, 
#even if there were no photos of a given species on that date (n_photos = 0) 
#or if the camera wasn't surveyed at all (surv = 0). 

#Columns are 'site', 'date', 'class' (or 'species'), 
#'n_photos', and 'surv' (0/1 if surveyed on this date)

dh_long <- fread('data/raw/detection_histories_all_long.csv')
head(dh_long)
sort(unique(dh_long$site)); length(unique(dh_long$site))


## Establish sampling periods (YEAR-SEASONS) -----------------------------------

#get start/end dates
(beg <- min(dh_long$date))
(end <- max(dh_long$date))
(n_days <- interval(beg, end)/days(1))

#set primary period start/end dates (seasons: wet=Dec-April, cool=May-Jul, dry=Aug-Nov)
(n_seasons <- 3) #3 per year
(n_yrs <- year(end) - year(beg))

dates_primary <- data.frame('season' = rep(c('wet','cool','dry'), length.out = (n_yrs+1) * n_seasons),
                            'season_year' = rep((year(beg)):year(end), each = n_seasons), 
                            #wet seasons will be called 'wet_2019' even tho it starts Dec 2018
                            'year' = c((year(beg)-1), rep((year(beg)):(year(end)-1), each = n_seasons), rep(year(end), n_seasons-1)), 
                            #hacky... but the dates need to start a year prior 
                            'start' = rep(c('12-01','05-01','08-01'), length.out = (n_yrs+1) * n_seasons))

dates_primary$primary <- paste(dates_primary$season_year, dates_primary$season, sep = '_')
dates_primary$start <- as.POSIXct(strptime(paste(dates_primary$year, dates_primary$start, sep = '-'),
                                           format = '%Y-%m-%d'), tz = 'Africa/Blantyre')
dates_primary <- dates_primary %>% arrange(start) %>% 
                                   select(-season_year, -year) %>%
                                   mutate(end = lead(start) - days(1)) 
                                   #subtract 1 to not include the next start date

#add last end date
dates_primary[dates_primary$primary == '2024_dry',]$end <- as.POSIXct('2024-11-30', tz = 'Africa/Blantyre') 

#list all dates
dates <- format(as.Date(beg) + days(0:(n_days)), format = '%Y-%m-%d') 
dates_df <- data.frame('date' = as.POSIXct(strptime(dates, '%Y-%m-%d'), tz = 'Africa/Blantyre'))

#and match up with primary periods
dates_df <- fuzzy_left_join(dates_df, dates_primary,
                            by = c('date' = 'start', 'date' = 'end'),
                            #match dates between start/end of each season
                            match_fun = list(`>=`, `<=`)) 
                            
dates_df <- dates_df %>% select(date, primary)

#set secondary periods (weeks)
secondary_length <- 7
dates_df <- dates_df %>% 
            group_by(primary) %>% 
            mutate(secondary = ceiling(row_number() / secondary_length)) %>% 
            ungroup()

#create a combined field ('year_season-week')
dates_df$occasion <- paste(dates_df$primary, 
                           formatC(as.numeric(dates_df$secondary), width = 2, flag = "0"), sep = '-')

#how many weeks in each season?
dates_df %>% group_by(primary) %>% summarise(max = max(secondary)) 

#wet: 22 weeks
#cool: 14 weeks
#dry: 18 weeks
#we'll remove or pad any incomplete seasons later

#save key to dates/time periods
write.csv(dates_df, 'data/cleaned/season_date_key.csv')


## Match detections w/ survey periods ------------------------------------------

head(dh_long)
head(dates_df)

#format dates
dh_long$date <- as.POSIXct(strptime(dh_long$date, format = '%Y-%m-%d'), tz = 'Africa/Blantyre')

#join
dh_long <- dh_long %>% left_join(dates_df, by = 'date')

#any NAs?
nrow(dh_long[is.na(dh_long$primary),])
nrow(dh_long[is.na(dh_long$secondary),])


## Aggregate by week -----------------------------------------------------------

#how many weeks will we need for each primary period?
(n_wks <- max(dh_long$secondary))
(n_primary <- length(unique(dh_long$primary)))

#loop thru for each species
dh_weekly <- list()
for (ss in unique(dh_long$class)){ 
  
  #sum n_photos for each site-year-season-week, preserving NAs
  dh_weekly_ss <- dh_long %>% filter(class == ss) %>% 
                          group_by(site, occasion) %>%
                          summarize(n_photos = if (all(is.na(n_photos))) NA_real_ else sum(n_photos, na.rm = TRUE), 
                                    surv = max(surv),
                                    start = min(date),
                                    .groups = 'drop',)
                          
  #now pad out so all seasons have the same number of weeks (fill with NAs)
  dh_weekly_ss_pad <- dh_weekly_ss %>% 
                          separate_wider_delim(occasion, delim = '-', names = c('year_season','week')) %>%
                          complete(site, year_season, week) %>%
                          arrange(site, year_season, week) 
  
  #store
  dh_weekly[[ss]] <- dh_weekly_ss_pad
  
}

names(dh_weekly)
head(dh_weekly$aardvark)
tail(dh_weekly$elephant)

#qc dimensions 
n_sites * n_wks * n_primary
nrow(dh_weekly$elephant) #good, should match

n_sites * n_wks
table(dh_weekly$elephant$year_season, useNA = 'a') #good, each season should match

#qc -- any photos outside surv period?
dh_weekly[[1]] %>% filter(surv == 0, !is.na(n_photos)) %>% nrow()
dh_weekly[[2]] %>% filter(surv == 0, !is.na(n_photos)) %>% nrow()
#etc.

#save
saveRDS(dh_weekly, 'data/cleaned/detections.RDS')


## Summaries to report ---------------------------------------------------------

head(dh_long)
head(dh_weekly)

#total photos across all species (but doesn't include empties and doesn't account for photos w/ multiple species)
sum(dh_long$n_photos, na.rm = TRUE) 

#read in 'fence' covariate (is each site north or south of the fence)
env_covar <- fread('data/cleaned/cleaned_covariates.csv')

#summarize by species (n sites, n photos, etc.)
dh_long_fence <- dh_long %>% left_join(env_covar %>% select(Site, fence), by = c('site' = 'Site'))
sp_summary <- dh_long_fence %>%
                  group_by(class, fence) %>%
                  filter(surv == 1) %>%
                  filter(!is.na(n_photos) & n_photos > 0) %>%
                  summarise(n_sites = length(unique(site)),
                            n_photos = sum(n_photos, na.rm = TRUE),
                            # n_sites_in_fence = length(unique(site[fence == 1])),
                            # n_photos_in_fence = sum(n_photos[fence == 1], na.rm = TRUE),
                            # n_sites_out_fence = length(unique(site[fence == 0])),
                            # n_photos_out_fence = sum(n_photos[fence == 0], na.rm = TRUE)
                            )
sp_summary

## For TABLE 1
sp_summary %>% pivot_wider(names_from = fence,
                           values_from = c(n_sites, n_photos),
                           names_prefix = 'fence_') %>%
  rename(n_sites_south = n_sites_fence_1,
         n_photos_south = n_photos_fence_1,
         n_sites_north = n_sites_fence_0,
         n_photos_north = n_photos_fence_0) %>%
  filter(class %in% c('buffalo','eland','elephant','impala','kudu',
                      'sable','warthog','waterbuck','zebra'))
  
length(unique(dh_long$site)) #200 sites
dh_long_fence %>% filter(surv == 1) %>%
  group_by(fence) %>%
  summarise(n_sites = length(unique(site))) #116 north, 84 south

#days surveyed at each camera
surv_dates <- dh_long %>% filter(surv == 1, class == 'elephant') %>%
  group_by(site) %>%
  summarise(n_days_surveyed = sum(surv),
            n_days_surveyed2 = length(unique(date)), #just to make sure they match
            n_weeks_surveyed = length(unique(occasion)))

surv_dates %>%
  arrange(n_days_surveyed) #min 1 day surveyed (site H06-B had only blank photos -- won't break the model but should have filtered it)

summary(surv_dates$n_days_surveyed) #mean 733 total days surveyed (range 1-2061)

surv_dates %>%
  arrange(n_weeks_surveyed) #min 1 week surveyed (H06-B), otherwise 5 weeks is the min

#per year?
surv_dates_year <- dh_long %>% filter(surv == 1, class == 'elephant') %>%
  mutate(year = year(date)) %>%
  group_by(site, year) %>%
  summarise(n_days_surveyed = sum(surv),
            n_days_surveyed2 = length(unique(date)), #just to make sure they match
            n_weeks_surveyed = length(unique(occasion)))

surv_dates_year %>%
  arrange(n_days_surveyed)

#summarize
surv_dates %>%
  ungroup() %>%
  summarise(mean_days = mean(n_days_surveyed),
            sd_days = sd(n_days_surveyed),
            min_days = min(n_days_surveyed),
            max_days = max(n_days_surveyed),
            mean_weeks = mean(n_weeks_surveyed),
            sd_weeks = sd(n_weeks_surveyed),
            min_weeks = min(n_weeks_surveyed),
            max_weeks = max(n_weeks_surveyed))

surv_dates_year %>%
  ungroup() %>%
  summarise(mean_days = mean(n_days_surveyed),
            sd_days = sd(n_days_surveyed),
            min_days = min(n_days_surveyed),
            max_days = max(n_days_surveyed),
            mean_weeks = mean(n_weeks_surveyed),
            sd_weeks = sd(n_weeks_surveyed),
            min_weeks = min(n_weeks_surveyed),
            max_weeks = max(n_weeks_surveyed))

summary(surv_dates$n_days_surveyed)


## Format into multi-season stacked design for spOcc ---------------------------

#read back in if necessary
# dh_weekly <- readRDS('data/cleaned/detections.RDS')

str(dh_weekly)
dim(dh_weekly$aardvark)
head(dh_weekly$aardvark)

#put the season_years in order
dh_weekly <- lapply(dh_weekly, function(x){
      x$year_season <- factor(x$year_season, 
                              levels = c('2018_dry',
                                         '2019_wet','2019_cool','2019_dry',
                                         '2020_wet','2020_cool','2020_dry',
                                         '2021_wet','2021_cool','2021_dry',
                                         '2022_wet','2022_cool','2022_dry',
                                         '2023_wet','2023_cool','2023_dry',
                                         '2024_wet','2024_cool','2024_dry'
                                         ))
      x
})
levels(dh_weekly$aardvark$year_season)
head(dh_weekly$aardvark)

#modify columns (use 'map' to apply to each dataframe in the list)
dh_weekly <- dh_weekly %>% map(~ .x %>%
                                 separate_wider_delim(year_season, delim = '_', names = c('year','season'), cols_remove = FALSE) %>%
                                 mutate(site_season = paste(site, season, sep = '_'))
                               )
head(dh_weekly$aardvark)

#for each species, we want an array with dimensions site_seasons*years*weeks, or sites*year_seasons*weeks
#values will be NA, 0, or 1

sapply(dh_weekly$aardvark, class) #check column classes

#transform to arrays using 'acast' and convert n_photos to 0/1 (preserving NAs)

#site_seasons * years
dh_arrays_season <- dh_weekly %>% map(~ .x %>%
                                       select(site_season, year, week, n_photos) %>%
                                       mutate(n_photos = ifelse(is.na(n_photos), NA, as.numeric(n_photos > 0))) %>%
                                       acast(site_season ~ year ~ week, value.var = 'n_photos'))

#sites * year_seasons
dh_arrays_year <- dh_weekly %>% map(~ .x %>%
                                     select(site, year_season, week, n_photos) %>%
                                     mutate(n_photos = ifelse(is.na(n_photos), NA, as.numeric(n_photos > 0))) %>%
                                     acast(site ~ year_season ~ week, value.var = 'n_photos'))


names(dh_arrays_season)
head(dh_arrays_season$elephant)
  dim(dh_arrays_season$elephant)
  n_sites * n_seasons #check sites*seasons
head(dh_arrays_year$elephant)
  dim(dh_arrays_year$elephant) #we have 1 extra season here (haven't deleted the first incomplete season yet, dry_2018)
  n_seasons * 6
  
#save
saveRDS(dh_arrays_season, 'data/cleaned/detection_arrays_seasonal.RDS')
saveRDS(dh_arrays_year, 'data/cleaned/detection_arrays_yearly.RDS')

## Move on now to 01b_prep_covariates.R

