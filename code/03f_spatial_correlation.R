# 5. Evaluate spatial parameters

library(data.table)

## read csvs
eval_files <- list.files('H_sp_occ/results/model17_MSversion/evaluation/', pattern = 'summaries.csv', full.names = TRUE)
eval_list <- lapply(eval_files, fread)

eval_df <- rbindlist(eval_list)
head(eval_df)
unique(eval_df$parameter)

## summarize sigma-sq and rho
spatial_summary <- eval_df %>% group_by(species) %>%
                      filter(parameter %in% c('sigma.sq','phi')) %>%
                      select(parameter, mean, sd, species) %>%
                      pivot_wider(names_from = parameter,
                                  values_from = c(mean, sd)) %>%
                      mutate(spatial_range = (3/mean_phi)/1000,
                             spatial_range_min = (3/(mean_phi+sd_phi))/1000,
                             spatial_range_max = (3/(mean_phi-sd_phi))/1000)

## REPORT IN TABLE 3
spatial_summary %>% select(mean_sigma.sq, sd_sigma.sq, mean_phi, sd_phi) %>%
  mutate(across(where(is.numeric), ~ round(.x, 4)))

spatial_summary %>% select(spatial_range, spatial_range_min, spatial_range_max) %>%
  mutate(across(where(is.numeric), ~ round(.x, 1))) %>%
  filter(species %in% c('buffalo','elephant','waterbuck','kudu'))


## CREATE TABLE OF PARAM EST, R-HAT, AND ESS
eval_df_summary <- eval_df %>% select(species, parameter_type, parameter, mean, `q2.5`, q50, `q97.5`, Rhat, ESS) %>%
                      filter(parameter_type != 'detection_RE') %>%
                      mutate(across(c(mean, `q2.5`, q50, `q97.5`, Rhat), ~ round(.x, 3))) %>%
                      mutate(across(c(ESS), ~ round(.x, 0)))
# arrange(species, parameter_type, parameter)

write.csv(eval_df_summary, 'H_sp_occ/results/model17_MSversion/evaluation/eval_df_summary.csv', row.names = FALSE)
