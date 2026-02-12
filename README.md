# Nkhotakota_spatial_occupancy

Code for spatial occupancy models with camera trap data from Nkhotakota Wildlife Reserve

Anonymized data are included in the `data` folder

Note that some of these scripts are meant to be run on a high performance computing cluster (HPC). If running multiple models (e.g., multiple species), the scripts can be called from the command line using arguments. For example:

`Rscript 02a_run_model.R zebra 20000 outputs/model_1`

The code is separated into 3 steps:

### 1. PREP:

prepares data for spOccupancy package

-   `01a_prep_detections.R`: reads in long-format detection data, specifies primary and secondary sampling periods, and performs aggregation to create detection histories for spOccupancy stacked format

-   `01b_prep_covariates.R`: formats covariate data for multi-season spOccupancy

-   `01c_prep_inputs.R`: takes detection histories, covariate data, and spatial coordinates to create input objects (RDS files) for spOccupancy

### 2. RUN

runs and evaluates the spatial occupancy model(s)

-   `02a_run_model.R`: runs the spatial occupancy model reported in the paper (single species spatiotemporal model: stPGOcc) [takes a long time! run on cluster]

-   `02b_run_ppc.R`: runs posterior predictive checks and calcualates Bayesian p-values

-   `02c_run_diagnostics.R`: reads in model output and generates trace plots and summary files for model diagnostics (coefficient estimates, Rhat, ESS) (used to be eval) (run the PPC within diagnostics now?)

### 3. REPORT

compute derived parameters, create plots, conduct post-hoc trend analysis

-   `03_prediction.R`: takes model outputs, covariate values, and spatial coordinates to generate predictions at new sites [takes a long time! run on cluster]

-   `03_prediction_maps.R`

-   `03_derive_parameters.R`

-   `03_marginal_plots.R`

-   `03_trend_analysis.R`
