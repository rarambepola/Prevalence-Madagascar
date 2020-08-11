# Prevalence-Madagascar
This repository contains code used for making monthly maps of prevalence of Plasmodium falciparum infection in individuals between 6 and 59 months of age in Madagascar between 2013 and 2016.

The routine case data (inlcuding health facility locatiosn) in this repository are **sample data**. Real data may be obtained by contacting the National Malaria Control Programme of Madagascar. Prevalence data is accurate and full Malaria Indicator Survey datasets are freely available from the DHS online data repository: https://dhsprogram.com/.

The outputs of each step are included so that any step can be run independently. Covariate rasters are not stored in this repository due to their size but can be found [here](https://drive.google.com/drive/folders/15KbwxvDxWnPD6yQcBY2QLA5JWI9Sl0BF?usp=sharing). Scripts that require these external data files are:
* `Incidence/1_incidence_surfaces.R` requires `raster_covs_static.tif`
* `Prevalence/1_prevalence_monthly_surfaces.R` requires the folder `Incidence_surfaces` and files `static_stack.tif`, `rain_stack.tif`, `lst_stack.tif`, `evi_stack.tif`
* `Prevalence/2_prevalence_monthly_surfaces_uncertainty.R` requires the folder `Incidence_surfaces` and files `static_stack.tif`, `rain_stack.tif`, `lst_stack.tif`, `evi_stack.tif`

Please note that the catchment model requires a large amount of information to be read into memory (around 30 GB) and therefore is unlikely to run on machines with less than 64GB of RAM. 

The causal inference scripts require the R packages [RCITcpp](https://github.com/rarambepola/RCITcpp) for efficient independence testing.
