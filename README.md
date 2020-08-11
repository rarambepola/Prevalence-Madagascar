# Prevalence-Madagascar
This repository contains code used for making monthly maps of prevalence of Plasmodium falciparum infection in individuals between 6 and 59 months of age in Madagascar between 2013 and 2016.

The routine case data (inlcuding health facility locatiosn) in this repository are **sample data**. Real data may be obtained by contacting the National Malaria Control Programme of Madagascar. Prevalence data is accurate and full Malaria Indicator Survey datasets are freely available from the DHS online data repository: https://dhsprogram.com/.

The outputs of each step are included so that any step can be run independently. Please note that the catchment model requires a large amount of information to be read into memory (around 30 GB) and therefore is unlikely to run on machines with less than 64GB of RAM. 

The causal inference scripts require two small R packages: [RCITcpp](https://github.com/rarambepola/RCITcpp), for efficient independence testing, and [causalInference](https://github.com/rarambepola/causalInference), for a two step, order-independent implementation of the PC algorithm.
