#set working directory
setwd(paste0(Sys.getenv("HOME"), "/Prevalence-Madagascar/Incidence"))

library(TMB)
library(raster)
library(malariaAtlas)
library(fasterize)
library(INLA)

#load data
load("sample_incidence_data.RData")

#load MDG raster to remove NA locations
raster_stack <- stack("../../Map_madagascar/Model_2020/raster_covs_static.tif")
r <- raster_stack[[1]]
raster_not_na <- !is.na(values(r))


model_folder <- NULL
model_name <- "fixed_hyperparameters"
model_path <- paste0(model_folder, model_name)
tryCatch(dyn.unload(dynlib(model_path)),
         error = function(e) print(e))
compile(paste0(model_path, ".cpp"))
dyn.load(dynlib(model_path))



mesh <- inla.mesh.2d(loc = hf_coords, cutoff = 0.1, max.edge = c(0.4, 1))

#model parameters from CV
log_rho <- 0.1
log_sigma <- -2

#make monthly incidence surfaces
for(month_use in 1:48){
  year <- 2012 + ceiling(month_use/12)

  
  ptm_start <- proc.time()
  print(month_use)
  N_total <- length(FIDs_complete)
  
  case_vector <- case_matrix[, month_use]
  load(paste0("../Catchment model/catchment_populations_", year, ".RData"))

  pops_vector <- population
  cases_na <- which(is.na(case_vector))
  cases_na_binary <- 1*!is.na(case_vector)
  rate_vector <- case_vector / pops_vector
  
  pops_vector_use <- pops_vector[!is.na(case_vector)]
  hf_coords_complete_use <- hf_coords_complete[!is.na(case_vector), ]
  rate_vector_use <- rate_vector[!is.na(case_vector)]
  case_vector_use <- case_vector[!is.na(case_vector)]
  n_obs <- length(case_vector)
  

  #create mesh
  alpha <- 2
  nu <- alpha - 1
  spde <- (inla.spde2.matern(mesh=mesh, alpha=alpha)$param.inla)[c("M0","M1","M2")]
  A <- inla.spde.make.A(mesh=mesh, loc=as.matrix(hf_coords_complete_use))
  n_s <- nrow(spde$M0)
  mesh_coords <- mesh$loc[, 1:2]
  
  holdout_set <- rep(0, length(case_vector_use))
  
  model_silent <- TRUE
    m <- MakeADFun(
      data = list(Y=case_vector_use,
                  pops=pops_vector_use,
                  A=A,
                  spde=spde,
                  prior_rho_min = 10,
                  prior_rho_prob = 0.001,
                  prior_sigma_max = 0.1,
                  prior_sigma_prob = 0.001,
                  holdout_set=holdout_set,
                  log_rho=log_rho,
                  log_sigma=log_sigma,
                  nu=nu
                  
      ),
      parameters = list(beta_0=runif(1, -1, 1),
                        S=rep(0, n_s)
      ),
      random = "S",
      DLL = model_name,
      silent=model_silent
    )
    ptm <- proc.time()
    fit <- nlminb(m$par, m$fn, m$gr, control=list(iter.max=300,eval.max=300))
    ptm2 <- proc.time()
    print(ptm2 - ptm)
    
    mesh_vals <- m$env$last.par.best[names(m$env$last.par.best) == "S"]
    beta_0 <- fit$par[names(fit$par) == "beta_0"]
    A_raster <- inla.spde.make.A(mesh=mesh, loc=as.matrix(coordinates(r)[raster_not_na, ]))
    field_vals <- A_raster %*% mesh_vals
    values(r)[raster_not_na] <- field_vals + beta_0

    writeRaster(r, filename=paste0("../Incidence_surfaces/log_inc_month_",
                                   month_use, ".tif"),
                overwrite=TRUE)

  
    ptm_end <- proc.time()
    print("total time")
    print(ptm_end - ptm_start)
  
}
