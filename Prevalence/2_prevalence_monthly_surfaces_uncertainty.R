#script to produce monthly prevalence uncertainty surfaces
setwd(paste0(Sys.getenv("HOME"), "/Prevalence-Madagascar/Prevalence"))

library(Rcpp)
library(TMB)
library(raster)
library(INLA)


load("prepared_pr_data.RData")
log_inc_vals_simple <- matrix(NA, dim(mdg_pr_new), 2)
months_used <- unique(mdg_pr_new$month_no)
for(month_used in months_used){
  month_index <- which(mdg_pr_new$month_no == month_used)
  log_inc_surface <- raster(paste0("../Incidence_surfaces/log_inc_month_", month_used, ".tif"))
  log_inc_surface1 <- raster(paste0("../Incidence_surfaces/log_inc_month_", month_used + 1, ".tif"))
  log_inc_vec <- extract(log_inc_surface, pr_coords[month_index, , drop=FALSE])
  log_inc_vec1 <- extract(log_inc_surface1, pr_coords[month_index, , drop=FALSE])
  
  if(sum(is.na(log_inc_vec)) > 0){
    which_NA <- which(is.na(log_inc_vec))
    log_inc_surface <- focal(log_inc_surface, matrix(1/9, 3, 3), sum, na.rm=T)
    log_inc_surface1 <- focal(log_inc_surface1, matrix(1/9, 3, 3), sum, na.rm=T)
    log_inc_vec[which_NA] <- extract(log_inc_surface, pr_coords[month_index[which_NA], , drop=FALSE])
    log_inc_vec1[which_NA] <- extract(log_inc_surface1, pr_coords[month_index[which_NA], , drop=FALSE])
  }
  
  log_inc_vals_simple[month_index, 1] <- log_inc_vec
  log_inc_vals_simple[month_index, 2] <- log_inc_vec1
}

log_inc_vals_simple[log_inc_vals_simple > -3.5] <- -3.5
log_inc_vals_simple[log_inc_vals_simple < -8] <- -8



#load latest points and just keep DHS
load("gbd2020_pts.RData")
pts_dhs <- which(mdg_pr_new$id %in% pf_pr_dhs$id)

mdg_pr_new <- mdg_pr_new[pts_dhs, ]
log_inc_vals_simple <- log_inc_vals_simple[pts_dhs, ]
pr_coords <- pr_coords[pts_dhs, ]


pr_cov_mat_norm_dynamic <- pr_cov_mat_norm_all
pr_cov_mat_norm_static <- pr_cov_mat_norm[, 13:20]
static_names <- c("Access", "AI", "Elevation", "PET", "Slope", "Night lights", "dtw", "TWI")
colnames(pr_cov_mat_norm_static) <- static_names
pr_cov_mat_use <- cbind(pr_cov_mat_norm_dynamic,
                        pr_cov_mat_norm_static)
covs_use <- c("Access", "AI", "dtw", "Rain0")
pr_cov_mat_use <- pr_cov_mat_use[, covs_use]

mesh <- inla.mesh.2d(loc = pr_coords, cutoff = 0.25, max.edge = c(1, 3))
plot(mesh, asp=1)
points(pr_coords, col="blue")
mesh_coords <- mesh$loc[, 1:2]
n_mesh <- mesh$n

alpha <- 2
nu <- alpha - 1
spde <- (inla.spde2.matern(mesh=mesh, alpha=alpha)$param.inla)[c("M0","M1","M2")]
A_pr <- inla.spde.make.A(mesh=mesh, loc=as.matrix(pr_coords))


month_offset <- 0
month_offset_i <- month_offset + 1
log_inc_vals_use <- log_inc_vals_simple[, month_offset_i]
log_inc_vals_use1 <- log_inc_vals_simple[, month_offset_i + 1]
plot(log_inc_vals_use , mdg_pr_new$pf_pr)


n_inc_mesh <- 15
log_inc_mesh <- seq(min(c(log_inc_vals_use,
                          log_inc_vals_use1)), 
                    max(c(log_inc_vals_use,
                          log_inc_vals_use1)), 
                    length.out = n_inc_mesh)

make_A <- function(mesh_points, obs_points){
  A <- c()
  n_mesh <- length(mesh_points)
  n_obs <- length(obs_points)
  for(i in 1:n_obs){
    obs_point <- obs_points[i]
    out_vec <- rep(0, n_mesh)
    if(sum(mesh_points < obs_point) == 0){
      out_vec[1] <- 1}
    else if(sum(mesh_points > obs_point) == 0){
      out_vec[n_mesh] <- 1
    }else{
      i_lower <- max(which(mesh_points < obs_point))
      i_upper <- i_lower + 1
      
      out_vec[i_upper] <- (obs_point - mesh_points[i_lower]) / (mesh_points[i_upper] - mesh_points[i_lower])
      out_vec[i_lower] <- (mesh_points[i_upper] - obs_point) / (mesh_points[i_upper] - mesh_points[i_lower])
    }
    
    A <- rbind(A, unname(out_vec))
  }
  return(A)
}

A_inc <- make_A(log_inc_mesh, log_inc_vals_use)
A_inc1 <- make_A(log_inc_mesh, log_inc_vals_use1)
inc_distances <- matrix(NA, n_inc_mesh, n_inc_mesh)
for(i in 1:n_inc_mesh){
  for(j in 1:n_inc_mesh){
    inc_distances[i, j] <- (log_inc_mesh[i] - log_inc_mesh[j])^2
  }
}

model_folder <- NULL
model_name <- "full_model_with_switches"
model_path <- paste0(model_folder, model_name)
tryCatch(dyn.unload(dynlib(model_path)),
         error = function(e) print(e))
compile(paste0(model_path, ".cpp"))
dyn.load(dynlib(model_path))



model_silent <- FALSE

nugget <- 0.001
holdout <- 0*(!mdg_pr_new$year_start == 2016)


log_inc_lambda_mean <- 2
log_inc_lambda_sd <- 1
prior_rho_min <- 5
prior_rho_prob <- 0.001
prior_sigma_max <- 0.1
prior_sigma_prob <- 0.05

##full model w covs
m <- MakeADFun(
  data = list(Y_pr=mdg_pr_new$pf_pos,
              N_pr=mdg_pr_new$examined,
              X_pr=pr_cov_mat_use,
              inc_distances=inc_distances,
              spde=spde,
              A_pr=A_pr,
              beta_sd=0.5,
              beta_inc_sd=0.5,
              log_inc_mesh=log_inc_mesh,
              A_inc=A_inc,
              A_inc1=A_inc1,
              nugget=nugget,
              holdout=holdout,
              nu=nu,
              log_rho_mean=2,
              log_rho_sd=0.1,
              log_sigma_mean=-1.5,
              log_sigma_sd=0.1,
              log_inc_lambda_mean=log_inc_lambda_mean,
              log_inc_lambda_sd=log_inc_lambda_sd,
              use_field=1,
              use_covs=1,
              use_inc=1,
              use_inc1=1,
              fix_hypers=0,
              log_rho_fixed=2.617,
              log_sigma_fixed=0,
              log_inc_lambda_fixed=2.968
  ),
  parameters = list(beta_0=runif(1, -1, 1),
                    inc_mesh_vals=runif(n_inc_mesh),
                    S=rep(0, n_mesh),
                    beta=rep(0, dim(pr_cov_mat_use)[2]),
                    log_rho=5,
                    log_sigma=0,
                    log_inc_lambda=log_inc_lambda_mean,
                    beta_inc=c(2, 0)
  ),
  random=c("S", "inc_mesh_vals"),
  DLL = model_name,
  silent=model_silent
)



ptm <- proc.time()
fit <- nlminb(m$par, m$fn, m$gr, control=list(iter.max=300,eval.max=300))
ptm2 <- proc.time()
print("Time to fit")
print(ptm2 - ptm)


mesh_vals_full <- m$report()[[2]]
plot(log_inc_mesh, mesh_vals_full)


rep <- sdreport(m, getJointPrecision = TRUE)
param_df <- data.frame(mean=fit$par,
                       names=names(fit$par))
Sigma <- solve(rep$jointPrecision)
names_all <- names(m$env$last.par.best)
which_beta <- which(names_all == "beta")
which_beta_inc <- which(names_all == "beta_inc")
which_param_get <- c(which(names_all == "beta_0"),
                     which_beta,
                     which(names_all == "log_rho"),
                     which(names_all == "log_sigma"),
                     which(names_all == "log_inc_lambda"),
                     which_beta_inc)

Sigma_sub <- Sigma[which_param_get, which_param_get]
library(MASS)
param_samples <- matrix(mvrnorm(1000, fit$par, Sigma_sub), ncol=length(which_param_get))
param_df$LI <- apply(param_samples, 2, function(s) quantile(s, 0.025))
param_df$UI <- apply(param_samples, 2, function(s) quantile(s, 0.975))

write.csv(param_df, file="para_df.csv")



which_param_get <- c(which(names_all == "beta_0"),
                     which(names_all == "beta"),
                     which(names_all == "beta_inc"),
                     which(names_all == "S"),
                     which(names_all == "inc_mesh_vals"))
param_mean <- m$env$last.par.best[which_param_get]
param_names <- names(param_mean)

Sigma_sub2 <- Sigma[which_param_get, which_param_get]
n_samples <- 200

param_samples <- matrix(mvrnorm(n_samples, param_mean, Sigma_sub2), ncol=length(which_param_get))


#make pixel predictions
library(Rcpp)
sourceCpp("get_GP_value_1D.cpp")
sourceCpp("get_quantiles.cpp")

log_incidence_surface <- raster("../Incidence_surfaces/log_inc_month_1.tif")
which_not_na <- which(!is.na(values(log_incidence_surface)))


#do spatial GP
A_pixel <- inla.spde.make.A(mesh=mesh, loc=as.matrix(coordinates(log_incidence_surface)[which_not_na, ]))
#do covariate bit
mesh_vals_full_mat <- param_samples[, param_names == "inc_mesh_vals"]
beta_incs_mat <- param_samples[, param_names == "beta_inc"]
beta_mat <- t(param_samples[, param_names == "beta"])
beta_0_vec <- param_samples[, param_names == "beta_0"]
S_mat <- param_samples[, param_names == "S"]

library(raster)
#load static covariates and normalise
load("normalised_hf_covariates.RData")
print(names(static_norm_covs))
static_stack <- raster::stack("static_stack.tif")
# raster::plot(static_stack)
static_covs <- values(static_stack)[which_not_na, ]
for(i in 1:dim(static_covs)[2]){
  static_covs[, i] <- (static_covs[, i] - static_cov_means[i]) / static_cov_sds[i]
}


rain_stack <- raster::stack("rain_stack.tif")
lst_stack <- raster::stack("lst_stack.tif")
evi_stack <- raster::stack("evi_stack.tif")


incs <- c()
covs <- c()
all <- c()
all_trans <- c()

cov_contr_mat <- c()

surface_names <- c("L95", "LQ", "median", "UQ", "U95")


for(month_i in 1:48){
    print(paste0("month ", month_i))
    rain0 <- rain_stack[[month_i]]
    lst0 <- lst_stack[[month_i]]
    
    
    raster_list <- list(rain0, lst0)
    cov_index <- c(1, 5)
    raster_mean <- pr_cov_mat_norm_all_mean[cov_index]
    raster_sd <- pr_cov_mat_norm_all_sd[cov_index]
    
    dynamic_covs <- c()
    for(j in 1:length(cov_index)){
      dynamic_covs <- cbind(dynamic_covs,
                            (values(raster_list[[j]])[which_not_na] - raster_mean[j]) / raster_sd[j]
      )
    }
    colnames(dynamic_covs) <- c("Rain0", "LST0")
    colnames(static_covs) <- static_names
    cov_mat <- cbind(dynamic_covs, static_covs)
    cov_mat <- cov_mat[, covs_use]

    cov_contr_mat <- cov_mat %*% beta_mat
    
    log_incidence_surface <- raster(paste0("../Incidence_surfaces/log_inc_month_",
                                           month_i, ".tif"))
    
    if(month_i == 48){
      log_incidence_surface1 <- raster(paste0("../Incidence_surfaces/log_inc_month_",
                                              month_i, ".tif"))
    }else{
      log_incidence_surface1 <- raster(paste0("../Incidence_surfaces/log_inc_month_",
                                              month_i + 1, ".tif"))
    }
    
    print("Sample")
    ptm_month_start <- proc.time()
    prevalence_samples <- matrix(nrow=n_samples, ncol=length(which_not_na))

    log_inc_vals <- values(log_incidence_surface)[which_not_na]
    log_inc_vals1 <- values(log_incidence_surface1)[which_not_na]
    
    for(i_sample in 1:n_samples){
      cat(i_sample)
      cat(" ")
      
      inc_gp_vals0 <- get_GP_val_faster(log_inc_mesh, log_inc_vals, mesh_vals_full_mat[i_sample, ])
      inc_gp_vals1 <- get_GP_val_faster(log_inc_mesh, log_inc_vals1, mesh_vals_full_mat[i_sample, ])
      spatial_GP_mesh <- S_mat[i_sample, ]
      
      pixel_GP <- as.vector(A_pixel %*% spatial_GP_mesh)
      inc_cont <-  beta_incs_mat[i_sample, 2] * inc_gp_vals1 + beta_incs_mat[i_sample, 1] * inc_gp_vals0
      logit_prev_i <- inc_cont+ pixel_GP + beta_0_vec[i_sample] + cov_contr_mat[, i_sample]
      prev_i <- 1 / (1 + exp(-logit_prev_i))
      
      prevalence_samples[i_sample, ] <- prev_i
    }
    cat("\n")
    print(proc.time() - ptm_month_start)
    
    prev_quartiles <- get_col_quantiles(prevalence_samples, c(0.025, 0.25, 0.5, 0.75, 0.975))

    prev_surface <- log_incidence_surface
    for(quartile_i in 1:length(surface_names)){
      values(prev_surface) <- NA
      values(prev_surface)[which_not_na] <- prev_quartiles[quartile_i, ]
      png(paste0("../Prevalence_surfaces_CI/plot_", surface_names[quartile_i], "_month", month_i, ".png"))
      plot(prev_surface, zlim=c(0, 0.6), main=month_i)
      dev.off()
      writeRaster(prev_surface, paste0("../Prevalence_surfaces_CI/", surface_names[quartile_i], 
                                       "_month", month_i, ".tif"), overwrite=TRUE)
    }
}


for(year_i in 1:4){
  print(year_i)
  prev_sum_mat <- matrix(0, nrow=length(surface_names), ncol=length(which_not_na))
  
  for(month_i in (1:12) + (year_i-1)*12){
    for(i in 1:length(surface_names)){
      surface_name <- surface_names[i]
      r <- raster(paste0("../Prevalence_surfaces_CI/", surface_name, "_month", month_i, ".tif"))
      prev_sum_mat[i, ] <- prev_sum_mat[i, ] + values(r)[which_not_na]
    }
  }
  
  for(i in 1:length(surface_names)){
    surface_name <- surface_names[i]
    values(prev_surface) <- NA
    values(prev_surface)[which_not_na] <- prev_sum_mat[i, ] / 12
    writeRaster(prev_surface, filename = paste0("../Prevalence_surfaces_CI/" ,surface_name, "_", year_i+2012, ".tif"))
  }
    
}
  


