#script to produce monthly prevalence surfaces
setwd(paste0(Sys.getenv("HOME"), "/Prevalence-Madagascar/Prevalence"))

library(TMB)
library(raster)
library(INLA)

load("../Causal_inference/causal_sets_all.RData")
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

#use correct feature set
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


log_inc_lambda_mean <- 3
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


#make pixel predictions

library(Rcpp)
sourceCpp("get_GP_value_1D.cpp")

log_incidence_surface <- raster("../Incidence_surfaces/log_inc_month_1.tif")
which_not_na <- which(!is.na(values(log_incidence_surface)))


#do spatial GP
A_pixel <- inla.spde.make.A(mesh=mesh, loc=as.matrix(coordinates(log_incidence_surface)[which_not_na, ]))
#do covariate bit
mesh_vals_full <- m$report()[[2]]
beta_incs <- fit$par[names(fit$par) == "beta_inc"]

library(raster)
#load static covariates and normalise
load("normalised_hf_covariates.RData")
print(names(static_norm_covs))
static_stack <- raster::stack("../../Map_madagascar/Model_2020/JustPrevalence/static_stack.tif")
static_covs <- values(static_stack)[which_not_na, ]
for(i in 1:dim(static_covs)[2]){
  static_covs[, i] <- (static_covs[, i] - static_cov_means[i]) / static_cov_sds[i]
}


rain_stack <- raster::stack("../../Map_madagascar/Model_2020/rain_stack.tif")
lst_stack <- raster::stack("../../Map_madagascar/Model_2020/lst_stack.tif")
evi_stack <- raster::stack("../../Map_madagascar/Model_2020/evi_stack.tif")


incs <- c()
covs <- c()
all <- c()
all_trans <- c()

cov_contr_mat <- c()


for(year_i in 1:4){
  year_prev <- rep(0, length(which_not_na))
  
  
  for(month_i in (1:12) + (year_i-1)*12){
    prev_sum <- 0
    rain0 <- rain_stack[[month_i]]
    rain1 <- rain_stack[[month_i + 1]]
    rain2 <- rain_stack[[month_i + 2]]
    rain3 <- rain_stack[[month_i + 3]]
    evi3 <- evi_stack[[month_i + 3]]
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
    beta <- fit$par[names(fit$par) == "beta"]
    
    cov_contr <- as.vector(cov_mat %*% beta)
    
    log_incidence_surface <- raster(paste0("../Incidence_surfaces/log_inc_month_",
                                           month_i, ".tif"))
    
    if(month_i == 48){
      log_incidence_surface1 <- raster(paste0("../Incidence_surfaces/log_inc_month_",
                                              month_i, ".tif"))
    }else{
      log_incidence_surface1 <- raster(paste0("../Incidence_surfaces/log_inc_month_",
                                              month_i + 1, ".tif"))
    }
    
    # print(i)
    inc_gp_vals0 <- get_GP_val(log_inc_mesh, values(log_incidence_surface)[which_not_na], mesh_vals_full)
    inc_gp_vals1 <- get_GP_val(log_inc_mesh, values(log_incidence_surface1)[which_not_na], mesh_vals_full)
    
    spatial_GP_mesh <- m$env$last.par.best[names(m$env$last.par.best) == "S"]
    
    pixel_GP <- as.vector(A_pixel %*% spatial_GP_mesh)
    inc_cont <-  beta_incs[2] * inc_gp_vals1 + beta_incs[1] * inc_gp_vals0
    logit_prev_no_cov <- inc_cont+ pixel_GP + fit$par["beta_0"] + cov_contr
    prev_no_cov <- 1 / (1 + exp(-logit_prev_no_cov))
    # 
    print(month_i)
    print(mean(pixel_GP, na.rm=T))
    print(mean(inc_cont, na.rm=T))
    print(mean(cov_contr, na.rm=T))
    print(fit$par["beta_0"])
    
    incs[month_i] <- mean(inc_cont, na.rm=T)
    covs[month_i] <- mean(cov_contr, na.rm=T)
    all[month_i] <- mean(logit_prev_no_cov, na.rm=T)
    all_trans[month_i] <- mean(prev_no_cov, na.rm=T)
    
    #do each covariate separately
    cov_each <- c()
    for(i in 1:dim(pr_cov_mat_use)[2]){
      cov_each[i] <- mean(beta[i] * cov_mat[, i], na.rm=T)
    }
    names(cov_each) <- colnames(pr_cov_mat_use)
    cov_contr_mat <- rbind(cov_contr_mat, cov_each)
    

    prev_surface <- log_incidence_surface
    values(prev_surface) <- NA
    values(prev_surface)[which_not_na] <- prev_no_cov
    png(paste0("../Prevalence_surfaces/plot_prevalence_month", month_i, ".png"))
    plot(prev_surface, zlim=c(0, 0.6), main=month_i)
    dev.off()
    
    # writeRaster(prev_surface, paste0("../Prevalence_surfaces/prevalence_month", month_i, ".tif"), overwrite=TRUE)
    year_prev <- year_prev + prev_no_cov
  }
  
  values(prev_surface)[which_not_na] <- year_prev / 12
  writeRaster(prev_surface, paste0("../Prevalence_surfaces/prevalence_", 2012+year_i, ".tif"), overwrite=TRUE)
  
}

