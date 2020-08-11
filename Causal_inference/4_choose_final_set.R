#Perform cross-validation on potential causal feature sets to choose final feature set
setwd(paste0(Sys.getenv("HOME"), "/Prevalence-Madagascar/Causal_inference"))
library(TMB)
library(raster)
library(INLA)

##prepare data
load("prepared_pr_data_parameter_tuning.RData")
log_inc_vals_simple[log_inc_vals_simple > -4] <- -4
load("causal_sets_all.RData")
pr_cov_mat_norm_static <- pr_cov_mat_norm[, 13:20]
#note that causal feature sets are matched by NAME later in this script
#so it doesn't matter than covariates that we don't use (e.g. 3 month timelags) are included
pr_cov_mat_norm_dynamic <- pr_cov_mat_norm_all 
static_names <- c("Access", "AI", "Elevation", "PET", "Slope", "Night lights", "dtw", "TWI")
colnames(pr_cov_mat_norm_static) <- static_names
pr_cov_mat_all <- cbind(pr_cov_mat_norm_static, pr_cov_mat_norm_dynamic)
month_offset <- 1
month_offset_i <- month_offset + 1
log_inc_vals_use <- log_inc_vals_simple[, month_offset_i]

#set up mesh for GMRF approximation of matern random field (see Lindgren 2011)
mesh <- inla.mesh.2d(loc = pr_coords, cutoff = 0.25, max.edge = c(1, 3))
# plot(mesh, asp=1)
# points(pr_coords, col="blue")
mesh_coords <- mesh$loc[, 1:2]
n_mesh <- mesh$n

alpha <- 2
nu <- alpha - 1
spde <- (inla.spde2.matern(mesh=mesh, alpha=alpha)$param.inla)[c("M0","M1","M2")]
A_pr <- inla.spde.make.A(mesh=mesh, loc=as.matrix(pr_coords)) #projection matrix

#prepare 1D mesh for indcidence-prevalence relationship
n_inc_mesh <- 100
log_inc_mesh <- seq(min(log_inc_vals_use), max(log_inc_vals_use), length.out = n_inc_mesh)
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
inc_distances <- matrix(NA, n_inc_mesh, n_inc_mesh)
for(i in 1:n_inc_mesh){
  for(j in 1:n_inc_mesh){
    inc_distances[i, j] <- (log_inc_mesh[i] - log_inc_mesh[j])^2
  }
}


#compile model code
model_folder <- "../Prevalence/"
model_name <- "full_model_no_pc"
model_path <- paste0(model_folder, model_name)
tryCatch(dyn.unload(dynlib(model_path)),
         error = function(e) print(e))
compile(paste0(model_path, ".cpp"))
dyn.load(dynlib(model_path))



#set up model parameters
model_silent <- TRUE
nugget <- 0.001
log_inc_lambda_mean <- 3
log_inc_lambda_sd <- 0.1

n_test_sets <- length(causal_sets_all) #how many potential feature sets
test_set_covs <- list()


##run model with different feature sets and holdout sets
#compute CV for each
for(test_set_i in 1:n_test_sets){
  print(paste("Test set: ", test_set_i, " out of ", n_test_sets))
  pr_cov_mat_use <- pr_cov_mat_all[, causal_sets_all[[test_set_i]]]
  
  
  count_i <- 0
  mesh_vals_full <- list()
  cov_set_cors <- c()
  for(holdout_set_i in 1:length(holdout_set_indicators)){
    print(holdout_set_i)
    holdout <- holdout_set_indicators[[holdout_set_i]]
    count_i <- count_i + 1
    
    ##full model w covs
    m <- MakeADFun(
      data = list(Y_pr=mdg_pr_new$pf_pos,
                  N_pr=mdg_pr_new$examined,
                  X_pr=pr_cov_mat_use,
                  inc_distances=inc_distances,
                  spde=spde,
                  A_pr=A_pr,
                  beta_sd=0.1,
                  log_inc_vals_simple=log_inc_vals_use,
                  log_inc_mesh=log_inc_mesh,
                  A_inc=A_inc,
                  nugget=nugget,
                  holdout=holdout,
                  nu=nu,
                  log_rho_mean=5,
                  log_rho_sd=0.1,
                  log_sigma_mean=0,
                  log_sigma_sd=0.1,
                  log_inc_lambda_mean=log_inc_lambda_mean,
                  log_inc_lambda_sd=log_inc_lambda_sd
      ),
      parameters = list(beta_0=runif(1, -1, 1),
                        inc_mesh_vals=runif(n_inc_mesh),
                        S=rep(0, n_mesh),
                        beta=rep(0, dim(pr_cov_mat_use)[2]),
                        log_rho=5,
                        log_sigma=0,
                        log_inc_lambda=log_inc_lambda_mean
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
    
    obs_vals <- mdg_pr_new$pf_pr
    pred_vals <- 1/(1 + exp( - m$report()[[1]]))
    which_holdout <- which(holdout == 1)
    which_fit <- which(holdout == 0)
    
    obs_vec <- obs_vals[which_holdout]
    pred_vec <- pred_vals[which_holdout]
    
    mesh_vals_full[[count_i]] <- m$report()[[2]]
    full_model_cor <- cor(obs_vec, pred_vec, use="complete")
    
    print(full_model_cor)
    cov_set_cors[holdout_set_i] <- full_model_cor
  }

 test_set_covs[[test_set_i]] <- cov_set_cors 
}


##plot results
cov_df <- data.frame(cov=unlist(test_set_covs),
                     holdout=rep(1:6, n_test_sets),
                     set=as.factor(rep(1:n_test_sets, each=6)))
library(ggplot2)
print(ggplot(cov_df, aes(x=holdout, y=cov, col=set)) + geom_point() + geom_line())


plot(sapply(test_set_covs, mean))
print(which.max(sapply(test_set_covs, mean)))
plot(sapply(test_set_covs, median))
print(which.max(sapply(test_set_covs, median)))




