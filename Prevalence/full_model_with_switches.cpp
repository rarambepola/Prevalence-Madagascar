#include <TMB.hpp>
#include <math.h>
#include <stdio.h>



template<class Type>
Type objective_function<Type>::operator() ()
{
  
  using namespace R_inla;
  using namespace density;
  using namespace Eigen;
  

  //prevalence stuff
  DATA_VECTOR(Y_pr);
  DATA_VECTOR(N_pr);
  DATA_MATRIX(X_pr);
  DATA_MATRIX(inc_distances);
  DATA_STRUCT(spde,spde_t);
  DATA_SPARSE_MATRIX(A_pr);
  DATA_SCALAR(beta_sd);
  DATA_SCALAR(beta_inc_sd);
  
  //(log) incidence data
  DATA_VECTOR(log_inc_mesh);
  DATA_MATRIX(A_inc);
  DATA_MATRIX(A_inc1);
  DATA_SCALAR(nugget);
  DATA_IVECTOR(holdout);
  DATA_SCALAR(nu);
  
  //prior settings
  DATA_SCALAR(log_rho_mean);
  DATA_SCALAR(log_rho_sd);
  DATA_SCALAR(log_sigma_mean);
  DATA_SCALAR(log_sigma_sd);
  DATA_SCALAR(log_inc_lambda_mean);
  DATA_SCALAR(log_inc_lambda_sd);
  DATA_INTEGER(use_field);
  DATA_INTEGER(use_covs);
  DATA_INTEGER(use_inc);
  DATA_INTEGER(use_inc1);
  DATA_INTEGER(fix_hypers);
  DATA_SCALAR(log_rho_fixed);
  DATA_SCALAR(log_sigma_fixed);
  DATA_SCALAR(log_inc_lambda_fixed);
  
  PARAMETER(beta_0);
  PARAMETER_VECTOR(inc_mesh_vals);
  PARAMETER_VECTOR(S);
  PARAMETER_VECTOR(beta);
  PARAMETER(log_rho);
  PARAMETER(log_sigma);
  PARAMETER(log_inc_lambda);
  PARAMETER_VECTOR(beta_inc);

  
  
  
  Type f=0;
  
  Type log_rho_use, log_sigma_use, log_inc_lambda_use;
  if(fix_hypers){
    log_rho_use=log_rho_fixed;
    log_sigma_use=log_sigma_fixed;
    log_inc_lambda_use=log_inc_lambda_fixed;
  }else{
    log_rho_use=log_rho;
    log_sigma_use=log_sigma;
    log_inc_lambda_use=log_inc_lambda;
  }
  
  Type sigma = exp(log_sigma_use);
  Type rho = exp(log_rho_use);
  Type kappa = sqrt(8.0) / rho;
  SparseMatrix<Type> Q = Q_spde(spde, kappa);
  Type scaling_factor = sqrt(exp(lgamma(nu)) / (exp(lgamma(nu + 1)) * 4 * M_PI * pow(kappa, 2*nu)));

  

  

  f -= dnorm(log_rho, log_rho_mean, log_rho_sd, true);
  f -= dnorm(log_sigma, log_sigma_mean, log_sigma_sd, true);
  
  
  //create incidence covariance matrix
  Type inc_lambda = exp(log_inc_lambda_use);
  f -= dnorm(log_inc_lambda, log_inc_lambda_mean, log_inc_lambda_sd, true);
  

  int n_mesh = log_inc_mesh.size();
  
  Type beta_mean=0.0;
  for(int i=0; i<beta.size(); i++){
    f -= dnorm(beta(i), beta_mean, beta_sd, true);
  }
  // f -= dnorm(beta_inc(0), beta_mean, beta_sd, true);
  // f -= dnorm(beta_inc(1), beta_mean, beta_sd, true);
  f -= dnorm(beta_inc(0), beta_mean, beta_inc_sd, true);
  f -= dnorm(beta_inc(1), beta_mean, beta_inc_sd, true);

  matrix<Type> K_inc(n_mesh, n_mesh);
  for(int i=0; i<n_mesh; i++){
    for(int j=0; j<(i+1); j++){
      K_inc(i, j) = exp(- inc_distances(i, j) / inc_lambda);
    }
  }
  
  for(int i=0; i<n_mesh; i++){
    K_inc(i, i) += nugget;
  }

  MVNORM_t<Type> cov_inc_mvnorm(K_inc);

  f += cov_inc_mvnorm(inc_mesh_vals);
  
  vector<Type> inc_GP_vals = A_inc * inc_mesh_vals;
  vector<Type> inc_GP_vals1 = A_inc1 * inc_mesh_vals;
  vector<Type> cov_contr = X_pr * beta;
  
  
  
  int n_obs_pr = Y_pr.size();
  vector<Type> logit_prev_out(n_obs_pr);
  vector<Type> field = A_pr * S;
  Type logit_prev_i;
  for(int i=0; i<n_obs_pr; i++){
    logit_prev_i=beta_0;
    if(use_covs) logit_prev_i += cov_contr(i);
    if(use_inc & use_inc1){
  
      logit_prev_i += beta_inc(0) * inc_GP_vals(i) + beta_inc(1) * inc_GP_vals1(i);
  
    }else{
      // if(use_inc) logit_prev_i += inc_GP_vals(i);
      // if(use_inc1) logit_prev_i += inc_GP_vals1(i);
      if(use_inc){;
        // f -= dnorm(beta_inc(0), beta_mean, beta_inc_sd, true);
        logit_prev_i += 2*beta_inc(0) * inc_GP_vals(i);
      }
      if(use_inc1){
        // f -= dnorm(beta_inc(1), beta_mean, beta_inc_sd, true);
        logit_prev_i += 2*beta_inc(1) * inc_GP_vals1(i);
      }
    }
    if(use_field) logit_prev_i += field(i);
    
    if(!holdout(i)){
      f -= dbinom_robust(Y_pr(i), N_pr(i), logit_prev_i, true);
    }
    logit_prev_out(i) = logit_prev_i;
  }
  
  f += SCALE(GMRF(Q), sigma / scaling_factor)(S);

  REPORT(logit_prev_out);
  REPORT(inc_mesh_vals);
  // REPORT(inc_mesh_vals1);
  return(f);
  
}