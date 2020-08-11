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
  
  //(log) incidence data
  DATA_VECTOR(log_inc_vals_simple);
  DATA_VECTOR(log_inc_mesh);
  DATA_MATRIX(A_inc);
  DATA_SCALAR(nugget);
  DATA_IVECTOR(holdout);
  DATA_SCALAR(nu);
  
  //prior settings
  // DATA_SCALAR(prior_rho_min);
  // DATA_SCALAR(prior_rho_prob);
  // DATA_SCALAR(prior_sigma_max);
  // DATA_SCALAR(prior_sigma_prob);
  DATA_SCALAR(log_rho_mean);
  DATA_SCALAR(log_rho_sd);
  DATA_SCALAR(log_sigma_mean);
  DATA_SCALAR(log_sigma_sd);
  DATA_SCALAR(log_inc_lambda_mean);
  DATA_SCALAR(log_inc_lambda_sd);
  
  PARAMETER(beta_0);
  PARAMETER_VECTOR(inc_mesh_vals);
  PARAMETER_VECTOR(S);
  PARAMETER_VECTOR(beta);
  PARAMETER(log_rho);
  PARAMETER(log_sigma);
  PARAMETER(log_inc_lambda);

  
  
  
  Type f=0;
  
  Type sigma = exp(log_sigma);
  Type rho = exp(log_rho);
  Type kappa = sqrt(8.0) / rho;
  SparseMatrix<Type> Q = Q_spde(spde, kappa);
  Type scaling_factor = sqrt(exp(lgamma(nu)) / (exp(lgamma(nu + 1)) * 4 * M_PI * pow(kappa, 2*nu)));
  // Type lambdatilde1 = -log(prior_rho_prob) * prior_rho_min;
  // Type lambdatilde2 = -log(prior_sigma_prob) / prior_sigma_max;
  // Type log_pcdensity = log(lambdatilde1) + log(lambdatilde2) - 2 * log_rho  
  //   -lambdatilde1 * pow(rho, -1) - lambdatilde2 * sigma;
  // f -= log_pcdensity + log_rho + log_sigma ;
  // 

  f -= dnorm(log_rho, log_rho_mean, log_rho_sd, true);
  f -= dnorm(log_sigma, log_sigma_mean, log_sigma_sd, true);
  
  
  //create incidence covariance matrix
  Type inc_lambda = exp(log_inc_lambda);
  f -= dnorm(log_inc_lambda, log_inc_lambda_mean, log_inc_lambda_sd, true);
  
  
  int n_mesh = log_inc_mesh.size();
  
  Type beta_mean=0.0;
  // Type beta_sd=0.5;
  for(int i=0; i<beta.size(); i++){
    f -= dnorm(beta(i), beta_mean, beta_sd, true);
  }

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
  vector<Type> cov_contr = X_pr * beta;
  
  // std::cout << inc_GP_vals.size() << "\n";
  
  
  int n_obs_pr = Y_pr.size();
  vector<Type> logit_prev_out(n_obs_pr);
  vector<Type> field = A_pr * S;
  Type logit_prev_i;
  for(int i=0; i<n_obs_pr; i++){
    logit_prev_i = beta_0 + inc_GP_vals(i) + cov_contr(i) + field(i);
    // logit_prev_i = beta_0 + inc_GP_vals(i) + field(i);
    if(!holdout(i)){
      f -= dbinom_robust(Y_pr(i), N_pr(i), logit_prev_i, true);
    }
    logit_prev_out(i) = logit_prev_i;
  }
  
  f += SCALE(GMRF(Q), sigma / scaling_factor)(S);

  REPORT(logit_prev_out);
  REPORT(inc_mesh_vals);
  return(f);
  
}