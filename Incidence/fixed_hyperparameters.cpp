#include <TMB.hpp>
#include <math.h>
#include <stdio.h>

const double pi = 3.141592653589793238462643383279502884;

template <class Type>
Type p2i(Type prev){
  return(2.616*prev - 3.596*pow(prev, 2) + 1.594*pow(prev, 3));
}

template <class Type>
Type invlogit_new(Type alpha){
  return(1 / (1 + exp(-alpha)));
}


template<class Type>
Type objective_function<Type>::operator() ()
{
  using namespace R_inla;
  using namespace density;
  using namespace Eigen;
  
  //incidence stuff
  DATA_VECTOR(Y);
  DATA_VECTOR(pops);
  DATA_SPARSE_MATRIX(A);
  DATA_STRUCT(spde,spde_t);
  DATA_SCALAR(prior_rho_min);
  DATA_SCALAR(prior_rho_prob);
  DATA_SCALAR(prior_sigma_max);
  DATA_SCALAR(prior_sigma_prob);
  DATA_IVECTOR(holdout_set);
  // DATA_IVECTOR(cases_not_na);
  // DATA_IVECTOR(holdout);
  DATA_SCALAR(log_rho);
  DATA_SCALAR(log_sigma);
  DATA_SCALAR(nu);
  // DATA_SCALAR(log_sigma_mean);
  // DATA_SCALAR(log_sigma_sd);
  // DATA_SCALAR(log_rho_mean);
  // DATA_SCALAR(log_rho_sd);
  // DATA_SCALAR(disp);
  
  
  PARAMETER(beta_0);
  PARAMETER_VECTOR(S);
  // PARAMETER(log_rho);
  // PARAMETER(log_sigma);
  
  Type f=0;
  
  Type sigma = exp(log_sigma);
  Type rho = exp(log_rho);
  Type kappa = sqrt(8.0) / rho;
  SparseMatrix<Type> Q = Q_spde(spde, kappa);
  
  // Type lambdatilde1 = -log(prior_rho_prob) * prior_rho_min;
  // Type lambdatilde2 = -log(prior_sigma_prob) / prior_sigma_max;
  // 
  // Type log_pcdensity = log(lambdatilde1) + log(lambdatilde2) - 2 * log_rho  
  //   -lambdatilde1 * pow(rho, -1) - lambdatilde2 * sigma;
  // f -= log_pcdensity + log_rho + log_sigma ;
  
  // f -= dnorm(log_sigma, log_sigma_mean, log_sigma_sd, true);
  // f -= dnorm(log_rho, log_rho_mean, log_rho_sd, true);  
  
  Type scaling_factor = sqrt(exp(lgamma(nu)) / (exp(lgamma(nu + 1)) * 4 * pi * pow(kappa, 2*nu)));
  
  //priors
  Type beta_mean = 0.0;
  Type beta_0_sd = 1.0;
  int n_obs = Y.size();
  int n_field = S.size();

  
  f -= dnorm(beta_0, beta_mean, beta_0_sd, true);
  
  
  vector<Type> field = A * S;
  
  Type log_rate;
  vector<Type> log_rate_vec(n_obs);
  
  for(int i=0; i<n_obs; i++){
    log_rate = beta_0 + field(i);
    log_rate_vec(i) = log_rate;
    if(!holdout_set(i)){
      f -= dpois(Y(i), pops(i) * exp(log_rate), true);
    }
  }
  

  
  f += SCALE(GMRF(Q), sigma / scaling_factor)(S);
  
  // REPORT(Q);
  // REPORT(log_tau);
  // REPORT(Sigma_beta);
  REPORT(log_rate_vec);
  return(f);
  
}