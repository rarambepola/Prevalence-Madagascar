#include <Rcpp.h>
#include <numeric> 
#include <algorithm>
// #include <boost/math/distributions/normal.hpp>
using namespace Rcpp;


// [[Rcpp::export]]
NumericVector get_quantiles_vec(NumericVector vals, NumericVector ps){
  NumericVector vals_copy=vals;
  std::sort(vals_copy.begin(), vals_copy.end());
  NumericVector quantiles_out(ps.size());
  
  int n_vals = vals.size();
  double p_i;
  int i_index;
  double i_diff;
  for(int i=0; i<ps.size(); i++){
    p_i = ps[i];
    i_index = 1 + (n_vals-1) * p_i;
    i_diff = 1 + ((double) n_vals-1)*p_i - (double)i_index;
    
    quantiles_out(i) = (vals_copy[i_index-1]*(1-i_diff)) + (vals_copy[i_index] * i_diff);
  }
  return(quantiles_out);
}

// // [[Rcpp::export]]
// NumericVector get_quantiles_boost_vec(NumericVector vals, NumericVector ps){
//   std::vector<double> vals_copy(vals.size());
//   for(int i=0; i<vals.size(); i++){
//     vals_copy[i] = vals(i);
//   }
//   NumericVector quantiles_out(ps.size());
//   double p_i;
//   for(int i=0; i<ps.size(); i++){
//     p_i = ps[i];
//     double quantile_out_i = boost::math::quantile(vals_copy, p_i);
//   }
//   return(quantiles_out);
// }


// [[Rcpp::export]]
NumericVector get_col_quantiles(NumericMatrix vals, NumericVector ps){
  NumericMatrix quantiles_out(ps.size(), vals.cols());
  
  for(int i=0; i<vals.cols(); i++){
    quantiles_out(_, i) = get_quantiles_vec(vals(_, i), ps);
  }
  
  return(quantiles_out);
}