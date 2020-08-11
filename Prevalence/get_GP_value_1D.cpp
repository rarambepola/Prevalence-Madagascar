#include <Rcpp.h>
#include <numeric> // for std::partial_sum used below
using namespace Rcpp;


// [[Rcpp::export]]
NumericVector get_GP_val(NumericVector mesh_points, 
                         NumericVector obs_points,
                         NumericVector mesh_gp_vals){
  
  int n_obs = obs_points.size();
  int n_mesh = mesh_points.size();
  // initialize an accumulator variable

  NumericVector obs_gp_vals(n_obs);
  
  
  for(int i=0; i<n_obs; i++){
    if(obs_points(i) < mesh_points(0)){
      obs_gp_vals(i) = mesh_gp_vals(0);
    }else{
      
      int below_index; 
      bool below_any = false;
      for(int j=0; j<n_mesh; j++){
        below_index=j;
        if(obs_points(i) < mesh_points(j)){
          below_any=true;
          break;
        }
      }
      
      if(below_any){
        double a=obs_points(i) - mesh_points(below_index-1);
        double x=mesh_gp_vals(below_index - 1);
        double y=mesh_gp_vals(below_index);
        
        obs_gp_vals(i) = x + a * (y-x) / (mesh_points(below_index) - mesh_points(below_index-1));
      }else{
        obs_gp_vals(i) = mesh_gp_vals(n_mesh-1);
      }

    }
    
  }
  
  return obs_gp_vals;
}


// [[Rcpp::export]]
NumericVector get_GP_val_faster(NumericVector mesh_points, 
                         NumericVector obs_points,
                         NumericVector mesh_gp_vals){
  
  int n_obs = obs_points.size();
  int n_mesh = mesh_points.size();
  // initialize an accumulator variable
  
  NumericVector obs_gp_vals(n_obs);
  
  //do it for i=0
  int below_index, below_index_old;
  bool below_any=false;
  bool below_same, above_same, same_as_previous;
  
  for(int i=0; i<n_obs; i++){
    if(obs_points(i) < mesh_points(0)){
      obs_gp_vals(i) = mesh_gp_vals(0);
    }else{
      same_as_previous = 0;
      //if previous was below any, check if the same one will work
      if(below_any){
        below_same = obs_points(i) < mesh_points(below_index_old);
        if(below_index_old == 0){
          same_as_previous = 1;
        }else{
          above_same = obs_points(i) > mesh_points(below_index_old - 1);
          same_as_previous = below_same & above_same;
        }
      }
      
      if(same_as_previous){
        below_index = below_index_old;
      }else{
        below_any = false;
        for(int j=0; j<n_mesh; j++){
          below_index=j;
          if(obs_points(i) < mesh_points(j)){
            below_any=true;
            break;
          }
        }
      }
      below_index_old = below_index;
      if(below_any){
        double a=obs_points(i) - mesh_points(below_index-1);
        double x=mesh_gp_vals(below_index - 1);
        double y=mesh_gp_vals(below_index);
        
        obs_gp_vals(i) = x + a * (y-x) / (mesh_points(below_index) - mesh_points(below_index-1));
      }else{
        obs_gp_vals(i) = mesh_gp_vals(n_mesh-1);
      }
      
    }
    
  }
  
  return obs_gp_vals;
}




// [[Rcpp::export]]
NumericVector get_GP_val_faster2(NumericVector mesh_points, 
                                NumericVector obs_points,
                                NumericVector mesh_gp_vals){
  
  int n_obs = obs_points.size();
  int n_mesh = mesh_points.size();
  // initialize an accumulator variable
  
  NumericVector obs_gp_vals(n_obs);
  
  //do it for i=0
  int below_index, below_index_old, index_start;
  bool below_any=false;
  bool below_same, above_same, same_as_previous;
  
  for(int i=0; i<n_obs; i++){
    if(obs_points(i) < mesh_points(0)){
      obs_gp_vals(i) = mesh_gp_vals(0);
    }else{
      same_as_previous = 0;
      //if previous was below any, check if the same one will work
      if(below_any){
        below_same = obs_points(i) < mesh_points(below_index_old);
        if(below_index_old == 0){
          same_as_previous = 1;
        }else{
          above_same = obs_points(i) > mesh_points(below_index_old - 1);
          same_as_previous = below_same & above_same;
        }
      }
      
      if(same_as_previous){
        below_index = below_index_old;
      }else{
        index_start=0;
        if(below_any & !below_same){
          index_start = below_index_old;
        }
        below_any = false;
        for(int j=index_start; j<n_mesh; j++){
          below_index=j;
          if(obs_points(i) < mesh_points(j)){
            below_any=true;
            break;
          }
        }
      }
      below_index_old = below_index;
      if(below_any){
        double a=obs_points(i) - mesh_points(below_index-1);
        double x=mesh_gp_vals(below_index - 1);
        double y=mesh_gp_vals(below_index);
        
        obs_gp_vals(i) = x + a * (y-x) / (mesh_points(below_index) - mesh_points(below_index-1));
      }else{
        obs_gp_vals(i) = mesh_gp_vals(n_mesh-1);
      }
      
    }
    
  }
  
  return obs_gp_vals;
}