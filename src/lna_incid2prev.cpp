// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace arma;
using namespace Rcpp;

//' Convert an LNA path from the counting process on transition events to the
//' compartment densities on their natural scale.
//'
//' @param path matrix containing the LNA path in terms of the counting
//'   processes on transition events (incidence)
//' @param flow_matrix stoichiometry matrix (the transpose of the flow matrix)
//' @param init_state initial compartment counts on the natural scale
//' @param forcing_inds logical vector of indicating at which times in the
//'   time-varying covariance matrix a forcing is applied.
//' @param forcing_matrix matrix containing the forcings.
//'
//' The process can be re-expressed by left-multiplying each row in the path
//' matrix by the stoichiometry matrix: \eqn{X_t = X_0 + A'N_t}.
//'
//' @export
// [[Rcpp::export]]
arma::mat lna_incid2prev(const arma::mat& path,
                         const arma::mat& flow_matrix,
                         const arma::rowvec& init_state,
                         const arma::mat& forcing_matrix,
                         const Rcpp::LogicalVector& forcing_inds,
                         const arma::uvec& forcing_tcov_inds,
                         const arma::mat& forcings_out,
                         const arma::cube& forcing_transfers) {
      
      int n_times = path.n_rows;
      int n_comps = flow_matrix.n_cols;
      int n_rates = flow_matrix.n_rows;
      int n_forcings = forcing_tcov_inds.n_elem;
      
      // for use with forcings
      double forcing_flow = 0;
      arma::rowvec forcing_distvec(n_comps, arma::fill::zeros);
      
      // initialize an object for the coverted path
      arma::mat conv_path(n_times, n_comps+1, arma::fill::zeros);
      conv_path.col(0) = path.col(0);
      
      // set the initial state
      arma::rowvec volumes = init_state + conv_path(0, arma::span(1, n_rates)) * flow_matrix;
      conv_path(0, arma::span(1,n_comps)) = volumes;
      
      if(forcing_inds[0]) {
            
            // distribute the forcings proportionally to the compartment counts in the applicable states
            for(int s=0; s < n_forcings; ++s) {
                  
                  forcing_flow     = forcing_matrix(0, forcing_tcov_inds[s]);
                  forcing_distvec  = (forcing_flow * normalise(forcings_out.col(s) % volumes.t(), 1)).t();
                  volumes         += (forcing_transfers.slice(s) * forcing_distvec.t()).t();
            }
      }
      
      // Loop through the path to compute the compartment counts
      for(int k = 1; k < n_times; ++k) {
            
            volumes += path(k, arma::span(1, n_rates)) * flow_matrix; 
            conv_path(k, arma::span(1, n_comps)) = volumes; 
            
            // apply forcings if called for - applied after censusing the path
            if(forcing_inds[k]) {
                  
                  // distribute the forcings proportionally to the compartment counts in the applicable states
                  for(int s=0; s < n_forcings; ++s) {
                        forcing_flow     = forcing_matrix(k, forcing_tcov_inds[s]);
                        forcing_distvec  = (forcing_flow * normalise(forcings_out.col(s) % volumes.t(), 1)).t();
                        volumes         += (forcing_transfers.slice(s) * forcing_distvec.t()).t();
                  }
            }
      }
      
      return conv_path;
}