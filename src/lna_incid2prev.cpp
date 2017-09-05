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
                         const Rcpp::LogicalVector& forcing_inds,
                         const arma::mat& forcing_matrix) {

        int n_times = path.n_rows;
        int n_comps = flow_matrix.n_cols;
        int n_rates = flow_matrix.n_rows;

        // initialize an object for the coverted path
        arma::mat conv_path(n_times, n_comps+1);
        conv_path.col(0) = path.col(0);

        // set the initial state
        conv_path(0, arma::span(1,n_comps)) = init_state;

        // Loop through the path to compute the compartment counts
        for(int k = 1; k < n_times; ++k) {
                conv_path(k, arma::span(1, n_comps)) =
                        conv_path(k-1, arma::span(1, n_comps)) + path(k, arma::span(1, n_rates)) * flow_matrix +
                        forcing_matrix.col(k-1).t();
        }

        return conv_path;
}