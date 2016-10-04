// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace arma;
using namespace Rcpp;

//' Convert an LNA path from the counting process on transition events to the
//' compartment densities on their natural scale.
//'
//' @param path matrix containing the LNA path in terms of the counting
//'   processes on transition events
//' @param flow_matrix stoichiometry matrix (the transpose of the flow matrix)
//' @param init_state initial compartment counts on the natural scale
//'
//' The process can be re-expressed by left-multiplying each row in the path
//' matrix by the stoichiometry matrix: \eqn{X_t = A'\phi(t, N_t)}.
//'
//' @export
// [[Rcpp::export]]
arma::mat convert_lna(const arma::mat& path, const arma::mat& flow_matrix, const arma::rowvec& init_state) {

        int n_times = path.n_rows;
        int n_comps = flow_matrix.n_cols;
        int n_rates = flow_matrix.n_rows;

        // initialize an object for the coverted path
        arma::mat conv_path(n_times, n_comps + 1, arma::fill::zeros);
        conv_path.col(0) = path.col(0);

        // fill out the path
        for(int j = 0; j < n_times; ++j) {
                conv_path(j, arma::span(1, n_comps)) = path(j, arma::span(1, n_rates)) * flow_matrix + init_state;
        }

        return conv_path;
}