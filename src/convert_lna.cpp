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
//' @param log_scale was the lna path computed on the log scale?
//'
//' The process can be re-expressed by left-multiplying each row in the path
//' matrix by the stoichiometry matrix: \eqn{X_t = A'\phi(t, N_t)}.
//'
//' @export
// [[Rcpp::export]]
arma::mat convert_lna(const arma::mat& path, const arma::mat& flow_matrix, const arma::rowvec& init_state, bool log_scale) {

        int n_times = path.n_rows;
        int n_comps = flow_matrix.n_cols;

        // initialize an object for the coverted path
        arma::mat conv_path(n_times, n_comps, arma::fill::zeros);

        // fill out the path
        if(log_scale) {
                for(int j = 0; j < n_times; ++j) {
                        conv_path.row(j) = (arma::exp(path.row(j)) - 1) * flow_matrix;
                }
        } else {
                for(int j = 0; j < n_times; ++j) {
                        conv_path.row(j) = path.row(j) * flow_matrix;
                }
        }

        conv_path.each_row() += init_state;
        return conv_path;
}