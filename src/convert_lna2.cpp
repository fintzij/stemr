// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace arma;
using namespace Rcpp;

//' Convert an LNA path from the counting process on transition events to the
//' compartment densities on their natural scale, making the conversion in place
//' for an existing census matrix.
//'
//' @param path matrix containing the LNA path in terms of the counting
//'   processes on transition events
//' @param flow_matrix stoichiometry matrix (the transpose of the flow matrix)
//' @param init_state initial compartment counts on the natural scale
//' @param statemat matrix where the compartment counts should be written
//'
//' The process can be re-expressed by left-multiplying each row in the path
//' matrix by the stoichiometry matrix: \eqn{X_t = A'\phi(t, N_t)}.
//'
//' @export
// [[Rcpp::export]]
void convert_lna2(const arma::mat& path, const arma::mat& flow_matrix, const arma::rowvec& init_state, arma::mat& statemat) {

        int n_times = path.n_rows;
        int n_comps = flow_matrix.n_cols;
        int n_rates = flow_matrix.n_rows;

        // transform the LNA path
        statemat.cols(1, n_comps) = arma::repmat(init_state, n_times, 1) + path.cols(1,n_rates) * flow_matrix;
}