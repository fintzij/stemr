#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

//' Identify which rates to update based on changes in the time-varying covariates.
//'
//' @param M time-varying covariate adjacency matrix
//' @param I logical vector indicating which covariates changed at a particular time.
//'
//' @return logical vector stating which rates need to be updated
//' @export
// [[Rcpp::export]]
Rcpp::LogicalVector rates_to_update(const arma::mat& M, const arma::rowvec I) {

        int n_rates = M.n_rows; // number of rate functions
        arma::uvec col_inds = arma::find(I == true);

        arma::mat M_sub = M.cols(col_inds);

        Rcpp::LogicalVector R(n_rates, FALSE); // logical vector of rates to update

        for(int k=0; k < n_rates; ++k) {
                R[k] = arma::any(M_sub.row(k));
        }

        return R;
}
