#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

//' Identify which rates to update based on changes in the time-varying covariates.
//'
//' @param rate_inds vector of rate indices to be modified
//' @param M time-varying covariate adjacency matrix
//' @param I logical vector indicating which covariates changed at a particular time.
//'
//' @return logical vector stating which rates need to be updated
//' @export
// [[Rcpp::export]]
void rate_update_tcovar(Rcpp::LogicalVector& rate_inds, const arma::mat& M, const arma::rowvec I) {

        int n_rates = rate_inds.size(); // number of rate functions
        arma::uvec col_inds = arma::find(I == true);

        arma::mat M_sub = M.cols(col_inds);

        for(int k=0; k < n_rates; ++k) {
                rate_inds[k] = arma::any(M_sub.row(k));
        }
}