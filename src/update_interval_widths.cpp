// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

//' Update factors and interval widths for automated factor slice sampling
//'
//' @param interval_widths vector of interval widths
//' @param slice_eigenvals eigenvalues of the posterior covariance
//' @param width_scaling scaling factor for the interval widths
//'
//' @return adapt interval widths in place
//' @export
// [[Rcpp::export]]
void update_interval_widths(arma::vec& interval_widths,
                            const arma::vec& slice_eigenvals,
                            const double width_scaling) {

      interval_widths = width_scaling * arma::sqrt(slice_eigenvals);
}