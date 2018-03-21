// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

//' Update factors and interval widths for automated factor slice sampling
//'
//' @param interval_widths vector of interval widths
//' @param contraction_rates rates of contractions along each slice direction
//' @param adaptation_factor weight for adaptation
//' @param target_contraction_rate expected number of contractions per iteration
//'
//' @return adapt interval widths in place
//' @export
// [[Rcpp::export]]
void update_interval_widths(arma::vec& interval_widths,
                            const arma::vec& contraction_rates,
                            const double adaptation_factor,
                            const double target_contraction_rate) {
      interval_widths = 
            arma::exp(arma::log(interval_widths) + adaptation_factor * (contraction_rates - target_contraction_rate)); 
}