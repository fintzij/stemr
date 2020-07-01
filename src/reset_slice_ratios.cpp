// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

//' Reset counters for interval expansions/contractions and slice ratios
//'
//' @param n_expansions vector with number of expansion
//' @param n_contractions vector with number of contractions
//' @param n_expansions_c cumulative numbers of expansions
//' @param n_contractions_c cumulative numbers of contractions
//' @param slice_ratios vector for storing ratio of cumulative number of 
//'   expansions over number of interval width changes
//'
//' @return reset objects in place
//' @export
// [[Rcpp::export]]
void reset_slice_ratios(arma::vec& n_expansions,
                        arma::vec& n_contractions,
                        arma::vec& n_expansions_c,
                        arma::vec& n_contractions_c,
                        arma::vec& slice_ratios) {
      
      // update the expansion-contraction ratios
      n_expansions.zeros();
      n_contractions.zeros();
      n_expansions_c.ones();
      n_contractions_c.ones();
      slice_ratios.fill(0.5);
}