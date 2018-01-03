// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

//' Update factors and interval widths for automated factor slice sampling
//'
//' @param interval_widths vector of interval widths
//' @param slice_singvals vector of singular values
//' @param slice_factors vector of singular vectors
//' @param slice_factors_t transpost matrix for singular vectors
//' @param kernel_cov empirical covariance matrix of model params
//' @param n_expansions vector with number of expansion
//' @param n_contractions vector with number of contractions
//' @param n_expansions_c cumulative numbers of expansions
//' @param n_contractions_c cumulative numbers of contractions
//' @param slice_ratios vector for storing ratio of cumulative number of 
//'   expansions over number of interval width changes
//' @param adaptation_factor 
//' @param adapt_factors should the factors and singular values be updated?
//'
//' @return update eigenvalues and eigenvectors in place
//' @export
// [[Rcpp::export]]
void update_factors(arma::vec& interval_widths,
                    arma::vec& slice_singvals,
                    arma::mat& slice_factors,
                    arma::mat& slice_factors_t,
                    const arma::mat& kernel_cov,
                    arma::vec& n_expansions,
                    arma::vec& n_contractions,
                    const arma::vec& n_expansions_c,
                    const arma::vec& n_contractions_c,
                    arma::vec& slice_ratios,
                    double adaptation_factor,
                    bool adapt_factors) {
      
      if(adapt_factors) {
            arma::svd_econ(slice_factors, slice_singvals, slice_factors_t, kernel_cov, "left"); 
      }
      
      // update the expansion-contraction ratios
      slice_ratios = n_expansions_c / (n_expansions_c + n_contractions_c);
      
      // substitute the slice ratios for where there are zero expansions or contractions
      arma::uvec exp_zeros = arma::find(n_expansions == 0);
      arma::uvec con_zeros = arma::find(n_contractions == 0);
      n_expansions.elem(exp_zeros) = slice_ratios.elem(exp_zeros);
      n_contractions.elem(con_zeros) = slice_ratios.elem(con_zeros);
      
      // update the interval widths -- Robbins-Monro recursion
      interval_widths = arma::exp(
            log(interval_widths) +
                  adaptation_factor * (n_expansions / (n_expansions + n_contractions) - 0.5));
      
      // reset the number of contractions and expansions
      n_expansions.zeros();
      n_contractions.zeros();      
}