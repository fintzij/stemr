// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

//' Update factors and interval widths for automated factor slice sampling
//'
//' @param 
//'
//' @return update eigenvalues and eigenvectors in place
//' @export
// [[Rcpp::export]]
void update_factors(arma::vec& interval_widths,
                    arma::vec& slice_eigenvals,
                    arma::mat& slice_factors,
                    const arma::mat& kernel_cov,
                    arma::vec& n_expansions,
                    arma::vec& n_contractions,
                    const arma::vec& n_expansions_c,
                    const arma::vec& n_contractions_c,
                    arma::vec& slice_ratios,
                    double adaptation_factor) {
      
      // update the factors
      arma::eig_sym(slice_eigenvals, slice_factors, kernel_cov); 
      
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