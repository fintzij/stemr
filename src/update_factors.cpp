// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

//' Update slice factor directions for automated factor slice sampling
//'
//' @param slice_singvals vector of singular values
//' @param slice_factors vector of singular vectors
//' @param kernel_cov empirical covariance matrix of model params
//' @param I_mat identity matrix of the same size as kernel_cov
//' 
//' @return update eigenvalues and eigenvectors in place
//' @export
// [[Rcpp::export]]
void update_factors(arma::vec& slice_singvals,
                    arma::mat& slice_factors,
                    arma::mat& slice_factors_t,
                    const arma::mat& kernel_cov) {
      
      arma::svd(slice_factors, slice_singvals, slice_factors_t, kernel_cov);
}