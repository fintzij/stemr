// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

//' Update slice factor directions for automated factor slice sampling
//'
//' @param slice_eigenvals vector of singular values
//' @param slice_eigenvecs vector of singular vectors
//' @param kernel_cov empirical covariance matrix of model params
//' 
//' @return update eigenvalues and eigenvectors in place
//' @export
// [[Rcpp::export]]
void update_factors(arma::vec& slice_eigenvals,
                    arma::mat& slice_eigenvecs,
                    const arma::mat& kernel_cov) {
      
      // compute the eigen decomposition
      arma::eig_sym(slice_eigenvals, slice_eigenvecs, kernel_cov);
      
      // double check that all eigenvalues are non-negative
      slice_eigenvals.elem(arma::find(slice_eigenvals < 0)).zeros();
}