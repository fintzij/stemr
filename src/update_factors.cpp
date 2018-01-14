// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

//' Update slice factor directions for automated factor slice sampling
//'
//' @param slice_singvals vector of singular values
//' @param slice_factors vector of singular vectors
//' @param slice_factors_t transpost matrix for singular vectors
//' @param kernel_cov empirical covariance matrix of model params
//' 
//' @return update eigenvalues and eigenvectors in place
//' @export
// [[Rcpp::export]]
void update_factors(arma::vec& slice_singvals,
                    arma::mat& slice_factors,
                    arma::mat& slice_factors_t,
                    const arma::mat& kernel_cov) {
      
      arma::svd_econ(slice_factors, slice_singvals, slice_factors_t, kernel_cov, "left"); 
}