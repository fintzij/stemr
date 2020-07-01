// [[Rcpp::depends(RcppArmadillo)]]
#define ARMA_DONT_PRINT_ERRORS
#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

//' Cholesky decomposition
//'
//' @param C matrix to be filled out with the cholesky of M
//' @param M symmetric positive definite matrix for which the upper triangle of the cholesky 
//'   is to be computed
//' @param nugget small positive constant to be added to the diagonal for numerical stability
//' 
//' @return set C equal to the matrix square root of M 
//' @export
// [[Rcpp::export]]
void comp_chol(arma::mat& C, arma::mat& M) {
      
      // attempt the cholesky
      bool success = arma::chol(C, M);
      
      // enforce diagonal dominance if necessary
      if(!success) {
            M.diag() = arma::max(M.diag(), arma::sum(arma::abs(M),1) - arma::abs(M.diag()));
            C = arma::chol(M);
      }
}