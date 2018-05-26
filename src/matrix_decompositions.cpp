// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

//' matrix square root
//'
//' @param S square root matrix to be filled out
//' @param M symmetric positive definite matrix for which square root is to be computed
//' @param nugget small positive constant to be added to the diagonal for numerical stability
//' 
//' @return set S equal to the matrix square root of M 
//' @export
// [[Rcpp::export]]
void comp_sqrtmat(arma::mat& S, const arma::mat& M, double nugget) {
      S = arma::sqrtmat_sympd(M + nugget * arma::eye(arma::size(M)));
}

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
void comp_chol(arma::mat& C, const arma::mat& M, double nugget) {
      C = arma::chol(M + nugget * arma::eye(arma::size(M)));
}