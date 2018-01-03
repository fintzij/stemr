// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

//' Update principal component Metropolis eigenvectors and eigenvalues
//'
//' @param eigenvalues vector of eigenvalues of the proposal covariance matrix
//' @param eigenvectors matrix of eigenvectors for the proposal covariance
//' @param comp_weights vector of component weights
//' @param kernel_cov proposal covariance matrix
//' @param update_weights should the weights for the direction sampling 
//'  distribution be recomputed. 
//'
//' @return update eigenvalues and eigenvectors in place
//' @export
// [[Rcpp::export]]
void update_princomps(arma::vec& eigenvalues,
                      arma::mat& eigenvectors,
                      arma::vec& comp_weights,
                      const arma::mat& kernel_cov,
                      bool update_weights) {
      
      // update the eigenvectors and eigenvalues
      arma::eig_sym(eigenvalues, eigenvectors, kernel_cov); 
      
      // compute the new direction weights
      if(update_weights) {
            comp_weights = arma::normalise(arma::pow(eigenvalues, 0.5), 1);
      }
}