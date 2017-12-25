// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

//' Adaptive principal component metropolis proposal.
//'
//' @param params_prop vector in which the proposed parameters should be stored
//' @param params_cur vector containing the current parameter vector
//' @param eigenvalues vector of eigenvalues of the proposal covariance matrix
//' @param eigenvectors matrix of eigenvectors for the proposal covariance
//' @param proposal_scaling vector of component scaling parameters
//' @param direction component to be used in the proposal
//'
//' @return propose new parameter values in place
//' @export
// [[Rcpp::export]]
void mvn_g_adaptive(arma::rowvec& params_prop,
                    const arma::rowvec& params_cur,
                    const arma::vec& eigenvalues,
                    const arma::mat& eigenvectors,
                    const arma::vec& proposal_scaling,
                    int direction) {
  
  int par_dim = params_cur.n_elem;
  
  params_prop = params_cur + 
    R::rnorm(0, sqrt(proposal_scaling[direction] * eigenvalues[direction])) * eigenvectors.col(direction).t();
}