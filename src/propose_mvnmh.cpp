// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

//' Multivariate normal Metropolis-Hastings proposal
//'
//' @param params_prop vector in which the proposed parameters should be stored
//' @param params_cur vector containing the current parameter vector
//' @param kernel_cov_chol cholesky of the kernel covariance
//' @param nugget zero if adaptation is not ongoing
//'
//' @return propose new parameter values in place
//' @export
// [[Rcpp::export]]
void propose_mvnmh(arma::rowvec& params_prop,
                    const arma::rowvec& params_cur,
                    const arma::mat& kernel_cov_chol,
                    double nugget) {

    int par_dim = params_cur.n_elem;

    if(nugget == 0) {
        params_prop = 
            params_cur + 
            Rcpp::as<arma::rowvec>(Rcpp::rnorm(par_dim)) * kernel_cov_chol;
    } else {
        params_prop = 
            params_cur + 
            nugget * Rcpp::as<arma::rowvec>(Rcpp::rnorm(par_dim)) +
            (1 - nugget) * Rcpp::as<arma::rowvec>(Rcpp::rnorm(par_dim)) * kernel_cov_chol;    
    }
}