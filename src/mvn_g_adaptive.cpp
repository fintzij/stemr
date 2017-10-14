// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

//' Global Metropolis random walk with global adaptive scaling
//'
//' @param params_prop vector in which the proposed parameters should be stored
//' @param params_cur vector containing the current parameter vector
//' @param kernel_cov vector of component proposal standard deviations
//' @param proposal_scaling scaling parameter for the proposal
//' @param nugget fixed covariance nugget contribution
//'
//' @return propose new parameter values in place
//' @export
// [[Rcpp::export]]
void mvn_g_adaptive(arma::rowvec& params_prop,
                    const arma::rowvec& params_cur,
                    const arma::mat& kernel_cov,
                    double proposal_scaling,
                    double nugget) {

        int par_dim = params_cur.n_elem;

        params_prop = params_cur + nugget * arma::randn(1, par_dim) +
                (1 - nugget) * proposal_scaling * arma::randn(1, par_dim) * arma::chol(kernel_cov, "upper");

}