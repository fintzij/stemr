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
void mvn_c_adaptive(arma::rowvec& params_prop,
                    const arma::rowvec& params_cur,
                    const arma::mat& kernel_cov,
                    const arma::vec& proposal_scaling,
                    arma::mat& sqrt_scalemat,
                    double nugget) {

        int par_dim = params_cur.n_elem;
        sqrt_scalemat.diag() = arma::sqrt(proposal_scaling);

        params_prop = params_cur + nugget * arma::randn(1, par_dim) +
                (1 - nugget) * arma::randn(1, par_dim) * arma::chol((sqrt_scalemat * kernel_cov * sqrt_scalemat), "upper");

}