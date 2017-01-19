// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

//' Adaptive random walk Metropolis-Hastings transition kernel based on Roberts
//' and Rosenthal (2009).
//'
//' @param params_prop vector in which the proposed parameters should be stored
//' @param params_cur vector containing the current parameter vector
//' @param covmat proposal covariance matrix
//' @param param_means sample means of each of the parameters
//' @param scaling scale factor
//' @param iteration MCMC iteration number
//' @param acceptances number of acceptances
//' @param nugget adaptive RWMH nugget distribution
//' @param nugget_weight adaptive RWMH nugget contribution
//' @param scale_start iteration at which to begin scaling the supplied covariance matrix
//' @param shape_start iteration at which to begin estimating the empirical covariance matrix
//' @param target target acceptance rate when varying the scale
//' @param scale_cooling scale cooling rate
//' @param max_scaling maximum scale factor
//'
//' @return propose new parameter values in place and modify scaling and/or the empirical covariance matrix in place
//' @export
// [[Rcpp::export]]
void mvn_rw(arma::rowvec& params_prop, const arma::rowvec& params_cur, arma::mat& covmat, arma::rowvec& param_means,
            double& scaling, const int& iteration, const int& acceptances, const arma::rowvec& nugget,
            const double& nugget_weight, const int& scale_start, const int& shape_start, const double& target,
            const double& scale_cooling, const double& max_scaling) {


}