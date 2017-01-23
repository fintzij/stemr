// [[Rcpp::depends(RcppArmadillo)]]
#include "stemr_types.h"
#include "stemr_utils.h"

using namespace Rcpp;
using namespace arma;

//' Adaptive random walk Metropolis-Hastings transition kernel based on Roberts
//' and Rosenthal (2009).
//'
//' @param params_prop vector in which the proposed parameters should be stored
//' @param params_cur vector containing the current parameter vector
//' @param covmat proposal covariance matrix
//' @param empirical_covmat empirical covariance matrix
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
//' @param opt_scaling 2.38^2/length(params)
//'
//' @return propose new parameter values in place and modify scaling and/or the empirical covariance matrix in place
//' @export
// [[Rcpp::export]]
void mvn_adaptive(arma::rowvec& params_prop, const arma::rowvec& params_cur, const arma::mat& covmat, arma::mat& empirical_covmat,
                  arma::rowvec& param_means, arma::vec& scaling, const double& iteration, const double& acceptances,
                  const arma::rowvec& nugget, const double& nugget_weight, const bool& adapt_scale, const bool& adapt_shape,
                  const int& scale_start, const int& shape_start, const double& target, const double& scale_cooling,
                  const double& max_scaling, const double& opt_scaling) {

        // if the shape should be adapted, update the empirical covariance matrix and sample means
        if(adapt_shape) {
                param_means      = ((iteration - 1) * param_means + params_cur) / iteration;
                arma::rowvec v   = params_cur - param_means;
                empirical_covmat = ((iteration - 1) * empirical_covmat + v.t() * v)/iteration;
        }

        if(adapt_scale && iteration >= scale_start && (!adapt_shape || acceptances < shape_start)) {

                // adapt the covariance scale toward the target acceptance rate.
                // N.B. 0 < scale_cooling < 1 so the finite adaptation criterion is satisfied
                double scale_factor = scaling[0] * exp(pow(scale_cooling, iteration - scale_start) * (acceptances / iteration - target));

                if(scale_factor < max_scaling) {
                        scaling[0] = scale_factor;
                } else {
                        scaling[0] = max_scaling;
                }

                // make the proposal using the scaled covariance matrix
                params_prop = params_cur + arma::randn(1,params_cur.n_elem) * arma::chol(pow(scaling[0],2.0)*covmat,"upper");

        } else if(adapt_shape && acceptances > shape_start) {

                // make the proposal using the empirical covariance matrix
                params_prop = nugget_weight*(params_cur + nugget % arma::randn(1, params_cur.n_elem)) +
                        (1-nugget_weight)*(params_cur+arma::randn(1,params_cur.n_elem)*arma::chol(opt_scaling*empirical_covmat,"upper"));

        } else {
                // make the proposal using the supplied covariance matrix
                params_prop = params_cur + arma::randn(1, params_cur.n_elem) * arma::chol(covmat, "upper");
        }
}
