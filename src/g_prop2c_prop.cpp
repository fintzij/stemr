// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

//' Get componentwise proposals from a global proposal.
//'
//' @param g2c_mat matrix in which to keep the componentwise proposals
//' @param params_cur vector containing the current parameter vector
//' @param params_prop vector in which the proposed parameters
//'
//' @return fill g2c_mat with componentwise proposals
//' @export
// [[Rcpp::export]]
void g_prop2c_prop(arma::mat& g2c_mat, const arma::rowvec& params_cur, const arma::rowvec& params_prop) {

        g2c_mat.each_row() = params_cur;
        g2c_mat.diag()     = params_prop;
}