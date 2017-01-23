#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

//' Identify which rates to update when a state transition event occurs.
//'
//' @param lnapars matrix of lna parameters, constants, and time-varying covars
//' @param parameters vector of parameters to be copied into the matrix
//'
//' @return modifies the lna parameter matrix in place
//' @export
// [[Rcpp::export]]
void pars2lnapars(arma::mat& lnapars, const arma::rowvec& parameters) {

        int n_pars  = parameters.n_elem;
        lnapars.cols(0, n_pars-1).each_row() = parameters;
}

//' Identify which rates to update when a state transition event occurs.
//'
//' @param dest destination row vector
//' @param orig origin row vector
//'
//' @return copy the elements of one row vector into another.
//' @export
// [[Rcpp::export]]
void copypars(arma::rowvec& dest, const arma::rowvec& orig) {

        dest = orig;
}

