// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include "stemr_types.h"

using namespace Rcpp;
using namespace arma;

//' Update rates by calling rate functions via Xptr.
//'
//' @param rates vector of rates to be modified
//' @param inds logical vector of indices of rates to be modified
//' @param state numeric vector of comaprtment counts
//' @param parameters numeric vector of parameter values
//' @param constants numeric vector of constants
//' @param tcovar numeric vector of time-varying covariate values
//' @param rate_ptr external pointer to rate function
//'
//' @export
// [[Rcpp::export]]
void CALL_RATE_FCN(Rcpp::NumericVector& rates, const Rcpp::LogicalVector& inds,
                   const arma::rowvec& state, const Rcpp::NumericVector& parameters,
                   const Rcpp::NumericVector& constants, const arma::rowvec& tcovar, SEXP rate_ptr) {
        Rcpp::XPtr<ratefcn_ptr> xpfun(rate_ptr);                // Receive the SEXP and put in Xptr
        ratefcn_ptr fun = *xpfun;                               // get function via pointer
        fun(rates, inds, state, parameters, constants, tcovar); // evaluate the funtion
}