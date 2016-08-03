// [[Rcpp::depends(Rcpp)]]

#include "stemr_types.h"

using namespace Rcpp;

//' Compute the Jacobian by calling the jacobian function via external Xptr.
//'
//' @param t time
//' @param state vector of compartment counts
//' @param parameters vector of model parameters, time-varying covariates, and constants
//' @param jacobian_ptr external pointer to function used to compute the jacobian
//'
//' @export
// [[Rcpp::export]]
arma::mat COMPUTE_JACOBIAN(double t, const arma::vec& state, const Rcpp::NumericVector& parameters, const Rcpp::NumericVector& constants, const Rcpp::NumericVector& tcovar, SEXP jacob_ptr) {

        Rcpp::XPtr<jacobian_ptr> xpfun(jacob_ptr);      // Receive the SEXP and put in Xptr
        jacobian_ptr fun = *xpfun;                      // get function via pointer

        return fun(t, state, parameters, constants, tcovar);               // Compute the Jacobian and return it
}