// [[Rcpp::depends(Rcpp)]]

#include "stemr_types.h"

using namespace Rcpp;

//' Compute the Jacobian by calling the jacobian function via external Xptr.
//'
//' @param t time
//' @param state vector of compartment counts and time-varying covariate values
//' @param parameters vector of model parameters and constants
//' @param jacobian_ptr external pointer to function used to compute the jacobian
//'
//' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix COMPUTE_JACOBIAN(double t, const Rcpp::NumericVector& state, const Rcpp::NumericVector& parameters, SEXP jacob_ptr) {

        Rcpp::XPtr<jacobian_ptr> xpfun(jacob_ptr);      // Receive the SEXP and put in Xptr
        jacobian_ptr fun = *xpfun;                      // get function via pointer

        return fun(t, state, parameters);               // Compute the Jacobian and return it
}