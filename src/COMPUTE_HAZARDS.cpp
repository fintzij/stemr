// [[Rcpp::depends(Rcpp)]]

#include "stemr_types.h"

using namespace Rcpp;

//' Compute the hazards by calling the hazard functions via external Xptr.
//'
//' @param t time
//' @param state vector of compartment counts and time-varying covariate values
//' @param parameters vector of model parameters and constants
//' @param hazard_ptr external pointer to function used to compute the hazards
//'
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector COMPUTE_HAZARD(double t, const Rcpp::NumericVector& state, const Rcpp::NumericVector& parameters, SEXP haz_ptr) {

        Rcpp::XPtr<hazard_ptr> xpfun(haz_ptr);      // Receive the SEXP and put in Xptr
        hazard_ptr fun = *xpfun;                      // get function via pointer

        return fun(t, state, parameters);               // Compute the Jacobian and return it
}