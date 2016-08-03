// [[Rcpp::depends(Rcpp)]]

#include "stemr_types.h"

using namespace Rcpp;

//' Compute the hazards by calling the hazard functions via external Xptr.
//'
//' @param t time
//' @param state vector of compartment counts
//' @param parameters vector of model parameters
//' @param constants vector of constants
//' @param vector of time-varying covariate values
//' @param hazard_ptr external pointer to function used to compute the hazards
//'
//' @export
// [[Rcpp::export]]
arma::vec COMPUTE_HAZARD(double t, const arma::vec& state, const Rcpp::NumericVector& parameters, const Rcpp::NumericVector& constants, const Rcpp::NumericVector& tcovar, SEXP haz_ptr) {

        Rcpp::XPtr<hazard_ptr> xpfun(haz_ptr);      // Receive the SEXP and put in Xptr
        hazard_ptr fun = *xpfun;                    // get function via pointer

        return fun(t, state, parameters, constants, tcovar);           // Compute the Jacobian and return it
}