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
void CALL_SET_ODE_PARAMS(Rcpp::NumericVector& p, SEXP set_ode_params_ptr) {

        Rcpp::XPtr<set_pars_ptr> xpfun(set_ode_params_ptr); // Receive the SEXP and put in Xptr
        set_pars_ptr fun = *xpfun;                   // get function via pointer

        fun(p);
}