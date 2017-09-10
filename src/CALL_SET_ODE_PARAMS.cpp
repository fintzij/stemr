// [[Rcpp::depends(Rcpp)]]
#include "stemr_types.h"

using namespace Rcpp;

//' Set the parameters for a system of ODEs via XPtr.
//'
//' @param p vector of parameters
//' @param set_ode_params_ptr external pointer to the ODE parameter setting function.
//'
//' @export
// [[Rcpp::export]]
void CALL_SET_ODE_PARAMS(Rcpp::NumericVector& p, SEXP set_ode_params_ptr) {

        Rcpp::XPtr<set_pars_ptr> xpfun(set_ode_params_ptr); // Receive the SEXP and put in Xptr
        set_pars_ptr fun = *xpfun;                          // get function via pointer

        // set the parameters
        fun(p);
}