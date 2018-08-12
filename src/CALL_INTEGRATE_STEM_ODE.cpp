// [[Rcpp::depends(Rcpp)]]

#include "stemr_types.h"

using namespace Rcpp;

//' Integrate a system of ODEs via external Xptr.
//'
//' @param init initial condition
//' @param start time at left endpoint of interval
//' @param end time at right endpoint
//' @param step_size set automatically by caller, required argument not specified by user
//' @param stem_ode_ptr external pointer for calling the ODE integrator
//'
//' @export
// [[Rcpp::export]]
void CALL_INTEGRATE_STEM_ODE(Rcpp::NumericVector& init, double start, double end, double step_size, SEXP stem_ode_ptr) {

        Rcpp::XPtr<ode_ptr> xpfun(stem_ode_ptr); // Receive the SEXP and put in Xptr
        ode_ptr fun = *xpfun;                   // get function via pointer

        return fun(init, start, end, step_size);
}