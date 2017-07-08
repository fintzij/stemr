// [[Rcpp::depends(Rcpp)]]

#include "stemr_types.h"

using namespace Rcpp;

//' Compute the hazards by calling the hazard functions via external Xptr.
//'
//' @param init initial condition
//' @param state time at left endpoint of interval
//' @param end time at right endpoint
//' @param step_size set automatically by caller, required argument not specified by user
//' @param lna_ode_ptr external pointer for calling the ODE integrator
//'
//' @export
// [[Rcpp::export]]
void CALL_INTEGRATE_STEM_LNA(Rcpp::NumericVector& init, double start, double end, double step_size, SEXP lna_ode_ptr) {

        Rcpp::XPtr<lna_ptr> xpfun(lna_ode_ptr); // Receive the SEXP and put in Xptr
        lna_ptr fun = *xpfun;                   // get function via pointer

        return fun(init, start, end, step_size);
}