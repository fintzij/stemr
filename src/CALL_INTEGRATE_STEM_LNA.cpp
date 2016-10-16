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
//' @param hazard_ptr external pointer to the LNA functions
//'
//' @export
// [[Rcpp::export]]
void CALL_INTEGRATE_STEM_LNA(Rcpp::NumericVector& init, double start, double end, double step_size, SEXP lna_ode_ptr) {

        Rcpp::XPtr<lna_ptr> xpfun(lna_ode_ptr); // Receive the SEXP and put in Xptr
        lna_ptr fun = *xpfun;                   // get function via pointer

        return fun(init, start, end, step_size);
}