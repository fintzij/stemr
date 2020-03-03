// [[Rcpp::depends(Rcpp)]]

#include "stemr_types.h"

using namespace Rcpp;

//' Simulate from the measurement process by calling measurement process
//' functions via external Xptr.
//'
//' @param obsmat observation matrix
//' @param emit_inds logical vector of measurement compartments to simulate
//' @param record_ind row in the observation matrix
//' @param state numeric vector of comaprtment counts
//' @param parameters numeric vector of parameter values
//' @param constants numeric vector of constants
//' @param tcovar numeric vector of time-varying covariate values
//' @param r_meas_ptr external pointer to measurement process simulation function
//'
//' @export
// [[Rcpp::export]]
void CALL_R_MEASURE(Rcpp::NumericMatrix& obsmat, const Rcpp::LogicalVector& emit_inds,
                    const int record_ind, const Rcpp::NumericVector& state, const Rcpp::NumericVector& parameters,
                    const Rcpp::NumericVector& constants, const Rcpp::NumericVector& tcovar, SEXP r_meas_ptr) {
        Rcpp::XPtr<r_measure_ptr> xpfun(r_meas_ptr);                              // Receive the SEXP and put in Xptr
        r_measure_ptr fun = *xpfun;                                               // get function via pointer
        fun(obsmat, emit_inds, record_ind, state, parameters, constants, tcovar); // evaluate the funtion
}