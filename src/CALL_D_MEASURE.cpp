// [[Rcpp::depends(Rcpp)]]

#include "stemr_types.h"

using namespace Rcpp;

//' Evaluate the log-density of the measurement process by calling measurement
//' process density functions via external Xptr.
//'
//' @param emitmat matrix of emission probabilities
//' @param emit_inds logical vector of measurement compartments to simulate
//' @param record_ind row in the observation and emission matrices
//' @param record vector of observed counts
//' @param state numeric vector of latent comaprtment counts
//' @param parameters numeric vector of parameter values
//' @param constants numeric vector of constants
//' @param tcovar numeric vector of time-varying covariate values
//' @param r_meas_ptr external pointer to measurement process simulation function
//'
//' @export
// [[Rcpp::export]]
void CALL_D_MEASURE(Rcpp::NumericMatrix& emitmat, const Rcpp::LogicalVector& emit_inds,
                    const int record_ind, const Rcpp::NumericVector& record, const Rcpp::NumericVector& state,
                    const Rcpp::NumericVector& parameters, const Rcpp::NumericVector& constants,
                    const Rcpp::NumericVector& tcovar, SEXP d_meas_ptr) {
        Rcpp::XPtr<d_measure_ptr> xpfun(d_meas_ptr);                              // Receive the SEXP and put in Xptr
        d_measure_ptr fun = *xpfun;                                               // get function via pointer
        fun(emitmat, emit_inds, record_ind, record, state, parameters, constants, tcovar); // evaluate the funtion
}