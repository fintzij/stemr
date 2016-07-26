// [[Rcpp::depends(Rcpp)]]

#include "stemr_types.h"
#include "stemr_utils.hpp"

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
void evaluate_d_measure(Rcpp::NumericMatrix& emitmat, const Rcpp::NumericMatrix& obsmat,
                        const Rcpp::NumericMatrix& statemat, const Rcpp::LogicalMatrix& measproc_indmat,
                        const Rcpp::NumericVector& parameters, const Rcpp::NumericVector& constants,
                        const Rcpp::NumericMatrix& tcovar_censusmat, SEXP d_meas_ptr) {

        // get dimensions
        Rcpp::IntegerVector emit_dims = emitmat.attr("dim");

        // evaluate the densities
        for(int j=0; j < emit_dims[0]; ++j) {
                // args: emitmat, emit_inds, record_ind, record, state, parameters, constants, tcovar, pointer
                CALL_D_MEASURE(emitmat, measproc_indmat.row(j), j, obsmat.row(j), statemat.row(j), parameters, constants, tcovar_censusmat.row(j), d_meas_ptr);
        }
}