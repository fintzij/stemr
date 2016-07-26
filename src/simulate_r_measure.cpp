// [[Rcpp::depends(Rcpp)]]
#include "stemr_types.h"
#include "stemr_utils.hpp"

using namespace Rcpp;

//' Simulate a data matrix from the measurement process of a stochastic epidemic
//' model.
//'
//' @param censusmat matrix of compartment counts at observation times
//' @param measproc_indmat logical matrix for which measure variables are
//'        observed at which times
//' @param parameters numeric vector of model parameters
//' @param constants numeric vector of constants
//' @param tcovar numeric matrix of time-varying covariate values at observation
//'        times.
//' @param r_measure_ptr external pointer to measurement process simulation fcn
//'
//' @return matrix with a simulated dataset from a stochastic epidemic model.
//' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix simulate_r_measure(Rcpp::NumericMatrix& censusmat, Rcpp::LogicalMatrix& measproc_indmat, Rcpp::NumericVector& parameters, Rcpp::NumericVector& constants, Rcpp::NumericMatrix& tcovar, SEXP r_measure_ptr) {

        // Get object dimensions
        Rcpp::IntegerVector obsmat_dims = measproc_indmat.attr("dim");

        // create observation matrix
        Rcpp::NumericMatrix obsmat(obsmat_dims[0], obsmat_dims[1] + 1);
        obsmat(_, 0) = censusmat(_, 0); // copy the observation times

        // simulate the dataset
        for(int j=0; j < obsmat_dims[0]; ++j) {
                // obsmat, emit_inds, record_ind, state, parameters, constants, tcovar, ptr
                CALL_R_MEASURE(obsmat, measproc_indmat(j, _), j, censusmat.row(j), parameters, constants, tcovar(j, _), r_measure_ptr);
        }

        return obsmat;
}