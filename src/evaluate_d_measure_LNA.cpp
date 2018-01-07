// [[Rcpp::depends(Rcpp)]]
#include "stemr_types.h"
#include "stemr_utils.h"

using namespace Rcpp;

//' Evaluate the log-density of a possibly time-verying measurement process
//' by calling measurement process density functions via external Xptr.
//'
//' @param emitmat matrix of emission probabilities
//' @param obsmat matrix containing the data
//' @param censusmat matrix containing the state of the latent process at
//'   observation times
//' @param measproc_indmat logical matrix indicating which compartments are
//'   observed at every observation time
//' @param lna_parameters matrix containing the LNA parameters, constants and
//'   time-varying coariates.
//' @param lna_param_vec container for storing the LNA parameters at each
//'   observation time.
//' @param lna_param_inds indices for the model parameters.
//' @param lna_const_inds indices for the constants.
//' @param lna_tcovar_inds indices for the time-varying covariates.
//' @param param_update_inds logical vector indicating when model parameters
//'   should be updated.
//' @param census_indices vector of indices when the LNA path has been censused.
//' @param lna_param_vec vector for keeping the current lna parameters
//' @param d_meas_ptr external pointer to measurement process density function
//'
//' @export
// [[Rcpp::export]]
void evaluate_d_measure_LNA(Rcpp::NumericMatrix& emitmat, const Rcpp::NumericMatrix& obsmat,
                        const Rcpp::NumericMatrix& censusmat, const Rcpp::LogicalMatrix& measproc_indmat,
                        const Rcpp::NumericMatrix& lna_parameters, const Rcpp::IntegerVector& lna_param_inds,
                        const Rcpp::IntegerVector& lna_const_inds, const Rcpp::IntegerVector& lna_tcovar_inds,
                        const Rcpp::LogicalVector& param_update_inds, const Rcpp::IntegerVector& census_indices,
                        Rcpp::NumericVector& lna_param_vec, SEXP d_meas_ptr) {

        // get constants
        int n_obstimes   = obsmat.nrow();
        int n_lna_params = lna_param_inds.size();
        int n_tcovar     = lna_tcovar_inds.size();

        // initialize parameters and time-varying covariates/parameters
        std::copy(lna_parameters.row(0).begin(), 
                  lna_parameters.row(0).begin() + n_lna_params, 
                  lna_param_vec.begin());
        
        std::copy(lna_parameters.row(0).end() - n_tcovar,
                  lna_parameters.row(0).end(),
                  lna_param_vec.end() - n_tcovar);
              
        // Rcpp::NumericVector lna_param_vec(lna_parameters.ncol());    // vector for storing the current parameter values
        // std::copy(lna_parameters.row(0).begin(), lna_parameters.row(0).end(), lna_param_vec.begin());
        // lna_param_vec = lna_parameters.row(0);                       // set the current parameter values

        // evaluate the densities
        for(int j=0; j < n_obstimes; ++j) {

                // update the model parameters if called for
                // measurement process is right continuous, hence indexing by j+1, not j
                if(param_update_inds[j]) {
                      std::copy(lna_parameters.row(j+1).begin(), 
                                lna_parameters.row(j+1).begin() + n_lna_params, 
                                lna_param_vec.begin());
                      
                      std::copy(lna_parameters.row(j+1).end() - n_tcovar,
                                lna_parameters.row(j+1).end(),
                                lna_param_vec.end() - n_tcovar);
                      
                         // std::copy(lna_parameters.row(j+1).begin(), lna_parameters.row(j+1).end(), lna_param_vec.begin());
                        // lna_param_vec = lna_parameters.row(census_indices[j+1]);
                }

                // args: emitmat, emit_inds, record_ind, record, state, parameters, constants, tcovar, pointer
                CALL_D_MEASURE(emitmat, measproc_indmat.row(j), j, obsmat.row(j), censusmat.row(j),
                               lna_param_vec[lna_param_inds], lna_param_vec[lna_const_inds], 
                               lna_param_vec[lna_tcovar_inds], d_meas_ptr);
        }
}