#ifndef stemr_UTILITIES_H
#define stemr_UTILITIES_H

#include <algorithm>

using namespace Rcpp;
using namespace arma;

// call rate functions
void CALL_RATE_FCN(Rcpp::NumericVector& rates, const Rcpp::LogicalVector& inds,
                   const arma::rowvec& state, const Rcpp::NumericVector& parameters,
                   const Rcpp::NumericVector& constants, const arma::rowvec& tcovar, SEXP rate_ptr);

// evaluate an emission density and update an existing emission matrix
void CALL_D_MEASURE(Rcpp::NumericMatrix& emitmat, const Rcpp::LogicalVector& emit_inds,
                    const int record_ind, const Rcpp::NumericVector& record, const Rcpp::NumericVector& state,
                    const Rcpp::NumericVector& parameters, const Rcpp::NumericVector& constants,
                    const Rcpp::NumericVector& tcovar, SEXP d_meas_ptr);

// simulate from the measurement process and update an existing observation matrix
void CALL_R_MEASURE(Rcpp::NumericMatrix& obsmat, const Rcpp::LogicalVector& emit_inds,
                    const int record_ind, const Rcpp::NumericVector& state, const Rcpp::NumericVector& parameters,
                    const Rcpp::NumericVector& constants, const Rcpp::NumericVector& tcovar, SEXP r_meas_ptr);

// integrate the LNA odes, call via XPtr
// Rcpp::NumericVector CALL_INTEGRATE_STEM_LNA(Rcpp::NumericVector& init, double start, double end, double step_size, SEXP lna_ode_ptr);
void CALL_INTEGRATE_STEM_LNA(Rcpp::NumericVector& init, double start, double end, double step_size, SEXP lna_ode_ptr);

// set the LNA parameters, call via XPtr
void CALL_SET_LNA_PARAMS(Rcpp::NumericVector& p, SEXP set_lna_params_ptr);

// update rates based on transition events or changes in time-varying covariates
void rate_update_tcovar(Rcpp::LogicalVector& rate_inds, const arma::mat& M, const arma::rowvec I);
void rate_update_event(Rcpp::LogicalVector& rate_inds, const Rcpp::LogicalMatrix& M, int event_code);

// gillespie simulation
arma::mat simulate_gillespie(const arma::mat& flow, const Rcpp::NumericVector& parameters, const Rcpp::NumericVector& constants, const arma::mat& tcovar, const arma::rowvec& init_states, const Rcpp::LogicalMatrix& rate_adjmat, const arma::mat& tcovar_adjmat, const arma::mat& tcovar_changemat, const Rcpp::IntegerVector init_dims, SEXP rate_ptr);

// simulation from the measurement process
Rcpp::NumericMatrix simulate_r_measure(const Rcpp::NumericMatrix& censusmat, const Rcpp::LogicalMatrix& measproc_indmat, const Rcpp::NumericVector& parameters, const Rcpp::NumericVector& constants, const arma::mat& tcovar, SEXP r_measure_ptr);

// build a census matrix with compartment counts at observation times
arma::mat build_census_path(Rcpp::NumericMatrix& path, Rcpp::NumericVector& census_times, Rcpp::IntegerVector& census_columns);

// update a census matrix with compartment counts at observation times
void retrieve_census_path(arma::mat& cencusmat, Rcpp::NumericMatrix& path, Rcpp::NumericVector& census_times, Rcpp::IntegerVector& census_columns);

// update the incidence in an existing census matrix
void compute_incidence(arma::mat& censusmat, arma::uvec col_inds);

// find the intervals for a vector
Rcpp::IntegerVector find_interval(Rcpp::NumericVector& x, Rcpp::NumericVector& breaks, bool rightmost_closed, bool all_inside);

// multivariate normal functions
arma::mat mvrn(int n, const arma::vec& mu, const arma::mat& sigma);

// copy functions to and from the lna ode statevec and the lna drift, residual, and diffusion process objects
void procs2vec(Rcpp::NumericVector& statevec, arma::vec& drift_proc, arma::vec& resid_proc, arma::mat& diffusion_proc);
void vec2procs(Rcpp::NumericVector& statevec, arma::vec& drift_proc, arma::vec& resid_proc, arma::mat& diffusion_proc);


// // compute and return the hazards
// arma::vec COMPUTE_HAZARD(double t, const arma::vec& state, const Rcpp::NumericVector& parameters, const Rcpp::NumericVector& constants, const Rcpp::NumericVector& tcovar, SEXP haz_ptr);
//
// // compute and return the Jacobian matrix
// arma::mat COMPUTE_JACOBIAN(double t, const arma::vec& state, const Rcpp::NumericVector& parameters, const Rcpp::NumericVector& constants, const Rcpp::NumericVector& tcovar, SEXP jacob_ptr);

#endif // stemr_UTILITIES_H