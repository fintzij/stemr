#ifndef stemr_UTILITIES_H
#define stemr_UTILITIES_H

#include <RcppArmadillo.h>
#include <algorithm>

using namespace Rcpp;
using namespace arma;

// call rate functions
void CALL_RATE_FCN(Rcpp::NumericVector& rates,
                   const Rcpp::LogicalVector& inds,
                   const arma::rowvec& state,
                   const Rcpp::NumericVector& parameters,
                   const Rcpp::NumericVector& constants,
                   const arma::rowvec& tcovar,
                   SEXP rate_ptr);

// evaluate an emission density and update an existing emission matrix
void CALL_D_MEASURE(Rcpp::NumericMatrix& emitmat,
                    const Rcpp::LogicalVector& emit_inds,
                    const int record_ind,
                    const Rcpp::NumericVector& record,
                    const Rcpp::NumericVector& state,
                    const Rcpp::NumericVector& parameters,
                    const Rcpp::NumericVector& constants,
                    const Rcpp::NumericVector& tcovar,
                    SEXP d_meas_ptr);

// simulate from the measurement process and update an existing observation matrix
void CALL_R_MEASURE(Rcpp::NumericMatrix& obsmat,
                    const Rcpp::LogicalVector& emit_inds,
                    const int record_ind,
                    const Rcpp::NumericVector& state,
                    const Rcpp::NumericVector& parameters,
                    const Rcpp::NumericVector& constants,
                    const Rcpp::NumericVector& tcovar,
                    SEXP r_meas_ptr);

// integrate the LNA odes, call via XPtr
void CALL_INTEGRATE_STEM_LNA(Rcpp::NumericVector& init,
                             double start,
                             double end,
                             double step_size,
                             SEXP lna_ode_ptr);

// set the LNA parameters, call via XPtr
void CALL_SET_LNA_PARAMS(Rcpp::NumericVector& p,
                         SEXP set_lna_params_ptr);

// integrate the ODEs, call via XPtr
void CALL_INTEGRATE_STEM_ODE(Rcpp::NumericVector& init,
                             double start,
                             double end,
                             double step_size,
                             SEXP stemr_ode_ptr);

// set the ODE parameters, call via XPtr
void CALL_SET_ODE_PARAMS(Rcpp::NumericVector& p, SEXP set_ode_params_ptr);

// update rates based on transition events or changes in time-varying covariates
void rate_update_tcovar(Rcpp::LogicalVector& rate_inds,
                        const arma::mat& M,
                        const arma::rowvec I);

void rate_update_event(Rcpp::LogicalVector& rate_inds,
                       const Rcpp::LogicalMatrix& M,
                       int event_code);

// gillespie simulation
arma::mat simulate_gillespie(const arma::mat& flow,
                             const Rcpp::NumericVector& parameters,
                             const Rcpp::NumericVector& constants,
                             const arma::mat& tcovar,
                             const arma::rowvec& init_states,
                             const Rcpp::LogicalMatrix& rate_adjmat,
                             const arma::mat& tcovar_adjmat,
                             const arma::mat& tcovar_changemat,
                             const Rcpp::IntegerVector init_dims,
                             const arma::vec& forcing_inds,
                             const arma::mat& forcing_matrix,
                             SEXP rate_ptr);

// simulation from the measurement process
Rcpp::NumericMatrix simulate_r_measure(const Rcpp::NumericMatrix& censusmat,
                                       const Rcpp::LogicalMatrix& measproc_indmat,
                                       const Rcpp::NumericVector& parameters,
                                       const Rcpp::NumericVector& constants,
                                       const arma::mat& tcovar,
                                       SEXP r_measure_ptr);

// build a census matrix with compartment counts at observation times
arma::mat build_census_path(Rcpp::NumericMatrix& path,
                            Rcpp::NumericVector& census_times,
                            Rcpp::IntegerVector& census_columns);

// census the lna path matrix, possibly computing prevalence and filling out cumulative incidence
void census_lna(const arma::mat& path,
                arma::mat& census_path,
                const arma::uvec& census_inds,
                const arma::uvec& lna_event_inds,
                const arma::mat& flow_matrix_lna,
                bool do_prevalence,
                const arma::rowvec& init_state,
                const arma::mat& forcing_matrix);

// update a census matrix with compartment counts at observation times
void retrieve_census_path(arma::mat& cencusmat,
                          Rcpp::NumericMatrix& path,
                          Rcpp::NumericVector& census_times,
                          Rcpp::IntegerVector& census_columns);

// update the incidence in an existing census matrix
void compute_incidence(arma::mat& censusmat,
                       arma::uvec& col_inds,
                       Rcpp::List& row_inds);

// find the intervals for a vector
Rcpp::IntegerVector find_interval(Rcpp::NumericVector& x,
                                  Rcpp::NumericVector& breaks,
                                  bool rightmost_closed,
                                  bool all_inside);

// multivariate normal functions
arma::mat rmvtn(int n,
                const arma::rowvec& mu,
                const arma::mat& sigma);

arma::vec dmvtn(const arma::mat& x,
                const arma::rowvec& mu,
                const arma::mat& sigma,
                bool logd = false);

// MCMC transition kernel functions
void c_rw(arma::rowvec& params_prop,
          const arma::rowvec& params_cur,
          int ind,
          const arma::vec& kernel_cov);

void c_rw_adaptive(arma::rowvec& params_prop,
                   const arma::rowvec& params_cur,
                   int ind,
                   const arma::vec& kernel_cov,
                   const arma::vec& proposal_scaling,
                   const arma::vec& nugget);

void mvn_rw(arma::rowvec& params_prop,
            const arma::rowvec& params_cur,
            const arma::mat& sigma_chol);

void mvn_g_adaptive(arma::rowvec& params_prop,
                    const arma::rowvec& params_cur,
                    const arma::mat& kernel_cov,
                    double proposal_scaling,
                    double nugget);

void mvn_c_adaptive(arma::rowvec& params_prop,
                    const arma::rowvec& params_cur,
                    const arma::mat& kernel_cov,
                    const arma::vec& proposal_scaling,
                    arma::mat& sqrt_scalemat,
                    double nugget);

// copy functions
void copy_elem(arma::rowvec& dest, const arma::rowvec& orig, int ind);
void copy_col(arma::mat& dest, const arma::mat& orig, int ind);
void copy_vec(arma::rowvec& dest, const arma::rowvec& orig);
void copy_mat(arma::mat& dest, const arma::mat& orig);
void pars2lnapars(arma::mat& lnapars, const arma::rowvec& parameters);
void g_prop2c_prop(arma::mat& g2c_mat, const arma::rowvec& params_cur, const arma::rowvec& params_prop);

// convert the lna from the counting process on transition events to its natural state space
// arma::mat convert_lna(const arma::mat& path, const arma::mat& flow_matrix, const arma::rowvec& init_state);
// void convert_lna2(const arma::mat& path, const arma::mat& flow_matrix, const arma::rowvec& init_state, arma::mat& statemat);

// compute the lna density
// retintegrating all LNA ODEs - required after updating parameters
// Rcpp::List lna_density(const Rcpp::List& path, const arma::colvec& lna_times, const Rcpp::NumericMatrix& lna_pars,
                       // const Rcpp::LogicalVector& param_update_inds, const arma::mat& flow_matrix,
                       // SEXP lna_pointer, SEXP set_pars_pointer);

// reintegrating just the drift and residual ODEs - sufficient after elliptical slice sampling
// Rcpp::List lna_density2(const Rcpp::List& path, const arma::colvec& lna_times, const Rcpp::NumericMatrix& lna_pars,
                        // const Rcpp::LogicalVector& param_update_inds, const arma::mat& flow_matrix,
                        // SEXP lna_pointer, SEXP set_pars_pointer);

#endif // stemr_UTILITIES_H