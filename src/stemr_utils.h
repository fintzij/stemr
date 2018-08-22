#ifndef stemr_UTILITIES_H
#define stemr_UTILITIES_H

#include <RcppArmadillo.h>
#include <algorithm>
#include <math.h>
#include <boost/numeric/odeint.hpp>

using namespace Rcpp;
using namespace arma;
namespace odeint = boost::numeric::odeint;

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
void CALL_INTEGRATE_STEM_ODE(Rcpp::NumericVector& init,
                             double start,
                             double end,
                             double step_size,
                             SEXP stem_ode_ptr);

// set the LNA parameters, call via XPtr
void CALL_SET_ODE_PARAMS(Rcpp::NumericVector& p,
                         SEXP set_ode_params_ptr);

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

void draw_normals(arma::vec& v);
void draw_normals2(arma::mat& M);
void sample_unit_sphere(arma::vec& v);

// MCMC transition kernel functions
void mvn_rw(arma::rowvec& params_prop,
            const arma::rowvec& params_cur,
            const arma::mat& sigma_chol);

void mvn_g_adaptive(arma::rowvec& params_prop,
                    const arma::rowvec& params_cur,
                    const arma::mat& kernel_cov_chol,
                    double nugget);

// copy functions
void add2vec(arma::rowvec& target, const arma::rowvec& increments, const arma::uvec& inds);
void copy_2_rows(arma::mat& dest, const arma::mat& orig, const arma::uvec& inds);
void copy_col(arma::mat& dest, const arma::mat& orig, int ind);
void copy_elem(arma::rowvec& dest, const arma::rowvec& orig, int ind);
void copy_elem2(arma::rowvec& dest, const arma::rowvec& orig, const arma::uvec& inds);
void g_prop2c_prop(arma::mat& g2c_mat, const arma::rowvec& params_cur, const arma::rowvec& params_prop);
void copy_pathmat(arma::mat& dest, const arma::mat& orig);
void copy_vec(arma::rowvec& dest, const arma::rowvec& orig);
void copy_vec2(arma::rowvec& dest, const arma::rowvec& orig, const arma::uvec& inds);
void copy_mat(arma::mat& dest, const arma::mat& orig);
void increment_elem(arma::vec& vec, int ind);
void insert_block(arma::mat& dest, const arma::mat& orig, const arma::uvec& rowinds, const arma::uvec& colinds);
void insert_tparam(arma::mat& tcovar, const arma::vec& values, int col_ind, const arma::uvec& tpar_inds);
void mat_2_arr(arma::cube& dest, const arma::mat& orig, int ind);
void pars2lnapars(arma::mat& lnapars, const arma::rowvec& parameters);
void pars2lnapars2(arma::mat& lnapars, const arma::rowvec& parameters, int c_start);
void reset_vec(arma::vec& v, double value = 0);

// comp_chol
void comp_chol(arma::mat& C, arma::mat& M);

// reset slice ratios
void reset_slice_ratios(arma::vec& n_expansions,
                        arma::vec& n_contractions,
                        arma::vec& n_expansions_c,
                        arma::vec& n_contractions_c,
                        arma::vec& slice_ratios);

// update factors
void update_factors(arma::vec& slice_eigenvals,
                    arma::mat& slice_eigenvecs,
                    const arma::mat& kernel_cov);

void update_interval_widths(arma::vec& interval_widths,
                            arma::vec& n_expansions_afss,
                            arma::vec& n_contractions_afss,
                            const arma::vec& c_expansions_afss,
                            const arma::vec& c_contractions_afss,
                            arma::vec& slice_ratios,
                            double adaptation_factor,
                            double target_ratio);

#endif // stemr_UTILITIES_H