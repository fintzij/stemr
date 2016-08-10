// [[Rcpp::depends(RcppArmadillo)]]
#include "stemr_types.h"
#include "stemr_utils.hpp"

using namespace Rcpp;
using namespace arma;

//' Simulate a data matrix from the measurement process of a stochastic epidemic
//' model.
//'
//' @param nsim integer for the number of simulations
//' @param lna_times vector of times at which the LNA must be evaluated
//' @param census_times vector of times at which the LNA should be censused
//' @param census_inds logical vector the same length as lna_times indicating
//'   whether the LNA path is censused at each of the lna times
//' @param restart_inds logical vector the same length as lna_times indicating
//'   whether the LNA is to be restarted at each of the lna times
//' @param lna_pars numeric matrix of parameters, constants, and time-varying
//'   covariates at each of the lna_times
//' @param init_states matrix of initial state vectors
//' @param incidence_codes vector of incidence codes, or 0 if no incidence
//'   compartments
//' @param lna_pointer external pointer to LNA integration fcn
//' @param set_pars_pointer external pointer to the function for setting the LNA
//'   parameters.
//'
//' @return array of paths simulated using the LNA.
//' @export
// [[Rcpp::export]]
arma::cube simulate_lna_paths(int nsim, const arma::colvec& lna_times, const arma::colvec& census_times, const Rcpp::LogicalVector& census_inds, const Rcpp::LogicalVector& restart_inds, const Rcpp::NumericMatrix& lna_pars, const arma::mat& init_states, const Rcpp::IntegerVector& incidence_codes, SEXP lna_pointer, SEXP set_pars_pointer) {

        // get the dimensions of various objects
        int n_comps = init_states.n_cols;         // number of model compartments
        int n_odes  = 2*n_comps + n_comps*n_comps;// number of ODEs
        int n_times = lna_times.n_elem;           // number of times at which the LNA must be evaluated
        int n_incidence = incidence_codes.size(); // number of incidence codes
        int n_census_times = census_times.n_elem; // number of census times

        // initialize the objects used in each time interval
        double t_L = 0;
        double t_R = 0;
        Rcpp::NumericVector current_params(lna_pars.ncol()); // vector for storing the current parameter values
        arma::mat lna_step(1, n_comps);                      // vector for storing the sampled lna values

        // create an armadillo cube to store the simulations
        arma::cube lna_paths(n_census_times, n_comps + 1, nsim);

        // initialize the LNA objects - the vector for storing the ODES, the state vector, and the Jacobian
        Rcpp::NumericVector lna_state_vec(n_odes);                        // vector to store the results of the ODEs
        arma::vec drift_process(n_comps, arma::fill::zeros);              // vector for the deterministic drift process
        arma::vec residual_process(n_comps,arma::fill::zeros);            // vector for the residual process
        arma::mat diffusion_process(n_comps, n_comps, arma::fill::zeros); // matrix for the diffusion

        // indices in the ODE vector for the drift, residual, and diffusion
        Rcpp::IntegerVector drift_inds = Rcpp::seq_len(n_comps) - 1;
        Rcpp::IntegerVector resid_inds = Rcpp::seq_len(n_comps) + n_comps - 1;
        Rcpp::IntegerVector diff_inds  = Rcpp::seq_len(n_comps * n_comps) + 2*n_comps - 1;

        // start simulating the LNA paths
        for(int k = 0; k < nsim; ++k) {

                // set the initial values of the LNA
                lna_paths.slice(k)(0, arma::span(1,n_comps)) = init_states.row(k); // assign the initial state to the path matrix
                lna_paths.slice(k).col(0) = census_times;                          // assign the census times to the current slice of the cube
                drift_process = init_states.row(k).t();                                // assign the initial state to the drift vector
                residual_process.zeros();                                          // zero out the residual vector
                diffusion_process.zeros();                                         // zero out the diffusion matrix
                procs2vec(lna_state_vec, drift_process, residual_process, diffusion_process);

                // iterate over the time sequence, solving the LNA over each interval
                for(int j=1; j < n_times; ++j) {

                        // set the times of the interval endpoints
                        t_L = lna_times[j-1];
                        t_R = lna_times[j];

                        // set the values of the
                        current_params = lna_pars.row(j-1);
                        CALL_SET_LNA_PARAMS(current_params, set_pars_pointer);

                        // integrate the LNA
                        // lna_state_vec = CALL_INTEGRATE_STEM_LNA(lna_state_vec, t_L, t_R, 1.0, lna_pointer);
                        CALL_INTEGRATE_STEM_LNA(lna_state_vec, t_L, t_R, 1.0, lna_pointer);

                        // transfer the elements of the lna_state_vec to the process objects
                        vec2procs(lna_state_vec, drift_process, residual_process, diffusion_process);

                        // sample the next value
                        lna_step = mvrn(1, drift_process + residual_process, diffusion_process);

                        // if any values are negative, resample them
                        while(arma::any(lna_step.row(0) < 0)) {
                                lna_step = mvrn(1, drift_process + residual_process, diffusion_process);
                        }

                        // insert the sampled value into the path if called for
                        if(census_inds[j]) {
                                lna_paths.slice(k)(j, arma::span(1, n_comps)) = lna_step;
                        }

                        // restart the LNA if called for
                        if(restart_inds[j]) {
                                drift_process = lna_step.row(0).t();    // restart the deterministic process at the sampled state
                                residual_process.zeros();               // set the residual process to zero
                                diffusion_process.zeros();              // the diffusion process is set to zero
                        } else {
                                residual_process = drift_process - lna_step.row(0).t(); // the residual process is the difference between the deterministic
                                                                                        // process and the sampled state
                                diffusion_process.zeros();                              // the diffusion process is set to zero
                        }

                        // copy the process objects back into the lna state vector
                        procs2vec(lna_state_vec, drift_process, residual_process, diffusion_process);
                }
        }

        // compute the incidence if necessary
        if(incidence_codes[0] != 0) {

                int incid_begin = incidence_codes[0];
                int incid_end   = incidence_codes[n_incidence -1];
                double maxval   = lna_paths.max();

                // iterate through the path array
                for(int k=0; k < nsim; ++k) {
                        lna_paths.subcube(arma::span(1,n_census_times-1), arma::span(incid_begin, incid_end), arma::span(k)) =
                                arma::diff(lna_paths.slice(k).cols(incid_begin, incid_end));

                        // clamp the slice to make sure there are no negative values
                        lna_paths.slice(k) = arma::clamp(lna_paths.slice(k), 0, maxval);
                }
        }

        // return the paths
        return lna_paths;
}