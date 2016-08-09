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
//' @param n_transitions number of possible transitions in the model
//' @param lna_pointer external pointer to LNA integration fcn
//'
//' @return array of paths simulated using the LNA.
//' @export
// [[Rcpp::export]]
arma::cube simulate_lna_paths(int nsim, const arma::colvec& lna_times, const arma::colvec& census_times, const Rcpp::LogicalVector& census_inds, const Rcpp::LogicalVector& restart_inds, const Rcpp::NumericMatrix& lna_pars, const arma::mat& init_states, const arma::uvec& incidence_codes, const int n_transitions, SEXP lna_pointer) {

        // get the dimensions of various objects
        int n_comps = init_states.n_cols;               // number of model compartments
        int n_odes  = n_comps + n_comps*n_transitions;  // number of ODEs = number of compartments + size of the Jacobian
        int n_times = lna_times.n_elem;                 // number of times at which the LNA must be evaluated
        int n_census_times = census_times.n_elem;       // number of census times

        // initialize the interval endpoints
        double t_L = 0;
        double t_R = 0;

        // create an armadillo cube to store the simulations
        arma::cube lna_paths(n_census_times, n_comps + 1, nsim);

        // initialize the LNA objects - the vector for storing the ODES, the state vector, and the Jacobian
        Rcpp::NumericVector x(n_comps * 2 + n_comps^2);                   // vector to store the results of the ODEs
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
                drift_process = init_states.row(k);                                // assign the initial state to the drift vector
                residual_process.zeros();                                          // zero out the residual vector
                diffusion_process.zeros();                                         // zero out the diffusion matrix
                procs2vec(x, drift_process, residual_process, diffusion_process);

                // iterate over the time sequence, solving the LNA over each interval
                for(int j=1; j < n_times; ++j) {

                        // set the times of the interval endpoints
                        t_L = lna_times[j-1];
                        t_R = lna_times[j];

                        // set the parameter values


                }

        }

        // return the paths
        return lna_paths;
}