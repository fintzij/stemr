// [[Rcpp::depends(RcppArmadillo)]]
#include "stemr_types.h"
#include "stemr_utils.hpp"

using namespace Rcpp;
using namespace arma;

//' Simulate a data matrix from the measurement process of a stochastic epidemic
//' model.
//'
//' @param ode_times vector of times at which the ODE must be evaluated
//' @param census_times vector of times at which the ODE should be censused
//' @param census_inds logical vector the same length as ode_times indicating
//'   whether the ODE path is censused at each of the ODE times
//' @param ode_pars numeric matrix of parameters, constants, and time-varying
//'   covariates at each of the ode_times
//' @param init_states matrix of initial state vectors
//' @param incidence_codes vector of incidence codes, or 0 if no incidence
//'   compartments
//' @param ode_pointer external pointer to ODE integration fcn
//' @param set_pars_pointer external pointer to the function for setting the ODE
//'   parameters.
//'
//' @return matrix containing the deterministic ODE path
//' @export
// [[Rcpp::export]]
arma::mat stem_ode_path(const arma::colvec& ode_times, const arma::colvec& census_times, const Rcpp::LogicalVector& census_inds, const Rcpp::NumericMatrix& ode_pars, const Rcpp::NumericVector& init_states, const Rcpp::IntegerVector& incidence_codes, SEXP ode_pointer, SEXP set_pars_pointer) {

        // get the dimensions of various objects
        int n_comps = init_states.size();         // number of model compartments
        int n_odes  = n_comps;                    // number of ODEs
        int n_times = ode_times.n_elem;           // number of times at which the ode must be evaluated
        int n_incidence = incidence_codes.size(); // number of incidence codes
        int n_census_times = census_times.n_elem; // number of census times

        // initialize the objects used in each time interval
        double t_L = 0;
        double t_R = 0;
        Rcpp::NumericVector current_params(ode_pars.ncol()); // vector for storing the current parameter values
        // arma::mat ode_step(1, n_comps);                      // vector for storing the sampled ode values

        // create an armadillo cube to store the simulations
        arma::mat ode_path(n_census_times, n_comps + 1);

        // vector to store the results of the ODEs
        Rcpp::NumericVector ode_state_vec(n_odes);
        ode_state_vec = init_states;

        // start simulating the ode path
        // set the initial values of the ode
        ode_path(0, arma::span(1,n_comps)) = Rcpp::as<arma::rowvec>(init_states); // assign the initial state to the path matrix
        ode_path.col(0) = census_times;                          // assign the census times to the current slice of the cube

        // iterate over the time sequence, solving the ode over each interval
        for(int j=1; j < n_times; ++j) {

                // set the times of the interval endpoints
                t_L = ode_times[j-1];
                t_R = ode_times[j];

                // set the values of the
                current_params = ode_pars.row(j-1);
                CALL_SET_ODE_PARAMS(current_params, set_pars_pointer);

                // integrate the ode
                // ode_state_vec = CALL_INTEGRATE_STEM_ode(ode_state_vec, t_L, t_R, 1.0, ode_pointer);
                CALL_INTEGRATE_STEM_ODE(ode_state_vec, t_L, t_R, 1.0, ode_pointer);

                // insert the sampled value into the path if called for
                if(census_inds[j]) {
                        ode_path(j, arma::span(1, n_comps)) = Rcpp::as<arma::rowvec>(ode_state_vec);
                }
        }


        // compute the incidence if necessary
        if(incidence_codes[0] != 0) {

                int incid_begin = incidence_codes[0];
                int incid_end   = incidence_codes[n_incidence -1];
                double maxval   = ode_path.max();

                // difference the incidence compartments
                ode_path(arma::span(1,n_census_times-1), arma::span(incid_begin, incid_end)) = arma::diff(ode_path.cols(incid_begin, incid_end));

                // clamp the slice to make sure there are no negative values
                ode_path = arma::clamp(ode_path, 0, maxval);
        }

        // return the paths
        return ode_path;
}