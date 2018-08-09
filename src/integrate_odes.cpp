// [[Rcpp::depends(RcppArmadillo)]]
#include "stemr_types.h"
#include "stemr_utils.h"

using namespace Rcpp;
using namespace arma;

//' Obtain the path of the deterministic mean of a stochastic epidemic model by
//' integrating the corresponding ODE functions.
//'
//' @param ode_times vector of interval endpoint times
//' @param ode_pars numeric matrix of parameters, constants, and time-varying
//'   covariates at each of the ode_times
//' @param init_start index in the parameter vector where the initial compartment
//'   volumes start
//' @param param_update_inds logical vector indicating at which of the times the
//'   ode parameters need to be updated.
//' @param stoich_matrix stoichiometry matrix giving the changes to compartments
//'   from each reaction
//' @param forcing_inds logical vector of indicating at which times in the
//'   time-varying covariance matrix a forcing is applied.
//' @param forcing_matrix matrix containing the forcings.
//' @param step_size initial step size for the ODE solver (adapted internally,
//' but too large of an initial step can lead to failure in stiff systems).
//' @param ode_pointer external pointer to ode integration function.
//' @param set_pars_pointer external pointer to the function for setting the ode
//'   parameters.
//'
//' @return List containing the ODE incidence and prevalence paths.
//'
//' @export
// [[Rcpp::export]]
Rcpp::List integrate_odes(const arma::rowvec& ode_times,
                       const Rcpp::NumericMatrix& ode_pars,
                       const Rcpp::IntegerVector& ode_param_inds,
                       const Rcpp::IntegerVector& ode_tcovar_inds,
                       const int init_start,
                       const Rcpp::LogicalVector& param_update_inds,
                       const arma::mat& stoich_matrix,
                       const Rcpp::LogicalVector& forcing_inds,
                       const arma::mat& forcing_matrix,
                       double step_size,
                       SEXP ode_pointer,
                       SEXP set_pars_pointer) {

        // get the dimensions of various objects
        int n_events = stoich_matrix.n_cols;         // number of transition events, e.g., S2I, I2R
        int n_comps  = stoich_matrix.n_rows;         // number of model compartments (all strata)
        int n_times  = ode_times.n_elem;             // number of times at which the ODEs must be evaluated
        int n_tcovar = ode_tcovar_inds.size();   // number of time-varying covariates or parameters

        // initialize the objects used in each time interval
        double t_L = 0;
        double t_R = 0;
        Rcpp::NumericVector current_params = ode_pars.row(0);   // vector for storing the current parameter values
        CALL_SET_ODE_PARAMS(current_params, set_pars_pointer);  // set the parameters in the odeintr namespace

        // initial state vector - copy elements from the current parameter vector
        arma::vec init_volumes(current_params.begin() + init_start, n_comps);

        // initialize the ODE objects - the vector for storing the current state
        Rcpp::NumericVector ode_state_vec(n_events);   // vector to store the ODEs

        // matrix in which to store the ODE path
        arma::mat incid_path(n_events+1, n_times, arma::fill::zeros); // incidence path
        arma::mat prev_path(n_comps+1, n_times, arma::fill::zeros);   // prevalence path (compartment volumes)
        incid_path.row(0) = ode_times;
        prev_path.row(0)  = ode_times;
        prev_path(arma::span(1,n_comps), 0) = init_volumes;

        // apply forcings if called for - applied after censusing at the first time
        if(forcing_inds[0]) {
                init_volumes += forcing_matrix.col(0);
        }

        // iterate over the time sequence, solving the ODEs over each interval
        for(int j=0; j < (n_times-1); ++j) {

                // set the times of the interval endpoints
                t_L = ode_times[j];
                t_R = ode_times[j+1];

                // Reset the ODE state vector and integrate the ODEs over the next interval
                std::fill(ode_state_vec.begin(), ode_state_vec.end(), 0.0);
                CALL_INTEGRATE_STEM_ODE(ode_state_vec, t_L, t_R, step_size, ode_pointer);

                // compute the compartment volumes
                init_volumes += stoich_matrix * Rcpp::as<arma::vec>(ode_state_vec);

                // Save the increment and volumes
                incid_path(arma::span(1, n_events), j+1) = Rcpp::as<arma::vec>(ode_state_vec);
                prev_path(arma::span(1, n_comps), j+1)  = init_volumes;

                // apply forcings if called for - applied after censusing the path
                if(forcing_inds[j+1]) {
                        init_volumes += forcing_matrix.col(j+1);
                }

                // ensure the initial volumes are non-negative
                try{
                        if(any(init_volumes < 0)) {
                                throw std::runtime_error("Negative compartment volumes.");
                        }

                } catch(std::runtime_error &err) {
                        forward_exception_to_r(err);

                } catch(...) {
                        ::Rf_error("c++ exception (unknown reason)");
                }

                // update the parameters if they need to be updated
                if(param_update_inds[j+1]) {
                      
                      // time-varying covariates and parameters
                      std::copy(ode_pars.row(j+1).end() - n_tcovar,
                                ode_pars.row(j+1).end(),
                                current_params .end() - n_tcovar);
                      
                }

                // copy the compartment volumes to the current parameters
                std::copy(init_volumes.begin(), init_volumes.end(), current_params.begin() + init_start);

                // set the ODE parameters and reset the ODE state vector
                CALL_SET_ODE_PARAMS(current_params, set_pars_pointer);
        }

        // return the paths
        return Rcpp::List::create(Rcpp::Named("incid_path") = incid_path.t(),
                                  Rcpp::Named("prev_path")  = prev_path.t());
}