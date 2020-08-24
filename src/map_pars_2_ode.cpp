// [[Rcpp::depends(RcppArmadillo)]]
#include "stemr_types.h"
#include "stemr_utils.h"

using namespace Rcpp;
using namespace arma;

//' Map parameters to the deterministic mean incidence increments for a stochastic
//' epidemic model.
//'
//' @param pathmat matrix where the ODE path should be stored
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
void map_pars_2_ode(arma::mat& pathmat,
                    const arma::rowvec& ode_times,
                    const Rcpp::NumericMatrix& ode_pars,
                    Rcpp::NumericVector& ode_param_vec,
                    const Rcpp::IntegerVector& ode_param_inds,
                    const Rcpp::IntegerVector& ode_tcovar_inds,
                    const int init_start,
                    const Rcpp::LogicalVector& param_update_inds,
                    const arma::mat& stoich_matrix,
                    const Rcpp::LogicalVector& forcing_inds,
                    const arma::uvec& forcing_tcov_inds,
                    const arma::mat& forcings_out,
                    const arma::cube& forcing_transfers,
                    double step_size,
                    SEXP ode_pointer,
                    SEXP set_pars_pointer) {

        // get the dimensions of various objects
        int n_events = stoich_matrix.n_cols;         // number of transition events, e.g., S2I, I2R
        int n_comps  = stoich_matrix.n_rows;         // number of model compartments (all strata)
        int n_times  = ode_times.n_elem;             // number of times at which the ODEs must be evaluated
        int n_tcovar = ode_tcovar_inds.size();       // number of time-varying covariates or parameters
        int n_forcings = forcing_tcov_inds.n_elem;   // number of forcings
        
        // for use with forcings
        double forcing_flow = 0;
        arma::vec forcing_distvec(n_comps, arma::fill::zeros);

        // initialize the objects used in each time interval
        double t_L = 0;
        double t_R = 0;

        // vector of parameters, initial compartment columes, constants, and time-varying covariates
        std::copy(ode_pars.row(0).begin(), ode_pars.row(0).end(), ode_param_vec.begin());
        
        // initial state vector - copy elements from the current parameter vector
        arma::vec init_volumes(ode_param_vec.begin() + init_start, n_comps);
        
        // set the parameters in the odeintr namespace
        CALL_SET_ODE_PARAMS(ode_param_vec, set_pars_pointer); 

        // initialize the ODE objects - the vector for storing the current state
        Rcpp::NumericVector ode_state_vec(n_events);   // vector to store the ODEs

        // apply forcings if called for - applied after censusing at the first time
        if(forcing_inds[0]) {
              
              // distribute the forcings proportionally to the compartment counts in the applicable states
              for(int j=0; j < n_forcings; ++j) {
                    
                    forcing_flow       = ode_pars(0, forcing_tcov_inds[j]);
                    forcing_distvec    = forcing_flow * normalise(forcings_out.col(j) % init_volumes, 1);
                    init_volumes      += forcing_transfers.slice(j) * forcing_distvec;
              }
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
                pathmat(j+1, arma::span(1, n_events)) = Rcpp::as<arma::rowvec>(ode_state_vec);

                // apply forcings if called for - applied after censusing the path
                if(forcing_inds[j+1]) {
                      
                      // distribute the forcings proportionally to the compartment counts in the applicable states
                      for(int s=0; s < n_forcings; ++s) {
                            
                            forcing_flow       = ode_pars(j+1, forcing_tcov_inds[s]);
                            forcing_distvec    = forcing_flow * normalise(forcings_out.col(s) % init_volumes, 1);
                            init_volumes      += forcing_transfers.slice(s) * forcing_distvec;
                      }
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
                                ode_param_vec.end() - n_tcovar);
                      
                }

                // copy the compartment volumes to the current parameters
                std::copy(init_volumes.begin(), init_volumes.end(), ode_param_vec.begin() + init_start);

                // set the ODE parameters and reset the ODE state vector
                CALL_SET_ODE_PARAMS(ode_param_vec, set_pars_pointer);
        }
}