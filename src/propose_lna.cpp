// [[Rcpp::depends(RcppArmadillo)]]
#include "stemr_types.h"
#include "stemr_utils.h"

using namespace Rcpp;
using namespace arma;

//' Simulate an LNA path using a non-centered parameterization for the
//' log-transformed counting process LNA.
//'
//' @param lna_times vector of interval endpoint times
//' @param lna_pars numeric matrix of parameters, constants, and time-varying
//'   covariates at each of the lna_times
//' @param init_start index in the parameter vector where the initial compartment
//'   volumes start
//' @param param_update_inds logical vector indicating at which of the times the
//'   LNA parameters need to be updated.
//' @param stoich_matrix stoichiometry matrix giving the changes to compartments
//'   from each reaction
//' @param forcing_inds logical vector of indicating at which times in the
//'   time-varying covariance matrix a forcing is applied.
//' @param forcing_matrix matrix containing the forcings.
//' @param step_size initial step size for the ODE solver (adapted internally,
//' but too large of an initial step can lead to failure in stiff systems).
//' @param reject_negatives logical for whether negative increments or compartment 
//' volumes should lead to rejection of the sampled path. If false, draws that lead 
//' to either are re-sampled in place instead of throwing an error. N.B. Resampling 
//' targets the WRONG distribution and should only be used when initializing an LNA 
//' path with the understanding that further warm-up under the correct distribution 
//' is required.
//' @param lna_pointer external pointer to the compiled LNA integration function.
//' @param set_pars_pointer external pointer to the function for setting LNA pars.
//' @return list containing the stochastic perturbations (i.i.d. N(0,1) draws) and
//' the LNA path on its natural scale which is determined by the perturbations.
//'
//' @export
// [[Rcpp::export]]
Rcpp::List propose_lna(const arma::rowvec& lna_times,
                       const Rcpp::NumericMatrix& lna_pars,
                       const int init_start,
                       const Rcpp::LogicalVector& param_update_inds,
                       const arma::mat& stoich_matrix,
                       const Rcpp::LogicalVector& forcing_inds,
                       const arma::mat& forcing_matrix,
                       bool reject_negatives,
                       double step_size,
                       SEXP lna_pointer,
                       SEXP set_pars_pointer) {
      
        // get the dimensions of various objects
        int n_events = stoich_matrix.n_cols;         // number of transition events, e.g., S2I, I2R
        int n_comps  = stoich_matrix.n_rows;         // number of model compartments (all strata)
        int n_odes   = n_events + n_events*n_events; // number of ODEs
        int n_times  = lna_times.n_elem;             // number of times at which the LNA must be evaluated
        int init_end = init_start + n_comps;         // index in the parameter vector to stop copying

        // initialize the objects used in each time interval
        double t_L = 0;
        double t_R = 0;
        Rcpp::NumericVector current_params = lna_pars.row(0);   // vector for storing the current parameter values
        CALL_SET_ODE_PARAMS(current_params, set_pars_pointer);  // set the parameters in the odeintr namespace

        // initial state vector - copy elements from the current parameter vector
        arma::vec init_volumes(current_params.begin() + init_start, n_comps);
        arma::vec init_volumes_prop(init_volumes.begin(), n_comps);

        // initialize the LNA objects - the vector for storing the current state
        Rcpp::NumericVector lna_state_vec(n_odes);   // vector to store the results of the ODEs

        arma::vec lna_drift(n_events, arma::fill::zeros); // incidence mean vector (natural scale)
        arma::mat lna_diffusion(n_events, n_events, arma::fill::zeros); // diffusion matrix

        arma::vec log_lna(n_events, arma::fill::zeros);  // LNA increment, log scale
        arma::vec nat_lna(n_events, arma::fill::zeros);  // LNA increment, natural scale

        // Objects for computing the eigen decomposition of the LNA covariance matrices
        arma::vec svd_d(n_events, arma::fill::zeros);
        arma::mat svd_U(n_events, n_events, arma::fill::zeros);
        arma::mat svd_V(n_events, n_events, arma::fill::zeros);
        bool good_svd = true;

        // matrix in which to store the LNA path
        arma::mat lna_path(n_events+1, n_times, arma::fill::zeros); // incidence path
        arma::mat prev_path(n_comps+1, n_times, arma::fill::zeros); // prevalence path (compartment volumes)
        lna_path.row(0) = lna_times;
        prev_path.row(0) = lna_times;
        prev_path(arma::span(1,n_comps), 0) = init_volumes;

        // apply forcings if called for - applied after censusing at the first time
        if(forcing_inds[0]) {
                init_volumes      += forcing_matrix.col(0);
                init_volumes_prop += forcing_matrix.col(0);
        }

        // indices at which the diffusion elements of lna_state vec start
        int diff_start = n_events;

        // sample the stochastic perturbations
        arma::mat draws(n_events, n_times-1, arma::fill::randn);

        // array for lna moments
        // arma::mat drift_vecs(n_events, n_times-1, arma::fill::zeros);
        // arma::cube diff_mats(n_events, n_events, n_times - 1, arma::fill::zeros);

        // iterate over the time sequence, solving the LNA over each interval
        for(int j=0; j < (n_times-1); ++j) {

                // set the times of the interval endpoints
                t_L = lna_times[j];
                t_R = lna_times[j+1];

                // Reset the LNA state vector and integrate the LNA ODEs over the next interval to 0
                std::fill(lna_state_vec.begin(), lna_state_vec.end(), 0.0);
                CALL_INTEGRATE_STEM_ODE(lna_state_vec, t_L, t_R, step_size, lna_pointer);

                // transfer the elements of the lna_state_vec to the process objects
                std::copy(lna_state_vec.begin(), lna_state_vec.begin() + n_events, lna_drift.begin());
                std::copy(lna_state_vec.begin() + n_events, lna_state_vec.end(), lna_diffusion.begin());

                // ensure symmetry of the diffusion matrix
                // lna_diffusion = arma::symmatu(lna_diffusion);

                // cache the moments
                // drift_vecs.col(j) = lna_drift;
                // diff_mats.slice(j) = lna_diffusion;

                // map the stochastic perturbation to the LNA path on its natural scale
                try{
                        if(lna_drift.has_nan() || lna_diffusion.has_nan()) {
                                good_svd = false;
                                throw std::runtime_error("Integration failed.");
                        } else {
                                good_svd = arma::svd(svd_U, svd_d, svd_V, lna_diffusion); // compute the SVD
                        }

                        if(!good_svd) {
                                throw std::runtime_error("SVD failed.");

                        } else {
                                svd_d.elem(arma::find(svd_d < 0)).zeros();          // zero out negative sing. vals
                                svd_V.each_row() %= arma::sqrt(svd_d).t();          // multiply rows of V by sqrt(d)
                                svd_U *= svd_V.t();                                 // complete svd_sqrt
                                svd_U.elem(arma::find(lna_diffusion == 0)).zeros(); // zero out numerical errors

                                log_lna = lna_drift + svd_U * draws.col(j);         // map the LNA draws
                        }

                } catch(std::exception & err) {

                        // reinstatiate the SVD objects
                        arma::vec svd_d(n_events, arma::fill::zeros);
                        arma::mat svd_U(n_events, n_events, arma::fill::zeros);
                        arma::mat svd_V(n_events, n_events, arma::fill::zeros);

                        // forward the exception
                        forward_exception_to_r(err);

                } catch(...) {
                        ::Rf_error("c++ exception (unknown reason)");
                }
                
                // If simulating, negative increments or volumes are rejected. 
                if(reject_negatives) {
                      
                      // compute the LNA increment
                      nat_lna = arma::exp(log_lna) - 1;
                      
                      // update the compartment volumes
                      init_volumes += stoich_matrix * nat_lna;
                      
                      // throw errors for negative increments or negative volumes
                      try{
                            if(any(nat_lna < 0)) {
                                  throw std::runtime_error("Negative increment.");
                            }
                            
                            if(any(init_volumes < 0)) {
                                  throw std::runtime_error("Negative compartment volumes.");
                            }
                            
                      } catch(std::exception &err) {
                            
                            forward_exception_to_r(err);
                            
                      } catch(...) {
                            ::Rf_error("c++ exception (unknown reason)");
                      }
                      
                      // save the increment and the prevalence
                      lna_path(arma::span(1,n_events), j+1) = nat_lna;
                      prev_path(arma::span(1,n_comps), j+1) = init_volumes;
                      
                      // apply forcings if called for - applied after censusing the path
                      if(forcing_inds[j+1]) {
                            init_volumes += forcing_matrix.col(j+1);
                            
                            // throw errors for negative increments or negative volumes
                            try{
                                  if(any(init_volumes < 0)) {
                                        throw std::runtime_error("Negative compartment volumes.");
                                  }
                                  
                            } catch(std::exception &err) {
                                  
                                  forward_exception_to_r(err);
                                  
                            } catch(...) {
                                  ::Rf_error("c++ exception (unknown reason)");
                            }
                      }
                      
                } else { 
                      // If Initializing, draws leading to negative compartments or volumes are resampled.
                      
                      // compute the LNA increment
                      nat_lna = arma::exp(log_lna) - 1;
                      
                      // update the compartment volumes
                      init_volumes_prop = init_volumes + stoich_matrix * nat_lna;
                      
                      // ensure monotonicity of the increment, compartment volumes, 
                      // and compartment volumes with forcings
                      while(any(nat_lna < 0) || any(init_volumes_prop < 0)) {
                            draws.col(j) = arma::randn(n_events);                          // draw a new vector of N(0,1)
                            log_lna      = lna_drift + svd_U * draws.col(j);               // map the new draws to
                            nat_lna      = arma::exp(log_lna) - 1;                         // compute the LNA increment
                            init_volumes_prop = init_volumes + stoich_matrix * nat_lna;    // compute new initial volumes
                      }
                      
                      // set the initial volumes to the proposed values
                      init_volumes = init_volumes_prop;
                      
                      // save the increment and the prevalence
                      lna_path(arma::span(1,n_events), j+1) = nat_lna;
                      prev_path(arma::span(1,n_comps), j+1) = init_volumes;
                      
                      // apply forcings if called for - applied after censusing the path
                      if(forcing_inds[j+1]) {
                            init_volumes += forcing_matrix.col(j+1);
                            
                            // throw errors for negative increments or negative volumes
                            try{
                                  if(any(init_volumes < 0)) {
                                        throw std::runtime_error("Negative compartment volumes after forcings.");
                                  }
                                  
                            } catch(std::exception &err) {
                                  forward_exception_to_r(err);
                                  
                            } catch(...) {
                                  ::Rf_error("c++ exception (unknown reason)");
                            }
                      }
                }

                // update the parameters if they need to be updated
                if(param_update_inds[j+1]) {
                        current_params = lna_pars.row(j+1);
                }

                // copy the compartment volumes to the current parameters
                std::copy(init_volumes.begin(), init_volumes.end(), current_params.begin() + init_start);

                // set the lna parameters and reset the LNA state vector
                CALL_SET_ODE_PARAMS(current_params, set_pars_pointer);
        }

        // return the paths
        return Rcpp::List::create(Rcpp::Named("draws")     = draws,
                                  Rcpp::Named("lna_path")  = lna_path.t(),
                                  Rcpp::Named("prev_path") = prev_path.t());

                                  // Rcpp::Named("drift_vecs") = drift_vecs,
                                  // Rcpp::Named("diff_mats") = diff_mats);
}