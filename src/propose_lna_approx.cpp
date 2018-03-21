// [[Rcpp::depends(RcppArmadillo)]]
#include "stemr_types.h"
#include "stemr_utils.h"

using namespace Rcpp;
using namespace arma;

//' Simulate an approximate LNA path using a non-centered parameterization for the
//' log-transformed counting process LNA. Resample the initial path in place, then
//' update with elliptical slice sampling.
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
//' @param max_attempts maximum number of tries to restart negative increments if 
//' reject_negatives is false.
//' @param nsim number of paths to simulate
//' @param ess_updates number of elliptical slice sampling updates between samples.
//' @param ess_warmup number of elliptical slice sampling updates before the first
//' sample.
//' @param step_size initial step size for the ODE solver (adapted internally,
//' but too large of an initial step can lead to failure in stiff systems).
//' @param lna_pointer external pointer to the compiled LNA integration function.
//' @param set_pars_pointer external pointer to the function for setting LNA pars.
//' @return list containing the stochastic perturbations (i.i.d. N(0,1) draws) and
//' the LNA path on its natural scale which is determined by the perturbations.
//'
//' @export
// [[Rcpp::export]]
Rcpp::List propose_lna_approx(const arma::rowvec& lna_times,
                              const Rcpp::NumericMatrix& lna_pars,
                              const int init_start,
                              const Rcpp::LogicalVector& param_update_inds,
                              const arma::mat& stoich_matrix,
                              const Rcpp::LogicalVector& forcing_inds,
                              const arma::mat& forcing_matrix,
                              int max_attempts,
                              int nsim,
                              int ess_updates,
                              int ess_warmup,
                              double step_size,
                              SEXP lna_pointer,
                              SEXP set_pars_pointer) {
      
        // get the dimensions of various objects
        int n_events = stoich_matrix.n_cols;         // number of transition events, e.g., S2I, I2R
        int n_comps  = stoich_matrix.n_rows;         // number of model compartments (all strata)
        int n_odes   = n_events + n_events*n_events; // number of ODEs
        int n_times  = lna_times.n_elem;             // number of times at which the LNA must be evaluated

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

        // arrays in which to store the LNA paths
        arma::cube incid_paths(n_times, n_events+1, nsim, arma::fill::zeros); // incidence path
        arma::cube prev_paths(n_times, n_comps+1, nsim, arma::fill::zeros); // prevalence path (compartment volumes)
        arma::cube lna_draws(n_times-1, n_events, nsim, arma::fill::zeros);     // N(0,1) draws that map to the path         
        
        // initialize the objects used in each time interval
        double t_L = 0;
        double t_R = 0;
        Rcpp::NumericVector current_params = lna_pars.row(0);   // vector for storing the current parameter values
        CALL_SET_ODE_PARAMS(current_params, set_pars_pointer);  // set the parameters in the odeintr namespace
        
        // initial state vector - copy elements from the current parameter vector
        arma::vec init_volumes(current_params.begin() + init_start, n_comps);
        arma::vec init_volumes_prop(init_volumes.begin(), n_comps);
        
        // matrices in which to store the current path
        arma::mat lna_path(n_events+1, n_times, arma::fill::zeros);
        arma::mat prev_path(n_comps+1, n_times, arma::fill::zeros);
        
        // fill out times and initial volumes
        lna_path.row(0) = lna_times;
        prev_path.row(0) = lna_times;
        prev_path(arma::span(1,n_comps), 0) = init_volumes;

        // apply forcings if called for - applied after censusing at the first time
        if(forcing_inds[0]) {
                init_volumes      += forcing_matrix.col(0);
                init_volumes_prop += forcing_matrix.col(0);
        }

        // sample the stochastic perturbations
        arma::mat draws_cur(n_events, n_times-1, arma::fill::randn);
        arma::mat draws_prop(n_events, n_times-1, arma::fill::randn);
        arma::mat draws_temp(n_events, n_times-1, arma::fill::randn);
        
        // integer for the attempt number
        int attempt = 0;

        // parameters for the elliptical slice sampler
        double theta = 0;
        double lower = 0; 
        double upper = 1;
        bool valid_path = false;

        // iterate over the time sequence, solving the LNA over each interval
        for(int j = 0; j < (n_times-1); ++j) {

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

                                log_lna = lna_drift + svd_U * draws_cur.col(j);         // map the LNA draws
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
                
                // when initializing, draws leading to negative compartments or volumes are resampled.
                // compute the LNA increment
                nat_lna = arma::exp(log_lna) - 1;
                
                // update the compartment volumes
                init_volumes_prop = init_volumes + stoich_matrix * nat_lna;
                
                // ensure monotonicity of the increment, compartment volumes, 
                // and compartment volumes with forcings
                attempt = 0;
                while((any(nat_lna < 0) || any(init_volumes_prop < 0)) && (attempt <= max_attempts)) {
                      attempt          += 1;
                      draws_cur.col(j)  = arma::randn(n_events);                          // draw a new vector of N(0,1)
                      log_lna           = lna_drift + svd_U * draws_cur.col(j);               // map the new draws to
                      nat_lna           = arma::exp(log_lna) - 1;                         // compute the LNA increment
                      init_volumes_prop = init_volumes + stoich_matrix * nat_lna;    // compute new initial volumes
                }
                
                try{
                      if(attempt > max_attempts) {
                            if(any(nat_lna < 0)) {
                                  throw std::runtime_error("Negative increment.");
                            }
                            
                            if(any(init_volumes < 0)) {
                                  throw std::runtime_error("Negative compartment volumes.");
                            }
                      }
                      
                } catch(std::exception &err) {
                      
                      forward_exception_to_r(err);
                      
                } catch(...) {
                      ::Rf_error("c++ exception (unknown reason)");
                }
                
                // set the initial volumes to the proposed values
                init_volumes = init_volumes_prop;
                
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

                // update the parameters if they need to be updated
                if(param_update_inds[j+1]) {
                        current_params = lna_pars.row(j+1);
                }

                // copy the compartment volumes to the current parameters
                std::copy(init_volumes.begin(), init_volumes.end(), current_params.begin() + init_start);

                // set the lna parameters and reset the LNA state vector
                CALL_SET_ODE_PARAMS(current_params, set_pars_pointer);
        }
        
        // now warm up the sampler
        for(int k = 0; k < ess_warmup; ++k) {
              
              // initialize the parameters and volumes
              current_params = lna_pars.row(0);   // vector for storing the current parameter values
              CALL_SET_ODE_PARAMS(current_params, set_pars_pointer);  // set the parameters in the odeintr namespace
              
              // initial state vector - copy elements from the current parameter vector
              std::copy(current_params.begin()+init_start, current_params.begin()+init_start+n_comps, init_volumes.begin());
              init_volumes_prop = init_volumes;
              
              // apply forcings if called for - applied after censusing at the first time
              if(forcing_inds[0]) {
                    init_volumes      += forcing_matrix.col(0);
                    init_volumes_prop += forcing_matrix.col(0);
              }
              
              // sample new perturbations
              draws_prop.randn();
              
              // center the bracket
              theta = runif(1, 0, 2*arma::datum::pi)[0];
              lower = theta - 2*arma::datum::pi;
              upper = theta;
              
              // construct the first proposal
              draws_temp = cos(theta) * draws_cur + sin(theta) * draws_prop;
              
              // initialize the log-likelihood (indicator for a valid proposal)
              valid_path = true;
              
              // get the LNA path
              for(int j = 0; j < (n_times-1); ++j) {
                    
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
                                
                                log_lna = lna_drift + svd_U * draws_temp.col(j);         // map the LNA draws
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
                    
                    // when initializing, draws leading to negative compartments or volumes are resampled.
                    // compute the LNA increment
                    nat_lna = arma::exp(log_lna) - 1;
                    
                    // update the compartment volumes
                    init_volumes_prop = init_volumes + stoich_matrix * nat_lna;
                    
                    // ensure monotonicity of the increment, compartment volumes, 
                    // and compartment volumes with forcings
                    if(any(nat_lna < 0) || any(init_volumes_prop < 0)) {
                          valid_path = false;
                          break;
                    }
                    
                    // set the initial volumes to the proposed values
                    init_volumes = init_volumes_prop;
                    
                    // apply forcings if called for - applied after censusing the path
                    if(forcing_inds[j+1]) {
                          init_volumes += forcing_matrix.col(j+1);
                          
                          // throw errors for negative increments or negative volumes
                          if(any(init_volumes < 0)) {
                                valid_path = false;
                                break;
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
              
              while((lower != upper) && !valid_path) {
                    
                    // initialize the parameters and volumes
                    current_params = lna_pars.row(0);   // vector for storing the current parameter values
                    CALL_SET_ODE_PARAMS(current_params, set_pars_pointer);  // set the parameters in the odeintr namespace
                    
                    // initial state vector - copy elements from the current parameter vector
                    std::copy(current_params.begin()+init_start, current_params.begin()+init_start+n_comps, init_volumes.begin());
                    init_volumes_prop = init_volumes;
                    
                    // apply forcings if called for - applied after censusing at the first time
                    if(forcing_inds[0]) {
                          init_volumes      += forcing_matrix.col(0);
                          init_volumes_prop += forcing_matrix.col(0);
                    }
                    
                    // construct the next proposal
                    if(theta<0) {
                          lower = theta;
                    } else {
                          upper = theta;
                    }
                    
                    theta = runif(1, lower, upper)[0];
                    draws_temp = cos(theta) * draws_cur + sin(theta) * draws_prop;
                    
                    // initialize the log-likelihood (indicator for a valid proposal)
                    valid_path = true;
                    
                    // apply forcings if called for - applied after censusing at the first time
                    if(forcing_inds[0]) {
                          init_volumes      += forcing_matrix.col(0);
                          init_volumes_prop += forcing_matrix.col(0);
                    }
                    
                    // get the LNA path
                    for(int j = 0; j < (n_times-1); ++j) {
                          
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
                                      
                                      log_lna = lna_drift + svd_U * draws_temp.col(j);         // map the LNA draws
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
                          
                          // when initializing, draws leading to negative compartments or volumes are resampled.
                          // compute the LNA increment
                          nat_lna = arma::exp(log_lna) - 1;
                          
                          // update the compartment volumes
                          init_volumes_prop = init_volumes + stoich_matrix * nat_lna;
                          
                          // ensure monotonicity of the increment, compartment volumes, 
                          // and compartment volumes with forcings
                          attempt = 0;
                          if(any(nat_lna < 0) || any(init_volumes_prop < 0)) {
                                valid_path = false;
                                break;
                          }
                          
                          // set the initial volumes to the proposed values
                          init_volumes = init_volumes_prop;
                          
                          // apply forcings if called for - applied after censusing the path
                          if(forcing_inds[j+1]) {
                                init_volumes += forcing_matrix.col(j+1);
                                
                                // throw errors for negative increments or negative volumes
                                if(any(init_volumes < 0)) {
                                      valid_path = false;
                                      break;
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
              }
              
              if(lower != upper) {
                    // update the draws
                    draws_cur = draws_temp;
              }
        }
        
        // now sample paths
        for(int k = 0; k < nsim; ++k) {
              for(int i = 0; i < ess_updates; ++i) {
                    
                    // initialize the parameters and volumes
                    current_params = lna_pars.row(0);   // vector for storing the current parameter values
                    CALL_SET_ODE_PARAMS(current_params, set_pars_pointer);  // set the parameters in the odeintr namespace
                    
                    // initial state vector - copy elements from the current parameter vector
                    std::copy(current_params.begin()+init_start, current_params.begin()+init_start+n_comps, init_volumes.begin());
                    init_volumes_prop = init_volumes;
                    
                    // apply forcings if called for - applied after censusing at the first time
                    if(forcing_inds[0]) {
                          init_volumes      += forcing_matrix.col(0);
                          init_volumes_prop += forcing_matrix.col(0);
                    }
                    
                    // sample new perturbations
                    draws_prop.randn();
                    
                    // center the bracket
                    theta = runif(1, 0, 2*arma::datum::pi)[0];
                    lower = theta - 2*arma::datum::pi;
                    upper = theta;
                    
                    // construct the first proposal
                    draws_temp = cos(theta) * draws_cur + sin(theta) * draws_prop;
                    
                    // initialize the log-likelihood (indicator for a valid proposal)
                    valid_path = true;
                    
                    // get the LNA path
                    for(int j = 0; j < (n_times-1); ++j) {
                          
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
                                      
                                      log_lna = lna_drift + svd_U * draws_temp.col(j);         // map the LNA draws
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
                          
                          // when initializing, draws leading to negative compartments or volumes are resampled.
                          // compute the LNA increment
                          nat_lna = arma::exp(log_lna) - 1;
                          
                          // update the compartment volumes
                          init_volumes_prop = init_volumes + stoich_matrix * nat_lna;
                          
                          // ensure monotonicity of the increment, compartment volumes, 
                          // and compartment volumes with forcings
                          attempt = 0;
                          if(any(nat_lna < 0) || any(init_volumes_prop < 0)) {
                                valid_path = false;
                                break;
                          }
                          
                          // set the initial volumes to the proposed values
                          init_volumes = init_volumes_prop;
                          
                          // save the increment and the prevalence
                          if(i == (ess_updates-1)) {
                                lna_path(arma::span(1,n_events), j+1) = nat_lna;
                                prev_path(arma::span(1,n_comps), j+1) = init_volumes;
                          }
                          
                          // apply forcings if called for - applied after censusing the path
                          if(forcing_inds[j+1]) {
                                init_volumes += forcing_matrix.col(j+1);
                                
                                // throw errors for negative increments or negative volumes
                                if(any(init_volumes < 0)) {
                                      valid_path = false;
                                      break;
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
                    
                    while((lower != upper) && !valid_path) {
                          
                          // initialize the parameters and volumes
                          current_params = lna_pars.row(0);   // vector for storing the current parameter values
                          CALL_SET_ODE_PARAMS(current_params, set_pars_pointer);  // set the parameters in the odeintr namespace
                          
                          // initial state vector - copy elements from the current parameter vector
                          std::copy(current_params.begin()+init_start, current_params.begin()+init_start+n_comps, init_volumes.begin());
                          init_volumes_prop = init_volumes;
                          
                          // apply forcings if called for - applied after censusing at the first time
                          if(forcing_inds[0]) {
                                init_volumes      += forcing_matrix.col(0);
                                init_volumes_prop += forcing_matrix.col(0);
                          }
                          
                          // construct the next proposal
                          if(theta<0) {
                                lower = theta;
                          } else {
                                upper = theta;
                          }
                          
                          theta = runif(1, lower, upper)[0];
                          draws_temp = cos(theta) * draws_cur + sin(theta) * draws_prop;
                          
                          // initialize the log-likelihood (indicator for a valid proposal)
                          valid_path = true;
                          
                          // get the LNA path
                          for(int j = 0; j < (n_times-1); ++j) {
                                
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
                                            
                                            log_lna = lna_drift + svd_U * draws_temp.col(j);         // map the LNA draws
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
                                
                                // when initializing, draws leading to negative compartments or volumes are resampled.
                                // compute the LNA increment
                                nat_lna = arma::exp(log_lna) - 1;
                                
                                // update the compartment volumes
                                init_volumes_prop = init_volumes + stoich_matrix * nat_lna;
                                
                                // ensure monotonicity of the increment, compartment volumes, 
                                // and compartment volumes with forcings
                                attempt = 0;
                                if(any(nat_lna < 0) || any(init_volumes_prop < 0)) {
                                      valid_path = false;
                                      break;
                                }
                                
                                // set the initial volumes to the proposed values
                                init_volumes = init_volumes_prop;
                                
                                // save the increment and the prevalence
                                if(i == (ess_updates-1)) {
                                      lna_path(arma::span(1,n_events), j+1) = nat_lna;
                                      prev_path(arma::span(1,n_comps), j+1) = init_volumes;
                                }
                                
                                // apply forcings if called for - applied after censusing the path
                                if(forcing_inds[j+1]) {
                                      init_volumes += forcing_matrix.col(j+1);
                                      
                                      // throw errors for negative increments or negative volumes
                                      if(any(init_volumes < 0)) {
                                            valid_path = false;
                                            break;
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
                    }
                    
                    if(lower != upper) {
                          // update the draws
                          draws_cur = draws_temp;
                    }
              }
              
              // save the increment and the prevalence
              incid_paths.slice(k) = lna_path.t();
              prev_paths.slice(k)  = prev_path.t();
              lna_draws.slice(k)   = draws_cur.t();
        }

        // return the paths
        return Rcpp::List::create(Rcpp::Named("draws")       = lna_draws,
                                  Rcpp::Named("incid_paths") = incid_paths,
                                  Rcpp::Named("prev_paths")  = prev_paths);
                                  // Rcpp::Named("drift_vecs") = drift_vecs,
                                  // Rcpp::Named("diff_mats") = diff_mats);
}