// [[Rcpp::depends(RcppArmadillo)]]
#define ARMA_DONT_PRINT_ERRORS
#include "stemr_types.h"
#include "stemr_utils.h"

using namespace Rcpp;
using namespace arma;

//' Map N(0,1) stochastic perturbations to an LNA path.
//'
//' @param pathmat matrix where the LNA path should be stored
//' @param draws matrix of N(0,1) draws to be mapped to an LNA path
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
//' @param svd_d vector in which to store SVD singular values
//' @param svd_U matrix in which to store the U matrix of the SVD
//' @param svd_V matrix in which to store the V matrix of the SVD
//' @param step_size initial step size for the ODE solver (adapted internally,
//' but too large of an initial step can lead to failure in stiff systems).
//' @param lna_pointer external pointer to LNA integration function.
//' @param set_pars_pointer external pointer to the function for setting the LNA
//'   parameters.
//'
//' @return fill out pathmat with the LNA path corresponding to the stochastic
//'   perturbations.
//'
//' @export
// [[Rcpp::export]]
void map_draws_2_lna(arma::mat& pathmat,
                     const arma::mat& draws,
                     const arma::rowvec& lna_times,
                     const Rcpp::NumericMatrix& lna_pars,
                     Rcpp::NumericVector& lna_param_vec,
                     const Rcpp::IntegerVector& lna_param_inds,
                     const Rcpp::IntegerVector& lna_tcovar_inds,
                     const int init_start,
                     const Rcpp::LogicalVector& param_update_inds,
                     const arma::mat& stoich_matrix,
                     const Rcpp::LogicalVector& forcing_inds,
                     const arma::mat& forcing_matrix,
                     arma::vec& svd_d,
                     arma::mat& svd_U,
                     arma::mat& svd_V,
                     double step_size,
                     SEXP lna_pointer,
                     SEXP set_pars_pointer) {

        // get the dimensions of various objects
        int n_events = stoich_matrix.n_cols;         // number of transition events, e.g., S2I, I2R
        int n_comps  = stoich_matrix.n_rows;         // number of model compartments (all strata)
        int n_odes   = n_events + n_events*n_events; // number of ODEs
        int n_times  = lna_times.n_elem;             // number of times at which the LNA must be evaluated
        int n_tcovar     = lna_tcovar_inds.size();   // number of time-varying covariates or parameters

        // initialize the objects used in each time interval
        double t_L = 0;
        double t_R = 0;

        // vector of parameters, initial compartment columes, constants, and time-varying covariates
        std::copy(lna_pars.row(0).begin(), lna_pars.row(0).end(), lna_param_vec.begin());

        CALL_SET_ODE_PARAMS(lna_param_vec, set_pars_pointer); // set the parameters in the odeintr namespace

        // initial state vector - copy elements from the current parameter vector
        arma::vec init_volumes(lna_param_vec.begin() + init_start, n_comps);

        // initialize the LNA objects
        bool good_svd = true;

        Rcpp::NumericVector lna_state_vec(n_odes); // the vector for storing the current state of the LNA ODEs

        arma::vec lna_drift(n_events, arma::fill::zeros);               // incidence mean vector (log scale)
        arma::mat lna_diffusion(n_events, n_events, arma::fill::zeros); // diffusion matrix

        arma::vec log_lna(n_events, arma::fill::zeros);  // LNA increment, log scale
        arma::vec nat_lna(n_events, arma::fill::zeros);  // LNA increment, natural scale

        // apply forcings if called for - applied after censusing at the first time
        if(forcing_inds[0]) {
                init_volumes += forcing_matrix.col(0);
        }

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
                lna_diffusion = arma::symmatu(lna_diffusion);

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

                // compute the LNA increment
                nat_lna = arma::vec(expm1(Rcpp::NumericVector(log_lna.begin(), log_lna.end())));

                // save the LNA increment
                pathmat(j+1, arma::span(1, n_events)) = nat_lna.t();

                // update the initial volumes
                init_volumes += stoich_matrix * nat_lna;

                // if any increments or volumes are negative, throw an error
                try{
                        if(any(nat_lna < 0)) {
                                throw std::runtime_error("Negative increment.");
                        }

                        if(any(init_volumes < 0)) {
                                throw std::runtime_error("Negative compartment volumes.");
                        }

                } catch(std::runtime_error &err) {

                        forward_exception_to_r(err);

                } catch(...) {
                        ::Rf_error("c++ exception (unknown reason)");
                }
                
                // apply forcings if called for - applied after censusing the path
                if(forcing_inds[j+1]) {
                      
                      init_volumes += forcing_matrix.col(j+1);
                      
                      // Check for positivity
                      try{
                            if(any(init_volumes < 0)) {
                                  throw std::runtime_error("Negative compartment volumes.");
                            }
                            
                      } catch(std::runtime_error &err) {
                            
                            forward_exception_to_r(err);
                            
                      } catch(...) {
                            ::Rf_error("c++ exception (unknown reason)");
                      }
                }

                // update the parameters if they need to be updated
                if(param_update_inds[j+1]) {
                      
                      // time-varying covariates and parameters
                      std::copy(lna_pars.row(j+1).end() - n_tcovar,
                                lna_pars.row(j+1).end(),
                                lna_param_vec.end() - n_tcovar);
                      
                }

                // copy the new initial volumes into the vector of parameters
                std::copy(init_volumes.begin(), init_volumes.end(), lna_param_vec.begin() + init_start);

                // set the lna parameters and reset the LNA state vector
                CALL_SET_ODE_PARAMS(lna_param_vec, set_pars_pointer);
        }
}