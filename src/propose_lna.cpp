// [[Rcpp::depends(RcppArmadillo)]]
#include "stemr_types.h"
#include "stemr_utils.h"

using namespace Rcpp;
using namespace arma;

//' Simulate a data matrix from the measurement process of a stochastic epidemic
//' model.
//'
//' @param lna_times vector of times at which the LNA must be evaluated
//' @param lna_pars numeric matrix of parameters, constants, and time-varying
//'   covariates at each of the lna_times
//' @param param_update_inds logical vector indicating at which of the times the
//'   LNA parameters need to be updated.
//' @param flow_matrix original flow matrix giving the changes to compartments
//'   from each reaction
//' @param lna_pointer external pointer to LNA integration fcn
//' @param set_pars_pointer external pointer to the function for setting the LNA
//'   parameters.
//'
//' @return array with slices containing the LNA path, the path of the residual
//'   process, the drift process, the residual process (conditional means), and
//'   diffusion process.
//' @export
// [[Rcpp::export]]
Rcpp::List propose_lna(const arma::colvec& lna_times, const Rcpp::NumericMatrix& lna_pars,
                       const Rcpp::LogicalVector& param_update_inds, const arma::mat& flow_matrix,
                       SEXP lna_pointer, SEXP set_pars_pointer) {

        // get the dimensions of various objects
        int n_comps = flow_matrix.n_rows;         // number of model compartments = number of rates
        int n_odes  = 2*n_comps + n_comps*n_comps;// number of ODEs
        int n_times = lna_times.n_elem;           // number of times at which the LNA must be evaluated

        // initialize the objects used in each time interval
        double t_L = 0;
        double t_R = 0;
        arma::mat lna_step(1, n_comps);                         // vector for storing the sampled lna values
        Rcpp::NumericVector current_params(lna_pars.ncol());    // vector for storing the current parameter values
        current_params = lna_pars.row(0);                       // set the current parameter values
        CALL_SET_LNA_PARAMS(current_params, set_pars_pointer);  // set the parameters in the odeintr namespace

        // initialize the LNA objects - the vector for storing the ODES, the state vector, and the Jacobian
        Rcpp::NumericVector lna_state_vec(n_odes);                        // vector to store the results of the ODEs
        arma::mat lna_path(n_times, n_comps + 1, arma::fill::zeros);      // matrix to store the path
        arma::mat residual_path(n_times, n_comps + 1, arma::fill::zeros); // matrix to store the residual path
        arma::mat drift_process(n_times, n_comps, arma::fill::zeros);     // vector for the deterministic drift process
        arma::mat residual_process(n_times, n_comps,arma::fill::zeros);   // vector for the residual process
        arma::cube diffusion_process(n_comps, n_comps, n_times, arma::fill::zeros); // matrix for the diffusion

        // indices at which the residual and diffusion elements of lna_state vec start
        int resid_start = n_comps;
        int diff_start  = 2*n_comps;

        // insert the census times into the path matrix
        lna_path.col(0)      = lna_times;
        residual_path.col(0) = lna_times;

        // iterate over the time sequence, solving the LNA over each interval
        for(int j=1; j < n_times; ++j) {

                // set the times of the interval endpoints
                t_L = lna_times[j-1];
                t_R = lna_times[j];

                // update the parameters if they need to be updated
                if(param_update_inds[j-1]) {
                        current_params = lna_pars.row(j-1);
                        CALL_SET_LNA_PARAMS(current_params, set_pars_pointer);
                }

                // integrate the LNA
                CALL_INTEGRATE_STEM_LNA(lna_state_vec, t_L, t_R, 1.0, lna_pointer);

                // transfer the elements of the lna_state_vec to the process objects
                drift_process.row(j) = Rcpp::as<arma::rowvec>(lna_state_vec).subvec(0,resid_start-1);
                residual_process.row(j) = Rcpp::as<arma::rowvec>(lna_state_vec).subvec(resid_start, diff_start-1);
                diffusion_process.slice(j) = arma::reshape(Rcpp::as<arma::rowvec>(lna_state_vec).subvec(diff_start, n_odes-1), n_comps, n_comps);

                // enforce diagonal dominance on the diffusion process and ensure its symmetry
                diffusion_process.slice(j) = arma::symmatu(diffusion_process.slice(j));

                // sample the next state and insert the sampled value into the path
                lna_step = rmvtn(1, residual_process.row(j), diffusion_process.slice(j));
                lna_path(j, arma::span(1, n_comps)) = drift_process.row(j) + lna_step.row(0);

                // // if any values are negative resample the next state
                // // n.b. unlike in the ESS we need to check that the path is monotonically increasing
                while(arma::any(lna_path(j, arma::span(1, n_comps)) < 0) ||
                      arma::any(lna_path(j, arma::span(1, n_comps)) < lna_path(j-1, arma::span(1,n_comps)))) {

                        lna_step = rmvtn(1, (residual_process.row(j)), diffusion_process.slice(j));
                        lna_path(j, arma::span(1, n_comps)) = drift_process.row(j) + lna_step.row(0);
                }

                // insert the sampled value into the path
                residual_path(j, arma::span(1, n_comps)) = lna_step.row(0);

                // copy the current value of the residual path to the LNA state vector
                std::copy(lna_step.begin(), lna_step.end(), lna_state_vec.begin() + resid_start);

                // zero out the diffusion part of the lna_state_vec
                std::fill(lna_state_vec.begin() + diff_start, lna_state_vec.end(), 0.0);
        }

        // return the paths
        return Rcpp::List::create(Rcpp::Named("lna_path")  = lna_path,
                                  Rcpp::Named("res_path")  = residual_path,
                                  Rcpp::Named("drift")     = drift_process,
                                  Rcpp::Named("residual")  = residual_process,
                                  Rcpp::Named("diffusion") = diffusion_process);
}