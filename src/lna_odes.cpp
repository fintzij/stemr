// [[Rcpp::depends(RcppArmadillo)]]
#include "stemr_types.h"
#include "stemr_utils.hpp"

using namespace arma;
using namespace Rcpp;

//' Construct a matrix containing the compartment counts at a sequence of census
//' times.
//'
//' @param t time
//' @param state vector containing the current values of the deterministic,
//'   stochastic, and variance parts of the LNA
//' @param stemr_lnamod list containing the stemr LNA model objects, including
//'   the external pointers to the hazard and jacobian functions (SEXPs), the
//'   vectors of model parameters, constants, and time-varying covariates, and
//'   the stoichiometry matrix (arma::mat).
//'
//' @return matrix containing the compartment counts at census times.
//' @export
// [[Rcpp::export]]
Rcpp::List lna_odes(double t, arma::vec& state, const Rcpp::List& stemr_lnamod) {

        // extract the model objects
        SEXP hazard_ptr   = stemr_lnamod["hazard_ptr"];                                         // external pointer to the hazard function
        SEXP jacobian_ptr = stemr_lnamod["jacobian_ptr"];                                       // external pointer to the jacobian function
        arma::mat flow_matrix = Rcpp::as<arma::mat>(stemr_lnamod["flow_matrix"]);               // stoichiometry matrix
        Rcpp::NumericVector params = Rcpp::as<Rcpp::NumericVector>(stemr_lnamod["parameters"]); // vector of parameters
        Rcpp::NumericVector consts = Rcpp::as<Rcpp::NumericVector>(stemr_lnamod["constants"]);  // vector of constants
        Rcpp::NumericVector tcovar = Rcpp::as<Rcpp::NumericVector>(stemr_lnamod["tcovar"]);     // vector of time-varying covariate values

        // get the size of various objects
        int n_comps = flow_matrix.n_rows;       // number of compartments
        int len_state = state.n_elem;           // length of the state vector

        // get the number of compartments and obtain the vectors for the odes
        arma::vec det_proc   = state.subvec(0, n_comps - 1);
        arma::vec innov_proc   = state.subvec(n_comps, 2*n_comps - 1);
        arma::mat resid_proc = arma::reshape(state.subvec(2*n_comps, len_state - 1), n_comps, n_comps);

        // vector to store the results
        arma::vec D_sys(len_state);

        // compute the hazard and the jacobian
        arma::vec hazards  = COMPUTE_HAZARD(t, det_proc, params, consts, tcovar, hazard_ptr);
        arma::mat jacobian = COMPUTE_JACOBIAN(t, det_proc, params, consts, tcovar, jacobian_ptr);

        arma::mat d_hazards = arma::diagmat(hazards);

        // construct the ODEs
        D_sys.subvec(0, n_comps - 1) = flow_matrix * hazards;
        D_sys.subvec(n_comps, 2*n_comps - 1) = flow_matrix * (jacobian * innov_proc);
        D_sys.subvec(2*n_comps, len_state - 1) = arma::vectorise(flow_matrix * jacobian * resid_proc + resid_proc * jacobian.t() * flow_matrix.t() + flow_matrix * arma::diagmat(hazards) * flow_matrix.t());

        return Rcpp::List::create(D_sys);
}