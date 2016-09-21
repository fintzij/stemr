// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include "stemr_types.h"

using namespace arma;
using namespace Rcpp;

//' Compute the LNA by calling the LNA functions via external Xptr.
//'
//' @param t time
//' @param state vector containing the drift, residual, and diffusion subvectors
//' @param parms list containing the LNA model parameters, stoichiometry matrix,
//'   and LNA pointer
//'
//' @export
// [[Rcpp::export]]
Rcpp::List CALL_COMPUTE_LNA(double t, arma::vec& state, Rcpp::List& parms) {

        SEXP lna_pointer = parms["lna_pointer"];
        Rcpp::XPtr<compute_lna_ptr> xpfun(lna_pointer); // Receive the SEXP and put in Xptr
        compute_lna_ptr fun = *xpfun;                           // get function via pointer

        Rcpp::NumericVector parameters = Rcpp::as<Rcpp::NumericVector>(parms["parameters"]); // vector of parameters
        arma::mat stoich = Rcpp::as<arma::mat>(parms["stoich"]);

        return fun(t, state, parameters, stoich);
}