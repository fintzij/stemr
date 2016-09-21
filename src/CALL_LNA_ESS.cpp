// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include "stemr_types.h"

using namespace arma;
using namespace Rcpp;

//' Compute the LNA by calling the LNA functions via external Xptr.
//'
//' @param t time
//' @param state vector of compartment counts
//' @param list containing the vector of model parameters, stoichiometry matrix,
//'   and elliptical slice sampling LNA pointer
//'
//' @export
// [[Rcpp::export]]
Rcpp::List CALL_LNA_ESS(double t, arma::vec& state, Rcpp::List& parms) {

        SEXP lna_ess_pointer = parms["lna_ess_pointer"];
        Rcpp::XPtr<lna_ess_ptr> xpfun(lna_ess_pointer); // Receive the SEXP and put in Xptr
        lna_ess_ptr fun = *xpfun;                       // get function via pointer

        Rcpp::NumericVector parameters = Rcpp::as<Rcpp::NumericVector>(parms["parameters"]); // vector of parameters
        arma::mat stoich = Rcpp::as<arma::mat>(parms["stoich"]);

        return fun(t, state, parameters, stoich);
}