#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

//' Identify which rates to update when a state transition event occurs.
//'
//' @param rate_inds vector of rate indices to be modified
//' @param M adjacency matrix for which rates need to be updated in response to a transition
//' @param event_code column in the rate adjacency matrix
//'
//' @return modifies logical vector stating which rates need to be updated
//' @export
// [[Rcpp::export]]
void rate_update_event(Rcpp::LogicalVector& rate_inds, const Rcpp::LogicalMatrix& M, int event_code) {

        rate_inds = M(_, event_code);
}