// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace arma;
using namespace Rcpp;

//' Construct a matrix containing the compartment counts and the incidence at a sequence of census times.
//'
//' @param path matrix containing the path to be censused (cumulative incidence).
//' @param census_path matrix to be filled out with the path.
//' @param census_inds vector of indices for census interval endpoints.
//' @param lna_event_inds vector of column indices in the path matrix for events that
//'   should be censused.
//' @param flow_matrix_lna matrix containing the flow matrix for the LNA (no incidence)
//' @param do_prevalence should the prevalence be computed
//' @param init_state the initial compartment counts
//' @param forcing_inds logical vector of indicating at which times in the
//'   time-varying covariance matrix a forcing is applied.
//' @param forcing_matrix matrix containing the forcings.
//'
//' @return matrix containing the compartment counts at census times.
//' @export
// [[Rcpp::export]]
void census_lna(const arma::mat& path,
                arma::mat& census_path,
                const arma::uvec& census_inds,
                const arma::uvec& lna_event_inds,
                const arma::mat& flow_matrix_lna,
                bool do_prevalence,
                const arma::rowvec& init_state,
                const arma::mat& forcing_matrix) {

        // get dimensions
        int n_census_times  = census_inds.n_elem;
        int n_census_events = lna_event_inds.n_elem;
        int n_comps         = flow_matrix_lna.n_cols;
        int n_rates         = flow_matrix_lna.n_rows;

        // get indices in the census_path matrix to keep incidence
        int incid_start = flow_matrix_lna.n_cols + 1;

        // census the incidence increments
        for(int k = 1; k < n_census_times; ++k) {

                for(int j = 0; j < n_census_events; ++j) {
                        census_path(k-1, incid_start + j) = arma::sum(path(arma::span(census_inds[k-1]+1, census_inds[k]),
                                                                           lna_event_inds[j]));
                }

        }

        // compute the prevalence if called for
        if(do_prevalence) {

                // initialize the state and increment vectors
                arma::rowvec state(init_state);
                arma::rowvec increment(n_rates, arma::fill::zeros);

                for(int k=1; k < n_census_times-1; ++k) {

                        increment = arma::sum(path(arma::span(census_inds[k-1] + 1, census_inds[k]),
                                                   arma::span(1, n_rates)), 0);

                        state += increment * flow_matrix_lna +
                                 arma::sum(forcing_matrix.cols(census_inds[k], census_inds[k+1] - 1), 1).t();

                        census_path(k-1, arma::span(1, n_comps)) = state;
                }
        }
}