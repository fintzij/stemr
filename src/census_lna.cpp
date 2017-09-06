// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace arma;
using namespace Rcpp;

//' Construct a matrix containing the compartment counts and the incidence at a sequence of census times.
//'
//' @param path matrix containing the path to be censused (cumulative incidence).
//' @param census_path matrix to be filled out with the path.
//' @param census_inds vector of indices for census interval endpoints
//' @param flow_matrix_lna matrix containing the flow matrix for the LNA (no incidence)
//' @param do_prevalence should the prevalence be computed
//' @param init_state the initial compartment counts
//' @param incidence_codes vector of column indices in the LNA path for which
//'   incidence will be computed.
//' @param forcing_inds logical vector of indicating at which times in the
//'   time-varying covariance matrix a forcing is applied.
//' @param forcing_matrix matrix containing the forcings.
//'
//' @return matrix containing the compartment counts at census times.
//' @export
// [[Rcpp::export]]
void census_lna(const arma::mat& path,
                arma::mat& census_path,
                arma::uvec& census_inds,
                const arma::mat& flow_matrix_lna,
                bool do_prevalence,
                const arma::rowvec& init_state,
                const arma::mat& forcing_matrix) {

        // get dimensions
        int n_census_times = census_inds.n_elem;
        int n_comps        = flow_matrix_lna.n_cols;
        int n_rates        = flow_matrix_lna.n_rows;
        int n_incidence    = census_path.n_cols - n_comps - 1;

        // get indices
        int incid_start = flow_matrix_lna.n_cols + 1;
        int incid_end   = flow_matrix_lna.n_rows + flow_matrix_lna.n_cols;

        // census the incidence increments
        for(int k = 1; k < n_census_times; ++k) {
                census_path(k-1, arma::span(incid_start, incid_end)) =
                        arma::sum(path(arma::span(census_inds[k-1]+1, census_inds[k]),
                                       arma::span(1, n_rates)),0);
        }

        if(do_prevalence) {

                arma::rowvec state(init_state);

                for(int k=0; k < n_census_times-1; ++k) {

                        state += census_path(k, arma::span(incid_start, incid_end)) * flow_matrix_lna +
                                 arma::sum(forcing_matrix.cols(census_inds[k], census_inds[k+1] - 1), 1).t();

                        census_path(k, arma::span(1, n_comps)) = state;
                }
        }
}