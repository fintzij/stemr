// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace arma;
using namespace Rcpp;

//' Construct a matrix containing the compartment counts and the incidence at a sequence of census times.
//'
//' @param path matrix containing the path to be censused.
//' @param census_path matrix to be filled out with the path.
//' @param flow_matrix_lna matrix containing the flow matrix for the LNA (no incidence)
//' @param do_prevalence should the prevalence be computed
//' @param init_state the initial compartment counts
//' @param incidence_codes_lna vector of column indices in the LNA path for which
//'   incidence will be computed. This function merely records the cumulative
//'   incidence at the census times from the LNA path.
//'
//' @return matrix containing the compartment counts at census times.
//' @export
// [[Rcpp::export]]
void census_lna(const arma::mat& path, arma::mat& census_path, const arma::uvec& census_inds, const arma::mat& flow_matrix_lna, bool do_prevalence,  const arma::rowvec& init_state, const arma::uvec& incidence_codes_lna) {

        // get dimensions
        int n_census_times = census_inds.n_elem;
        int n_comps        = flow_matrix_lna.n_cols;
        int n_rates        = flow_matrix_lna.n_rows;
        int n_incidence    = census_path.n_cols - n_comps - 1;

        if(do_prevalence) {
                Rcpp::IntegerVector prev_col_inds = Rcpp::seq_len(n_rates);
                census_path.cols(1, n_comps) = arma::repmat(init_state, n_census_times, 1) + path(census_inds, Rcpp::as<arma::uvec>(prev_col_inds)) * flow_matrix_lna;
        }

        if(n_incidence != 0) {
                for(int k=0; k < n_incidence; ++k) {
                        for(int j=0; j < n_census_times; ++j) {
                                census_path(j, k + n_comps + 1) = path(census_inds[j], incidence_codes_lna[k]);
                        }
                }
        }
}