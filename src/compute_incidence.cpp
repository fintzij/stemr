// [[Rcpp::depends(RcppArmadillo)]]
#include "stemr_types.h"
#include "stemr_utils.h"

using namespace arma;
using namespace Rcpp;

//' Difference an incidence variable in a census matrix.
//'
//' @param censusmat matrix of compartment counts at census times, to be updated
//' @param col_inds column indices for which incidence should be computed
//' @param row_inds list of vectors of row indices for each of the incidence variables
//'
//'
//' @return update the census matrix in place
//' @export
// [[Rcpp::export]]
void compute_incidence(arma::mat& censusmat, arma::uvec& col_inds, Rcpp::List& row_inds) {

        int n_vars = col_inds.n_elem;           // number of incidence variables
        int last_row = censusmat.n_rows - 1;    // index of the final row in the census matrix

        for(int k=0; k < n_vars; ++k) {

                arma::uvec census_inds = row_inds[k];
                int col_ind = col_inds[k];
                int n_times = census_inds.n_elem;

                if(census_inds[0] != 0) {
                        censusmat(arma::span(0, census_inds[0]-1), col_ind).fill(0);
                }

                for(int j=0; j < n_times - 1; ++j) {
                        censusmat(arma::span(census_inds[j], census_inds[j+1]-1), col_ind).fill(censusmat(census_inds[j], col_ind));
                }

                censusmat(arma::span(census_inds[n_times - 1], last_row), col_ind).fill(censusmat(census_inds[n_times - 1], col_ind));

                censusmat(arma::span(1, last_row), col_ind) = arma::diff(censusmat.col(col_ind));
        }
}