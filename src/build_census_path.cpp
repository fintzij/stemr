// [[Rcpp::depends(RcppArmadillo)]]
#include "stemr_types.h"
#include "stemr_utils.h"

using namespace arma;
using namespace Rcpp;

//' Construct a matrix containing the compartment counts at a sequence of census times.
//'
//' @param path matrix containing the path to be censused.
//' @param census_times vector of census times.
//' @param census_columns vector of column indices to be censused (C++ indexing
//'   beginning at 0).
//'
//' @return matrix containing the compartment counts at census times.
//' @export
// [[Rcpp::export]]
arma::mat build_census_path(Rcpp::NumericMatrix& path, Rcpp::NumericVector& census_times, Rcpp::IntegerVector& census_columns) {

        // get dimensions
        int n_census_times = census_times.size();
        int n_comps = census_columns.size();

        Rcpp::IntegerVector path_dims = path.attr("dim");

        // get path matrix pointers
        arma::mat path_mat(path.begin(), path_dims[0], path_dims[1], false);

        // initialize census matrix
        arma::mat census_matrix(n_census_times, n_comps + 1);
        census_matrix.col(0) = Rcpp::as<arma::colvec>(census_times);

        // get vector of times in the path matrix and the census time indices
        Rcpp::NumericVector timevec = path(_,0);
        Rcpp::IntegerVector census_inds(n_census_times);
        census_inds = find_interval(census_times, timevec, true, true) - 1; // rightmost_closed = true, all.inside = true

        // fill out the census matrix
        census_matrix.cols(1, n_comps) = path_mat.submat(Rcpp::as<arma::uvec>(census_inds), Rcpp::as<arma::uvec>(census_columns));

        return census_matrix;
}