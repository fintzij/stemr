// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

//' Draw new N(0,1) values and fill a vector.
//'
//' @param v vector to fill with new N(0,1) draws
//'
//' @return draw new values in place
//' @export
// [[Rcpp::export]]
void draw_normals(arma::vec& v) {
      std::copy_n(Rcpp::rnorm(v.n_elem).begin(), v.n_elem, v.begin());
}

//' Draw new N(0,1) values and fill a matrix.
//'
//' @param M matrix to fill with new N(0,1) draws
//'
//' @return draw new values in place
//' @export
// [[Rcpp::export]]
void draw_normals2(arma::mat& M) {
      std::copy_n(Rcpp::rnorm(M.n_elem).begin(), M.n_elem, M.begin());
}

//' Sample the unit sphere.
//'
//' @param v vector to fill with a vector of draws on the unit sphere
//'
//' @return draw new values in place
//' @export
// [[Rcpp::export]]
void sample_unit_sphere(arma::vec& v) {
      std::copy_n(Rcpp::rnorm(v.n_elem).begin(), v.n_elem, v.begin());
      v = arma::normalise(v, 2);
}
