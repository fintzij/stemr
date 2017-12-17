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
      
      v.randn();
}