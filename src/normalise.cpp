// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

//' normalise a vector in place
//'
//' @param v vector to be normalised
//' @param p norm
//' 
//' @return normalise vector in place
//' @export
// [[Rcpp::export]]
void normalise(arma::vec& v, int p) {
      v = arma::normalise(v, p);
}

