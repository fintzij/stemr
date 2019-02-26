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

//' return a normalised vector
//'
//' @param v vector to be normalised
//' @param p norm
//' 
//' @return normalised vector 
//' @export
// [[Rcpp::export]]
arma::vec normalise2(arma::vec& v, int p) {
      return arma::normalise(v, p);
}