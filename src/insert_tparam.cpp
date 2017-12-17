// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

//' Insert time-varying parameters into a tcovar matrix.
//'
//' @param tcovar matrix into which the parameter values should be copied
//' @param values vector of values that should be copied in
//' @param col_ind C++ index for the column where the values should go
//' @param tpar_inds C++ indices for the vector elements that go in each row
//'
//' @return copy the values in place
//' @export
// [[Rcpp::export]]
void insert_tparam(arma::mat& tcovar, const arma::vec& values, int col_ind, const arma::uvec& tpar_inds) {
      
      tcovar.col(col_ind) = values.elem(tpar_inds);
}