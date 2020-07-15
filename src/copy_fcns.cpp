#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

//' Insert parameters into each row of a parameter matrix
//'
//' @param lnapars matrix of lna parameters, constants, and time-varying covars
//' @param parameters vector of parameters to be copied into the matrix
//'
//' @return modifies the lna parameter matrix in place
//' @export
// [[Rcpp::export]]
void pars2lnapars(arma::mat& lnapars, const arma::rowvec& parameters) {

        int n_pars  = parameters.n_elem;
        lnapars.cols(0, n_pars-1).each_row() = parameters;
}

//' Insert parameters into the first row of a parameter matrix
//'
//' @param lnapars matrix of lna/ode parameters, constants, and time-varying covars
//' @param parameters vector of parameters to be copied into the matrix
//' @param c_start index of the initial column
//'
//' @return modifies the lna parameter matrix in place
//' @export
// [[Rcpp::export]]
void pars2lnapars2(arma::mat& lnapars, const arma::rowvec& parameters, int c_start) {
      
      int c_end = parameters.n_elem - 1;
      lnapars(0, arma::span(c_start, c_start + c_end)) = parameters;
}

//' Insert parameters into the first row of a parameter matrix
//' 
//' @param parmat parameter matrics
//' @param pars vector of parameters to insert
//' @param colinds vector of column indices
//' @param rowinds vector of row indices, just the first row by default.
void pars2parmat(arma::mat& parmat, 
                 const arma::rowvec pars,
                 const arma::uvec colinds,
                 const arma::uvec rowinds = 0) {
  parmat.submat(rowinds, colinds) = pars; 
}

//' Copy an element from one vector into another
//'
//' @param dest destination row vector
//' @param orig origin row vector
//' @param ind C++ style index for the element to be copied
//'
//' @return copy an element of one row vector into another.
//' @export
// [[Rcpp::export]]
void copy_elem(arma::rowvec& dest, const arma::rowvec& orig, int ind) {

        dest[ind] = orig[ind];
}

//' Copy an multiple elements from one vector into another
//'
//' @param dest destination row vector
//' @param orig origin row vector
//' @param ind C++ style index for the element to be copied
//'
//' @return copy an element of one row vector into another.
//' @export
// [[Rcpp::export]]
void copy_elem2(arma::rowvec& dest, const arma::rowvec& orig, const arma::uvec& inds) {
      
      dest.elem(inds) = orig.elem(inds);
}

//' Increment an element of a vector by 1
//'
//' @param vec destination row vector
//' @param ind C++ style index for the element to be copied
//'
//' @return Add 1 to an element of a vector
//' @export
// [[Rcpp::export]]
void increment_elem(arma::vec& vec, int ind) {
      
      vec[ind] += 1;
}

//' Copy the contents of one vector into another
//'
//' @param dest destination row vector
//' @param orig origin row vector
//'
//' @return copy the elements of one row vector into another.
//' @export
// [[Rcpp::export]]
void copy_vec(arma::rowvec& dest, const arma::rowvec& orig) {

        dest = orig;
}

//' Copy the contents of one vector into another
//'
//' @param dest destination row vector
//' @param orig origin row vector
//' @param inds vector of indices in the destination
//'
//' @return copy the elements of one row vector into another.
//' @export
// [[Rcpp::export]]
void copy_vec2(arma::rowvec& dest, const arma::rowvec& orig, const arma::uvec& inds) {
      
      dest.elem(inds) = orig;
}

//' Copy the contents of one matrix into another
//'
//' @param dest destination matrix
//' @param orig origin matrix
//'
//' @return copy the elements of one matrix into another.
//' @export
// [[Rcpp::export]]
void copy_mat(arma::mat& dest, const arma::mat& orig) {
      
      dest = orig;
}

//' Insert one matrix into another
//'
//' @param dest destination matrix
//' @param orig origin matrix
//' @param rowinds vector of row indices
//' @param colinds vector of column indices
//'
//' @return copy the elements of one matrix into another.
//' @export
// [[Rcpp::export]]
void insert_block(arma::mat& dest, const arma::mat& orig, const arma::uvec& rowinds, const arma::uvec& colinds) {
      
      dest(rowinds, colinds) = orig;
}

//' Copy the contents of one matrix into another
//'
//' @param dest destination matrix
//' @param orig origin matrix
//' @param ind column index
//'
//' @return copy the elements of one matrix into another.
//' @export
// [[Rcpp::export]]
void copy_col(arma::mat& dest, const arma::mat& orig, int ind) {
      
      dest.col(ind) = orig.col(ind);
}

//' Copy the columns of one matrix into another
//'
//' @param dest destination matrix
//' @param orig origin matrix
//' @param ind column index
//'
//' @return copy the elements of one matrix into another.
//' @export
// [[Rcpp::export]]
void copy_pathmat(arma::mat& dest, const arma::mat& orig) {
      
      dest.cols(1, orig.n_cols-1) = orig.cols(1, orig.n_cols-1);
}

//' Copy some of the rows of one matrix into another
//'
//' @param dest destination matrix
//' @param orig origin matrix
//' @param inds row indices
//'
//' @return copy the elements of one matrix into another.
//' @export
// [[Rcpp::export]]
void copy_2_rows(arma::mat& dest, const arma::mat& orig, const arma::uvec& inds) {
      
      dest.rows(inds) = orig;
}

//' Copy a matrix into a slice of an array
//'
//' @param dest array into which to copy
//' @param orig matrix to copy
//' @param ind slice index (C++)
//'
//' @return copy a matrix into an array.
//' @export
// [[Rcpp::export]]
void mat_2_arr(arma::cube& dest, const arma::mat& orig, int ind) {
      
      dest.slice(ind) = orig;
}

//' Reset a vector by filling it with an element
//'
//' @param v vector to fill with zeros
//' @param value to insert
//'
//' @return reset vector in place
//' @export
// [[Rcpp::export]]
void reset_vec(arma::vec& v, double value = 0) {
      v.fill(value);
}

//' Add the contents of one vector to another vector
//'
//' @param dest target vector
//' @param orig vector to be added
//' @param indices in the target
//'
//' @return add the elements of one row vector to another.
//' @export
// [[Rcpp::export]]
void add2vec(arma::rowvec& target, const arma::rowvec& increments, const arma::uvec& inds) {
      target.elem(inds) += increments;
}