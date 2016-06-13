#ifndef stemr_UTILITIES_H
#define stemr_UTILITIES_H

#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

void CALL_RATE_FCN(Rcpp::NumericVector& rates, const Rcpp::LogicalVector& inds,
                   const arma::rowvec& state, const Rcpp::NumericVector& parameters,
                   const Rcpp::NumericVector& constants, const arma::rowvec& tcovar, SEXP rate_ptr);

Rcpp::LogicalVector rates_to_update(const arma::mat& M, const arma::rowvec I);
#endif // stemr_UTILITIES_H