#ifndef stemr_UTILITIES_H
#define stemr_UTILITIES_H

#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

void CALL_RATE_FCN(Rcpp::NumericVector& rates, const Rcpp::LogicalVector& inds,
                   const arma::rowvec& state, const Rcpp::NumericVector& parameters,
                   const Rcpp::NumericVector& constants, const arma::rowvec& tcovar, SEXP rate_ptr);

void rate_update_tcovar(Rcpp::LogicalVector& rate_inds, const arma::mat& M, const arma::rowvec I);
void rate_update_event(Rcpp::LogicalVector& rate_inds, const Rcpp::LogicalMatrix& M, int event_code);

#endif // stemr_UTILITIES_H