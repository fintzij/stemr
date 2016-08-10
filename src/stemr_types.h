#ifndef stemr_types_h
#define stemr_types_h

#include <RcppArmadillo.h>

using namespace arma;
using namespace Rcpp;

typedef void(*ratefcn_ptr)(Rcpp::NumericVector& rates, const Rcpp::LogicalVector& inds,
             const arma::rowvec& state, const Rcpp::NumericVector& parameters,
             const Rcpp::NumericVector& constants, const arma::rowvec& tcovar);

typedef void(*d_measure_ptr)(Rcpp::NumericMatrix& emitmat, const Rcpp::LogicalVector& emit_inds,
             const int record_ind, const Rcpp::NumericVector& record, const Rcpp::NumericVector& state,
             const Rcpp::NumericVector& parameters, const Rcpp::NumericVector& constants,
             const Rcpp::NumericVector& tcovar);

typedef void(*r_measure_ptr)(Rcpp::NumericMatrix& obsmat, const Rcpp::LogicalVector& emit_inds,
             const int record_ind, const Rcpp::NumericVector& state, const Rcpp::NumericVector& parameters,
             const Rcpp::NumericVector& constants, const Rcpp::NumericVector& tcovar);

// typedef Rcpp::NumericVector(*lna_ptr)(Rcpp::NumericVector& init, double start, double end, double step_size);
typedef void(*lna_ptr)(Rcpp::NumericVector& init, double start, double end, double step_size);
typedef void(*set_pars_ptr)(Rcpp::NumericVector& p);

// typedef arma::vec(*hazard_ptr)(double t, const arma::vec& state, const Rcpp::NumericVector& parameters, const Rcpp::NumericVector& constants, const Rcpp::NumericVector& tcovar);
// typedef arma::mat(*jacobian_ptr)(double t, const arma::vec& state, const Rcpp::NumericVector& parameters, const Rcpp::NumericVector& constants, const Rcpp::NumericVector& tcovar);

#endif