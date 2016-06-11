// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace arma;
using namespace Rcpp;
void RATES_LUMPED(Rcpp::NumericVector& rates, const Rcpp::LogicalVector& inds, const arma::rowvec& state, const Rcpp::NumericVector& parameters, const Rcpp::NumericVector& constants, const arma::rowvec& tcovar) {

 if(inds[0]) rates[0] = (parameters[0] * state[4] + parameters[2] * tcovar[2] + parameters[1] * (log(state[5]) + log(state[6]) + log(state[7])) + parameters[5]*sin(2*M_PI*tcovar[2]/52) + parameters[6]*cos(2*M_PI*tcovar[2]/52)) * state[0];
 if(inds[1]) rates[1] = (parameters[0] * state[5] + parameters[2] * tcovar[2] + parameters[1] * (log(state[4]) + log(state[6]) + log(state[7])) + parameters[7]*sin(2*M_PI*tcovar[2]/52) + parameters[8]*cos(2*M_PI*tcovar[2]/52)) * state[1];
 if(inds[2]) rates[2] = (parameters[0] * state[6] + parameters[2] * tcovar[2] + parameters[1] * (log(state[4]) + log(state[5]) + log(state[7])) + parameters[9]*sin(2*M_PI*tcovar[2]/52) + parameters[10]*cos(2*M_PI*tcovar[2]/52)) * state[2];
 if(inds[3]) rates[3] = (parameters[0] * state[7] + parameters[2] * tcovar[2] + parameters[1] * (log(state[4]) + log(state[5]) + log(state[6])) + parameters[11]*sin(2*M_PI*tcovar[2]/52) + parameters[12]*cos(2*M_PI*tcovar[2]/52)) * state[3];
 if(inds[4]) rates[4] = (parameters[3]) * state[4];
 if(inds[5]) rates[5] = (parameters[3]) * state[5];
 if(inds[6]) rates[6] = (parameters[3]) * state[6];
 if(inds[7]) rates[7] = (parameters[3]) * state[7];
 if(inds[8]) rates[8] = (parameters[4]) * state[10];
 if(inds[9]) rates[9] = (parameters[4]) * state[11];
}
typedef void(*ratefcn_ptr)(Rcpp::NumericVector& rates, const Rcpp::LogicalVector& inds, const arma::rowvec& state, const Rcpp::NumericVector& parameters, const Rcpp::NumericVector& constants, const arma::rowvec& tcovar);
// [[Rcpp::export]]
Rcpp::XPtr<ratefcn_ptr> putRatePtrInXPtr(void) {
return(Rcpp::XPtr<ratefcn_ptr>(new ratefcn_ptr(&RATES_LUMPED)));
}
