#include <Rcpp.h>
using namespace Rcpp;

//' Transform the parameters from their natural scale to the estimation scale
//'
//' @param natural_params numeric vector of parameter values on their natural
//'   scale
//' @param scaled_params numeric vector of scaled parameter values, to be
//'   modified in pleace
//' @param scales character vector of estimation scales, either 'linear', 'log',
//'   or 'logit'
//'
//' @return updates the scaled_params vector in place with the scaled parameter
//'   values
//' @export
// [[Rcpp::export]]
void to_estimation_scale(Rcpp::NumericVector& natural_params, Rcpp::NumericVector& scaled_params, Rcpp::CharacterVector& scales) {

        // get number of parameters
        int n_pars = natural_params.size();

        for(int r=0; r < n_pars; ++r) {

                if(scales[r] == "linear") {
                        scaled_params[r] = natural_params[r];

                } else if(scales[r] == "log") {
                        scaled_params[r] = log(natural_params[r]);

                } else if(scales[r] == "logit") {
                        scaled_params[r] = -log(1/natural_params[r] - 1);

                }
        }
}

//' Transform the parameters from estimation scales to their natural scales
//'
//' @param natural_params numeric vector of parameter values on their natural
//'   scale, to be modified in place
//' @param scaled_params numeric vector of scaled parameter values
//' @param scales character vector of estimation scales, either 'linear', 'log',
//'   or 'logit'
//'
//' @return updates the natural_params vector in place with the inverse of the
//'   scale transformation functions
//' @export
// [[Rcpp::export]]
void from_estimation_scale(Rcpp::NumericVector& natural_params, Rcpp::NumericVector& scaled_params, Rcpp::CharacterVector& scales) {

        // get number of parameters
        int n_pars = natural_params.size();

        for(int r=0; r < n_pars; ++r) {

                if(scales[r] == "linear") {
                        natural_params[r] = scaled_params[r];

                } else if(scales[r] == "log") {
                        natural_params[r] = exp(scaled_params[r]);

                } else if(scales[r] == "logit") {
                        natural_params[r] = 1 / (1 + exp(-scaled_params[r]));

                }
        }
}
