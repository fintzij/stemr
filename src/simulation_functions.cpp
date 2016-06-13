// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadilloExtensions/sample.h>
#include "stemr_types.h"
#include "stemr_utils.hpp"

using namespace arma;
using namespace Rcpp;

//' Simulate a stochastic epidemic model path via Gillespie's direct method and
//' returns a matrix containing a simulated path from a stochastic epidemic model.
//'
//' @param flow Flow matrix
//' @param parameters Vector of parameters
//' @param constants vector of constants
//' @param tcovar matrix of time-varying covariates
//' @param init_states vector of initial compartment counts
//' @param rate_adjmat adjacency matrix for updating rates after each event
//' @param tcovar_adjmat adjacency matrix for updating rates after each time a
//'   covariate changes
//' @param tcovar_changemat indicator matrix identifying which covariates change
//'   at each time
//' @param init_dims initial estimate for dimensions of the bookkeeping matrix,
//'   calculated as sum_strata(stratum size x number states x 3), rounded to the
//'   next greatest power of 2.
//' @param rate_ptr external function pointer to the lumped rate functions.
//'
//' @return matrix with a simulated path from a stochastic epidemic model.
//' @export
// [[Rcpp::export]]

arma::mat simulate_gillespie(const Rcpp::IntegerMatrix& flow, const Rcpp::NumericVector& parameters, const Rcpp::NumericVector& constants, const arma::mat& tcovar, const arma::rowvec& init_states, const Rcpp::LogicalMatrix& rate_adjmat, const arma::mat& tcovar_adjmat, const Rcpp::LogicalMatrix& tcovar_changemat, const Rcpp::IntegerVector init_dims, SEXP rate_ptr) {

        // Get dimensions of various objects
        int n_comps = init_states.size();
        Rcpp::IntegerVector flow_dims = flow.attr("dim");
        Rcpp::IntegerVector tcovar_dims(2);
        tcovar_dims[0] = tcovar.n_rows;
        tcovar_dims[1] = tcovar.n_cols;

        // initialize bookkeeping matrix
        arma::mat path(init_dims[0], init_dims[1]);

        // insert the initial compartment counts
        Rcpp::IntegerVector comp_inds   = Rcpp::seq_len(n_comps) + 1;
        path(0, span(2,init_dims[1]-1)) = init_states;

        // initialize the time varying covariates and the left and right
        // endpoints of the first piecewise homogeneous interval
        int tcov_ind(1);                                // row index in the time-varying covariate matrix
        arma::rowvec tcovs = tcovar.row(tcov_ind);      // initialize the time-varying covariates
        double t_L = tcovar(0,0);                       // left-endpoint of first interval
        double t_R = tcovar(1,0);                       // right-endpoint of first interval
        double t_cur = t_L;                             // current time
        double t_max = tcovar(tcovar_dims[0] - 1, 0);   // maximum time
        Rcpp::NumericVector dt(1);                      // time increment

        // initialize the rates
        Rcpp::LogicalVector rate_inds(flow_dims[0], true);
        Rcpp::NumericVector rates(flow_dims[0]);
        CALL_RATE_FCN(rates, rate_inds, init_states, parameters, constants, tcovs, rate_ptr);

        // set keep_going and the row index
        bool keep_going = true;
        int ind = 1;     // row from which to begin inserting
        int ind_inc = 0; // variable for counting insertions

        // start simulating
        while(keep_going) {

                // sample the next event time
                dt = Rcpp::rexp(1, sum(rates));
                t_cur += dt[0];

                if(t_cur > t_R) {

                        if(t_cur > t_max) {        // stop simulating
                                keep_going = false;

                        } else {                   // increment the time-homogeneous interval and the rates
                                t_cur = t_R;                    // left-endpoint
                                t_L   = t_R;                    // right-endpoint
                                tcov_ind += 1;                  // index in the time-varying covariate matrix
                                tcovs = tcovar.row(tcov_ind);   // time-varying covariate vector

                                rate_inds =
                        }

                }
        }

        return path;
}
