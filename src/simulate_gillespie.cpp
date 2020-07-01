// [[Rcpp::depends(RcppArmadillo)]]
#include "stemr_types.h"
#include "stemr_utils.h"
#include <RcppArmadilloExtensions/sample.h>

using namespace arma;
using namespace Rcpp;

//' Simulate a stochastic epidemic model path via Gillespie's direct method and
//' returns a matrix containing a simulated path from a stochastic epidemic
//' model.
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
//' @param forcing_inds logical vector of indicating at which times in the
//'   time-varying covariance matrix a forcing is applied.
//' @param forcing_matrix matrix containing the forcings.
//' @param rate_ptr external function pointer to the lumped rate functions.
//'
//' @return matrix with a simulated path from a stochastic epidemic model.
//' @export
// [[Rcpp::export]]
arma::mat simulate_gillespie(const arma::mat& flow,
                             const Rcpp::NumericVector& parameters,
                             const Rcpp::NumericVector& constants,
                             const arma::mat& tcovar,
                             double t_max,
                             const arma::rowvec& init_states,
                             const Rcpp::LogicalMatrix& rate_adjmat,
                             const arma::mat& tcovar_adjmat,
                             const arma::mat& tcovar_changemat,
                             const Rcpp::IntegerVector init_dims,
                             const Rcpp::LogicalVector& forcing_inds,
                             const arma::uvec& forcing_tcov_inds,
                             const arma::mat& forcings_out,
                             const arma::cube& forcing_transfers,
                             SEXP rate_ptr) {
      
      // Get dimensions of various objects
      Rcpp::IntegerVector flow_dims(2);       // size of flow matrix
      Rcpp::IntegerVector tcovar_dims(2);     // size of tcovar matrix
      flow_dims[0]    = flow.n_rows;
      flow_dims[1]    = flow.n_cols;
      tcovar_dims[0]  = tcovar.n_rows;
      tcovar_dims[1]  = tcovar.n_cols;
      int n_forcings  = forcing_tcov_inds.n_elem;
      
      // for use with forcings
      double forcing_flow = 0;
      arma::vec forcing_distvec(init_dims[1]-2, arma::fill::zeros);
      
      // initialize bookkeeping matrix
      arma::mat path(init_dims[0], init_dims[1]);
      path(0, 1) = -1;
      int path_nrows = path.n_rows;
      
      // vector of event codes
      Rcpp::IntegerVector events = Rcpp::seq_len(flow_dims[0]) - 1; // vector of event codes
      Rcpp::IntegerVector next_event(1);                            // next event
      
      // initialize the time varying covariates and the left and right
      // endpoints of the first piecewise homogeneous interval
      int tcov_ind = 0;                                // row index in the time-varying covariate matrix
      arma::rowvec tcovs = tcovar.row(tcov_ind);      // initialize the time-varying covariates
      double t_L = tcovar(tcov_ind,0);                // left-endpoint of first interval
      double t_R = tcovar(tcov_ind + 1,0);            // right-endpoint of first interval
      double t_cur = t_L;                             // current time
      Rcpp::NumericVector dt(1);                      // time increment
      
      // insert the initial compartment counts
      path(0, arma::span(2,init_dims[1]-1)) = init_states;
      
      // initialize a state vector
      arma::rowvec state = init_states;
      
      // apply forcings if necessary
      if(forcing_inds[tcov_ind]) {
            for(int j=0; j < n_forcings; ++j) {
                  
                  // distribute the forcings proportionally to the compartment counts in the applicable states
                  forcing_flow = tcovar(tcov_ind, forcing_tcov_inds[j]);
                  forcing_distvec = arma::round(forcing_flow * normalise(forcings_out.col(j) % state.t(), 1));
                  state += (forcing_transfers.slice(j) * forcing_distvec).t();
            }
      }
      
      // insert the initial time into the path matrix
      path(0, 0) = t_cur;
      
      // initialize the rates
      Rcpp::LogicalVector rate_inds(flow_dims[0], true); // logical vector of rates to update
      Rcpp::NumericVector rates(flow_dims[0]);           // initialize vector of rates
      CALL_RATE_FCN(rates, rate_inds, state, parameters, constants, tcovs, rate_ptr); // compute rates
      
      // vector for storing the event probabilities - Rcpp's sample function modifies probs (and hence rates) in place
      Rcpp::NumericVector event_probs(flow_dims[0]);
      event_probs = rates / sum(rates);
      
      // set keep_going and the row index
      bool keep_going = true;
      int ind_start = 1;       // row from which to begin inserting, updated throughout
      int ind_cur = ind_start; // row index into which to actually insert
      
      // start simulating
      while(keep_going) {
            
            // sample the next event time
            dt = Rcpp::rexp(1, sum(rates));
            t_cur += dt[0];
            
            if(t_cur > t_R) {
                  
                  if((t_R == t_max) && (t_cur > t_max)) {
                        
                        // stop simulating
                        keep_going = false;
                        
                  } else {
                        
                        // increment the time-homogeneous interval and the rates
                        tcov_ind += 1;                       // increment the index in the time-varying covariate matrix
                        tcovs     = tcovar.row(tcov_ind);    // time-varying covariate vector
                        t_L       = t_R;                     // set left endpoint to right endpoint
                        t_cur     = t_R;                     // set current time to right endpoint
                        t_R       = tcovar(tcov_ind + 1, 0); // increment the right endpoint
                        
                        // identify rates that need to be updated
                        rate_update_tcovar(rate_inds, tcovar_adjmat, tcovar_changemat.row(tcov_ind));
                        
                        // apply forcings if necessary
                        if(forcing_inds[tcov_ind]) {
                              
                              for(int j=0; j < n_forcings; ++j) {
                                    
                                    // distribute the forcings proportionally to the compartment counts in the applicable states
                                    forcing_flow = tcovar(tcov_ind, forcing_tcov_inds[j]);
                                    forcing_distvec = arma::round(forcing_flow * normalise(forcings_out.col(j) % state.t(), 1));
                                    state += (forcing_transfers.slice(j) * forcing_distvec).t();
                              }
                              
                              // throw errors for negative volumes
                              try{
                                    if(any(state < 0)) {
                                          throw std::runtime_error("Negative compartment volumes.");
                                    }
                                    
                              } catch(std::exception &err) {
                                    
                                    forward_exception_to_r(err);
                                    
                              } catch(...) {
                                    ::Rf_error("c++ exception (unknown reason)");
                              }
                        }
                        
                        // insert the time, event, and new state vector into the path matrix
                        path(ind_cur, 0) = t_cur;                             // insert new time
                        path(ind_cur, 1) = -1;                                // event code
                        path(ind_cur, arma::span(2,init_dims[1]-1)) = state;  // state
                        
                        // increment the index
                        ind_cur += 1;
                        
                        // if there are no empty rows in the path matrix, add some
                        if(ind_cur == path_nrows) {
                              path.insert_rows(ind_cur, init_dims[0]);
                              path_nrows = path.n_rows;
                        }
                        
                        // update the rate functions
                        CALL_RATE_FCN(rates, rate_inds, state, parameters, constants, tcovs, rate_ptr);
                        
                        // update event probabilities
                        event_probs = rates / sum(rates);
                  }
                  
            } else {
                  
                  // sample the next event
                  next_event = Rcpp::RcppArmadillo::sample(events, 1, false, event_probs);
                  
                  // update the state vector
                  state += flow.row(next_event[0]);
                  
                  // insert the time, event, and new state vector into the path matrix
                  path(ind_cur, 0) = t_cur;                       // insert new time
                  path(ind_cur, 1) = next_event[0];               // event code
                  path(ind_cur, arma::span(2,init_dims[1]-1)) = state;  // state
                  
                  // increment the index
                  ind_cur += 1;
                  
                  // if there are no empty rows in the path matrix, add some
                  if(ind_cur == path_nrows) {
                        path.insert_rows(ind_cur, init_dims[0]);
                        path_nrows = path.n_rows;
                  }
                  
                  // update the rates
                  rate_update_event(rate_inds, rate_adjmat, next_event[0]);                       // identify rates that need to be updated
                  CALL_RATE_FCN(rates, rate_inds, state, parameters, constants, tcovs, rate_ptr); // update the rate functions
                  
                  event_probs = rates / sum(rates); // update the event probabilities
                  
                  // if all rates equal zero, stop simulating
                  // keep_going = Rcpp::is_true(any(rates != 0));
            }
      }
      
      path.shed_rows(ind_cur, path.n_rows - 1);
      
      // ensure that t_max is the time of the last row in path. if not, add it
      if(path(path.n_rows-1, 0) != t_max) {
            arma::rowvec last_row = path.row(path.n_rows - 1);
            last_row(0) = t_max;
            last_row(1) = -1;
            path.insert_rows(path.n_rows, last_row);
      }
      
      return path;
}