// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace arma;
using namespace Rcpp;

//' Construct a matrix containing the compartment counts and the incidence at a sequence of census times.
//'
//' @param path matrix containing the path to be censused (cumulative incidence).
//' @param census_path matrix to be filled out with the path.
//' @param census_inds vector of indices for census interval endpoints.
//' @param event_inds vector of column indices in the path matrix for events that
//'   should be censused.
//' @param flow_matrix matrix containing the flow matrix for the LNA (no incidence)
//' @param do_prevalence should the prevalence be computed
//' @param init_state the initial compartment counts
//' @param forcing_inds logical vector of indicating at which times in the
//'   time-varying covariance matrix a forcing is applied.
//' @param parmat matrix with parameters, constants, and time varying
//'   covariates and parameters.
//'
//' @return matrix containing the compartment counts at census times.
//' @export
// [[Rcpp::export]]
void census_latent_path(
                const arma::mat& path,
                arma::mat& census_path,
                const arma::uvec& census_inds,
                const Rcpp::Nullable<Rcpp::IntegerVector>& event_inds,
                const arma::mat& flow_matrix,
                bool do_prevalence,
                const arma::mat& parmat,
                const arma::uvec& initdist_inds,
                const Rcpp::LogicalVector& forcing_inds,
                const arma::uvec& forcing_tcov_inds,
                const arma::mat& forcings_out,
                const arma::cube& forcing_transfers,
                arma::uvec row0 = 0) {

        // get dimensions
        int n_census_times  = census_inds.n_elem;
        int n_comps         = flow_matrix.n_cols;
        int n_rates         = flow_matrix.n_rows;
        int n_forcings      = forcing_tcov_inds.n_elem;

        // for use with forcings
        double forcing_flow = 0;
        arma::vec forcing_distvec(n_comps, arma::fill::zeros);

        // get incidence indices
        if (event_inds.isNotNull()) {

                // convert incidence indices to an arma::uvec
                arma::uvec incid_inds = Rcpp::as<arma::uvec>(event_inds);

                // get indices in the census_path matrix to keep incidence
                int n_census_events = incid_inds.n_elem;
                int incid_start = flow_matrix.n_cols + 1;

                // census the incidence increments
                for(int k = 1; k < n_census_times; ++k) {

                        for(int j = 0; j < n_census_events; ++j) {
                                census_path(k-1, incid_start + j) =
                                        arma::sum(path(arma::span(census_inds[k-1]+1, census_inds[k]), incid_inds[j]));
                        }

                }
        }

        // compute the prevalence if called for
        if(do_prevalence) {

              // initialize the state and increment vectors
              arma::rowvec state(n_comps);
              arma::rowvec increment(n_rates, arma::fill::zeros);

              state = parmat.submat(row0, initdist_inds);

              for(int k=1; k < n_census_times-1; ++k) {

                    // numbers of transitions
                    increment = arma::sum(path(arma::span(census_inds[k-1] + 1, census_inds[k]),
                                               arma::span(1, n_rates)), 0);

                    // new state
                    state += increment * flow_matrix;

                    // save state
                    census_path(k-1, arma::span(1, n_comps)) = state;

                    // apply forcings if called for - applied after censusing the path
                    if(forcing_inds[k]) {

                          // distribute the forcings proportionally to the compartment counts in the applicable states
                          for(int s=0; s < n_forcings; ++s) {
                                forcing_flow     = parmat(k, forcing_tcov_inds[s]);
                                forcing_distvec  = forcing_flow * normalise(forcings_out.col(s) % state, 1);
                                state           += forcing_transfers.slice(s) * forcing_distvec;
                          }
                    }
              }
        }
}
