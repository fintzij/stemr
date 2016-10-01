// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

//' Copy drift, residual, and diffusion processes into the lna state vector
//'
//' @param statevec vector to be copied into, i.e. the vectorised ODE state
//' @param driftvec armadillo vector for the drift process
//' @param residvec armadillo vector for the residual process
//' @param diffmat armadillo matrix for the diffusion process
//'
//' @return modifies the statevec in place
//' @export
// [[Rcpp::export]]
void procs2vec(Rcpp::NumericVector& statevec, arma::mat& drift, arma::mat& resid, Rcpp::NumericVector& diff, int ind) {

        int n_comps = drift.n_cols;
        Rcpp::IntegerVector arr_dims = diff.attr("dim");

        // instatiate the diffusion array
        arma::cube diffusion(diff.begin(), arr_dims[0], arr_dims[1], arr_dims[2], false);

        // make the copies
        for(int j=0; j < n_comps; ++j) {
                statevec[j]         = drift(ind, j);
                statevec[j+n_comps] = resid(ind, j);
        }

        std::copy(diffusion.slice(ind).begin(), diffusion.slice(ind).end(), statevec.begin() + 2*n_comps);
}

//' Copy the LNA ODE state into the drift, residual, and diffusion processes
//' vectors
//'
//' @param statevec vector to be copied into, i.e. the vectorised ODE state
//' @param driftvec armadillo vector for the drift process
//' @param residvec armadillo vector for the residual process
//' @param diffmat armadillo matrix for the diffusion process
//'
//' @return modifies the ode objects in place
//' @export
// [[Rcpp::export]]
void vec2procs(Rcpp::NumericVector& statevec, arma::mat& drift, arma::mat& resid, Rcpp::NumericVector& diff, int ind) {

        int n_comps = drift.n_cols;
        Rcpp::IntegerVector arr_dims = diff.attr("dim");

        // instatiate the diffusion array
        arma::cube diffusion(diff.begin(), arr_dims[0], arr_dims[1], arr_dims[2], false);

        // make the copies
        for(int j=0; j < n_comps; ++j) {
                drift(ind, j) = statevec[j];
                resid(ind, j) = statevec[j+n_comps];
        }

        std::copy(statevec.begin() + 2*n_comps, statevec.end(), diffusion.slice(ind).begin());
}

// //' Copy the LNA ODE state vector into the drift, residual, and diffusion
// //' processes matrices, imposing diagonal dominance for the diffusion matrix
// //'
// //' @param statevec vector to be copied into, i.e. the vectorised ODE state
// //' @param driftmat armadillo matrix for the drift process
// //' @param residmat armadillo matrix for the residual process
// //' @param diffarr armadillo cube for the diffusion process
// //' @param ind row index in the ODE matrices (slice in diffusion cube)
// //'
// //' @return modifies the ode objects in place
// //' @export
// // [[Rcpp::export]]
// void vec2procmats(Rcpp::NumericVector& statevec, arma::mat& driftmat, arma::mat& residmat, Rcpp::NumericVector& diffarr, int ind) {
//
//         int row_ind = ind - 1;
//         int n_comps = driftmat.n_cols;
//         Rcpp::IntegerVector arr_dims = diffarr.attr("dim");
//
//         // instatiate the diffusion array
//         arma::cube diffusion_array(diffarr.begin(), arr_dims[0], arr_dims[1], arr_dims[2], false);
//
//         // make the copies
//         for(int j=0; j < n_comps; ++j) {
//                 driftmat(row_ind, j) = statevec[j];
//                 residmat(row_ind, j) = statevec[j+n_comps];
//         }
//
//         std::copy(statevec.begin() + 2*n_comps, statevec.end(), diffusion_array.slice(row_ind).begin());
//
//         // enforce diagonal dominance by adding 1e-7
//         diffusion_array.slice(row_ind).diag() += 0.0000001;
// }
//
// //' Insert the lna step into the path matrix
// //'
// //' @param path matrix containing the full LNA path
// //' @param lna_step vector containing the sampled states
// //' @param ind row index
// //'
// //' @return modifies the path matrix in place
// //' @export
// // [[Rcpp::export]]
// void insert_lna_step(arma::mat& path, arma::rowvec& lna_step, int ind) {
//
//         int row_ind = ind-1;
//         path(row_ind, arma::span(1, path.n_cols-1)) = lna_step;
// }

// //' insert the next step of the residual process into the residual process matrix
// //'
// //' @param residual_process matrix containing the residual process
// //' @param lna_step vector containing the lna step
// //' @param drift_process matrix containing the drift process
// //' @param ind row index
// //'
// //' @return modifies the residual process matrix in place
// //' @export
// // [[Rcpp::export]]
// void insert_resids(arma::mat& residual_process, const arma::rowvec& lna_step, const arma::mat& drift_process, int ind) {
//
//         int row_ind = ind-1;
//         residual_process.row(row_ind) = lna_step - drift_process.row(row_ind);
//
// }

// //' Copy drift, residual, and diffusion processes into the lna state vector
// //'
// //' @param statevec vector to be copied into, i.e. the vectorised ODE state
// //' @param driftmat armadillo matrix for the drift process
// //' @param diffarr armadillo cube for the diffusion process
// //' @param ind row index in the ODE matrices (slice in diffusion cube)
// //' @param zero_diffusion should the diffusion process be zeroed out
// //'
// //' @return modifies the statevec in place
// //' @export
// // [[Rcpp::export]]
// void procmats2vec(Rcpp::NumericVector& statevec, arma::mat& driftmat, Rcpp::NumericVector& diffarr, int ind, bool zero_diffusion) {
//
//         int row_ind = ind - 1;
//         int n_comps = driftmat.n_cols;
//         Rcpp::IntegerVector arr_dims = diffarr.attr("dim");
//
//         // instatiate the diffusion array
//         arma::cube diffusion_array(diffarr.begin(), arr_dims[0], arr_dims[1], arr_dims[2], false);
//
//         // make the copies
//         for(int j=0; j < n_comps; ++j) {
//                 statevec[j] = driftmat(row_ind, j);
//         }
//
//         if(zero_diffusion) {
//                 std::fill(statevec.begin() + 2*n_comps, statevec.end(), 0.0);
//         } else {
//                 std::copy(diffusion_array.slice(row_ind).begin(), diffusion_array.slice(row_ind).end(), statevec.begin() + 2*n_comps);
//         }
// }
