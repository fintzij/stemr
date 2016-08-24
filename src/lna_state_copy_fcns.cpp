// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

//' Copy drift, residual, and diffusion processes into std::vector
//'
//' @param statevec vector to be copied into, i.e. the vectorised ODE state
//' @param drift_proc armadillo vector for the drift process
//' @param resid_proc armadillo vector for the residual process
//' @param diffusion_proc armadillo matrix for the diffusion process
//'
//' @return modifies the statevec in place
//' @export
// [[Rcpp::export]]
void procs2vec(Rcpp::NumericVector& statevec, arma::vec& driftvec, arma::vec& residvec, arma::mat& diffmat) {

        int n_comps  = driftvec.n_elem;

        // make the copies
        std::copy(driftvec.begin(), driftvec.end(), statevec.begin());
        std::copy(residvec.begin(), residvec.end(), statevec.begin() + n_comps);
        std::copy(diffmat.begin(), diffmat.end(), statevec.begin() + 2*n_comps);

}

//' Copy the LNA ODE state into the drift, residual, and diffusion processes
//' vectors
//'
//' @param statevec vector to be copied into, i.e. the vectorised ODE state
//' @param driftvec armadillo vector for the drift process
//' @param residvec armadillo vector for the residual process
//' @param diffmat armadillo matrix for the diffusion process
//'
//' @return modifies the statevec in place
//' @export
// [[Rcpp::export]]
void vec2procs(Rcpp::NumericVector& statevec, arma::vec& driftvec, arma::vec& residvec, arma::mat& diffmat) {

        int n_comps  = driftvec.n_elem;

        // make the copies
        std::copy(statevec.begin(), statevec.begin() + n_comps, driftvec.begin());
        std::copy(statevec.begin() + n_comps, statevec.begin() + 2*n_comps, residvec.begin());
        std::copy(statevec.begin() + 2*n_comps, statevec.end(), diffmat.begin());

}