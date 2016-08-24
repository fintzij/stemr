// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace arma;
using namespace Rcpp;

//' Produce samples from a multivariate normal density using the singular value
//' decomposition as in the mvtnorm R package.
//'
//' @param n number of samples
//' @param mu mean vector
//' @param sigma covariance matrix
//'
//' @export
// [[Rcpp::export]]
arma::mat mvrn(int n, const arma::vec& mu, const arma::mat& sigma) {

        // dimension of the system and output vectors
        int p = mu.n_elem;

        // objects for the SVD
        arma::vec s(p);
        arma::mat U(p,p);
        arma::mat V(p,p);

        // compute the SVD
        arma::svd(U, s, V, sigma);

        // check that the matrix is positive definite
        double maxval = s.max();
        if(!arma::all(s >= (-1.490116e-08 * sqrt(maxval)))) {
                Rcpp::stop("sigma is not numerically positive definite");
        }

        // compute the matrix square root
        arma::mat R = U * arma::diagmat(arma::sqrt(arma::clamp(s, 0.0, maxval))) * V;

        // simulate the data
        mat X = arma::randn(n, p) * R;
        X.each_row() += mu.t();

        return X;
}
