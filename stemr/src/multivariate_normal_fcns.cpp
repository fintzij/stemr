// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace arma;
using namespace Rcpp;

const double log2pi = std::log(2.0 * M_PI);

//' Produce samples from a multivariate normal density using the Cholesky
//' decomposition
//'
//' @param n number of samples
//' @param mu mean vector
//' @param sigma covariance matrix
//'
//' (source: http://gallery.rcpp.org/articles/simulate-multivariate-normal/)
//'
//' @export
// [[Rcpp::export]]
arma::mat rmvtn(int n, const arma::rowvec& mu, const arma::mat& sigma) {

        // dimension of the system and output vectors
        int p = mu.n_elem;

        // generate independent standard normal RVs
        Rcpp::NumericVector draws = Rcpp::rnorm(n*p);
        arma::mat X(draws.begin(), n, p, false, true);

        // add the mean and multiply by the upper triangular portion of the cholesky
        return arma::repmat(mu, n, 1) + X * arma::chol(sigma, "upper");
}

//' Multivariate normal density
//'
//' @param x matrix of draws for which to evaluate the density
//' @param mu mean vector of the distribution
//' @param sigma covariance matrix
//' @param logd should the log be returned
//'
//' (source: http://gallery.rcpp.org/articles/dmvnorm_arma/)
//'
//' @export
// [[Rcpp::export]]
arma::vec dmvtn(const arma::mat& x, const arma::rowvec& mu, const arma::mat& sigma, bool logd = false) {

        int n = x.n_rows;
        int xdim = x.n_cols;
        arma::vec out(n);
        arma::mat rooti = arma::trans(arma::inv(arma::chol(sigma, "upper")));
        double rootisum = arma::sum(log(rooti.diag()));
        double constants = -(static_cast<double>(xdim)/2.0) * log2pi;

        for (int i=0; i < n; i++) {
                arma::vec z = rooti * arma::trans(x.row(i) - mu) ;
                out(i)      = constants - 0.5 * arma::sum(z%z) + rootisum;
        }

        if (logd == false) {
                out = exp(out);
        }

        return out;
}
