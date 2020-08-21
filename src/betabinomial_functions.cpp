#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::plugins(cpp11)]]

// Beta-binomial functions, pilfered from the extraDistr package until that 
// package exposes C++ code to make it callable from C++.
// 
// Citation: 
//
// "Tymoteusz Wolodzko (2019). extraDistr: Additional Univariate and Multivariate
// Distributions. R package version 1.8.11.
// https://CRAN.R-project.org/package=extraDistr"

// [[Rcpp::export]]
Rcpp::NumericVector dbbinom(
        const Rcpp::NumericVector& x,
        const Rcpp::NumericVector& size,
        const Rcpp::NumericVector& alpha,
        const Rcpp::NumericVector& beta,
        const bool& log_prob = false) {
    
    // length of the returned pmf vector
    int Nmax = std::max({
        x.length(),
        size.length(),
        alpha.length(),
        beta.length()
    });

    // initialize vector
    Rcpp::NumericVector p(Nmax);
    
    // compute lpmfs
    for (int i = 0; i < Nmax; ++i) {
        p[i] = 
            R::lchoose(size[i % size.length()], x[i % x.length()]) + 
            R::lbeta(x[i % x.length()] + alpha[i % alpha.length()], 
                     size[i % size.length()] - x[i % x.length()] + beta[i % beta.length()]) - 
            R::lbeta(alpha[i % alpha.length()], 
                     beta[i % beta.length()]);
    }
    
    // exponentiate if needed
    if (!log_prob)
        p = Rcpp::exp(p);
    
    return p;
}

// [[Rcpp::export]]
Rcpp::NumericVector rbbinom(
        const int& n,
        const NumericVector& size,
        const NumericVector& alpha,
        const NumericVector& beta) {
    
    Rcpp::NumericVector x(n);
    double prob;
    
    for (int i = 0; i < n; i++) {
        // sample probability
        prob = R::rbeta(alpha[i % alpha.length()], beta[i % beta.length()]);
        
        // sample count
        x[i] = R::rbinom(size[i % size.length()], prob);   
    }
    
    return x;
}
