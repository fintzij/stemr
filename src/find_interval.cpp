// [[Rcpp::depends(Rcpp)]]
#include "stemr_types.h"
#include "stemr_utils.hpp"
#include <Rcpp.h>
using namespace Rcpp;

//' Given a vector of interval endpoints \code{breaks}, determine in which
//' intervals the elements of a vector \code{x} fall.
//'
//' @param x vector for whose elements the corresponding intervals are
//'   identified
//' @param breaks vector containing the elements
//' @param rightmost_closed logical; if true, the results for x[j]=breaks[N] is
//'   N-1.
//' @param all_inside logical; if true, 0 is mapped to 1, and N is mapped to N-1
//'
//' The rightmost interval is assumed to be closed. Compares to the behavior of
//' the \code{findInterval} function in \code{R}, when \code{rightmost.closed =
//'  TRUE}
//'
//' @return matrix containing the compartment counts at census times.
//' @export
// [[Rcpp::export]]
Rcpp::IntegerVector find_interval(Rcpp::NumericVector& x, Rcpp::NumericVector& breaks, bool rightmost_closed, bool all_inside) {

        int x_size = x.size();
        Rcpp::IntegerVector out(x_size);

        Rcpp::NumericVector::iterator x_it = x.begin(), x_end = x.end(), breaks_it = breaks.begin(), breaks_end = breaks.end();
        Rcpp::IntegerVector::iterator out_it = out.begin(), out_end = out.end();
        Rcpp::NumericVector::iterator ubound;

        if(!all_inside & !rightmost_closed) {
                for(; x_it != x_end; x_it++, out_it++) {
                        ubound = std::upper_bound(breaks_it, breaks_end, *x_it);
                        *out_it = std::distance(breaks_it,ubound);
                }

        } else if(all_inside) {
                int n_breaks = breaks.size();

                for(; x_it != x_end; x_it++, out_it++) {
                        ubound = std::upper_bound(breaks_it, breaks_end, *x_it);
                        *out_it = std::distance(breaks_it,ubound);

                        if(*out_it == n_breaks) {
                                *out_it -= 1;
                        } else if(*out_it == 0) {
                                *out_it += 1;
                        }
                }

        } else if(!all_inside & rightmost_closed) {
                double max_breaks = Rcpp::max(breaks);
                int n_ints = breaks.size() - 1;

                for(; x_it != x_end; x_it++, out_it++) {

                        if(*x_it == max_breaks) {
                                *out_it = n_ints;
                        } else {
                                ubound = std::upper_bound(breaks_it, breaks_end, *x_it);
                                *out_it = std::distance(breaks_it,ubound);
                        }
                }
        }

        return out;
}