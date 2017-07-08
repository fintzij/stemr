#include <Rcpp.h>

using namespace Rcpp;

// [[Rcpp::export]]
double takeLog(double val) {
        try {
                if (val <= 0.0) {         	// log() not defined here
                        throw std::range_error("Inadmissible value");
                }
                return log(val);
        } catch(std::exception &ex) {
                forward_exception_to_r(ex);
        } catch(...) {
                ::Rf_error("c++ exception (unknown reason)");
        }
        return NA_REAL;             // not reached
}