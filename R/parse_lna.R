#' Instatiate the C++ rate functions for the LNA.
#'
#' @param rates list of rate functions
#' @param rate_derivs character vector of partial derivatives of the rate
#'   functions.
#' @param messages logical; print a message that the rates are being compiled
#'
#' @return Vector of strings that serve as function pointers.
#' @export
parse_lna <- function(rates, hazard_fcns, rate_derivs, messages = TRUE) {

        # get the number of rates and the number of compartments
        n_rates         <- length(hazard_fcns)
        n_compartments  <- length(rate_derivs) / n_rates

        # first compile the function to compute the Jacobian
        jacobian_args <- hazard_args <- "double t, const arma::vec& state, const Rcpp::NumericVector& parameters, const Rcpp::NumericVector& constants, const Rcpp::NumericVector& tcovar"
        jacobian_inds <- matrix(c(rep(seq(0, n_rates - 1), each = n_compartments), rep(seq(0, n_compartments - 1), n_rates)), ncol = 2)
        hazard_fcn    <- paste0("arma::vec hazards(", n_rates,");")
        jacobian_fcn  <- paste0("arma::mat jacobian(", n_rates, ", ", n_compartments, ");")

        for(j in seq_along(hazard_fcns)) {
                hazard_fcn <- paste(hazard_fcn, paste0("hazards[", j-1, "] = ", hazard_fcns[j],";"), sep = "\n")
        }

        for(i in seq_along(rate_derivs)) {
                if(!identical(rate_derivs[i], "0")) {
                        jacobian_fcn <- paste(jacobian_fcn, paste0("jacobian(", jacobian_inds[i,1],",",jacobian_inds[i,2],") = ",
                                                                   rate_derivs[i],";"), sep = "\n")
                }
        }

        hazard_fcn   <- paste(hazard_fcn, "return hazards;", sep = "\n")
        jacobian_fcn <- paste(jacobian_fcn, "return jacobian;", sep = "\n")

        # compile function for computing the Jacobian
        hazard_code <- paste("// [[Rcpp::depends(RcppArmadillo)]]",
                             "#include <RcppArmadillo.h>",
                             "using namespace arma;",
                             "using namespace Rcpp;",
                             paste0("arma::vec HAZARD(",hazard_args,") {"),
                             hazard_fcn,
                             "}",
                             paste0("typedef arma::vec(*hazard_ptr)(", hazard_args,");"),
                             "// [[Rcpp::export]]",
                             "Rcpp::XPtr<hazard_ptr> HAZARD_XPtr() {",
                             "return(Rcpp::XPtr<hazard_ptr>(new hazard_ptr(&HAZARD)));",
                             "}", sep = "\n")

        jacobian_code <- paste("// [[Rcpp::depends(RcppArmadillo)]]",
                               "#include <RcppArmadillo.h>",
                               "using namespace arma;",
                               "using namespace Rcpp;",
                               paste0("arma::mat JACOBIAN(",jacobian_args,") {"),
                               jacobian_fcn,
                               "}",
                               paste0("typedef arma::mat(*jacobian_ptr)(", jacobian_args,");"),
                               "// [[Rcpp::export]]",
                               "Rcpp::XPtr<jacobian_ptr> JACOBIAN_XPtr() {",
                               "return(Rcpp::XPtr<jacobian_ptr>(new jacobian_ptr(&JACOBIAN)));",
                               "}", sep = "\n")

        if(messages) {
                print("Compiling LNA functions.")
        }

        Rcpp::sourceCpp(code = hazard_code, env = globalenv())
        Rcpp::sourceCpp(code = jacobian_code, env = globalenv())

        lna_pointers <- c(hazard_ptr = HAZARD_XPtr(), jacobian_ptr = JACOBIAN_XPtr())

        return(lna_pointers)
}