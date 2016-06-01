#' Instatiate the C++ rate functions for a stochastic epidemic model and return
#' a vector of function pointers.
#'
#' @param rates list of rate functions
#' @param param_codes vector of parameter codes
#' @param compartment_codes vector of compartment codes
#' @param const_codes vector of constant names
#' @param tcovar_codes vector of time-varying covariate codes
#'
#' @return Two vector of strings that serve as function pointers. Also loads
#'   functions for forward simulation and a function for solving the system of
#'   ODEs into the global environment.
#' @export
parse_rates <- function(rates, param_codes, compartment_codes, const_codes, tcovar_codes) {

        rates_lumped   <- paste0("RATE", 1:length(rates),"_LUMPED")
        rates_unlumped <- paste0("RATE", 1:length(rates),"_UNLUMPED")

        arg_strings <- "const Rcpp::IntegerVector& state, const Rcpp::NumericVector& parameters"
        if(!is.null(const_codes))  arg_strings <- paste(arg_strings, "const Rcpp::NumericVector& constants", sep = ", ")
        if(!is.null(tcovar_codes)) arg_strings <- paste(arg_strings, "const Rcpp::NumericVector& tcovar", sep = ", ")

        cat("Compiling rate functions:")

        for(s in seq_along(rates)) {

                cat("RATE",s,sep="\n")

                code_lumped <- paste("// [[Rcpp::depends(RcppArmadillo)]]",
                                 "#include <RcppArmadillo.h>",
                                 "using namespace arma;",
                                 "using namespace Rcpp;",
                                 "// [[Rcpp::export]]",
                                 paste0("double ", rates_lumped[s],"(", arg_strings,") {"),
                                 paste0("double rate = ",rates[[s]]$lumped,";"),
                                 "return rate;}",
                                 sep = "\n"
                )

                code_unlumped <- paste("// [[Rcpp::depends(RcppArmadillo)]]",
                                     "#include <RcppArmadillo.h>",
                                     "using namespace arma;",
                                     "using namespace Rcpp;",
                                     "// [[Rcpp::export]]",
                                     paste0("double ", rates_unlumped[s],"(", arg_strings,") {"),
                                     paste0("double rate = ",rates[[s]]$unlumped,";"),
                                     "return rate;}",
                                     sep = "\n"
                )

                Rcpp::sourceCpp(code = code_lumped, env = globalenv())
                Rcpp::sourceCpp(code = code_unlumped, env = globalenv())
        }
}