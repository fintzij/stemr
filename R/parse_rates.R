#' Instatiate the C++ rate functions for a stochastic epidemic model and return
#' a vector of function pointers.
#'
#' @param rates list of rate functions
#' @param messages logical; print a message that the rates are being compiled
#'
#' @return Two vector of strings that serve as function pointers.
#' @export
parse_rates <- function(rates, messages = TRUE) {

        arg_strings <- "Rcpp::NumericVector& rates, const Rcpp::LogicalVector& inds, const arma::rowvec& state, const Rcpp::NumericVector& parameters, const Rcpp::NumericVector& constants, const arma::rowvec& tcovar"

        fcns_lumped <- fcns_unlumped <- character(0)

        for(i in seq_along(rates)) {
                fcns_lumped <- paste(fcns_lumped,paste0("if(inds[",i-1,"]) rates[",i-1,"] = ", rates[[i]]$lumped,";"), sep = "\n ")
                fcns_unlumped <- paste(fcns_unlumped,paste0("if(inds[",i-1,"]) rates[",i-1,"] = ", rates[[i]]$unlumped,";"), sep = "\n ")
        }

        # compile function for updating elements a rate vector
        code_lumped <- paste("// [[Rcpp::depends(RcppArmadillo)]]",
                           "#include <RcppArmadillo.h>",
                           "using namespace arma;",
                           "using namespace Rcpp;",
                           paste0("void RATES_LUMPED(",arg_strings,") {"),
                           fcns_lumped,
                           "}",
                           paste0("typedef void(*ratefcn_ptr)(", arg_strings,");"),
                           "// [[Rcpp::export]]",
                           "Rcpp::XPtr<ratefcn_ptr> LUMPED_XPtr() {",
                           "return(Rcpp::XPtr<ratefcn_ptr>(new ratefcn_ptr(&RATES_LUMPED)));",
                           "}", sep = "\n")

        code_unlumped <- paste("// [[Rcpp::depends(RcppArmadillo)]]",
                             "#include <RcppArmadillo.h>",
                             "using namespace arma;",
                             "using namespace Rcpp;",
                             paste0("void RATES_UNLUMPED(",arg_strings,") {"),
                             fcns_unlumped,
                             "}",
                             paste0("typedef void(*ratefcn_ptr)(", arg_strings,");"),
                             "// [[Rcpp::export]]",
                             "Rcpp::XPtr<ratefcn_ptr> UNLUMPED_XPtr() {",
                             "return(Rcpp::XPtr<ratefcn_ptr>(new ratefcn_ptr(&RATES_UNLUMPED)));",
                             "}", sep = "\n")

        if(messages) {
                print("Compiling rate functions.")
        }

        Rcpp::sourceCpp(code = code_lumped, env = globalenv())
        Rcpp::sourceCpp(code = code_unlumped, env = globalenv())

        rate_pointers <- c(lumped_ptr = LUMPED_XPtr(), unlumped_ptr = UNLUMPED_XPtr())

        return(rate_pointers)
}