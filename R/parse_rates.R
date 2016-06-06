#' Instatiate the C++ rate functions for a stochastic epidemic model and return
#' a vector of function pointers.
#'
#' @param rates list of rate functions
#' @param param_codes vector of parameter codes
#' @param compartment_codes vector of compartment codes
#' @param const_codes vector of constant names
#' @param tcovar_codes vector of time-varying covariate codes
#' @param allrates compile each individual rate function, defaults to FALSE
#'
#' @return Two vector of strings that serve as function pointers. Also loads
#'   functions for forward simulation and a function for solving the system of
#'   ODEs into the global environment.
#' @export
parse_rates <- function(rates, param_codes, compartment_codes, const_codes, tcovar_codes, allrates = FALSE) {

        rates_lumped   <- paste0("RATE", 1:length(rates),"_LUMPED")
        rates_unlumped <- paste0("RATE", 1:length(rates),"_UNLUMPED")

        arg_strings <- "const Rcpp::IntegerVector& state, const Rcpp::NumericVector& parameters"
        if(!is.null(const_codes))  arg_strings <- paste(arg_strings, "const Rcpp::NumericVector& constants", sep = ", ")
        if(!is.null(tcovar_codes)) arg_strings <- paste(arg_strings, "const Rcpp::NumericVector& tcovar", sep = ", ")

        if(allrates) {
                for(s in seq_along(rates)) {
                        cat("Compiling rate functions:")
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


        # compile function for updating elements a rate vector
        code_lumped <- code_lumped <- paste("// [[Rcpp::depends(RcppArmadillo)]]",
                             "#include <RcppArmadillo.h>",
                             "using namespace arma;",
                             "using namespace Rcpp;",
                             "// [[Rcpp::export(name = \".RATES_LUMPED\")]]",
                             paste0("void RATES_LUMPED(Rcpp::NumericVector& rates, const Rcpp::LogicalVector inds,", arg_strings, ") {"), sep = "\n ")
        code_unlumped <- paste("// [[Rcpp::depends(RcppArmadillo)]]",
                               "#include <RcppArmadillo.h>",
                               "using namespace arma;",
                               "using namespace Rcpp;",
                               "// [[Rcpp::export(name = \".RATES_UNLUMPED\")]]",
                               paste0("void RATES_UNLUMPED(Rcpp::NumericVector& rates, const Rcpp::LogicalVector inds,", arg_strings, ") {"), sep = "\n ")

        for(i in seq_along(rates_lumped)) {
                code_lumped <- paste(code_lumped,
                                     paste0("if(inds[",i-1,"]) rates[",i-1,"] = ", rates[[i]]$lumped,";"), sep = "\n ")
                code_unlumped <- paste(code_unlumped,
                                       paste0("if(inds[",i-1,"]) rates[",i-1,"] = ", rates[[i]]$unlumped,";"), sep = "\n ")
        }
        code_lumped <- paste(code_lumped, "}" ,sep = "\n "); code_unlumped <- paste(code_unlumped, "}", sep = "\n ")

        Rcpp::sourceCpp(code = code_lumped, env = globalenv())
        Rcpp::sourceCpp(code = code_unlumped, env = globalenv())

}