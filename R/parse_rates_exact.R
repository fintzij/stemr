#' Instatiate the C++ rate functions for a stochastic epidemic model and return
#' a vector of function pointers.
#'
#' @param rates list of rate functions
#' @param compile_rates if TRUE, code will be generated and compiled. If a
#'   character string for the name of a file that does not yet exist, code will
#'   be generated but not compiled. If the name of a file that exists in the
#'   current working directory, the code in the file will be compiled.
#' @param messages logical; print a message that the rates are being compiled
#'
#' @return Two vector of strings that serve as function pointers.
#' @export
parse_rates_exact <- function(rates, compile_rates, messages = TRUE) {

        if(is.logical(compile_rates) && compile_rates) {
                generate_code <- TRUE
                compile_code  <- TRUE
                load_file     <- FALSE
        } else {
                files <- list.files()
                if(compile_rates %in% list.files()) {
                        generate_code <- FALSE
                        compile_code  <- TRUE
                        load_file     <- TRUE
                } else {
                        generate_code <- TRUE
                        compile_code  <- FALSE
                        load_file     <- FALSE
                }
        }

        if(generate_code) {
                if(messages) {
                      exp_warning <- FALSE
                      for(r in seq_along(rates)) {
                            exp_warning <- exp_warning || grepl("\\^", rates[[r]]$lumped) || grepl("\\^", rates[[r]]$unlumped)
                            if(exp_warning) warning("If there is an exponentiated term in a rate, check that the base and exponent are both enclosed in parentheses, e.g., (base)^(exponent). Ignore this warning if rates are correctly specified.")
                            if(exp_warning) break
                      }
                }

                # ensure powers are converted
                for(r in seq_along(rates)) {
                        if(!identical(rates[[r]]$unlumped, character(0))) rates[[r]]$unlumped <- sub_powers(rates[[r]]$unlumped)
                        if(!identical(rates[[r]]$lumped, character(0)))   rates[[r]]$lumped   <- sub_powers(rates[[r]]$lumped)
                }

                arg_strings <- "Rcpp::NumericVector& rates, const Rcpp::LogicalVector& inds, const arma::rowvec& state, const Rcpp::NumericVector& parameters, const Rcpp::NumericVector& constants, const arma::rowvec& tcovar"

                fcns_lumped <- vector("list", length = length(rates))
                fcns_unlumped <- vector("list", length = length(rates))

                for(i in seq_along(rates)) {
                        fcns_lumped[[i]] <- paste(paste0("if(inds[",i-1,"]) rates[",i-1,"] = ", rates[[i]]$lumped,";"), sep = "\n ")
                        fcns_unlumped[[i]] <- paste(paste0("if(inds[",i-1,"]) rates[",i-1,"] = ", rates[[i]]$unlumped,";"), sep = "\n ")
                }

                # generate lumped code
                fcns_lumped <- paste(unlist(fcns_lumped), collapse = "\n")
                code_lumped <- paste("// [[Rcpp::depends(RcppArmadillo)]]",
                                     "#include <RcppArmadillo.h>",
                                     "using namespace arma;",
                                     "using namespace Rcpp;",
                                     paste0("void RATES_LUMPED(",arg_strings,") {"),
                                     fcns_lumped,
                                     "}\n",
                                     paste0("typedef void(*ratefcn_ptr)(", arg_strings,");"),
                                     "// [[Rcpp::export]]",
                                     "Rcpp::XPtr<ratefcn_ptr> LUMPED_XPtr() {",
                                     "return(Rcpp::XPtr<ratefcn_ptr>(new ratefcn_ptr(&RATES_LUMPED)));",
                                     "}", sep = "\n")

                exact_code <- code_lumped

                # generate unlumped code
                unlumped_inds <- sapply(rates, function(x) !identical(x$unlumped, character(0)))

                if(sum(unlumped_inds) == length(rates)) {

                        fcns_unlumped <- paste(unlist(fcns_unlumped), collapse = "\n")

                        code_unlumped <- paste("// [[Rcpp::depends(RcppArmadillo)]]",
                                               "#include <RcppArmadillo.h>",
                                               "using namespace arma;",
                                               "using namespace Rcpp;",
                                               paste0("void RATES_UNLUMPED(",arg_strings,") {"),
                                               fcns_unlumped,
                                               "}\n",
                                               paste0("typedef void(*ratefcn_ptr)(", arg_strings,");"),
                                               "// [[Rcpp::export]]",
                                               "Rcpp::XPtr<ratefcn_ptr> UNLUMPED_XPtr() {",
                                               "return(Rcpp::XPtr<ratefcn_ptr>(new ratefcn_ptr(&RATES_UNLUMPED)));",
                                               "}", sep = "\n")

                        exact_code <- paste(exact_code, code_unlumped, sep = "\n\n")
                }

                if(is.character(compile_rates)) {
                        filename <- ifelse(substr(compile_rates, nchar(compile_rates)-3, nchar(compile_rates)) != ".txt",
                                           paste0(compile_rates,".txt"), compile_rates)
                        cat(exact_code, file = filename)
                }
        }

        if(load_file) {
                exact_code <- readChar(compile_rates, nchar = 1e6)
        }

        if(compile_code) {
                if(messages) {
                        print("Compiling rate functions.")
                }

                Rcpp::sourceCpp(code = exact_code, env = globalenv(), verbose = FALSE)

                rate_pointers <- c(lumped_ptr = LUMPED_XPtr())

                if(sum(unlumped_inds) == length(rates)) rate_pointers <- c(rate_pointers, unlumped_ptr = UNLUMPED_XPtr())

                return(rate_pointers)
        }
}