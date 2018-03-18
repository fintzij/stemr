#' Construct and compile the functions for proposing an LNA path, with
#' integration of the LNA ODEs accomplished using the Boost odeint library.
#'
#' @param lna_rates list containing the LNA rate functions, derivatives, and
#'   parameter codes
#' @param compile_lna if TRUE, code will be generated and compiled. If a
#'   character string for the name of a file that does not yet exist, code will
#'   be generated but not compiled. If the name of a file that exists in the
#'   current working directory, the code in the file will be compiled.
#' @param messages should messages be printed
#'
#' @return list containing the LNA pointers and calling code
#' @export
load_lna <- function(lna_rates, compile_lna, messages, atol, rtol, stepper) {

        if(is.logical(compile_lna) && compile_lna) {
                generate_code <- TRUE
                compile_code  <- TRUE
                load_file     <- FALSE
        } else {
                files <- list.files()
                if(compile_lna %in% list.files()) {
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
                # get the number of rates and the number of compartments
                n_rates         <- length(lna_rates$lna_rates)
                n_params        <- length(lna_rates$lna_param_codes)
                n_odes          <- n_rates + n_rates^2
                init_state_ind  <- grep(pattern = "_0", names(lna_rates$lna_param_codes))[1] - 1 # num. params before initial comp counts

                # construct the body of the lna ODEs.
                # The first n_rates compartments are the odes for the hazard functions.
                diff_start      <- n_rates
                jacobian_inds   <- matrix(c(rep(seq(0, n_rates - 1), each = n_rates),
                                            rep(seq(0, n_rates - 1), n_rates)), ncol = 2)
                drift_inds      <- seq_len(n_rates)-1
                diffusion_inds  <- seq(n_rates, n_odes-1, by = 1)

                # strings to construct the drift and diffusion vectors, and to exponentiate the current state
                # exponentiate the current state
                exp_Z_terms     <- paste(paste0("odeintr::Z = arma::vec(x).subvec(0,",n_rates-1,");"),
                                         "Z.elem(arma::find(Z<0)).zeros();", # ensures compartment counts are nonnegative
                                         "odeintr::exp_Z = arma::exp(odeintr::Z);",
                                         "odeintr::exp_neg_Z = arma::exp(-odeintr::Z);",
                                         "odeintr::exp_neg_2Z = arma::square(odeintr::exp_neg_Z);", sep = "\n")

                # strings to compute the ito terms, hazards, drift, and jacobian
                haz_terms      <- paste(paste("odeintr::hazards[",0:(n_rates-1),"]", " = ",
                                              lna_rates$lna_rates,";", sep = ""), collapse = "\n")

                non_zero_inds  <- which(lna_rates$derivatives != "0")
                jacobian_terms <- paste(paste("odeintr::jacobian(",
                                              jacobian_inds[non_zero_inds,1], ", ", jacobian_inds[non_zero_inds,2], ") = ",
                                              lna_rates$derivatives[non_zero_inds],";", sep = ""), collapse = "\n")

                # diffusion_ode  <- paste0("odeintr::diffusion_ode = arma::vectorise(odeintr::diffusion * odeintr::jacobian.t() + ",
                #                          "arma::diagmat(odeintr::exp_neg_2Z % odeintr::hazards) + ",
                #                          "odeintr::jacobian * odeintr::diffusion, 0);")

                diffusion_terms <- paste(paste0("odeintr::diffusion = arma::reshape(arma::vec(x).subvec(",
                                                n_rates, ",", n_odes-1, "),", n_rates,",", n_rates,");"),
                                         paste0("odeintr::diffusion = odeintr::diffusion * odeintr::jacobian.t() + ",
                                              "odeintr::jacobian * odeintr::diffusion;"),
                                         "odeintr::diffusion.diag() += odeintr::exp_neg_2Z % odeintr::hazards;",sep = "\n")

                # dxdt strings
                dxdt_drift     <- paste("dxdt[", drift_inds, "] = ",
                                        paste0(lna_rates$ito_coefs,"*odeintr::hazards[", seq_along(drift_inds) - 1, "];"),
                                        collapse = "\n", sep = "")
                dxdt_diffusion <- paste("dxdt[", diffusion_inds, "] = odeintr::diffusion[", seq_along(diffusion_inds)-1, "];",
                                        collapse = "\n", sep = "")

                # concatenate everything
                LNA_odes <- paste(exp_Z_terms, haz_terms, jacobian_terms,
                                  diffusion_terms, dxdt_drift, dxdt_diffusion, sep = "\n\n")

                # generate the stemr_lna functions that will actually be called
                LNA_integrator <- paste("void INTEGRATE_STEM_LNA(Rcpp::NumericVector& init, double start, double end, double step_size = 0.001) {",
                                        "std::copy(init.begin(), init.end(), odeintr::state.begin());",
                                        "odeint::integrate_adaptive(odeintr::stepper, odeintr::sys, odeintr::state, start, end, step_size);",
                                        "std::copy(odeintr::state.begin(), odeintr::state.end(), init.begin());",
                                        "}\n",
                                        "typedef void(*ode_ptr)(Rcpp::NumericVector& init, double start, double end, double step_size);",
                                        "// [[Rcpp::export]]",
                                        "Rcpp::XPtr<ode_ptr> LNA_XPtr() {",
                                        "return(Rcpp::XPtr<ode_ptr>(new ode_ptr(&INTEGRATE_STEM_LNA)));",
                                        "}", sep = "\n")

                # function to set the LNA parameters,
                # differs from the odeintr function in that it takes a numeric vector, not a std::vector
                param_setter   <- paste("void SET_LNA_PARAMS(Rcpp::NumericVector& p) {",
                                        "std::copy(p.begin(), p.end(), odeintr::pars.begin());",
                                        "}\n",
                                        "typedef void(*set_pars_ptr)(Rcpp::NumericVector& p);",
                                        "// [[Rcpp::export]]",
                                        "Rcpp::XPtr<set_pars_ptr> LNA_set_params_XPtr() {",
                                        "return(Rcpp::XPtr<set_pars_ptr>(new set_pars_ptr(&SET_LNA_PARAMS)));",
                                        "}",sep = "\n")

                # paste the LNA integrator and parameter setting functions together
                stemr_LNA_code <- paste(LNA_integrator, param_setter, sep = "\n \n")

                # get the code for the LNA ODEs
                LNA_code <- odeintr::compile_sys(name = "INTEGRATE_LNA",
                                                 sys = LNA_odes,
                                                 sys_dim = n_rates + n_rates^2,
                                                 pars = n_params,
                                                 method = stepper,
                                                 rtol = rtol,
                                                 atol = atol,
                                                 globals = paste(paste(
                                                         "\n",
                                                         paste0("static arma::vec Z(", n_rates,",arma::fill::zeros);"),
                                                         paste0("static arma::vec exp_Z(", n_rates,",arma::fill::zeros);"),
                                                         paste0("static arma::vec exp_neg_Z(", n_rates,",arma::fill::zeros);"),
                                                         paste0("static arma::vec exp_neg_2Z(", n_rates,",arma::fill::zeros);"),
                                                         paste0("static arma::vec hazards(",n_rates,",arma::fill::zeros);"),
                                                         paste0("static arma::mat jacobian(", n_rates,",",n_rates, ",arma::fill::zeros);"), sep = "\n"),
                                                         paste0("static arma::mat diffusion(", n_rates,",",n_rates,",arma::fill::zeros);"),
                                                         sep = "\n"),
                                                 headers = paste("// [[Rcpp::depends(RcppArmadillo)]]",
                                                                 "#include <RcppArmadillo.h>",
                                                                 "using namespace arma;",
                                                                 sep = "\n"),
                                                 compile = F) # get the C++ code

                # RcppArmadillo is included, so remove the include tag for Rcpp
                LNA_code <- gsub("#include <Rcpp.h>", "", LNA_code)

                # we don't need the odeintr functions exported to the R global environment,
                # so remove the export attributes
                LNA_code <- gsub("// \\[\\[Rcpp::export\\]\\]", "", LNA_code)

                # paste the stemr LNA ode code to the end of the odeintr generated code
                LNA_code <- paste(LNA_code, stemr_LNA_code, sep = "\n")

                if(is.character(compile_lna)) {
                        filename <- ifelse(substr(compile_lna, nchar(compile_lna)-3, nchar(compile_lna)) != ".txt",
                                           paste0(compile_lna,".txt"), compile_lna)
                        cat(LNA_code, file = filename)
                }
        }

        if(load_file) {
                LNA_code <- readChar(compile_lna, nchars = 1e6)
        }

        if(compile_code) {
                # compile the LNA code
                if(messages) print("Compiling LNA functions.")
                Rcpp::sourceCpp(code = LNA_code, env = globalenv(), verbose = FALSE)

                # get the LNA function pointers
                lna_pointer <- c(lna_ptr = LNA_XPtr(),
                                 set_lna_params_ptr = LNA_set_params_XPtr(),
                                 LNA_code = LNA_code)

                return(lna_pointer)
        }
}