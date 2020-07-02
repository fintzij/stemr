#' Construct and compile the functions for proposing an ODE path, with
#' integration of the deterministic mean ODEs accomplished using the Boost
#' odeint library.
#'
#' @param ode_rates list containing the ODE rate functions, derivatives, and
#'   parameter codes
#' @param compile_ode if TRUE, code will be generated and compiled. If a
#'   character string for the name of a file that does not yet exist, code will
#'   be generated but not compiled. If the name of a file that exists in the
#'   current working directory, the code in the file will be compiled.
#' @param messages should messages be printed
#'
#' @return list containing the ODE pointers and calling code
#' @export
load_ode <- function(ode_rates, compile_ode, messages, atol, rtol, stepper) {

        ODE_XPtr = NULL
        ODE_set_params_XPtr = NULL
      
        if(is.logical(compile_ode) && compile_ode) {
                generate_code <- TRUE
                compile_code  <- TRUE
                load_file     <- FALSE
        } else {
                if(compile_ode %in% list.files()) {
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
                n_rates         <- length(ode_rates$hazards)
                n_params        <- length(ode_rates$ode_param_codes)

                # construct the body of the ode ODEs.
                # The first n_rates compartments are the odes for the hazard functions.
                drift_inds      <- seq_len(n_rates)-1

                # dxdt strings
                ODE_odes     <- paste("dxdt[", drift_inds, "] = ", ode_rates$hazards, ";",
                                        collapse = "\n", sep = "")

                # generate the stemr_ode functions that will actually be called
                ODE_integrator <- paste("void INTEGRATE_STEM_ODE(Rcpp::NumericVector& init, double start, double end, double step_size = 0.001) {",
                                        "std::copy(init.begin(), init.end(), odeintr::state.begin());",
                                        "odeint::integrate_adaptive(odeintr::stepper, odeintr::sys, odeintr::state, start, end, step_size);",
                                        "std::copy(odeintr::state.begin(), odeintr::state.end(), init.begin());",
                                        "}\n",
                                        "typedef void(*ode_ptr)(Rcpp::NumericVector& init, double start, double end, double step_size);",
                                        "// [[Rcpp::export]]",
                                        "Rcpp::XPtr<ode_ptr> ODE_XPtr() {",
                                        "return(Rcpp::XPtr<ode_ptr>(new ode_ptr(&INTEGRATE_STEM_ODE)));",
                                        "}", sep = "\n")

                # function to set the ODE parameters,
                # differs from the odeintr function in that it takes a numeric vector, not a std::vector
                param_setter   <- paste("void SET_ODE_PARAMS(Rcpp::NumericVector& p) {",
                                        "std::copy(p.begin(), p.end(), odeintr::pars.begin());",
                                        "}\n",
                                        "typedef void(*set_pars_ptr)(Rcpp::NumericVector& p);",
                                        "// [[Rcpp::export]]",
                                        "Rcpp::XPtr<set_pars_ptr> ODE_set_params_XPtr() {",
                                        "return(Rcpp::XPtr<set_pars_ptr>(new set_pars_ptr(&SET_ODE_PARAMS)));",
                                        "}",sep = "\n")

                # paste the ODE integrator and parameter setting functions together
                stemr_ODE_code <- paste(ODE_integrator, param_setter, sep = "\n \n")

                # get the code for the ODE ODEs
                ODE_code <- odeintr::compile_sys(name = "INTEGRATE_ODE",
                                                 sys = ODE_odes,
                                                 sys_dim = n_rates,
                                                 pars = n_params,
                                                 method = stepper,
                                                 rtol = rtol,
                                                 atol = atol,
                                                 headers = paste("// [[Rcpp::depends(RcppArmadillo)]]",
                                                                 "#include <RcppArmadillo.h>",
                                                                 "using namespace arma;",
                                                                 sep = "\n"),
                                                 compile = F) # get the C++ code

                # RcppArmadillo is included, so remove the include tag for Rcpp
                ODE_code <- gsub("#include <Rcpp.h>", "", ODE_code)

                # we don't need the odeintr functions exported to the R global environment,
                # so remove the export attributes
                ODE_code <- gsub("// \\[\\[Rcpp::export\\]\\]", "", ODE_code)

                # paste the stemr ODE ode code to the end of the odeintr generated code
                ODE_code <- paste(ODE_code, stemr_ODE_code, sep = "\n")

                if(is.character(compile_ode)) {
                        filename <- ifelse(substr(compile_ode, nchar(compile_ode)-3, nchar(compile_ode)) != ".txt",
                                           paste0(compile_ode,".txt"), compile_ode)
                        cat(ODE_code, file = filename)
                }
        }

        if(load_file) {
                ODE_code <- readChar(compile_ode, nchars = 1e6)
        }

        if(compile_code) {
                # compile the ODE code
                if(messages) print("Compiling ODE functions.")
                Rcpp::sourceCpp(code = ODE_code, 
                                # env = globalenv(),
                                rebuild = TRUE, 
                                verbose = FALSE)

                # get the ODE function pointers
                ode_pointer <- c(ode_ptr = ODE_XPtr(),
                                 set_ode_params_ptr = ODE_set_params_XPtr(),
                                 ODE_code = ODE_code)

                return(ode_pointer)
        }
}