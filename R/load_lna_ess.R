#' Construct and compile the functions for proposing an LNA path within the
#' elliptical slice sampling algorithm, with integration of the LNA ODEs
#' accomplished using the Boost odeint library.
#'
#' @param lna_rates list containing the LNA rate functions, derivatives, and
#'   parameter codes
#' @param messages should messages be printed
#'
#' @return list containing the LNA pointers and calling code
#' @export
load_lna_ess <- function(lna_rates, messages, atol, rtol) {

        # get the number of rates and the number of compartments
        n_rates         <- length(lna_rates$rates)
        n_compartments  <- length(lna_rates$rates)
        n_params        <- length(lna_rates$lna_param_codes)
        n_odes          <- 2*n_compartments

        # construct the body of the lna ODEs. The first n_rates compartments are
        # the odes driven by the hazard functions.
        resid_start   <- n_compartments
        resid_end     <- 2*n_compartments - 1
        jacobian_inds <- matrix(c(rep(seq(0, n_compartments - 1), each = n_compartments),
                                  rep(seq(0, n_compartments - 1), n_compartments)), ncol = 2)

        drift_inds    <- seq_len(n_compartments)-1
        resid_inds    <- drift_inds + n_compartments

        # string to construct the drift, residual, and diffusion vectors
        resid_terms     <- paste0("odeintr::resid = arma::vec(x).subvec(", resid_start,",", resid_end,");")

        # strings to compute the rates, hazards (cumulative hazards = A'phi) and jacobian
        rate_terms     <- paste(paste("odeintr::rates[",0:(n_rates-1),"]", " = ", lna_rates$rates,";", sep = ""),
                                collapse = "\n")
        jacobian_terms <- paste(paste("odeintr::jacobian(", jacobian_inds[,1], ", ", jacobian_inds[,2], ") = ",
                                      lna_rates$derivatives,";", sep = ""), collapse = "\n")

        drift_ode      <- "odeintr::drift_ode = odeintr::rates;"
        resid_ode      <- "odeintr::resid_ode = odeintr::jacobian * odeintr::resid;"

        # dxdt strings
        dxdt_drift     <- paste("dxdt[", drift_inds, "] = odeintr::drift_ode[", seq_along(drift_inds) - 1, "];",
                                collapse = "\n", sep = "")
        dxdt_resid     <- paste("dxdt[", resid_inds, "] = odeintr::resid_ode[", seq_along(resid_inds) - 1, "];",
                                collapse = "\n", sep = "")

        # concatenate everything
        LNA <- paste(resid_terms, rate_terms, jacobian_terms,
                     drift_ode, resid_ode, dxdt_drift, dxdt_resid, sep = "\n")

        # generate the stemr_lna functions that will actually be called
        LNA_integrator <- paste("void INTEGRATE_LNA_ESS(Rcpp::NumericVector& init, double start, double end, double step_size=1.0) {",
                                "LNA_ESS_set_state(init);",
                                "odeint::integrate_adaptive(odeintr::stepper, odeintr::sys, odeintr::state, start, end, step_size);",
                                "init = Rcpp::wrap(LNA_ESS_get_state());",
                                "}\n",
                                "typedef void(*lna_ptr)(Rcpp::NumericVector& init, double start, double end, double step_size);",
                                "// [[Rcpp::export]]",
                                "Rcpp::XPtr<lna_ptr> LNA_XPtr() {",
                                "return(Rcpp::XPtr<lna_ptr>(new lna_ptr(&INTEGRATE_LNA_ESS)));",
                                "}", sep = "\n")

        param_setter   <- paste("void SET_LNA_ESS_PARAMS(Rcpp::NumericVector& p) {",
                                "std::copy(p.begin(), p.end(), odeintr::pars.begin());",
                                "}\n",
                                "typedef void(*set_pars_ptr)(Rcpp::NumericVector& p);",
                                "// [[Rcpp::export]]",
                                "Rcpp::XPtr<set_pars_ptr> LNA_ESS_set_params_XPtr() {",
                                "return(Rcpp::XPtr<set_pars_ptr>(new set_pars_ptr(&SET_LNA_ESS_PARAMS)));",
                                "}",sep = "\n")

        # paste the LNA integrator and parameter setting functions together
        stemr_LNA_ESS_code <- paste(LNA_integrator, param_setter, sep = "\n \n")

        # get the code for the LNA ODEs
        LNA_ESS_code <- odeintr::compile_sys(name = "LNA_ESS",
                                             sys = LNA,
                                             sys_dim = 2*n_compartments,
                                             pars = n_params,
                                             method = "rk78_a",
                                             rtol = rtol,
                                             atol = atol,
                                             globals = paste(
                                                     paste(paste0("static arma::vec resid(", n_compartments,");"),
                                                           paste0("static arma::vec drift_ode(", n_compartments,");"),
                                                           paste0("static arma::vec resid_ode(", n_compartments,");"),
                                                           paste0("static arma::vec rates(",n_rates,");"),
                                                           paste0("static arma::mat jacobian(", n_compartments,",",n_compartments, ");"),
                                                           sep = "\n"), sep = "\n"),
                                             headers = paste("// [[Rcpp::depends(RcppArmadillo)]]",
                                                             "#include <RcppArmadillo.h>",
                                                             "using namespace arma;",
                                                             sep = "\n"),
                                             compile = F) # get the C++ code

        # RcppArmadillo is included, so remove the include tag for Rcpp
        LNA_ESS_code <- gsub("#include <Rcpp.h>", "", LNA_ESS_code)

        # we don't need the odeintr functions exported to the R global environment,
        # so remove the export attributes
        LNA_ESS_code <- gsub("// \\[\\[Rcpp::export\\]\\]", "", LNA_ESS_code)

        # paste the stemr LNA ode code to the end of the odeintr generated code
        LNA_ESS_code <- paste(LNA_ESS_code, stemr_LNA_ESS_code, sep = "\n")

        # compile the LNA code
        if(messages) print("Compiling LNA inference functions.")
        Rcpp::sourceCpp(code = LNA_ESS_code, env = globalenv())

        # get the LNA function pointers
        lna_pointer <- c(lna_ess_ptr = LNA_XPtr(), set_lna_ess_params_ptr = LNA_ESS_set_params_XPtr(), LNA_ess_code = LNA_ESS_code)

        return(lna_pointer)

}