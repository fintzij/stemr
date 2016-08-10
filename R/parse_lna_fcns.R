#' Instatiate the C++ rate functions for the LNA using the odeintr package.
#'
#' @param lna_rates list containing the strings vectors of ODES for the hazards
#'   and jacobian matrix, along with the lna parameter codes.
#' @param flow_matrix
#' @param messages logical; print a message that the rates are being compiled
#'
#' @return External function pointer to the _no_record LNA integrator
#' @export
parse_lna_fcns <- function(lna_rates, flow_matrix, messages = TRUE) {

        # get the number of rates and the number of compartments
        n_rates         <- nrow(flow_matrix)
        n_compartments  <- ncol(flow_matrix)
        n_params        <- length(lna_rates$lna_param_codes)

        # construct the stochiometry matrix constructor
        stoich_mtx <- t(flow_matrix)
        stoich_mtx_code <- paste("static arma::mat stoich = {",paste("{", apply(stoich_mtx, 1, paste, collapse = ", "), "}", collapse = ",\n "),"};")

        # construct the body of the lna ODEs. The first n_rates compartments are
        # the odes driven by the hazard functions. The next n_rates x n_compartment
        # rates are the jacobian odes
        jacobian_inds <- matrix(c(rep(seq(0, n_rates - 1), each = n_compartments), rep(seq(0, n_compartments - 1), n_rates)), ncol = 2)
        drift_inds    <- 0:(n_compartments-1)
        resid_inds    <- n_compartments:(2*n_compartments-1)
        diffusion_inds<- (2*n_compartments):(2*n_compartments + n_compartments^2 - 1)
        diff_mat_inds <- matrix(c(rep(seq(0, n_compartments - 1), n_compartments), rep(seq(0, n_compartments - 1), each = n_compartments)), ncol = 2)

        # string to construct the drift, residual, and diffusion vectors
        drift_terms     <- paste0("odeintr::drift[", drift_inds,"] = x[",drift_inds,"];", collapse = "\n")
        resid_terms     <- paste0("odeintr::resid[", resid_inds,"] = x[",resid_inds,"];", collapse = "\n")
        diffusion_terms <- paste0("odeintr::diffusion(", diff_mat_inds[,1],",",diff_mat_inds[,2],") = x[",diffusion_inds,"];", collapse = "\n")

        # strings to compute the hazards and jacobian
        hazard_terms    <- paste(paste("odeintr::hazard[",0:(n_rates-1),"]", " = ", lna_rates$hazards, ";", sep = ""), collapse = "\n")
        jacobian_terms  <- paste(paste("odeintr::jacobian(", jacobian_inds[,1], ", ", jacobian_inds[,2], ") = ",lna_rates$derivatives,";", sep = ""),
                                 collapse = "\n")

        # strings to compute the LNA odes
        drift_ode       <- "odeintr::drift_ode = odeintr::stoich * odeintr::hazard;"
        resid_ode       <- "odeintr::resid_ode = odeintr::stoich * odeintr::jacobian * odeintr::resid;"
        diffusion_ode   <- paste0("odeintr::diffusion_ode = arma::vectorise(odeintr::diffusion * odeintr::jacobian.t() * odeintr::stoich.t() + ",
                                  "odeintr::stoich * arma::diagmat(odeintr::hazard) * odeintr::stoich.t() + ",
                                  "odeintr::stoich * odeintr::jacobian * odeintr::diffusion);")

        # dxdt strings
        dxdt_drift      <- paste("dxdt[", drift_inds, "] = odeintr::drift_ode[", drift_inds, "];", collapse = "\n", sep = "")
        dxdt_resid      <- paste("dxdt[", resid_inds, "] = odeintr::resid_ode[", resid_inds, "];", collapse = "\n", sep = "")
        dxdt_diffusion  <- paste("dxdt[", diffusion_inds, "] = odeintr::diffusion_ode[", diffusion_inds, "];", collapse = "\n", sep = "")

        # concatenate everything
        LNA_odes <- paste(drift_terms, resid_terms, diffusion_terms,
                          hazard_terms, jacobian_terms,
                          drift_ode, resid_ode, diffusion_ode,
                          dxdt_drift, dxdt_resid, dxdt_diffusion, sep = "\n")

        # generate the stemr_lna functions that will actually be called
        # stemr_LNA_integrator is equivalent to the _no_record odeintr function, but takes the endpoints of the interval and returns a numeric vector.
        LNA_integrator <- paste("Rcpp::NumericVector INTEGRATE_STEM_LNA(Rcpp::NumericVector init, double start, double end, double step_size=1.0) {",
                                      "INTEGRATE_LNA_set_state(init);",
                                      "odeint::integrate_adaptive(odeintr::stepper, odeintr::sys, odeintr::state, start, end, step_size);",
                                      "return Rcpp::wrap(INTEGRATE_LNA_get_state());",
                                      "}\n",
                                "typedef Rcpp::NumericVector(*lna_ptr)(Rcpp::NumericVector init, double start, double end, double step_size);",
                                "// [[Rcpp::export]]",
                                "Rcpp::XPtr<lna_ptr> LNA_XPtr() {",
                                "return(Rcpp::XPtr<lna_ptr>(new lna_ptr(&INTEGRATE_STEM_LNA)));",
                                "}", sep = "\n")

        # function to set the LNA parameters, differs from the odeintr function in that it takes a numeric vector, not a std::vector
        param_setter   <- paste("void SET_LNA_PARAMS(Rcpp::NumericVector p) {",
                                "std::copy(p.begin(), p.end(), odeintr::pars.begin());",
                                "}\n",
                                "typedef void(*set_pars_ptr)(Rcpp::NumericVector p);",
                                "// [[Rcpp::export]]",
                                "Rcpp::XPtr<set_pars_ptr> LNA_set_params_XPtr() {",
                                "return(Rcpp::XPtr<set_pars_ptr>(new set_pars_ptr(&SET_LNA_PARAMS)));",
                                "}",sep = "\n")

        # paste the LNA integrator and parameter setting functions together
        stemr_LNA_code <- paste(LNA_integrator, param_setter, sep = "\n \n")


        # get the code for the LNA ODEs
        LNA_ODE_code <- odeintr::compile_sys(name = "INTEGRATE_LNA",
                                             sys = LNA_odes,
                                             sys_dim = n_compartments*2 + n_compartments^2,
                                             pars = n_params,
                                             globals = paste(paste(paste0("static arma::vec drift(", n_compartments,");"),
                                                                   paste0("static arma::vec resid(", n_compartments,");"),
                                                                   paste0("static arma::mat diffusion(", n_compartments,",",n_compartments,");"),
                                                                   paste0("static arma::vec drift_ode(", n_compartments,");"),
                                                                   paste0("static arma::vec resid_ode(", n_compartments,");"),
                                                                   paste0("static arma::vec diffusion_ode(", n_compartments,");"),
                                                                   paste0("static arma::vec hazard(", n_rates,");"),
                                                                   paste0("static arma::mat jacobian(", n_rates,",",n_compartments, ");"), sep = "\n"),
                                                             stoich_mtx_code, sep = "\n"),
                                             headers = paste("// [[Rcpp::depends(RcppArmadillo)]]",
                                                             "#include <RcppArmadillo.h>",
                                                             "using namespace arma;",
                                                             sep = "\n"),
                                             compile = F) # get the C++ code

        # RcppArmadillo is included, so remove the include tag for Rcpp
        LNA_ODE_code <- gsub("#include <Rcpp.h>", "", LNA_ODE_code)

        # we don't need the odeintr functions exported to the R global environment,
        # so remove the export attributes
        LNA_ODE_code <- gsub("// \\[\\[Rcpp::export\\]\\]", "", LNA_ODE_code)

        # paste the stemr LNA ode code to the end of the odeintr generated code
        LNA_ODE_code <- paste(LNA_ODE_code, stemr_LNA_code, sep = "\n")

        # compile the LNA code
        if(messages) print("Compiling LNA functions.")
        Rcpp::sourceCpp(code = LNA_ODE_code, env = globalenv())

        # get the LNA function pointers
        lna_pointer <- c(lna_ptr = LNA_XPtr(), set_lna_params_ptr = LNA_set_params_XPtr())

        return(lna_pointer)
}