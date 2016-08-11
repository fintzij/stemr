#' Compile the ODE functions for a stochastic epidemic model
#'
#' @param ode_rates character vector of hazard rates, as outputted by
#'   \code{\link{build_ode_rates}}
#' @inheritParams parse_lna_fcns
#'
#' @return external function pointers for integrating the ODE
#' @export
parse_ode_fcns <- function(ode_rates, flow_matrix, messages) {

        # get the number of rates and the number of compartments
        n_rates         <- nrow(flow_matrix)
        n_compartments  <- ncol(flow_matrix)
        n_params        <- length(ode_rates$ode_param_codes)

        # construct the stochiometry matrix constructor
        stoich_mtx <- t(flow_matrix)
        stoich_mtx_code <- paste("static arma::mat ode_stoich = {",
                                 paste("{", apply(stoich_mtx, 1, paste, collapse = ", "),
                                       "}", collapse = ",\n"),"};")

        # construct the body of the ODEs.
        dxdt_inds <- 0:(n_compartments-1)

        # strings to compute the hazards and jacobian
        hazard_terms <- paste(paste("odeintr::ode_hazard[",0:(n_rates-1),"]", " = ", ode_rates$hazards, ";", sep = ""), collapse = "\n ")

        # strings to compute the LNA odes
        ode_vec <- "odeintr::ode_vec = odeintr::ode_stoich * odeintr::ode_hazard;"

        # dxdt strings
        dxdt_ode <- paste("dxdt[", dxdt_inds, "] = odeintr::ode_vec[", dxdt_inds, "];", collapse = "\n", sep = "")

        # concatenate everything
        ODEs <- paste(hazard_terms, ode_vec, dxdt_ode, sep = "\n")

        # generate the stemr_lna functions that will actually be called
        # stemr_LNA_integrator is equivalent to the _no_record odeintr function, but takes the endpoints of the interval and returns a numeric vector.
        ODE_integrator <- paste("void INTEGRATE_STEM_ODE(Rcpp::NumericVector& init, double start, double end, double step_size=1.0) {",
                                "INTEGRATE_ODE_set_state(init);",
                                "odeint::integrate_adaptive(odeintr::stepper, odeintr::sys, odeintr::state, start, end, step_size);",
                                "init = Rcpp::wrap(INTEGRATE_ODE_get_state());",
                                "}\n",
                                "typedef void(*lna_ptr)(Rcpp::NumericVector& init, double start, double end, double step_size);",
                                "// [[Rcpp::export]]",
                                "Rcpp::XPtr<lna_ptr> ODE_XPtr() {",
                                "return(Rcpp::XPtr<lna_ptr>(new lna_ptr(&INTEGRATE_STEM_ODE)));",
                                "}", sep = "\n")

        # function to set the LNA parameters, differs from the odeintr function in that it takes a numeric vector, not a std::vector
        ODE_param_setter   <- paste("void SET_ODE_PARAMS(Rcpp::NumericVector& p) {",
                                "std::copy(p.begin(), p.end(), odeintr::pars.begin());",
                                "}\n",
                                "typedef void(*set_pars_ptr)(Rcpp::NumericVector& p);",
                                "// [[Rcpp::export]]",
                                "Rcpp::XPtr<set_pars_ptr> ODE_set_params_XPtr() {",
                                "return(Rcpp::XPtr<set_pars_ptr>(new set_pars_ptr(&SET_ODE_PARAMS)));",
                                "}",sep = "\n")

        # paste the ODE integrator and parameter setting functions together
        stemr_ODE_code <- paste(ODE_integrator, ODE_param_setter, sep = "\n \n")


        # get the code for the LNA ODEs
        ODE_code <- odeintr::compile_sys(name = "INTEGRATE_ODE",
                                         sys = ODEs,
                                         sys_dim = n_compartments,
                                         pars = n_params,
                                         globals = paste(paste(paste0("static arma::vec ode_vec(", n_compartments,");"),
                                         paste0("static arma::vec ode_hazard(", n_rates,");"), sep = "\n"),
                                         stoich_mtx_code, sep = "\n"),
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

        # paste the stemr LNA ode code to the end of the odeintr generated code
        ODE_code <- paste(ODE_code, stemr_ODE_code, sep = "\n")

        # compile the LNA code
        if(messages) print("Compiling ODE functions.")
        Rcpp::sourceCpp(code = ODE_code, env = globalenv())

        # get the LNA function pointers
        ode_pointers <- c(ode_ptr = ODE_XPtr(), set_ode_params_ptr = ODE_set_params_XPtr())

        return(ode_pointers)

}