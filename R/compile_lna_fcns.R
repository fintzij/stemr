#' Instatiate the C++ rate functions for the LNA using the odeintr package.
#'
#' @param lna_rates list containing the strings vectors of ODES for the hazards
#'   and jacobian matrix, along with the lna parameter codes.
#' @param flow_matrix matrix indicating the changes to each compartment per
#'   event
#' @param lna_scale either "log" or "linear"
#' @param messages logical; print a message that the rates are being compiled
#'
#' @return External function pointers for integrating the LNA
#' @export
compile_lna_fcns <- function(lna_rates, flow_matrix, lna_scale, messages = TRUE) {

        # get the number of rates and the number of compartments
        n_rates         <- nrow(flow_matrix)
        n_compartments  <- ncol(flow_matrix)
        n_params        <- length(lna_rates$lna_param_codes)

        # header code
        headers <- paste("// [[Rcpp::depends(RcppArmadillo)]]",
                         "// [[Rcpp::plugins(cpp11)]]",
                        "#include <RcppArmadillo.h>",
                        "using namespace arma;",
                        "using namespace Rcpp;",
                        sep = "\n")

        # construct the stochiometry matrix constructor
        stoich_mtx      <- t(flow_matrix)
        stoich_mtx_code <- paste("static arma::mat stoich = {",paste("{", apply(stoich_mtx, 1, paste, collapse = ", "), "}", collapse = ",\n "),"};")

        # code for the global variables
        globals <- paste("namespace stemr_lna {",
                         paste0("static arma::vec lna_out_vec(",2*n_compartments + n_compartments^2,");"),
                         paste0("static arma::vec drift(", n_compartments,");"),
                         paste0("static arma::vec resid(", n_compartments,");"),
                         paste0("static arma::mat diffusion(",n_compartments, ",", n_compartments,");"),
                         paste0("static arma::vec rates(", n_rates,");"),
                         paste0("static arma::vec hazard(", n_compartments,");"),
                         paste0("static arma::mat jacobian(", n_compartments, ",", n_compartments,");"),
                         paste0("static int n_compartments = ", n_compartments,";"),
                         paste0("static int n_odes = ", 2*n_compartments + n_compartments^2,";"),
                         stoich_mtx_code,
                         "}", sep = "\n")

        # function title
        first_line <- paste0("Rcpp::List compute_lna(double t, arma::vec& state, Rcpp::NumericVector& parms) {")

        # construct the body of the lna ODEs. The first n_rates compartments are
        # the odes driven by the hazard functions. The next n_rates x n_compartment
        # rates are the jacobian odes
        jacobian_inds <- matrix(c(rep(seq(0, n_compartments - 1), each = n_compartments), rep(seq(0, n_compartments - 1), n_compartments)), ncol = 2)
        drift_inds    <- 0:(n_compartments-1)
        resid_inds    <- n_compartments:(2*n_compartments-1)
        diffusion_inds<- (2*n_compartments):(2*n_compartments + n_compartments^2 - 1)
        diff_mat_inds <- matrix(c(rep(seq(0, n_compartments - 1), n_compartments), rep(seq(0, n_compartments - 1), each = n_compartments)), ncol = 2)

        # string to construct the drift, residual, and diffusion vectors
        # drift_terms     <- paste0("stemr_lna::drift[", drift_inds,"] = state[",drift_inds,"];", collapse = "\n")
        # resid_terms     <- paste0("stemr_lna::resid[", resid_inds,"] = state[",resid_inds,"];", collapse = "\n")
        # diffusion_terms <- paste0("stemr_lna::diffusion(", diff_mat_inds[,1],",",diff_mat_inds[,2],") = x[",diffusion_inds,"];", collapse = "\n")
        vec2procs_step  <- paste("std::copy(state.begin(), state.begin() + stemr_lna::n_compartments, stemr_lna::drift.begin());",
                                 "std::copy(state.begin() + stemr_lna::n_compartments, state.begin() + 2*stemr_lna::n_compartments, stemr_lna::resid.begin());",
                                 "std::copy(state.begin() + 2*stemr_lna::n_compartments, state.end(), stemr_lna::diffusion.begin());", sep = "\n")

        # strings to compute the rates, hazards (cumulative hazards = A'phi) and jacobian
        rate_terms      <- paste(paste("stemr_lna::rates[",0:(n_rates-1),"]", " = ", lna_rates$rates,";", sep = ""), collapse = "\n")
        hazard_terms    <- paste(paste("stemr_lna::hazard[",0:(n_compartments-1),"]", " = ", lna_rates$hazards, ";", sep = ""), collapse = "\n")
        jacobian_terms  <- paste(paste("stemr_lna::jacobian(", jacobian_inds[,1], ", ", jacobian_inds[,2], ") = ",lna_rates$derivatives,";",
                                       sep = ""),collapse = "\n")

        # strings to compute the LNA odes
        drift_ode       <- "stemr_lna::lna_out_vec(arma::span(0, stemr_lna::n_compartments-1)) = stemr_lna::hazard;"
        resid_ode       <- "stemr_lna::lna_out_vec(arma::span(stemr_lna::n_compartments, 2*stemr_lna::n_compartments-1)) = stemr_lna::jacobian * stemr_lna::resid;"

        if(lna_scale == "log") {
                diffusion_ode   <- paste0("stemr_lna::lna_out_vec(arma::span(2*stemr_lna::n_compartments, stemr_lna::n_odes-1)) = arma::vectorise(stemr_lna::diffusion * stemr_lna::jacobian.t() + ",
                                          " arma::diagmat(arma::exp(-stemr_lna::drift)) * stemr_lna::stoich * arma::diagmat(stemr_lna::rates) * stemr_lna::stoich.t() * arma::diagmat(arma::exp(-stemr_lna::drift)) + ",
                                          "stemr_lna::jacobian * stemr_lna::diffusion);")
        } else if(lna_scale == "linear") {
                diffusion_ode   <- paste0("stemr_lna::lna_out_vec(arma::span(2*stemr_lna::n_compartments, stemr_lna::n_odes-1)) = arma::vectorise(stemr_lna::diffusion * stemr_lna::jacobian.t() + ",
                                          "stemr_lna::stoich * arma::diagmat(stemr_lna::rates) * stemr_lna::stoich.t() + ",
                                          "stemr_lna::jacobian * stemr_lna::diffusion);")
        }

        # return command
        return_step    <- paste("return Rcpp::List::create(stemr_lna::lna_out_vec);","}", sep = "\n")

        # generate the stemr_lna functions that will actually be called
        LNA_caller     <- paste("typedef Rcpp::List(*compute_lna_ptr)(double t, arma::vec& state, Rcpp::NumericVector& parms);",
                                "// [[Rcpp::export]]",
                                "Rcpp::XPtr<compute_lna_ptr> COMPUTE_LNA_XPtr() {",
                                "return(Rcpp::XPtr<compute_lna_ptr>(new compute_lna_ptr(&compute_lna)));",
                                "}", sep = "\n")

        # concatenate everything
        LNA_code <- paste(headers, globals, first_line, vec2procs_step,
                          rate_terms, hazard_terms, jacobian_terms,
                          drift_ode, resid_ode, diffusion_ode,
                          return_step, LNA_caller,sep = "\n")

        # compile the LNA code
        if(messages) print("Compiling LNA functions.")
        Rcpp::sourceCpp(code = LNA_code, env = globalenv())

        # get the LNA function pointers
        lna_pointer <- c(lna_ptr = COMPUTE_LNA_XPtr(), LNA_code = LNA_code)

        return(lna_pointer)
}