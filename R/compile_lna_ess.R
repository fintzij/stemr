#' Instatiate the C++ rate functions for the LNA elliptical slice sampling in
#' which only the residual and drift processes need to be integrated.
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
compile_lna_ess <- function(lna_rates, flow_matrix, lna_scale, messages = TRUE) {

        # get the number of rates and the number of compartments
        n_rates         <- nrow(flow_matrix)
        n_compartments  <- ncol(flow_matrix)
        n_params        <- length(lna_rates$lna_param_codes)

        # header code
        headers <- paste("// [[Rcpp::depends(RcppArmadillo)]]",
                        "#include <RcppArmadillo.h>",
                        "using namespace arma;",
                        "using namespace Rcpp;",
                        sep = "\n")

        # function title
        first_line <- paste("Rcpp::List lna_ess_fcn(double t, arma::vec& state, Rcpp::NumericVector& parms, arma::mat& stoich) {",
                            "int n_compartments = stoich.n_rows;",
                            "int n_odes = state.n_elem;",
                            "arma::vec drift = state.subvec(0,n_compartments-1);",
                            "arma::vec resid = state.subvec(n_compartments, 2*n_compartments - 1);",
                            "arma::vec hazard(n_compartments, arma::fill::zeros);",
                            "arma::mat jacobian(n_compartments, n_compartments, arma::fill::zeros);",
                            "arma::vec lna_ess_vec(n_odes, arma::fill::zeros);",sep = "\n")

        # construct the body of the lna ODEs. The first n_rates compartments are
        # the odes driven by the hazard functions. The next n_rates x n_compartment
        # rates are the jacobian odes
        jacobian_inds <- matrix(c(rep(seq(0, n_compartments - 1), each = n_compartments),
                                  rep(seq(0, n_compartments - 1), n_compartments)), ncol = 2)
        drift_inds    <- 0:(n_compartments-1)
        resid_inds    <- n_compartments:(2*n_compartments-1)

        # strings to compute the rates, hazards (cumulative hazards = A'phi) and jacobian
        hazard_terms    <- paste(paste("hazard[",0:(n_compartments-1),"]", " = ", lna_rates$hazards, ";", sep = ""), collapse = "\n")
        jacobian_terms  <- paste(paste("jacobian(", jacobian_inds[,1], ", ", jacobian_inds[,2], ") = ",lna_rates$derivatives,";", sep = ""),collapse = "\n")

        # strings to compute the LNA odes
        drift_ode       <- "lna_ess_vec(arma::span(0, n_compartments-1)) = hazard;"
        resid_ode       <- "lna_ess_vec(arma::span(n_compartments, 2*n_compartments-1)) = jacobian * resid;"

        # return command
        return_step    <- paste("return Rcpp::List::create(lna_ess_vec);","}", sep = "\n")

        # generate the stemr_lna functions that will actually be called
        LNA_ess_caller <- paste("typedef Rcpp::List(*lna_ess_ptr)(double t, arma::vec& state, Rcpp::NumericVector& parms, arma::mat& stoich);",
                                "// [[Rcpp::export]]",
                                "Rcpp::XPtr<lna_ess_ptr> LNA_ESS_XPtr() {",
                                "return(Rcpp::XPtr<lna_ess_ptr>(new lna_ess_ptr(&lna_ess_fcn)));",
                                "}", sep = "\n")

        # concatenate everything
        LNA_ess_code <- paste(headers, first_line,
                          hazard_terms, jacobian_terms,
                          drift_ode, resid_ode,
                          return_step, LNA_ess_caller,sep = "\n")

        # compile the LNA code
        if(messages) print("Compiling LNA elliptical slice sampling functions.")
        Rcpp::sourceCpp(code = LNA_ess_code, env = globalenv())

        # get the LNA function pointers
        lna_pointer <- c(lna_ess_ptr = LNA_ESS_XPtr(), LNA_ess_code = LNA_ess_code)

        return(lna_pointer)
}