#' Instatiate the C++ rate functions for the LNA with integration done using
#' deSolve.
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
compile_lna <- function(lna_rates, flow_matrix, lna_scale, messages = TRUE) {

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
        first_line <- paste("Rcpp::List compute_lna(double t, arma::vec& state, Rcpp::NumericVector& parms, arma::mat& stoich) {",
                            "int n_compartments = stoich.n_rows;",
                            "int n_rates = stoich.n_cols;",
                            "int n_odes = state.n_elem;",
                            "arma::vec drift = state.subvec(0,n_compartments-1);",
                            "arma::vec resid = state.subvec(n_compartments, 2*n_compartments - 1);",
                            "arma::mat diffusion = arma::reshape(state.subvec(2*n_compartments, n_odes - 1), n_compartments, n_compartments);",
                            "arma::vec rates(n_rates, arma::fill::zeros);",
                            "arma::vec hazard(n_compartments, arma::fill::zeros);",
                            "arma::mat jacobian(n_compartments, n_compartments, arma::fill::zeros);",
                            "arma::vec lna_out_vec(n_odes, arma::fill::zeros);",sep = "\n")

        # construct the body of the lna ODEs. The first n_rates compartments are
        # the odes driven by the hazard functions. The next n_rates x n_compartment
        # rates are the jacobian odes
        jacobian_inds <- matrix(c(rep(seq(0, n_compartments - 1), each = n_compartments), rep(seq(0, n_compartments - 1), n_compartments)), ncol = 2)
        drift_inds    <- 0:(n_compartments-1)
        resid_inds    <- n_compartments:(2*n_compartments-1)
        diffusion_inds<- (2*n_compartments):(2*n_compartments + n_compartments^2 - 1)
        diff_mat_inds <- matrix(c(rep(seq(0, n_compartments - 1), n_compartments), rep(seq(0, n_compartments - 1), each = n_compartments)), ncol = 2)


        # strings to compute the rates, hazards (cumulative hazards = A'phi) and jacobian
        rate_terms      <- paste(paste("rates[",0:(n_rates-1),"]", " = ", lna_rates$rates,";", sep = ""), collapse = "\n")
        hazard_terms    <- paste(paste("hazard[",0:(n_compartments-1),"]", " = ", lna_rates$hazards, ";", sep = ""), collapse = "\n")
        jacobian_terms  <- paste(paste("jacobian(", jacobian_inds[,1], ", ", jacobian_inds[,2], ") = ",lna_rates$derivatives,";", sep = ""),collapse = "\n")

        # strings to compute the LNA odes
        drift_ode       <- "lna_out_vec(arma::span(0, n_compartments-1)) = hazard;"
        resid_ode       <- "lna_out_vec(arma::span(n_compartments, 2*n_compartments-1)) = jacobian * resid;"

        if(lna_scale == "log") {
                diffusion_ode   <- paste0("lna_out_vec(arma::span(2*n_compartments, n_odes-1)) = arma::vectorise(diffusion * jacobian.t() + ",
                                          " arma::diagmat(arma::exp(-drift)) * stoich * arma::diagmat(rates) * stoich.t() * arma::diagmat(arma::exp(-drift)) + ",
                                          "jacobian * diffusion);")
        } else if(lna_scale == "linear") {
                diffusion_ode   <- paste0("lna_out_vec(arma::span(2*n_compartments, n_odes-1)) = arma::vectorise(diffusion * jacobian.t() + ",
                                          "stoich * arma::diagmat(rates) * stoich.t() + ",
                                          "jacobian * diffusion);")
        }

        # set numerical error to zero
        # zero_error <- paste("lna_out_vec(arma::find(abs(lna_out_vec) < 1e-6)).zeros();")

        # return command
        return_step    <- paste("return Rcpp::List::create(lna_out_vec);","}", sep = "\n")

        # generate the stemr_lna functions that will actually be called
        LNA_caller     <- paste("typedef Rcpp::List(*compute_lna_ptr)(double t, arma::vec& state, Rcpp::NumericVector& parms, arma::mat& stoich);",
                                "// [[Rcpp::export]]",
                                "Rcpp::XPtr<compute_lna_ptr> COMPUTE_LNA_XPtr() {",
                                "return(Rcpp::XPtr<compute_lna_ptr>(new compute_lna_ptr(&compute_lna)));",
                                "}", sep = "\n")

        # concatenate everything
        LNA_code <- paste(headers, first_line,
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