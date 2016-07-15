#' Instatiate the C++ emission probability functions for simulation and density
#' evaluation for a stochastic epidemic model measurement process and return a
#' vector of function pointers.
#'
#' @param meas_procs list of measurement process functions
#' @param messages logical; print a message that the rates are being compiled?
#'
#' @return Two vector of strings that serve as function pointers.
#' @export
parse_meas_procs <- function(meas_procs, messages = TRUE) {

        # emitmat is the matrix of emission probabilities
        d_measure_args <- "Rcpp::NumericMatrix& emitmat, const Rcpp::LogicalVector& emit_inds, const int record_ind, const arma::rowvec& record, const arma::rowvec& state, const Rcpp::NumericVector& parameters, const Rcpp::NumericVector& constants, const arma::rowvec& tcovar"

        # obsmat is a matrix of observations
        r_measure_args <- "Rcpp::NumericMatrix& obsmat, const Rcpp::LogicalVector& emit_inds, const int record_ind, const arma::rowvec& state, const Rcpp::NumericVector& parameters, const Rcpp::NumericVector& constants, const arma::rowvec& tcovar"

        emits_rmeas <- emits_dmeas <- character(0)

        for(i in seq_along(emissions_rmeas)) {
                emits_rmeas <- paste(emits_rmeas,paste0("if(meas_vars[",i-1,"]) emitmat(record_ind,",i-1,") = ", meas_procs[[i]]$r_measure,";"), sep = "\n ")
                emits_dmeas <- paste(emits_dmeas,paste0("if(meas_vars[",i-1,"]) obsmat(record_ind,",i-1,") = ", meas_procs[[i]]$d_measure,";"), sep = "\n ")
        }

        # compile function for updating elements a rate vector
        code_r_measure <- paste("// [[Rcpp::depends(RcppArmadillo)]]",
                             "#include <RcppArmadillo.h>",
                             "using namespace arma;",
                             "using namespace Rcpp;",
                             paste0("void R_MEASURE(",r_measure_args,") {"),
                             emits_rmeas,
                             "}",
                             paste0("typedef void(*r_measure_ptr)(", r_measure_args,");"),
                             "// [[Rcpp::export]]",
                             "Rcpp::XPtr<r_measure_ptr> D_MEASURE_XPtr() {",
                             "return(Rcpp::XPtr<r_measure_ptr>(new r_measure_ptr(&R_MEASURE)));",
                             "}", sep = "\n")

        code_d_measure <- paste("// [[Rcpp::depends(RcppArmadillo)]]",
                               "#include <RcppArmadillo.h>",
                               "using namespace arma;",
                               "using namespace Rcpp;",
                               paste0("void D_MEASURE(",d_measure_args,") {"),
                               emits_dmeas,
                               "}",
                               paste0("typedef void(*d_measure_ptr)(", d_measure_args,");"),
                               "// [[Rcpp::export]]",
                               "Rcpp::XPtr<d_measure_ptr> D_MEASURE_XPtr() {",
                               "return(Rcpp::XPtr<d_measure_ptr>(new d_measure_ptr(&D_MEASURE)));",
                               "}", sep = "\n")


        if(messages) {
                print("Compiling rate functions.")
        }

        Rcpp::sourceCpp(code = code_r_measure, env = globalenv())
        Rcpp::sourceCpp(code = code_d_measure, env = globalenv())

        rate_pointers <- c(r_measure_ptr = D_MEASURE_XPtr(), d_measure_ptr = R_MEASURE_XPtr())

        return(rate_pointers)
}