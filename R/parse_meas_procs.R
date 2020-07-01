#' Instatiate the C++ emission probability functions for simulation and density
#' evaluation for a stochastic epidemic model measurement process and return a
#' vector of function pointers.
#'
#' @param meas_procs list of measurement process functions
#' @param messages logical; print a message that the rates are being compiled?
#'
#' @return Two vector of strings that serve as function pointers.
#' @export
parse_meas_procs <- function(meas_procs, compile_moments = FALSE, messages = TRUE) {
      
        # set to null to avoid warnings
      R_MEASURE_XPtr <- NULL
      D_MEASURE_XPtr <- NULL
      MEAS_MEAN_XPtr <- NULL
      MEAS_VAR_XPtr  <- NULL

        # emitmat is the matrix of emission probabilities
        d_measure_args <- "Rcpp::NumericMatrix& emitmat, const Rcpp::LogicalVector& emit_inds, const int record_ind, const Rcpp::NumericVector& record, const Rcpp::NumericVector& state, const Rcpp::NumericVector& parameters, const Rcpp::NumericVector& constants, const Rcpp::NumericVector& tcovar"

        # obsmat is a matrix of observations
        r_measure_args <- "Rcpp::NumericMatrix& obsmat, const Rcpp::LogicalVector& emit_inds, const int record_ind, const Rcpp::NumericVector& state, const Rcpp::NumericVector& parameters, const Rcpp::NumericVector& constants, const Rcpp::NumericVector& tcovar"

        r_meas <- d_meas <- m_meas <- v_meas <- character(0)

        for(i in seq_along(meas_procs)) {

                d_meas <- paste(d_meas,paste(paste0("if(emit_inds[",i-1,"]) {"),
                                             paste0("Rcpp::NumericVector obs(1,",meas_procs[[i]]$meas_var,");"),
                                             paste0("emitmat(record_ind,",i,") = ", meas_procs[[i]]$dmeasure,"[0];"),
                                             "}", sep = "\n"), sep = "\n ")

                r_meas <- paste(r_meas,paste0("if(emit_inds[",i-1,"] && (",
                                              meas_procs[[i]]$emission_params[1]," != 0)) obsmat(record_ind,",i,") = ",
                                              meas_procs[[i]]$rmeasure,"[0];"), sep = "\n ")

                m_meas <- paste(m_meas,paste0("if(emit_inds[",i-1,"]) obsmat(record_ind,",i,") = ",
                                              meas_procs[[i]]$mmeasure,";"), sep = "\n ")

                v_meas <- paste(v_meas,paste0("if(emit_inds[",i-1,"]) obsmat(record_ind,",i,") = ",
                                              meas_procs[[i]]$vmeasure,";"), sep = "\n ")
        }

        # compile function for updating elements a rate vector
        code_r_measure <- paste("// [[Rcpp::depends(Rcpp)]]",
                             "#include <Rcpp.h>",
                             "using namespace Rcpp;",
                             paste0("void R_MEASURE(",r_measure_args,") {"),
                             r_meas,
                             "}",
                             paste0("typedef void(*r_measure_ptr)(", r_measure_args,");"),
                             "// [[Rcpp::export]]",
                             "Rcpp::XPtr<r_measure_ptr> R_MEASURE_XPtr() {",
                             "return(Rcpp::XPtr<r_measure_ptr>(new r_measure_ptr(&R_MEASURE)));",
                             "}", sep = "\n")

        code_d_measure <- paste("// [[Rcpp::depends(Rcpp)]]",
                               "#include <Rcpp.h>",
                               "using namespace Rcpp;",
                               paste0("void D_MEASURE(",d_measure_args,") {"),
                               d_meas,
                               "}",
                               paste0("typedef void(*d_measure_ptr)(", d_measure_args,");"),
                               "// [[Rcpp::export]]",
                               "Rcpp::XPtr<d_measure_ptr> D_MEASURE_XPtr() {",
                               "return(Rcpp::XPtr<d_measure_ptr>(new d_measure_ptr(&D_MEASURE)));",
                               "}", sep = "\n")

        code_m_measure <- paste("// [[Rcpp::depends(Rcpp)]]",
                                "#include <Rcpp.h>",
                                "using namespace Rcpp;",
                                paste0("void MEAS_MEAN(",r_measure_args,") {"),
                                m_meas,
                                "}",
                                paste0("typedef void(*r_measure_ptr)(", r_measure_args,");"),
                                "// [[Rcpp::export]]",
                                "Rcpp::XPtr<r_measure_ptr> MEAS_MEAN_XPtr() {",
                                "return(Rcpp::XPtr<r_measure_ptr>(new r_measure_ptr(&MEAS_MEAN)));",
                                "}", sep = "\n")

        code_v_measure <- paste("// [[Rcpp::depends(Rcpp)]]",
                                "#include <Rcpp.h>",
                                "using namespace Rcpp;",
                                paste0("void MEAS_VAR(",r_measure_args,") {"),
                                v_meas,
                                "}",
                                paste0("typedef void(*r_measure_ptr)(", r_measure_args,");"),
                                "// [[Rcpp::export]]",
                                "Rcpp::XPtr<r_measure_ptr> MEAS_VAR_XPtr() {",
                                "return(Rcpp::XPtr<r_measure_ptr>(new r_measure_ptr(&MEAS_VAR)));",
                                "}", sep = "\n")

        if(messages) {
                print("Compiling measurement process functions.")
        }

        Rcpp::sourceCpp(code = code_r_measure, env = globalenv(), rebuild = TRUE)
        Rcpp::sourceCpp(code = code_d_measure, env = globalenv(), rebuild = TRUE)

        measproc_pointers <- c(r_measure_ptr = R_MEASURE_XPtr(),
                               d_measure_ptr = D_MEASURE_XPtr(),
                               meas_proc_code = paste(code_r_measure, code_d_measure, sep = "\n\n"))

        if(compile_moments) {
                Rcpp::sourceCpp(code = code_m_measure, env = globalenv(), rebuild = TRUE)
                Rcpp::sourceCpp(code = code_v_measure, env = globalenv(), rebuild = TRUE)

                measproc_pointers <- c(measproc_pointers,
                                       m_measure_ptr = MEAS_MEAN_XPtr(),
                                       v_measure_ptr = MEAS_VAR_XPtr())

                measproc_pointers$meas_proc_code = paste(measproc_pointers$meas_proc_code,
                                                         code_m_measure, code_d_measure,
                                                         collapse = "\n\n")
        }

        return(measproc_pointers)
}