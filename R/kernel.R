#' Construct a list for building a transition kernel to be used in proposing new
#' parameter values.
#'
#' @param parameters character vector with names of parameters
#' @param method algorithm to be used in generating parameter updates -
#'   currently only multivariate normal random walk is supported.
#' @param cov_mtx variance-covariance matrix for random walk proposals with
#'   column and row names corresponding to the parameters.
#'
#' @return list containing objects to be used in the proposals.
#' @export
kernel <- function(parameters, method = "mvn_rw", cov_mtx) {

        if(!all.equal(length(parameters), nrow(cov_mtx), ncol(cov_mtx))) {
                stop("The covariance matrix must have number of rows and number of columns equal to the number of parameters")
        }

        if(!all(parameters %in% rownames(cov_mtx)) || !all(parameters %in% colnames(cov_mtx))) {
                stop("The row and column names of the covariance matrix must match the names of the parameters")
        }

        if(!identical(rownames(cov_mtx), colnames(cov_mtx))) {
                stop("The variables in the covariance matrix rows and columns must be identically ordered")
        }

        if(method != "mvn_rw") {
                stop("Only multivariate normal random walk is currently supported")
        }

        return(list(parameters = parameters, method = method, cov_mtx = cov_mtx))
}