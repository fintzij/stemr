#' Construct a list for building a transition kernel to be used in proposing new
#' parameter values.
#'
#' @param parameters character vector with names of parameters
#' @param method algorithm to be used in generating parameter updates -
#'   currently only multivariate normal random walk is supported.
#' @param cov_mtx variance-covariance matrix for random walk proposals, with
#'   rows and columns corresponding to the parameters argument
#'
#' @return list containing objects to be used in the proposals.
#' @export
kernel <- function(parameters, method = "mvn_rw", cov_mtx) {

        if(!all.equal(length(parameters), nrow(cov_mtx), ncol(cov_mtx))) {
                stop("The covariance matrix must have number of rows and number of columns equal to the number of parameters")
        }

        if(method != "mvn_rw") {
                stop("Only multivariate normal random walk is currently supported")
        }

        colnames(cov_mtx) <- rownames(cov_mtx) <- parameters

        return(list(parameters = parameters, method = method, cov_mtx = cov_mtx))
}