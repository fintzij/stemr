#' Specify an MCMC transition kernel for updating the model parameters.
#'
#' @description Model parameters can be updated either using random walk
#'   Metropolis-Hastings (RWMH) with a fixed proposal covariance matrix, or
#'   using a modified version of the adaptive RWMH of Roberts and Rosenthal
#'   (2009).
#'
#' @param method algorithm to be used in generating parameter updates -
#'   currently only multivariate normal random walk and an adaptive version are
#'   supported.
#' @param cov_mtx variance-covariance matrix for random walk proposals, with
#'   rows and columns corresponding to the parameters argument.
#' @param scale_start iteration at which to begin scaling the proposal
#'   covariance matrix, defaults to 100.
#' @param shape_start number of acceptances required before estimating the
#'   scaled empirical covariance matrix in the proposals, defaults to 100.
#' @param scale_cooling rate at which to cool the adaptation, defaults to 0.999.
#' @param max_scaling maximum scale factor, defaults to 50.
#' @param target target acceptance rate when adapting the scale, defaults to
#'   0.234.
#' @param nugget either a scalar or vector for the nugget covariance in the
#'   adaptive RWMH algorithm (defaults to 0.01 * diag(cov_mtx)) if not
#'   specified).
#' @param nugget_weight weight that the nugget contribution receives in the
#'   adaptive RWMH proposal (defaults to 0.05).
#'
#' @details The adaptive algorithm will adapt the scale of the covariance matrix
#'   in order to achieve the target acceptance rate. Starting at iteration
#'   \code{scale_start}, the scaling factor is calculated as \eqn{\lambda =
#'   min(scaling * e^(scale_cooling^(iteration - scale_start) *
#'   (acceptances/iteration-target)), max_scaling)}, and the initial covariance
#'   matrix is scaled by \eqn{\lambda^2}. Once \code{shape_start} acceptances
#'   have been obtained, the scaling factor is set to 2.38^2/d, where d is the
#'   number of parameters. The proposal distribution at iteration n is
#'   constructed as a mixture of Gaussian kernels, as \deqn{Q_n(x,\cdot) = (1 -
#'   nugget_weight) * N(x, 2.38^2 \Sigma_n / d) + nugget_weight * N(x, diag(nugget) / d).}
#'
#'   References: G.O. Roberts and J.S. Rosenthal. "Examples of adaptive MCMC."
#'   *Journal of Computational and Graphical Statistics* 18.2 (2009): 349-367.
#'
#' @return list containing the method and covariance matrix for the MCMC kernel.
#' @export
kernel <- function(method = "mvn_rw", cov_mtx, scale_start=100, shape_start=100, scale_cooling=0.999, max_scaling=50, target=0.234, nugget = NULL, nugget_weight = 0.05) {

        if(!all.equal(length(parameters), nrow(cov_mtx), ncol(cov_mtx))) {
                stop("The covariance matrix must have number of rows and number of columns equal to the number of parameters")
        }

        if(!identical(rownames(cov_mtx), colnames(cov_mtx))) {
                stop("The variables in the covariance matrix rows and columns must be identically ordered")
        }

        if(!method %in% c("mvn_rw", "mvn_adaptive")) {
                stop("Only multivariate normal random walk or adaptive multivariate normal random walk are currently supported.")
        }

        if(target <= 0 || target >=1) {
                stop("The target acceptance rate must be between 0 and 1.")
        }

        if(scale_cooling <= 0 || scale_cooling >=1) {
                stop("The scale cooling rate must be between 0 and 1.")
        }

        if(is.null(nugget)) {
                nugget <- 0.01 * diag(cov_mtx) / nrow(cov_mtx)
        } else if(length(nugget) == 1) {
                nugget <- rep(nugget, nrow(cov_mtx)) / nrow(cov_mtx)
        } else {
                nugget <- nugget / nrow(cov_mtx)
        }

        names(nugget) <- rownames(cov_mtx)

        kernel_settings <- list(scale_start   = scale_start,
                                scale_cooling = scale_cooling,
                                max_scaling   = max_scaling,
                                shape_start   = shape_start,
                                target        = target,
                                nugget        = nugget,
                                nugget_weight = nugget_weight)

        return(list(method = method, cov_mtx = cov_mtx, kernel_settings = kernel_settings))
}