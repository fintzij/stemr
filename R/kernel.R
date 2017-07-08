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
#' @param covmat variance-covariance matrix for random walk proposals, with rows
#'   and columns corresponding to the parameters argument.
#' @param scale_start iteration at which to begin scaling the proposal
#'   covariance matrix, defaults to 100.
#' @param shape_start number of acceptances required before estimating the
#'   scaled empirical covariance matrix in the proposals, defaults to 100.
#' @param scale_cooling rate at which to cool the adaptation, defaults to 0.999.
#' @param max_scaling maximum scale factor, defaults to 50.
#' @param target target acceptance rate when adapting the scale, defaults to
#'   0.234.
#' @param nugget a scalar giving the ratio of the nugget standard deviation to
#'   the proposal standard deviation, defaults to 0.01.
#' @param nugget_weight weight that the nugget contribution receives in the
#'   adaptive RWMH proposal (defaults to 0.05).
#' @param reset_nugget By default, the nugget is reset to nugget *
#'   diag(empirical covmat) when shape_start acceptances have been obtained. The
#'   nugget can also be reset at whenever acceptances%%reset_nugget == 0 in
#'   addition to this.
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
#'   nugget_weight) * N(x, 2.38^2 \Sigma_n / d) + nugget_weight * N(x,
#'   diag(nugget) / d).}
#'
#'   References: G.O. Roberts and J.S. Rosenthal. "Examples of adaptive MCMC."
#'   *Journal of Computational and Graphical Statistics* 18.2 (2009): 349-367.
#'
#'   Code for the adaptive algorithm is also based on a slightly different
#'   algorithm implemented in the \code{pomp} packages
#'   (https://github.com/kingaa/pomp/blob/master/R/proposals.R).
#'
#' @return list containing the method and covariance matrix for the MCMC kernel.
#' @export
kernel <- function(method = "mvn_rw", covmat, scale_start=NULL, shape_start=NULL, scale_cooling=0.999, max_scaling=50, target=0.234, nugget = NULL, nugget_weight = 0.05) {

        if(!all.equal(nrow(covmat), ncol(covmat))) {
                stop("The covariance matrix must be square.")
        }

        if(!identical(rownames(covmat), colnames(covmat))) {
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

        if(method == "mvn_adaptive" && !is.null(nugget) && (nugget < 0 || nugget > 1)) {
                stop("The nugget must be between 0 and 1.")
        }

        # we take square roots b/c we will just multiply the nugget contribution by the standard deviation
        if(is.null(nugget)) {
                nugget <- 0.01
        }

        kernel_settings <- list(scale_start   = as.integer(scale_start),
                                scale_cooling = as.numeric(scale_cooling),
                                max_scaling   = as.numeric(max_scaling),
                                shape_start   = as.integer(shape_start),
                                target        = as.numeric(target),
                                nugget        = nugget,
                                nugget_weight = as.numeric(nugget_weight))

        return(list(method = method, covmat = covmat, kernel_settings = kernel_settings))
}