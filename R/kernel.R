#' Specify an MCMC transition kernel for updating the model parameters.
#'
#' @description Model parameters can be updated either using either
#'   componentwise or global Metropolis, with options for global or
#'   componentwise adaptive scaling.
#'
#' @param method algorithm to be used in generating parameter updates, see
#'   details below.
#' @param sigma vector of variances or variance-covariance matrix for random
#'   walk proposals, with names or rows and columns corresponding to the names
#'   of parameters on their estimation scales.
#' @param scale_constant constant multiple of the adaptations determined by
#'   \code{scale_cooling}.
#' @param scale_cooling rate at which to cool the adaptation, defaults to 0.5.
#'   Adaptation contributions are governed by a harmonic sequence:
#'   scale_constant/(iteration/step_size+1)^scale_cooling. The
#'   \code{plot_adaptations} function may be used to plot the adaptation
#'   factors.
#' @param step_size adaptation increment for each iteration, defaults to 1.
#' @param max_scaling maximum scale factor, defaults to Inf.
#' @param target_g target acceptance rate for global Metropolis proposals.
#' @param target_c target acceptance rate for componentwise Metropolis proposals
#'   or for principal components Metropolis proposals.
#' @param nugget fixed nugget contribution, defaults to 0.01.
#' @param stop_adaptation iteration at which to stop adapting the proposal
#'   distribution.
#' @param weight_update_interval number of iterations between updates to the
#'   direction sampling distribution if using adaptive principal components
#'   Metropolis. Defaults to 100 times the number of parameters if not
#'   specified. The weights are the n^th root of the eigenvalues, where n is the
#'   number of model parameters, normalized to sum to 1.
#' @param messages should messages be printed?
#'
#' @details Specifies a Metropolis transition kernel wtih symmetric Gaussian
#'   proposals. The options for the method are as follows: 1) c_rw:
#'   componentwise random walk, updating each parameter in turn. 2) mvn_rw:
#'   global random walk, updating all parameters jointly. 3) c_rw_adaptive:
#'   componentwise adaptive Metropolis with componentwise adaptive scaling (Alg.
#'   5 in Andrieu and Thoms). 4) mvn_rw_c_adaptive: global adaptive Metropolis
#'   with componentwise adaptive scaling (Alg. 6 in Andrieu and Thoms). 5)
#'   mvn_rw_g_adaptive: global adaptive Metropolis with global adaptive scaling
#'   (Alg. 4 in Andrieu and Thoms). 6) pca_adaptive: adaptive principal
#'   components Metropolis (Alg. 8 in Andrieu and Thoms).
#'
#'   References: Andrieu, Christophe, and Johannes Thoms. "A tutorial on
#'   adaptive MCMC." Statistics and computing 18.4 (2008): 343-373.
#'
#' @return list containing the method and covariance matrix for the MCMC kernel.
#' @export
kernel <-
        function(method = "c_rw",
                 sigma,
                 scale_constant = 0.5,
                 scale_cooling = 0.5+1e-5,
                 step_size = 1,
                 max_scaling = Inf,
                 target_g = 0.234,
                 target_c = 0.44,
                 nugget = NULL,
                 stop_adaptation = Inf,
                 weight_update_interval = NULL,
                 messages = TRUE) {

        if(!method %in% c("c_rw", "mvn_rw", "c_rw_adaptive", "mvn_c_adaptive", "mvn_g_adaptive", "pca_adaptive")) {
                stop("The method for the MCMC kernel is not correctly specified.")
        }
          
        if(scale_cooling <=0.5 | scale_cooling > 1) {
                stop("The cooling rate must be between 0.5 and 1.")
        }

        if(target_g < 0 || target_g >1) {
                stop("The target acceptance rate must be between 0 and 1.")
        }

        if(target_c < 0 || target_c >1) {
                stop("The target acceptance rate must be between 0 and 1.")
        }

        if(is.null(nugget)) {
                nugget <- 0.001 * min(diag(sigma))
        }

        if(method == "c_rw_adaptive" && length(nugget) != nrow(sigma)) {
                if(messages) warning("The same nugget contribution will be used for all components, which may not be appropriate.")

                nugget <- rep(nugget, nrow(sigma))
        }

        if(method == "pca_adaptive" && is.null(weight_update_interval)) {
              weight_update_interval <- 100 * nrow(sigma)
        }              

        kernel_settings <- list(scale_constant = scale_constant,
                                scale_cooling = scale_cooling,
                                step_size     = step_size,
                                max_scaling   = as.numeric(max_scaling),
                                target_g      = target_g,
                                target_c      = target_c,
                                nugget        = nugget,
                                stop_adaptation = stop_adaptation,
                                weight_update_interval = weight_update_interval)

        return(list(method = method, sigma = sigma, kernel_settings = kernel_settings))
}