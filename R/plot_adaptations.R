#' Plot the adaptation schedule for an adaptive MCMC algorithm for a given
#' number of iterations, given as \code{scale_constant/(1 + iteration *
#' step_size)^scale_cooling}
#'
#' @param scale_constant maximum scaling, must be less than 1.
#' @param scale_cooling rate at which to cool the adaptation, must be between
#'   0.5 and 1.
#' @param iterations number of MCMC iterations
#' @param step_size increment per iteration
#'
#' @return produce a plot of the cooling schedule
#' @export
plot_adaptations <- function(scale_constant, scale_cooling, iterations, step_size) {
    if(scale_cooling <= 0.5 | scale_cooling >1) warning("The cooling rate must be between 0 and 1 in order to ensure that the chain is irreducible and converges.")
        
  plot(x = seq_len(iterations+1),  y = scale_constant * (seq(0,iterations) * step_size+1)^-scale_cooling,
             "l", ylim = c(0,1), main = paste0("Scale cooling = ", scale_cooling),
             xlab = "Iteration", ylab = "Adaptation factor")
}