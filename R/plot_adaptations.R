#' Plot the adaptation schedule for an adaptive MCMC algorithm for a given
#' number of iterations, given as \code{max(1, scale_constant/(1 + iteration *
#' step_size)^scale_cooling)}
#'
#' @param scale_constant maximum scaling.
#' @param scale_cooling rate at which to cool the adaptation, must be between
#'   0.5 and 1.
#' @param iterations number of MCMC iterations
#' @param step_size increment per iteration
#'
#' @return produce a plot of the cooling schedule
#' @export
plot_adaptations <- function(iterations, scale_cooling, scale_constant = 1, step_size = 1, adaptation_offset = 0, plot_limits = c(0,1)) {
    if(scale_cooling <= 0.5 | scale_cooling >1) warning("The cooling rate must be in (0.5,1] in order to ensure that the chain is irreducible and converges.")
        
  plot(x = seq_len(iterations+1), 
       y = pmin(1, scale_constant * (seq(0,iterations) * step_size + adaptation_offset + 1)^-scale_cooling),
             "l", ylim = plot_limits, main = paste0("Scale cooling = ", scale_cooling),
             xlab = "Iteration", ylab = "Adaptation factor")
  abline(h = 0)
  abline(h = scale_constant * (adaptation_offset + 1)^-scale_cooling, lty = 2)
  abline(h = scale_constant * (iterations * step_size + adaptation_offset + 1)^-scale_cooling, lty = 2)
}