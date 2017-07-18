#' Plot the adaptation schedule for an adaptive MCMC algorithm for a given
#' number of iterations, given as 1/iteration^scale_cooling.
#'
#' @param scale_cooling rate at which to cool the adaptation, must be between 0 and 1.
#' @param iterations number of MCMC iterations
#'
#' @return produce a plot of the cooling schedule
#' @export
plot_adaptations <- function(scale_cooling, iterations) {
        plot(x = seq_len(iterations),  y = (seq_len(iterations)+1)^-scale_cooling, "l",
             main = paste0("Scale cooling = ", scale_cooling),xlab = "Iteration", ylab = "Adaptation factor")
}