#' Specify the prior for t0 and the standard deviation of the random--walk
#' proposal distribution.
#'
#' @param rw_sd standard deviation of the random walk proposal.
#' @param mean prior mean of the truncated normal distribution
#' @param sd prior standard deviation of the truncated normal distribution
#' @param lower lower bound, defaults to -Inf if not specified.
#' @param upper upper bound, defaults to the first observation time if not
#'   specified.
#'
#' @return arguments for specifying prior and proposal distribution for t0
#' @export
t0_kernel <- function(rw_sd, mean, sd, lower = NULL, upper = NULL) {

        return(list(rw_sd = rw_sd, mean = mean, sd = sd, lower = lower, upper = upper))

}