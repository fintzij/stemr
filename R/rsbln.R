#' Stick-Breaking Log-Normal
#'
#' @param n number of observations
#' @param logit_stick_means vector of stick means on the logit scale
#' @param stick_sds vector of stick sds on the logit scale
#' @param stick_size size to scale the samples by
#'
#' @return samples from a stick breaking process
#' @export
#'
#' @examples
rsbln <- function(n, logit_stick_means, stick_sds, stick_size) {
  n_compartments <- length(logit_stick_means) + 1

  partial_stick_proportions <- cbind(
    matrix(
      expit(rnorm(n = n * n_half_parameters,
                  mean = logit_stick_means,
                  sd = stick_sds)),
      nrow = n, ncol = n_compartments - 1, byrow = T),
    1)

  full_stick_proportions <- matrix(nrow = n, ncol = n_compartments)
  full_stick_proportions[, 1] = partial_stick_proportions[, 1]
  for (i in 2:n_compartments) full_stick_proportions[, i] <- (1 - rowSums(full_stick_proportions[ , 1:(i-1), drop = F])) * partial_stick_proportions[, i, drop = F]
  full_stick_proportions * stick_size
}
