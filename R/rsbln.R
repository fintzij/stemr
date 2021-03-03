#' Stick-Breaking Log-Normal
#'
#' @param n number of observations
#' @param stick_means vector of stick means on the logit scale
#' @param stick_sds vector of stick sds on the logit scale
#' @param stick_size size to scale the samples by
#'
#' @return samples from a stick breaking process
#' @export
#'
#' @examples
rsbln <- function(n, stick_means, stick_sds, stick_size) {
  n_compartments <- length(stick_means) + 1

  normal_draws <-
    matrix(rnorm(n = n * (n_compartments - 1)),
           nrow = n,
           ncol = n_compartments - 1)

  sbln_normal_to_volume(normal_draws, stick_means, stick_sds, stick_size)
}

#' Stick-Breaking Log-Normal
#'
#' @param normal_draws a matrix of normal draws to convert to compartment volumes
#' @param stick_means vector of stick means on the logit scale
#' @param stick_sds vector of stick sds on the logit scale
#' @param stick_size size to scale the samples by
#'
#' @return samples from a stick breaking process
#' @export
#'
#' @examples
sbln_normal_to_volume <- function(normal_draws, stick_means, stick_sds, stick_size) {

  return_vector <- is.null(dim(normal_draws))

  if (return_vector) normal_draws <- matrix(normal_draws, nrow = 1)

  partial_stick_proportions <- cbind(
    expit(normal_draws *
            rep(stick_sds, each = dim(normal_draws)[1]) +
            rep(stick_means, each = dim(normal_draws)[1])),
    1)

  full_stick_proportions <- matrix(nrow = dim(partial_stick_proportions)[1], ncol = dim(partial_stick_proportions)[2])
  full_stick_proportions[, 1] = round(partial_stick_proportions[, 1] * stick_size)
  for (i in 2:dim(full_stick_proportions)[2]) full_stick_proportions[, i] <- round((stick_size - rowSums(full_stick_proportions[ , 1:(i-1), drop = F])) * partial_stick_proportions[, i, drop = F])

  if (return_vector) full_stick_proportions <- as.vector(full_stick_proportions)

  full_stick_proportions
}
