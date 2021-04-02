#' Stick-Breaking Logit-Normal
#'
#' @param n number of observations
#' @param stick_means vector of stick means on the logit scale
#' @param stick_sds vector of stick sds on the logit scale
#' @param stick_size size by which to scale the observations (population size)
#'
#'
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

#' Convert Normal Draws to Stick-Breaking Logit-Normal Draws
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


#' Stick-Breaking Logit-Normal Explorer
#'
#' This function facilitates the creation of Stick-Breaking Logit-Normal priors
#' with user-sepcified properties (median, lower, and upper quantiles for each
#' compartment).
#'
#' @param n number of samples to simulate when fitting parameters
#' @param target_median Named vector containing target median values of each
#' compartment. The sum of target_median is used as the population size.
#' @param target_lower Named vector containing target lower quantile value of
#' each compartment. One compartment (the same as in target_upper) should be
#' specified as NA.
#' @param target_upper Named vector containing target upper quantile value of
#' each compartment. One compartment (the same as in target_lower) should be
#' specified as NA.
#' @param width float indicating the width of the interval used to determine the
#'  quantiles for target_lower and target_upper
#' @param samples_and_summary logical indicating if the samples from the
#' stick-breaking logit-normal should be returns and compared to the desired
#' median, lower, and upper quantiles
#'
#' @return R code that can be used to call stem_initializer to describe an sbln
#' with the desired properties. If samples_and_summary is TRUE, samples from
#' this prior and a summary of these samples are also returned.
#' @export
#'
#' @examples
#' sbln_explorer(n = 10000, target_median = c(S = 980, I = 20, R = 1), target_lower = c(S = 950, I = 5, R = NA), target_upper = c(S = 990, I = 40, R = NA), samples_and_summary = T)
sbln_explorer <- function(n, target_median, target_lower, target_upper, width = 0.9, samples_and_summary = F) {
  popsize <- sum(target_median)
  lower_p <- (1 - width) / 2
  upper_p <- 1 - ((1 - width) / 2)

  compartments_original_order <- names(target_median)
  target_median_original_order <- target_median
  target_lower_original_order <- target_lower
  target_upper_original_order <- target_upper

  n_compartments <- length(compartments_original_order)
  missing_compartment_index <- which(is.na(target_lower))

  compartments <- c(names(sort(target_median[-missing_compartment_index], decreasing = T)), compartments_original_order[missing_compartment_index])
  target_median <- target_median[compartments]
  target_lower <- target_lower[compartments]
  target_upper <- target_upper[compartments]

  stick_means <- head(target_median / rev(cumsum(rev(target_median))), -1)
  stick_means <- logit(stick_means)

  stick_sd_upper_limit <- sapply(1:(n_compartments - 1), function(i) optimize(f = function(x) norm(as.matrix(expit(qnorm(p = c(lower_p, upper_p), mean = stick_means[[i]], sd = x)) * popsize - c(target_lower[i], target_upper[i])), type = "f"), interval = c(0, 5))$minimum)
  names(stick_sd_upper_limit) <- names(stick_means)

  stick_sds <- numeric(n_compartments - 1)
  names(stick_sds) <- names(stick_means)
  stick_sds[1] <- stick_sd_upper_limit[1]
  to_allocate <- (1 - expit(rnorm(n, stick_means[1], stick_sd_upper_limit[1]))) * popsize

  for (i in 2:(n_compartments - 1)) {
    stick_sds[i] <- optimize(f = function(x) norm(as.matrix(quantile(to_allocate * expit(rnorm(n, stick_means[i], x)), c(lower_p, upper_p)) - c(target_lower[i], target_upper[i]))), lower = 0, upper = stick_sd_upper_limit[i])$minimum
    to_allocate <- (1 - expit(rnorm(n, stick_means[i], stick_sds[i]))) * to_allocate
  }
  capture.output(init_states_chr <- deparse1(dput(target_median)), file = "/dev/null")
  capture.output(prior_chr <- deparse1(dput(unname(c(stick_means, stick_sds)))))

  stem_initializer_code <-
    paste0("stem_initializer(",
           "init_states = ", init_states_chr, ", ",
           "prior = ", prior_chr, ", ",
           'dist = "sbln"', ", ",
           "fixed = F",
           ")")

  message("Use this code to create stem_initializer:")
  cat(stem_initializer_code)
  cat("\n")

  sbln_summary <-
    data.frame(compartment_name = compartments,
               target_median,
               target_lower,
               target_upper)

  if (samples_and_summary) {
    samples <- rsbln(n = n, stick_means, stick_sds, stick_size = popsize)
    colnames(samples) <- compartments

    sbln_summary <-
      data.frame(compartment_name = compartments,
                 target_median,
                 observed_median = apply(samples, 2, median),
                 target_lower,
                 observed_lower = apply(samples, 2, function(x) quantile(x, probs = lower_p)),
                 target_upper,
                 observed_upper = apply(samples, 2, function(x) quantile(x, probs = upper_p)))
    return(
      list(
        summary = sbln_summary,
         samples = samples))
  } else {
    return(NULL)
  }
}
