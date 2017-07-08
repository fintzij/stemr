#' Specify the prior distribution for a model parameter
#'
#' @param parameter character string giving the name of the model parameter
#' @param distribution character string giving the prior distribution for the
#'   parameter, may be specified according to the table below
#' @param hyperpars named vector of prior distribution hyperparameters with
#'   names matching the arguments of the appropriate distribution in the table
#'   below
#' @param scale on which the prior is specified. Defaults to "linear", may be
#'   specified as "log" or "logit". The appropriate transformation will be
#'   applied automatically for the proposals.
#'
#' @section Specification of prior distributions:
#'
#'   The distributions and hyperparameters should be specified according to the
#'   table below
#'
#'   \tabular{lll}{ distribution \tab abbreviation \tab hyperparameters \cr Beta
#'   \tab beta \tab shape1, shape2, ncp \cr Cauchy \tab cauchy \tab location,
#'   scale \cr Gamma \tab gamma \tab shape, rate = 1, scale = 1/rate \cr Normal
#'   \tab norm \tab mean, sd \cr Truncated normal (a,b) \tab tnorm \tab mean,
#'   sd, a, b \cr Student t \tab t \tab df, ncp \cr Uniform \tab unif \tab min,
#'   max}
#'
#'   N.B. The prior for t0 should be specified through the
#'   \code{\link{t0_prior}} function.
#'
#' @return list to be used in evaluating the prior densities of model
#'   parameters.
#'
#' @export
#'
#' @examples prior("theta", "norm", c(mean = 1, sd = 2), scale = "log")
prior <- function(parameter, distribution, hyperpars, scale = "linear") {

        if(!scale %in% c("linear", "log", "logit")) {
                stop("The parameter estimation scale must be one of 'linear', 'log', or 'logit'.")
        }

        if(is.null(names(hyperpars))) {
                stop("The vector of hyperparameters must have named elements.")
        }

        return(list(parameter = parameter, distribution = distribution, hyperpars = hyperpars, scale = scale))
}