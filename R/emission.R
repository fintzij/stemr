#' Generates an emission list to be supplied to the \code{\link{stem_measure}}
#' function.
#'
#' @param meas_var character vector supplying the name of the measurement
#'   variable, either the name of a column in a dataset, or the name of a
#'   variable in a dataset to be simulated.
#' @param distribution emission distribution, one of "binomial", "poisson",
#'   "negbinomial", "gaussian"
#' @param emission_params character vector of emission distribution parameters,
#'   generally a function of model compartments and parameters.
#' @param incidence do the data represent incidence counts (as opposed to
#'   prevalence counts)? defaults to FALSE. If true, the incidence compartments
#'   will be used.
#' @param obstimes numeric vector of observation times, required if the
#'   measurement process is to be simulated from.
#' @param strata strata to which the emission distribution applies, may be
#'   supplied as "ALL"
#'
#' @section Specifying emission distributions:
#'
#'   Each emission list specifies, at a minimum, the names of the measurement
#'   variables, the emission distribution, and the parameters of the emission
#'   distribution. Optionally, the user may supply a character vector of strata,
#'   possibly using the reserved word "ALL" to indicate all model strata. The
#'   \code{meas_var} arguments are then appended with the strata names and a
#'   separate emission distribution is generated for each stratum. All emission
#'   distributions in the stochastic epidemic model are assumed to be
#'   conditionally independent.
#'
#'   The available emission distributions are the binomial, poisson, negative
#'   binomial, and gaussian distributions. The user specifies the parameters as
#'   strings in the canonical order they are presented. Thus, \enumerate{\item
#'   poisson: lambda \item binomial: size, prob \item negative binomial: size,
#'   prob \item gaussian: mean, sd}. The \code{strata} argument may either be
#'   supplied as a character vector of model strata, or may be specified as
#'   "ALL" to indicate all model strata. The case-sensitive key word, "SELF",
#'   may be used in the name of the measurement variable and the parameter
#'   strings in conjunction with the \code{strata} argument to facilitate
#'   specification of multiple measurement processes. Thus, if there are
#'   multiple strata whose observations are accrued at the same set of times,
#'   the user would only need to supply a single emission function since "SELF"
#'   will be parsed and replaced with the names of the strata.
#'
#'   The parameters for each distribution are supplied as character vectors, and
#'   will typically be functions of model compartments, parameters, time-varying
#'   covariates, and constants. As in specification of the rate functions, the
#'   strings must be valid C++ code, although if a string is a single line, a
#'   semi-colon need not be included. Some examples of emission lists follow
#'   here: \enumerate{\item \code{emission("prevalence", "binomial", c("I",
#'   "rho"))}: the observed variable, names "prevalence" is a binomial sample of
#'   the number of individuals in compartment "I", with sampling probability
#'   "rho". \item \code{emission("I_SELF", "binomial", c("I_SELF", "rho"),
#'   incidence = TRUE, strata = "ALL")}: the observed incidence for each stratum
#'   is a binomial sample of the true incidence in that stratum.}
#'
#' @return an emission list to be parsed in the \code{\link{stem_measure}}
#'   function.
#' @export
emission <- function(meas_var, distribution, emission_params, incidence = FALSE, obstimes = NULL, strata = NULL) {

        if(!distribution %in% c("binomial", "poisson", "negbinomial", "gaussian")) {
                stop("The emission distribution must be one of 'binomial', 'poisson', 'negbinomial', or 'gaussian'.")
        }

        # check that the correct number of parameters are supplied for each distribution
        if(distribution == "poisson" && length(emission_params) != 1) {
                stop("The mean parameter must be specified for poisson emission probabilities.")
        } else if((distribution == "binomial" || distribution == "negbinomial") && length(emission_params) != 2) {
                stop("The size and sampling probability parameters must be specified for binomial and negative binomial emission probabilities.")
        } else if(distribution == "normal" && length(emission_params) != 2) {
                stop("The mean and standard deviation must be specified for gaussian emission probabilities.")
        }

        # make sure that neither the measurement variable or emission parameters use "SELF" if strata are not specified
        if(is.null(strata)) {
                if(grepl("SELF", meas_var) || sapply(emission_params, grepl, pattern = "SELF")) {
                        stop("The 'SELF' key word may only be used when the strata to which the emission distribution applies are specified.")
                }
        }

        # generate the list
        list(meas_var = meas_var, distribution = distribution, emission_params = emission_params, incidence = incidence, obstimes = obstimes, strata = strata)
}