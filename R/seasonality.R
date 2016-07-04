#' Produces a string for seasonal terms to be added to a rate function, supplied
#' as the last argument in a rate list.
#'
#' @param period period of the season, defaults to 52 for weeks. May also be
#'   supplied as a vector if S>1 in order to specify multiple terms with
#'   different periods.
#' @param s_param_values seasonality parameter values (or initial values), a
#'   2*length(period) numeric vector, with all sine parameters, then all cosine
#'   parameters, and supplied in the same order as period. If no values are
#'   supplied, the parameter values default to 0. Thus, if \code{period = c(12,
#'   52)}, the parameter values would be, for example, \code{c(SIN_12 = 1,
#'   SIN_52 = 2, COS_12 = 3, COS_52 = 4)}.
#' @param log should the seasonal covariates enter into the unlumped rates on
#'   the log scale? If true (default), each rate, \eqn{\lambda}, is log-linear
#'   in the seasonal terms. For example, if \eqn{\lambda_{S,I} = (\beta I +
#'   \alpha(t))S}, then \deqn{\alpha(t) = exp(\sum \gamma_p \times sin(2\pi
#'   p\times t) + \delta_p \times cos(2\pi p\times t))}
#' @param common_seasonlity logical for whether the seasonal terms should be
#'   shared among strata (TRUE) or whether each stratum should have its own
#'   seasonality
#'
#' @description Supplies a string for the seasonal terms \deqn{\sum \gamma_p
#'   \times sin(2\pi p\times t) + \delta_p \times cos(2\pi p\times t)} to be
#'   added to a rate function. The parameters \eqn{\gamma period} and
#'   \eqn{\delta period} will automatically be named to reflect whether the
#'   seasonal terms are common to all strata, and then appended to the parameter
#'   vector. If \code{log=TRUE}, the seasonal terms are exponentiated. Thus,
#'   \deqn{\sum \gamma_p \times sin(2\pi p\times t) + \delta_p \times cos(2\pi
#'   p\times t)}
#'
#' @return string formula for the seasonal terms to be added to a rate.
#' @export
#'
#' @examples # add two sine and cosine terms, one with a period of 12 (for months), the other with a period of 52 (for weeks)
#' list("beta * I", S, I, strata = NULL, seasonality(period = c(12, 52), s_params = c(0.01, 0.05)))
seasonality <- function(period = 52, s_params = NULL, log = TRUE, common_seasonality = TRUE) {

        form <- vector(mode = "character")
        nterms <- length(period)

        # generate param names
        if(is.null(s_params)) s_params <- rep(0, 2*length(period))
        s_param_names <- names(s_params) <- paste0(rep(c("SIN", "COS"), each = length(period)), period, ifelse(common_seasonality, "", "_SELF"))

        # generate sine and cosine terms
        sin_terms <- paste0("sin(", 2*seq(1,nterms), "*M_PI*TIME/", period, ")")
        cos_terms <- paste0("cos(", 2*seq(1,nterms), "*M_PI*TIME/", period, ")")

        # generate the formula
        f <- paste(paste(s_param_names, c(sin_terms, cos_terms), sep = "*"), collapse = " + ")

        if(log) f <- paste0("exp(",f,")")

        return(list(seasonality = f, s_params = s_params, common_seasonality = common_seasonality, period = period))
}