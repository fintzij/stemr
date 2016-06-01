#' Produces a string for seasonal terms to be added to a rate function, supplied
#' as the last argument in a rate list.
#'
#' @param period period of the season, defaults to 52 for weeks. May also be
#'   supplied as a vector if S>1 in order to specify multiple terms with
#'   different periods.
#' @param common_seasonlity logical for whether the seasonal terms should be
#'   shared among strata (TRUE) or whether each stratum should have its own
#'   seasonality
#' @description Supplies a string for the seasonal terms \deqn{\sum \gamma_p
#'   \times sin(2\pi p\times t) + \delta_p \times cos(2\pi p\times t)} to be
#'   added to a rate function. The parameters \eqn{\gamma period} and
#'   \eqn{\delta period} will automatically be named to reflect whether the
#'   seasonal terms are common to all strata, and then appended to the parameter
#'   vector.
#'
#' @return string formula for the seasonal terms to be added to a rate.
#' @export
#'
#' @examples # add two sine and cosine terms, one with a period of 12 (for months), the other with a period of 52 (for weeks)
#' list("beta * I", S, I, strata = NULL, seasonality(period = 12))
seasonality <- function(period = 52, common_seasonality = TRUE) {

        form <- vector(mode = "character")
        nterms <- length(period)

        # generate param names
        s_params <- paste0(rep(c("SIN", "COS"), each = length(period)), period, ifelse(common_seasonality, "", "_SELF"))

        # generate sine and cosine terms
        sin_terms <- paste0("sin(", 2*seq(1,nterms), "*M_PI*tcovar[0]/", period, ")")
        cos_terms <- paste0("cos(", 2*seq(1,nterms), "*M_PI*tcovar[0]/", period, ")")


        # generate the formula
        f <- paste(paste(s_params, c(sin_terms, cos_terms), sep = "*"), collapse = " + ")

        return(list(seasonality = f, s_params = s_params, common_seasonality = common_seasonality))
}