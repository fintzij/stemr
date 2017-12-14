#' Expit transformation, i.e. inverse logit
#'
#' @param x
#'
#' @return expit(x)
#' @export
expit <- function(x) {
        1/(1+exp(-x))
}

#' Logit tranformation
#'
#' @param y
#'
#' @return logit(y)
#' @export
logit <- function(y) {
        -log(1/y - 1)
}