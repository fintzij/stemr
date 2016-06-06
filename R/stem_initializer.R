#' Determines the initial state probabilities at the first observation time.
#'
#' The state of each stratum in the system is initialized at the first
#' observation time, either with a vector of compartment counts, or with a
#' vector of state probabilities, where it is assumed that \eqn{X(t_1)\sim
#' Categorical(\cdot, \mathbf{p}_{t_1})}. The initializer argument in the
#' \code{\link{stem_dynamics}} function is supplied as a list of calls to
#' \code{stem_initializer}.
#'
#' To specify the initial state probabilities, a vector of compartment counts is
#' supplied in the \code{init_states} argument and \code{fixed} is set to
#' \code{FALSE}. The vector is then normalized to give the initial state
#' probabilities, though the size of the stratum population is first read off. A
#' vector of initial state probabilities will be generated. If a non-random
#' vector is initial counts is supplied, the \code{fixed} argument is set to
#' \code{TRUE}.
#'
#' @param init_states named vector of initial compartment counts, normalized to
#'   initial state probabilities.
#' @param fixed should the vector of compartment counts be treated as random?
#' @param strata character vector strata to which the vector of initial
#'   compartment counts or probabilities applies, possibly "ALL".
#' @param shared_params should the compartments that are specified together in
#'   the strata argument also share parameters? Defaults to FALSE.
#'
#' @return list of settings used to initialize the initial state at the first
#'   observation time.
#' @export
#'
#' @examples # All strata are of size 1000, with \eqn{P(X_j=S)=0.5,P(X_j=I)=0.2,
#' # P(X_j=R=0.3)}
#' stem_initializer(c(S=500, I = 200, R = 300), fixed = FALSE, strata = "ALL")
#'
#' # All strata are of size 1000. At the first observation time, there are 200
#' # susceptible, 100 infected, and 700 recovered adults, and there are 900
#' # susceptibles children, 50 infected children, and 50 recovered children. The
#' # distribution among the elderly is the same as among children.
#' list(stem_initializer(c(S=200,I=100,R=700), fixed = TRUE, strata = "adults"),
#'      stem_initializer(c(S=900,I=50,R=50), fixed=TRUE, strata = c("children",
#'      "elderly")))
stem_initializer <- function(init_states, fixed, strata = NULL, shared_params = FALSE) {

        if(is.null(names(init_states))) {
                stop("Compartment names must be specified in each initial probability vector.")
        }

        list(init_states = init_states, fixed = fixed, strata = strata, shared_params = shared_params)
}