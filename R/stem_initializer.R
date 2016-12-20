#' Determines the initial state probabilities at the first observation time.
#' This function is applied internally in different ways depending on whether
#' inference (or simulation) is accomplished using the LNA or using Bayesian
#' data augmentation (Gillespie's direct algorithm for simulation).
#'
#' The initial state of the system is initialized at the first observation time,
#' either with a fixed vector of compartment counts, or with a vector of
#' unnormalized state probabilities, \eqn{\mathbf{p}_{t_0}}, the sum of which
#' determine the size of the population or stratum. \eqn{X(t_0)} is assigned a
#' prior distribution \eqn{X(t_0)\sim Categorical(\cdot, \mathbf{p}_{t_0}),} if
#' doing exact inference, or \eqn{X(t_0)\sim Multinomial(\cdot,
#' \mathbf{p}_{t_0}),} if performing approximate inference/simulation using the
#' LNA. The distinction is important as the exact data augmentation method
#' implemented in this package requires augmentation on the unlumped state space
#' of subject-level state labels. The \code{state_initializer} argument in the
#' \code{\link{stem_dynamics}} function is supplied as a list of calls to
#' \code{stem_initializer}.
#'
#' There are three possible ways of initializing the state at time \eqn{t_0},
#' which are determined by the combination of the \code{fixed} argument, and the
#' \code{prior} argument: 1) fixed initial state - \code{fixed = TRUE}; 2)
#' random initial state with initial state probabilities that correspond to the
#' normalized compartment counts in the \code{init_state} argument - \code{fixed
#' = TRUE, prior = NULL}; 3) random initial state with initial state  -
#' \code{fixed = TRUE, prior = vector of hyperparameters}. When the parameters
#' for the prior distribution are specified, the compartment counts are still
#' initialized in the first iteration using the vector supplied in the
#' \code{init_states} argument. The vector of compartment counts is also used
#' for determining the population size and the size of the various strata.
#'
#' When approximate inference is carried out using the LNA and the initial state
#' is random, the compartment counts at time \eqn{t_0} are updated in a
#' Metropolis-Hastings step using the multinomial prior as a proposal
#' distribution.
#'
#' @param init_states named vector of initial compartment counts.
#' @param fixed should the vector of compartment counts be treated as random?
#' @param strata character vector strata to which the vector of initial
#'   compartment counts or probabilities applies, possibly "ALL".
#' @param prior multinomial distribution state probability parameters.
#'
#' @return list of settings used to initialize the initial state at the first
#'   observation time.
#' @export
#'
#' @examples # All strata are of size 1000, with \eqn{P(X_j=S)=0.5,P(X_j=I)=0.2,
#' # P(X_j=R=0.3)}, and one set of initial state probabilities govern the
#' distribution of initial compartment counts in all strata:
#' stem_initializer(c(S=500, I = 200, R = 300), fixed = FALSE, strata = "ALL")
#'
#' # All strata are of size 1000. At the first observation time, there are 200
#' # susceptible, 100 infected, and 700 recovered adults, and there are 900
#' # susceptibles children, 50 infected children, and 50 recovered children. The
#' # distribution among the elderly is the same as among children:
#' list(stem_initializer(c(S=200,I=100,R=700), fixed = TRUE, strata = "adults"),
#'      stem_initializer(c(S=900,I=50,R=50), fixed=TRUE, strata = c("children",
#'      "elderly")))
stem_initializer <- function(init_states, fixed, strata = NULL, prior = NULL) {

        if(is.null(names(init_states))) {
                stop("Compartment names must be specified in each initial probability vector.")
        }

        # if a prior was not supplied but the initial state is random, compute the prior probabilities
        if(is.null(prior) && !fixed) {
                prior <- init_states / sum(init_states)
        }

        # make sure prior is normalized
        if(!is.null(prior)) {
                prior <- prior / sum(prior)
        }

        list(init_states = init_states, fixed = fixed, strata = strata, prior = prior)
}