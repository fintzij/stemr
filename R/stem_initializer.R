#' Determines the initial state probabilities or concentrations at the first
#' observation time. This function is applied internally in different ways
#' depending on whether inference (or simulation) is accomplished using the LNA
#' or using Bayesian data augmentation (Gillespie's direct algorithm for
#' simulation).
#'
#' The initial state of the system is initialized at the first observation time,
#' either with a fixed vector of compartment counts, or is sampled by drawing
#' the initial compartment counts from a dirichlet-multinomial distribution for
#' a Markov jump process (MJP), or by the normal approximation to the
#' dirichlet-multinomial (LNA/ODE). If the initial state is not fixed and no
#' prior is specified, then the hyperparameters are handled internally as
#' follows: For the MJP, if the multinomial hyperprior parameters are not
#' specified, the normalized probabilities are computed from the init_states
#' vector and used in the prior slot. For the LNA or ODE, the init_state vector
#' is used to construct the multivariate normal approximation of the
#' multinomial.
#'
#' @param init_states named vector of initial compartment counts (MJP) or
#'   volumes (LNA).
#' @param fixed should the vector of compartment counts be treated as random?
#' @param strata character vector strata to which the vector of initial
#'   compartment counts or probabilities applies, possibly "ALL".
#' @param prior hyperparameters for state probabilities at time t0.
#' @param dist one of "multinom" if the initial distribution (or its normal
#'   approximation) is multinomial or "dirmultinom" if dirichlet-multinomial.
#'
#' @return list of settings used to initialize the initial state at the first
#'   observation time.
#' @export
#'
#' @examples Fixed initial state, all strata are of size 1000:
#' stem_initializer(c(S=500, I=200, R=300), fixed=TRUE, strata = "ALL")
#'
#' Random initial state where for a single subject P(X_j=S)=0.5,P(X_j=I)=0.2,
#' P(X_j=R=0.3), and one set of initial state probabilities govern the
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
stem_initializer <- function(init_states, fixed, strata = NULL, prior = NULL, dist = "multinom") {

        if(is.null(names(init_states))) {
                stop("Compartment names must be specified in each initial probability vector.")
        }

        list(init_states = init_states, fixed = fixed, strata = strata, prior = prior, dist = dist)
}