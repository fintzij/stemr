#' Construct an adjacency matrix for specifying which rates must be updated when
#' a transition occurs, with adjacency determined at the lumped population
#' level.
#'
#' Constructs an adjacency matrix with rates given in the rows, and reactions on
#' which they each depend given in columns. Thus, when a reaction occurs, all of
#' the non-zero entries in the corresponding column must be updated.
#'
#' @param rates intermediate list of rate functions created within the
#'   stem_dynamics function
#' @param compartment_codes vector of compartment codes.
#'
#' @return adjacency matrix
#' @export
build_rate_adjmat <- function(rates, compartment_codes) {

        depends_on  <- matrix(0, nrow = length(rates), ncol = length(compartment_codes))
        affects     <- matrix(0, nrow = length(rates), ncol = length(compartment_codes))

        colnames(depends_on) <- colnames(affects) <- names(compartment_codes)
        rownames(depends_on) <- rownames(affects) <- paste0("RATE",1:length(rates))

        code_strings <- paste0("state\\[",compartment_codes,"\\]")

        # determine dependencies
        for(r in seq_along(rates)) {
                depends <- sapply(code_strings, grepl, rates[[r]]$lumped)
                depends_on[r, depends] <- 1

                affects[r, c(rates[[r]]$from, rates[[r]]$to)] <- 1
        }

        # construct the adjacency matrix
        rate_adjmat <- ifelse((depends_on %*% t(affects)) > 0, TRUE, FALSE)

        return(rate_adjmat)
}