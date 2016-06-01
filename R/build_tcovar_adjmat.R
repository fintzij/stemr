#' Construct an adjacency matrix for which rates need to be update when there is
#' a change in the time-varying covariates.
#'
#' Constructs an adjacency matrix with rates given in the rows, time, and
#' time-varying covariates on which they each depend given in columns. Thus,
#' when time or a time-varying covariate changes, all of the non-zero entries in
#' the corresponding column must be updated.
#'
#' @param rates list of rates in ste_dynamics
#' @param tcovar_codes character vector of time-varying covariate names
#'
#' @return adjacency matrix for which rates need to be updated when a
#'   time-varying covariate changes.
#' @export
build_tcovar_adjmat <- function(rates, tcovar_codes = NULL) {
        tcovar_adjmat <- matrix(0, ncol = length(tcovar_codes), nrow = length(rates))
        colnames(tcovar_adjmat) <- adjmat_names <- names(tcovar_codes)
        rownames(tcovar_adjmat) <- paste0("RATE", 1:length(rates))

        tcovar_names <- paste0("tcovar\\[",tcovar_codes,"\\]")

        for(t in seq_along(rates)) {
                tcovar_adjmat[t,] <- sapply(tcovar_names, grepl, rates[[t]]$lumped)
        }

        return(tcovar_adjmat)
}