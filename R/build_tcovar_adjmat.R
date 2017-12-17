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
#' @param forcings list of forcings
#'
#' @return adjacency matrix for which rates need to be updated when a
#'   time-varying covariate changes.
#' @export
build_tcovar_adjmat <- function(rates, tcovar_codes = NULL, forcings = NULL) {
        tcovar_adjmat <- matrix(FALSE, ncol = length(tcovar_codes), nrow = length(rates))
        colnames(tcovar_adjmat) <- adjmat_names <- names(tcovar_codes)
        rownames(tcovar_adjmat) <- unlist(lapply(rates, function(x) paste0(x$from,"2",x$to)))

        tcovar_names <- paste0("tcovar\\[",tcovar_codes,"\\]")

        for(t in seq_along(rates)) {
                tcovar_adjmat[t,] <- sapply(tcovar_names, grepl, rates[[t]]$lumped)
        }
        
        # identify rates changed by forcings
        if(!is.null(forcings)) {
              for(s in seq_along(forcings)) {
                    adj_col  <- which(colnames(tcovar_adjmat) == forcings[[s]]$tcovar_name)
                    affected <- grepl(forcings[[s]]$from, rownames(tcovar_adjmat)) |
                          grepl(forcings[[s]]$to, rownames(tcovar_adjmat))
                    tcovar_adjmat[affected, adj_col] <- TRUE
              }
        }

        return(tcovar_adjmat)
}