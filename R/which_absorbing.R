#' Detects which states in a model are absorbing states
#'
#' @param flow_matrix
#'
#' @return logical vector indicating which model compartments are absorbing.
#' @export
which_absorbing <- function(flow_matrix) {

        comp_names <- colnames(flow_matrix); incidence_comps <- grepl("INCIDENCE", comp_names)
        absorbing_states <- rep(FALSE, ncol(flow_matrix) - sum(incidence_comps))
        names(absorbing_states) <- comp_names[!incidence_comps]

        for(j in which(!incidence_comps)) {
                absorbing_states[j] <-  ifelse((1 %in% flow_matrix[,j]) & (!-1 %in% flow_matrix[,j]), TRUE, FALSE)
        }

        return(absorbing_states)
}