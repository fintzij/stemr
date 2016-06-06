#' Detects which states in a model are absorbing states
#'
#' @param flow_matrix
#'
#' @return logical vector indicating which model compartments are absorbing.
#' @export
which_absorbing <- function(flow_matrix) {

        absorbing_states <- rep(FALSE, ncol(flow_matrix))
        names(absorbing_states) <- colnames(flow_matrix)

        for(j in 1:ncol(flow_matrix)) {
                absorbing_states[j] <-  ifelse((1 %in% flow_matrix[,j]) & (!-1 %in% flow_matrix[,j]), TRUE, FALSE)
        }

        return(absorbing_states)
}