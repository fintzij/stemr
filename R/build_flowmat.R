#' Given a list of (unparsed) rate functions, construct the matrix of
#' compartment flows.
#'
#' @param rates Unparsed list of rate functions (intermediate object within the
#'   \code{stem_dynamics} function)
#' @param compartment_names vector of compartment_strata names
#'
#' @return matrix whose elements indicate the change in each compartment
#'   (columns) associated with each reaction event (rows).
#' @export
#'
build_flowmat <- function(rates, compartment_names) {

        # compartments are columns, reactions are rows
        flow_matrix             <- matrix(0L, ncol = length(compartment_names), nrow = length(rates))
        colnames(flow_matrix)   <- compartment_names

        # fill out the flow matrix
        for(k in 1:length(rates)) {
                flow_matrix[k, rates[[k]][[3]]] <- -1
                flow_matrix[k, rates[[k]][[4]]] <- 1
        }

        return(flow_matrix)
}