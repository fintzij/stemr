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
        flow_matrix           <- matrix(0L, ncol = length(compartment_names), nrow = length(rates))
        colnames(flow_matrix) <- compartment_names
        rownames(flow_matrix) <- rate_names <- paste0(sapply(rates, "[[", "from"), "2", sapply(rates, "[[", "to"))

        # fill out the flow matrix
        for(k in 1:length(rates)) {
                flow_matrix[k, rates[[k]][["from"]]] <- -1
                flow_matrix[k, rates[[k]][["to"]]] <- 1
        }

        incidence_rates <- sapply(rates, "[[", "incidence")

        if(any(incidence_rates)) {
                incidence_matrix           <- diag(1, nrow = length(rates), ncol = sum(incidence_rates))
                rownames(incidence_matrix) <- rate_names
                colnames(incidence_matrix) <- rate_names[incidence_rates]
                flow_matrix                <- cbind(flow_matrix, incidence_matrix)
        }

        return(flow_matrix)
}