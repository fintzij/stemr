#' Determines whether a model is progressive based on its flow matrix
#'
#' @param flow_matrix
#'
#' @return TRUE/FALSE indicating whether a model is progressive
#' @export
is_progressive <- function(flow_matrix) {

        dag <- diag(0, ncol(flow_matrix))

        for(j in 1:nrow(flow_matrix)) {
                dag[which(flow_matrix[j,, drop = FALSE] == -1), which(flow_matrix[j,, drop = FALSE] == 1)] <- 1
        }

        progressive <- ifelse(any(dag[lower.tri(dag)] != 0), FALSE, TRUE)

        return(progressive)
}