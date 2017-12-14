#' Interact character vectors to get a character vector of concatenated
#' combinations.
#'
#' @param ... character vectors
#'
#' @return character vector whose elements are the concatenated elements of the
#'   input character vectors
#' @export
#'
#' @examples interact(sex=c("male", "female"), age=c("young", "adult", "old"))
interact <- function(...) {
        # get vectors to be combined
        args <- list(...)

        # generate matrix of combinations
        combns <- expand.grid(args)

        # concatenate rows of the matrix
        interactions <- apply(combns, 1, paste, collapse = "_")

        return(interactions)
}