#' Return the log prior density over all parameter blocks
#'
#' @param param_blocks list of parameter blocks
#'
#' @return log prior density
#' @export
calc_params_logprior <- function(param_blocks) {
    return(sum(sapply(param_blocks, "[[", "log_pd")))
}