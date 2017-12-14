#' Reset the nugget based based on a sample
#'
#' @param nugget current value, a scalar
#' @param parameter_samples parameter samples on their estimation
#' @param start beginning of range to include
#' @param end end of range to include
#' @param lna_initdist_inds indices for which parameters are initial distribution parameters
#'
#' @return vector of nugget standard deviations
#' @export
reset_nugget <- function(nugget, parameter_samples_est, start, end, lna_initdist_inds) {

        if(start == end) {
                new_nugget <- rep(nugget, sum(!(colnames(parameter_samples_est) %in% names(lna_initdist_inds))))
        } else {
                # transform the samples
                transformed_samples <- parameter_samples_est[start:end, -c(lna_initdist_inds+1), drop = FALSE]

                new_nugget <- nugget / sqrt(ncol(parameter_samples_est) - length(lna_initdist_inds)) * apply(transformed_samples, 2, sd)
        }
        return(new_nugget)
}