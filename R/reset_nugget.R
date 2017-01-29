#' Reset the nugget based based on a sample
#'
#' @param nugget current value, a scalar
#' @param parameter_samples parameter samples on their natural scale
#' @param start beginning of range to include
#' @param end end of range to include
#' @param estimation_scales character vector of estimation scales
#'
#' @return vector of nugget standard deviations
#' @export
reset_nugget <- function(nugget, parameter_samples, start, end, estimation_scales) {

        if(start == end) {
                new_nugget <- rep(nugget, length(estimation_scales))
        } else {
                # transform the samples
                transformed_samples <- parameter_samples[start:end, , drop = FALSE]

                for(s in seq_len(ncol(transformed_samples))) {
                        if(estimation_scales[s] == "log") {
                                transformed_samples[,s] <- log(transformed_samples[,s])
                        } else if(estimation_scales[s] == "logit") {
                                transformed_samples[,s] <- logit(transformed_samples[,s])
                        }
                }

                new_nugget <- nugget / sqrt(length(estimation_scales)) * apply(transformed_samples, 2, sd)
        }
        return(new_nugget)
}