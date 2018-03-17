#' Generate a list of settings for automated factor slice sampling
#'
#' @param target_ratio target ratio of expansions/(expansions + contractions),
#'   defaults to 0.5. Smaller values overrelax the interval widths.
#' @param n_harss_updates number of hit-and-run updates.
#' @param initial_width initial bracket width, defaults to 1.
#'
#' @return list with additional settings for hit and run slice sampling
#' @export
harss_settings <-
      function(n_harss_updates = 1,
               initial_width = 1,
               target_ratio = 0.5) {
      
      list(n_harss_updates          = n_harss_updates,
           initial_width            = initial_width,
           target_ratio             = target_ratio)
}
