#' Generate a list of settings for automated factor slice sampling
#'
#' @param n_harss_updates number of hit-and-run updates.
#' @param bracket_update_interval how often the bracket update interval should be updated
#'
#' @return list with additional settings for hit and run slice sampling
#' @export
harss_settings <-
      function(n_harss_updates = 1,
               bracket_update_interval = 100) {
      
      list(n_harss_updates          = n_harss_updates,
           bracket_update_interval  = bracket_update_interval)
}
