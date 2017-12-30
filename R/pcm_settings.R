#' Specify settings for principal component metropolis
#'
#'@param weight_update_interval number of iterations between updates to the
#'  direction sampling distribution if using adaptive principal components
#'  Metropolis. Defaults to 100 times the number of parameters if not specified.
#'  The weights are the sqrt(n)^th root of the eigenvalues, where n is the number of
#'  model parameters, normalized to sum to 1.
#'@param pcm_update_schedule either "random" (default) or "sequential". If
#'  "random", each iteration consists of one parameter update in a direction
#'  chosen at random according to the weight distribution. If "sequential", each
#'  MCMC iteration consists of an update in each direction.
#'@param n_pcm_updates If the principal component metropolis schedule is random,
#'  the number of directions, chosen at random, that should be updated per
#'  iteration.
#'
#' @return list containing settings for principal component metropolis
#' @export
pcm_settings <-
      function(weight_update_interval = NULL,
               pcm_update_schedule = "random",
               n_pcm_updates = 1) {
            
      return(list(weight_update_interval = weight_update_interval,
                  pcm_update_schedule = pcm_update_schedule, 
                  n_pcm_updates = n_pcm_updates))
}