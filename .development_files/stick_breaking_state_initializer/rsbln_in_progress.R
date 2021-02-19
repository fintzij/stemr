dynamics$initdist_priors


sample_inits <- function(n, compartments, stick_means, stick_sds) {
  p <- length(compartments)
  partial_stick_proportions <- map2_dfr(logit(stick_means), stick_sds, ~expit(rnorm(n, .x, .y))) %>% as.matrix() %>% cbind(1)
  full_stick_proportions <- matrix(nrow = n, ncol = p)
  full_stick_proportions[, 1] = partial_stick_proportions[, 1]
  for (i in 2:p) full_stick_proportions[, i] <- (1 - rowSums(full_stick_proportions[ , 1:(i-1), drop = F])) * partial_stick_proportions[, i, drop = F]
  colnames(full_stick_proportions) <- compartments
  as_tibble(full_stick_proportions * popsize)
}

initdist_params <- stem_object$dynamics$initdist_params
initdist_priors <- stem_object$dynamics$initdist_priors

initdist_params <- c(S_0 = 2739459, E_0 = 17488, I_0 = 52463, R_0 = 365205, D_0 = 1077)
initdist_priors <- c(S = NA, E = 0.500344086093721, I = 0.515752688227904, R = 0.0869139093238033, D = 0.0434436515566288)

rsbln <- function(n, initdist_params, initdist_priors) {
  names(initdist_params) <- substr(names(initdist_params), start = 0, stop = nchar(names(initdist_params)) - 2)
  rev(initdist_params[!is.na(initdist_priors)]
  stick_means <- head(initdist_params / rev(cumsum(rev(initdist_params))), -1)
}



# Here are the names that have priors
names_with_priors <- names(initdist_priors[!is.na(initdist_priors)])
sorted_names <- names(sort(initdist_params, decreasing = T))
sorted_names_with_priors <- sorted_names[sorted_names %in% names_with_priors]

n_compartments <- length(sorted_names)

# Q: Why is this called stick_means?
# A: It's the mean length of the compartment's stick compared to the stick it is broken from
stick_means <- (initdist_params[sorted_names] / rev(cumsum(rev(initdist_params[sorted_names]))))[-n_compartments]
logit_stick_means <- logit(stick_means)

# tring to replicate this:
partial_stick_proportions <- map2_dfr(logit_stick_means, stick_sds, ~expit(rnorm(n, .x, .y))) %>% as.matrix() %>% cbind(1)
n <- 100
sapply(seq_along(logit_stick_means), expit(rnorm(n, logit_stick_means, initdist_priors[names(logit_stick_means)]))
