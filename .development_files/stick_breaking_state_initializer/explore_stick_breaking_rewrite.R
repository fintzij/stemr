library(stemr)
library(tidyverse)






stemr::rsbln()


rsbln <- function(n, compartments, logit_stick_means, stick_sds) {
  n_compartments <- length(compartments)

  partial_stick_proportions <- cbind(matrix(expit(rnorm(n * (n_compartments - 1), mean = logit_stick_means, sd = stick_sds)), nrow = n, ncol = n_compartments - 1, byrow = T), 1)
  colnames(partial_stick_proportions) <- compartments

  full_stick_proportions <- matrix(nrow = n, ncol = n_compartments, dimnames = dimnames(partial_stick_proportions))
  full_stick_proportions[, 1] = partial_stick_proportions[, 1]
  for (i in 2:n_compartments) full_stick_proportions[, i] <- (1 - rowSums(full_stick_proportions[ , 1:(i-1), drop = F])) * partial_stick_proportions[, i, drop = F]
  full_stick_proportions
}

explore_stick_breaking <- function(n, target_median, target_lower, target_upper, width = 0.9) {
  popsize <- sum(target_median)
  lower_p <- (1 - width) / 2
  upper_p <- 1 - ((1 - width) / 2)

  compartments_original_order <- names(target_median)
  target_median_original_order <- target_median
  target_lower_original_order <- target_lower
  target_upper_original_order <- target_upper

  n_compartments <- length(compartments_original_order)
  missing_compartment_index <- which(is.na(target_lower))

  compartments <- c(names(sort(target_median[-missing_compartment_index], decreasing = T)), compartments_original_order[missing_compartment_index])
  target_median <- target_median[compartments]
  target_lower <- target_lower[compartments]
  target_upper <- target_upper[compartments]

  stick_means <- head(target_median / rev(cumsum(rev(target_median))), -1)
  logit_stick_means <- logit(stick_means)

  stick_sd_upper_limit <- sapply(1:(n_compartments - 1), function(i) optimize(f = function(x) norm(as.matrix(expit(qnorm(p = c(lower_p, upper_p), mean = logit_stick_means[[i]], sd = x)) * popsize - c(target_lower_original_order[-missing_compartment_index][[i]], target_upper_original_order[-missing_compartment_index][[i]])), type = "f"), interval = c(0, 5))$minimum)
  names(stick_sd_upper_limit) <- names(stick_means)

  stick_sds <- numeric(n_compartments - 1)
  names(stick_sds) <- names(stick_means)
  stick_sds[1] <- stick_sd_upper_limit[1]
  to_allocate <- (1 - expit(rnorm(n, logit_stick_means[1], stick_sd_upper_limit[1]))) * popsize

  for (i in 2:(n_compartments - 1)) {
    stick_sds[i] <- optimize(f = function(x) norm(as.matrix(quantile(to_allocate * expit(rnorm(n, logit_stick_means[i], x)), c(lower_p, upper_p)) - c(target_lower[i], target_upper[i]))), lower = 0, upper = stick_sd_upper_limit[i])$minimum
    to_allocate <- (1 - expit(rnorm(n, logit_stick_means[i], stick_sds[i]))) * to_allocate
  }

  samples <- rsbln(n = n, compartments, logit_stick_means, stick_sds) * popsize

  reduce(.x =
           list(enframe(compartments_original_order, name = NULL, value = "name"),
                enframe(logit_stick_means, value = "stick_mean"),
                enframe(stick_sds, value = "stick_sd"),
                enframe(target_median, value = "target_median"),
                samples %>% as_tibble() %>% pivot_longer(everything()) %>% group_by(name) %>% summarize(observed_median = quantile(value, 0.5), .groups = "drop"),
                enframe(target_lower, value = "target_lower"),
                samples %>% as_tibble() %>% pivot_longer(everything()) %>% group_by(name) %>% summarize(observed_lower = quantile(value, lower_p), .groups = "drop"),
                enframe(target_upper, value = "target_upper"),
                samples %>% as_tibble() %>% pivot_longer(everything()) %>% group_by(name) %>% summarize(observed_upper = quantile(value, upper_p), .groups = "drop")),
         .f = full_join, by = "name")
}

################################################################################




summarize_stick_breaking <- function(n, target_median, target_lower, target_upper, width = 0.9) {
  explore_stick_breaking(n, target_median, target_lower, target_upper, width) %>%
    select(compartment = name, starts_with("target"), starts_with("observed")) %>%
    pivot_longer(-compartment) %>%
    separate(col = name, into = c("a", "target_type")) %>%
    pivot_wider(names_from = a) %>%
    mutate(rel_diff = abs(observed - target) / target) %>%
    arrange(desc(rel_diff))
}


summarize_stick_breaking_bad <- function(n, target_median, target_lower, target_upper, width = 0.9) {
  explore_stick_breaking_bad(n, target_median, target_lower, target_upper, width) %>%
    select(compartment = name, starts_with("target"), starts_with("observed")) %>%
    pivot_longer(-compartment) %>%
    separate(col = name, into = c("a", "target_type")) %>%
    pivot_wider(names_from = a) %>%
    mutate(rel_diff = abs(observed - target) / target) %>%
    arrange(desc(rel_diff))
}

target_median_full = c(S = 2739459, E = 17488, I = 52463, R = 365205, D = 1077)
target_lower_full = c(S = 2680000, E = 10000, I = 20000, R = 3e+05, D = 1000)
target_upper_full = c(S = 2770000, E = 40000, I = 120000, R = 4e+05, D = 1200)


target_lower <- target_lower_full
target_lower[1] <- NA

target_upper <- target_upper_full
target_upper[1] <- NA

explore_stick_breaking(n = 1e5, target_median_full, target_lower, target_upper)
summarize_stick_breaking(n = 1e5, target_median = target_median_full, target_lower, target_upper)

tmp <- map_dfr(seq_along(target_median), ~{
  target_median <- target_median_full

  target_lower <- target_lower_full
  target_lower[.] <- NA

  target_upper <- target_upper_full
  target_upper[.] <- NA
  missing_var <- names(target_median)[.]

  summarize_stick_breaking(n = 1e6,
                           target_median_full,
                           target_lower_full,
                           target_upper_full,
                           width = 0.9) %>%
    mutate(missing_var = missing_var)
})

tmp %>% arrange(desc(rel_diff))

tmp_bad <- map_dfr(seq_along(target_median), ~{
  target_median <- target_median_full

  target_lower <- target_lower_full
  target_lower[.] <- NA

  target_upper <- target_upper_full
  target_upper[.] <- NA
  missing_var <- names(target_median)[.]

  summarize_stick_breaking_bad(n = 1e6,
                           target_median,
                           target_lower,
                           target_upper,
                           width = 0.9) %>%
    mutate(missing_var = missing_var)
})

tmp %>% arrange(desc(rel_diff))
tmp_bad %>% arrange(desc(rel_diff))

summarize_stick_breaking(n = 1e6,
                         target_median = c(S = 2739459, E = 17488, I = 52463, R = 365205, D = 1077),
                         target_lower = c(S = NA, E = 10000, I = 20000, R = 3e+05, D = 1000),
                         target_upper = c(S = NA, E = 40000, I = 120000, R = 4e+05, D = 1200),
                         width = 0.9)

# Way different results
summarize_stick_breaking(n = 1e6,
                         target_median = c(S = 2739459, E = 17488, I = 52463, R = 365205, D = 1077),
                         target_lower = c(S = 2659325, E = 10000, I = NA, R = 3e+05, D = 1000),
                         target_upper = c(S = 2790651, E = 40000, I = NA, R = 4e+05, D = 1200),
                         width = 0.9)

