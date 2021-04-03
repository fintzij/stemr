library(stemr)
library(tidyverse)


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
  stick_means <- logit(stick_means)


  stick_sd_upper_limit <- sapply(1:(n_compartments - 1), function(i) optimize(f = function(x) norm(as.matrix(expit(qnorm(p = c(lower_p, upper_p), mean = stick_means[[i]], sd = x)) * popsize - c(target_lower[i], target_upper[i])), type = "f"), interval = c(0, 5))$minimum)
  names(stick_sd_upper_limit) <- names(stick_means)

  stick_sds <- numeric(n_compartments - 1)
  names(stick_sds) <- names(stick_means)
  stick_sds[1] <- stick_sd_upper_limit[1]
  to_allocate <- (1 - expit(rnorm(n, stick_means[1], stick_sd_upper_limit[1]))) * popsize

  for (i in 2:(n_compartments - 1)) {
    stick_sds[i] <- optimize(f = function(x) norm(as.matrix(quantile(to_allocate * expit(rnorm(n, stick_means[i], x)), c(lower_p, upper_p)) - c(target_lower[i], target_upper[i]))), lower = 0, upper = stick_sd_upper_limit[i])$minimum
    to_allocate <- (1 - expit(rnorm(n, stick_means[i], stick_sds[i]))) * to_allocate
  }

  samples <- rsbln(n = n, stick_means, stick_sds, stick_size = popsize) %>%
    `colnames<-`(compartments)

  reduce(.x =
           list(enframe(compartments, name = NULL, value = "name"),
                enframe(stick_means, value = "stick_mean"),
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


target_median_full = c(S = 2739459, E = 17488, I = 52463, R = 365205, D = 1077)
target_lower_full = c(S = 2680000, E = 10000, I = 20000, R = 3e+05, D = 1000)
target_upper_full = c(S = 2770000, E = 40000, I = 120000, R = 4e+05, D = 1200)


target_lower <- target_lower_full
target_lower[1] <- NA

target_upper <- target_upper_full
target_upper[1] <- NA
popsize <- sum(target_median_full)

a <- explore_stick_breaking(n = 1e6, target_median_full, target_lower, target_upper)
dput(a$name)
dput(target_median_full[a$name])
dput(c(head(a$stick_mean, -1), head(a$stick_sd, -1)))

summarize_stick_breaking(n = 1e6, target_median = target_median_full, target_lower, target_upper)

tmp <-
  map_dfr(seq_along(target_median_full), ~{
    target_lower <- target_lower_full
    target_lower[.] <- NA

    target_upper <- target_upper_full
    target_upper[.] <- NA
    missing_var <- names(target_median_full)[.]

    summarize_stick_breaking(n = 1e6,
                             target_median_full,
                             target_lower,
                             target_upper,
                             width = 0.9) %>%
      mutate(missing_var = missing_var)
})

tmp %>%
  group_by(missing_var) %>%
  arrange(desc(rel_diff))


target_median_full_02 = c(S = 2795419, E = 3498, I = 10493, R = 365205, D = 1077)
target_lower_02       = c(S = NA,     E = 2000,  I = 4000, R = 3e+05, D = 1000)
target_upper_02       = c(S = NA,     E = 8000,  I = 24000, R = 4e+05, D = 1200)


a <- explore_stick_breaking(n = 1e6, target_median_full_02, target_lower_02, target_upper_02)
dput(a$name)
dput(target_median_full[a$name])
dput(c(head(a$stick_mean, -1), head(a$stick_sd, -1)))



# March Conditions --------------------------------------------------------
target_median_full_03 = c(S = 3168000, E = 2690, I = 5000, R = 1, D = 1)
target_lower_03       = c(S = 3157000, E = NA,  I = 1500, R = 0, D = 0)
target_upper_03       = c(S = 3174000, E = NA,  I = 13000, R = 2, D = 2)


# popsize <- 3175692
# S <- rbeta(n = 2000, 983,2.7) * popsize
# E_plus_I <- popsize - S
# # I
# quantile(rbeta(n = 2000, 41.3, 17.26) * E_plus_I, probs = c(0.05, 0.5, 0.95))
# # E
# quantile(E_plus_I - (rbeta(n = 2000, 41.3, 17.26) * E_plus_I), probs = c(0.05, 0.5, 0.95))

a <- explore_stick_breaking(n = 1e6, target_median_full_03, target_lower_03, target_upper_03)

summarize_stick_breaking(n = 1e6, target_median = target_median_full_03, target_lower_03, target_upper_03)
dput(a$name)
dput(target_median_full_03[a$name])
dput(c(head(a$stick_mean, -1), head(a$stick_sd, -1)))
