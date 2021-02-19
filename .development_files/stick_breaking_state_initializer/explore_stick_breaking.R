library(stemr)
library(tidyverse)
library(tidybayes)

sample_inits <- function(n, compartments, stick_means, stick_sds) {
  p <- length(compartments)
  partial_stick_proportions <- map2_dfr(logit(stick_means), stick_sds, ~expit(rnorm(n, .x, .y))) %>% as.matrix() %>% cbind(1)
  full_stick_proportions <- matrix(nrow = n, ncol = p)
  full_stick_proportions[, 1] = partial_stick_proportions[, 1]
  for (i in 2:p) full_stick_proportions[, i] <- (1 - rowSums(full_stick_proportions[ , 1:(i-1), drop = F])) * partial_stick_proportions[, i, drop = F]
  colnames(full_stick_proportions) <- compartments
  as_tibble(full_stick_proportions * popsize)
}

stick_breaking_explorer <- function(target_median, target_lower, target_upper, n = 1000, width = 0.9) {
  popsize <- sum(target_median)
  lower_p <- (1 - width) / 2
  upper_p <- 1 - ((1 - width) / 2)

  compartments_original_order <- names(target_median)
  target_median_original_order <- target_median

  compartments <- unique(c(names(target_lower), compartments))
  target_median <- target_median[compartments]

  stick_means <- head(target_median / rev(cumsum(rev(target_median))), -1)
  logit_stick_means <- logit(stick_means)
  optimized_sd <- sapply(1:length(target_lower), function(i) optimize(f = function(x) norm(as.matrix(expit(qnorm(p = c(lower_p, upper_p), mean = logit_stick_means[[i]], sd = x)) * popsize - c(target_lower[[i]], target_upper[[i]])), type = "f"), interval = c(0, 2))$minimum) %>%
    `names<-`(names(stick_means))

  samples <- sample_inits(n, compartments, stick_means, optimized_sd) %>% select(!!(compartments_original_order))

  samples_plot <-
    samples %>%
    pivot_longer(everything()) %>%
    mutate(name = fct_inorder(name)) %>%
    ggplot(aes(value, fill = name)) +
    facet_wrap(. ~ name, scale = "free_x") +
    stat_halfeye(normalize = "xy", show.legend = F) +
    cowplot::theme_minimal_grid() +
    scale_x_continuous(name = "People", labels = scales::comma) +
    ylab(NULL)

  print(samples_plot)

  reduce(.x =
           list(enframe(compartments_original_order, name = NULL, value = "name"),
                enframe(logit_stick_means, value = "stick_mean"),
                enframe(optimized_sd, value = "stick_sd"),
                enframe(target_lower, value = "target_lower"),
                samples %>% pivot_longer(everything()) %>% group_by(name) %>% summarize(observed_lower = quantile(value, lower_p), .groups = "drop"),
                enframe(target_upper, value = "target_upper"),
                samples %>% pivot_longer(everything()) %>% group_by(name) %>% summarize(observed_upper = quantile(value, upper_p), .groups = "drop")),
         .f = full_join, by = "name")
}


target_median <- c(S = 2739459, E = 17488, I = 52463, R = 365205, D = 1077)
target_lower <- c(E = 10000, I = 20000, R = 3e+05, D = 1000)
target_upper <- c(E = 40000, I = 120000, R = 4e+05, D = 1200)
width <- 0.9

tmp <- stick_breaking_explorer(target_median, target_lower, target_upper, n = 5000, width = 0.9)
tmp$stick_sd %>% dput()
