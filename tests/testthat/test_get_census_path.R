library(stemr)

context("Obtaining the compartment counts at census times")

test_that("Compartment counts at census times are correctly obtained", {

        set.seed(12511)
        compartments <- c("S","I","R")
        rates <- list(rate("beta * I", "S", "I"),
                      rate("mu", "I", "R"))
        state_initializer <- stem_initializer(c(S = 15, I = 2, R = 0), fixed = TRUE)
        parameters <- c(beta = 0.02, mu = 1/7, rho = 0.5)
        tcovar <- NULL
        constants <- NULL
        strata <- NULL
        t0 <- 0; tmax <- 52

        dynamics <- stem_dynamics(rates = rates, parameters = parameters, state_initializer = state_initializer, tmax = tmax, compartments=compartments, strata = strata, tcovar = tcovar, messages = FALSE)

        stem_object <- stem(dynamics = dynamics)

        path <- simulate_stem(stem_object, paths = TRUE, tmax = 10, nsim = 1, messages = FALSE)$paths[[1]]

        cens_path <-matrix(c(0,15,2,0,
                             1,15,2,0,
                             2,15,2,0,
                             3,15,2,0,
                             4,15,2,0,
                             5,13,3,1,
                             6,13,2,2,
                             7,12,3,2,
                             8,12,3,2,
                             9,12,2,3,
                             10,10,3,4), ncol = 4, byrow = TRUE)

        expect_identical(build_census_path(path, 0:10, 2:4), cens_path)
})