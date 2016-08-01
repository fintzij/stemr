library(stemr)

context("Obtaining the compartment counts at census times")

test_that("Compartment counts at census times are correctly obtained", {

        set.seed(12511)
        compartments <- c("S","I","R")
        rates <- list(rate("beta * I", "S", "I"),
                      rate("mu", "I", "R"))
        state_initializer <- stem_initializer(c(S = 20, I = 2, R = 0), fixed = FALSE)
        parameters <- c(beta = 0.02, mu = 1/7, rho = 0.5)
        tcovar <- NULL
        constants <- NULL
        strata <- NULL
        t0 <- 0; tmax <- 52

        dynamics <- stem_dynamics(rates = rates, parameters = parameters, state_initializer = state_initializer, tmax = tmax, compartments=compartments, strata = strata, tcovar = tcovar, messages = FALSE)

        stem_object <- stem(dynamics = dynamics)

        path <- simulate_stem(stem_object, paths = TRUE, tmax = 10, messages = FALSE)$paths[[1]]

        cens_path <-matrix(c(0,20,2,0,
                             1,18,4,0,
                             2,16,5,1,
                             3,14,6,2,
                             4,14,5,3,
                             5,14,5,3,
                             6,13,6,3,
                             7,11,7,4,
                             8,10,8,4,
                             9,9,7,6,
                             10,8,7,7), ncol = 4, byrow = TRUE)

        expect_identical(build_census_path(path, 0:10, 2:4), cens_path)
})