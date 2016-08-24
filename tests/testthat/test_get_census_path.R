library(stemr)

context("Obtaining the compartment counts at census times")

test_that("Compartment counts at census times are correctly obtained", {
        skip_on_cran()

        set.seed(12511)
        compartments <- c("S","I","R")
        rates <- list(rate("beta * I", "S", "I"),
                      rate("mu", "I", "R"))
        state_initializer <- stem_initializer(c(S = 10, I = 2, R = 0), fixed = TRUE)
        parameters <- c(beta = 0.02, mu = 1/7, rho = 0.5)
        tcovar <- NULL
        constants <- NULL
        strata <- NULL
        t0 <- 0; tmax <- 52

        dynamics <- stem_dynamics(rates = rates, parameters = parameters, state_initializer = state_initializer, tmax = tmax, compartments=compartments, strata = strata, tcovar = tcovar, compile_lna = F, messages = FALSE)

        stem_object <- stem(dynamics = dynamics)

        path <- simulate_stem(stem_object, paths = TRUE, tmax = 10, nsim = 1, messages = FALSE)$paths[[1]]

        cens_path <-matrix(c(0,10,2,0,
                             1,10,1,1,
                             2,9,2,1,
                             3,8,1,3,
                             4,7,1,4,
                             5,7,1,4,
                             6,7,1,4,
                             7,7,1,4,
                             8,6,2,4,
                             9,5,3,4,
                             10,5,3,4), ncol = 4, byrow = TRUE)

        expect_identical(build_census_path(path, 0:10, 2:4), cens_path)
})