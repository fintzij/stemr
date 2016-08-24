library(stemr)

context("Calling rate functions")

# tests relate to CALL_RATE_FCN
test_that("Rates are computed properly for a simple system", {
        skip_on_cran()

        # SIR example
        compartments <- c("S","I","R")
        rates <- list(rate("beta * I", "S", "I"),
                      rate("mu", "I", "R"))
        state_initializer <- stem_initializer(c(S = 10, I = 1, R = 5), fixed = FALSE)
        parameters <- c(beta = 1, mu = 1/7, rho = 0.5)
        tcovar <- NULL
        constants <- NULL
        strata <- NULL
        t0 <- 0; tmax <- 52

        # set dynamics
        dynamics <- stem_dynamics(rates = rates, parameters = parameters, tmax = tmax, state_initializer = state_initializer, compartments=compartments, strata = strata, tcovar = tcovar, compile_lna = F, compile_ode = F, messages = FALSE)

        # initialize stem object
        stem_object <- stem(dynamics = dynamics)

        # initialize model objects
        timestep <- tmax - t0

        stem_object$dynamics$tcovar <- build_tcovar_matrix(tcovar = stem_object$dynamics$tcovar, timestep = timestep, t0 = t0, tmax = tmax)
        stem_object$dynamics$tcovar_codes <- seq_len(ncol(stem_object$dynamics$tcovar) - 1)
        names(stem_object$dynamics$tcovar_codes) <- colnames(stem_object$dynamics$tcovar)[2:ncol(stem_object$dynamics$tcovar)]

        # compile rate functions and get pointers
        rate_ptrs <- parse_rates(rates = stem_object$dynamics$rates, messages = FALSE)

        # objects with which to evaluate the rates
        rates_lumped   <- c(0,0) # lumped rate vector to be updated
        rates_unlumped <- c(0,0) # rate vector to be updated
        inds           <- c(T,T) # rates to update
        states         <- matrix(c(1:10, 10:1, 20:29), ncol = 3)
        tcovar         <- c(0,0)
        constants      <- 16

        for(j in 1:nrow(states)) {
                state <- states[j,]
                CALL_RATE_FCN(rates_lumped, inds, state, parameters, constants, tcovar, rate_ptrs$lumped_ptr) # update lumped rate
                CALL_RATE_FCN(rates_unlumped, inds, state, parameters, constants, tcovar, rate_ptrs$unlumped_ptr) # update unlumped rate

                testthat::expect_equivalent(rates_lumped, c(state[1] * state[2], state[2]/7))
                testthat::expect_equivalent(rates_unlumped, c(state[2], 1/7))
        }
})

test_that("Rates are computed properly for a simple system with two time-varying covariates", {
        skip_on_cran()

        set.seed(12511)
        compartments <- c("S","I","R")
        rates <- list(rate("beta * I", "S", "I", seasonality = seasonality(period = 52, intercept = 0.005, trend = 0.01, s_params = c(SIN_52 = 0.2, COS_52 = 0.3), log = TRUE)),
                      rate("mu * alpha", "I", "R"))
        state_initializer <- stem_initializer(c(S = 10, I = 1, R = 5), fixed = FALSE)
        parameters <- c(beta = 1, mu = 1/7)
        tcovar <- matrix(c(0,10,1, 0.5), ncol = 2); colnames(tcovar) <- c("t", "alpha")
        constants <- NULL
        strata <- NULL
        t0 <- 0; tmax <- 52

        # set dynamics
        dynamics <- stem_dynamics(rates = rates, parameters = parameters, state_initializer = state_initializer, compartments=compartments, strata = strata, tcovar = tcovar, t0 = t0, tmax = tmax, compile_lna = F, messages = FALSE)

        # initialize stem object
        stem_object <- stem(dynamics = dynamics)

        # initialize model objects
        timestep <- (tmax - t0)/50

        stem_object$dynamics$tcovar <- build_tcovar_matrix(tcovar = stem_object$dynamics$tcovar, timestep = timestep, t0 = t0, tmax = tmax)
        stem_object$dynamics$tcovar_codes <- seq_len(ncol(stem_object$dynamics$tcovar) - 1)
        names(stem_object$dynamics$tcovar_codes) <- colnames(stem_object$dynamics$tcovar)[2:ncol(stem_object$dynamics$tcovar)]
        stem_object$dynamics$n_tcovar <- ncol(stem_object$dynamics$tcovar) - 1
        stem_object$dynamics$tcovar_changemat <- build_tcovar_changemat(stem_object$dynamics$tcovar)

        timecode <- which(names(stem_object$dynamics$tcovar_codes) == "TIME")

        # rebuild the time-varying covariate adjacency matrix
        stem_object$dynamics$tcovar_adjmat <- build_tcovar_adjmat(stem_object$dynamics$rates, stem_object$dynamics$tcovar_codes)

        rate_ptrs <- parse_rates(rates = stem_object$dynamics$rates, messages = FALSE)

        # objects with which to evaluate the rates
        rates_lumped   <- c(0,0) # lumped rate vector to be updated
        rates_unlumped <- c(0,0) # rate vector to be updated
        inds           <- c(T,T) # rates to update
        states         <- matrix(c(1:10, 10:1, 20:29), ncol = 3)
        eval_times     <- sort(runif(10, 0, 52))
        time_inds      <- findInterval(eval_times, stem_object$dynamics$tcovar[,1], rightmost.closed = TRUE)
        tcovars        <- stem_object$dynamics$tcovar[time_inds,]
        constants      <- stem_object$dynamics$constants
        parameters     <- stem_object$dynamics$parameters

        for(j in 1:nrow(states)) {
                state <- states[j,]
                tcov <- tcovars[j,]
                CALL_RATE_FCN(rates_lumped, inds, state, parameters, constants, tcov, rate_ptrs$lumped_ptr) # update lumped rate
                CALL_RATE_FCN(rates_unlumped, inds, state, parameters, constants, tcov, rate_ptrs$unlumped_ptr) # update unlumped rate

                testthat::expect_equivalent(rates_lumped,
                                            unname(c((state[2] + exp(0.005 + 0.01 * tcov[3] + 0.2*sin(2*pi*tcov[3]/52) + 0.3*cos(2*pi*tcov[3]/52)))*state[1], state[2]*tcov[2]/7)))
                testthat::expect_equivalent(rates_unlumped,
                                            unname(c(state[2] + exp(0.005 + 0.01 * tcov[3] + 0.2*sin(2*pi*tcov[3]/52) + 0.3*cos(2*pi*tcov[3]/52)), tcov[2]/7)))
        }
})

test_that("Rates are computed properly for a complex system with time varying covariates", {
        skip_on_cran()

        set.seed(12511)
        strata <- interact(sex=c("male", "female"), age=c("young", "old"))
        compartments <- list(S = "ALL", I = "ALL", R = "ALL", D = c("male_old", "female_old"))
        rates <- list(rate("beta * I_SELF + gamma * comp_fcn(I_ADJ, sum)", "S", "I", "ALL",
                           seasonality(period = 52, intercept = 1e-5, trend = 5e-5, common_seasonality = FALSE)),
                      rate("mu", "I", "R", "ALL"),
                      rate("delta", "R", "D", c("male_old", "female_old")))
        state_initializer <- list(stem_initializer(c(S = 900, I = 10, R = 90), fixed = FALSE, strata = c("male_young", "female_young"), shared_params = FALSE),
                                  stem_initializer(c(S = 400, I = 5, R = 500, D = 0), fixed = FALSE, strata = c("male_old", "female_old"), shared_params = FALSE))
        adjacency <- matrix(1, nrow = length(strata), ncol = length(strata)); diag(adjacency) <- 0
        colnames(adjacency) = rownames(adjacency) = strata
        parameters = c(gamma = 1e-4, mu = 1/7, delta = 0.5)
        tcovar <- matrix(c(0,10,1e-3,9e-4), ncol = 2); colnames(tcovar) <- c("t", "beta")
        constants <- NULL
        t0 <- 0; tmax <- 52

        # set dynamics
        dynamics <- stem_dynamics(rates = rates, parameters = parameters,state_initializer = state_initializer, compartments=compartments, strata = strata, adjacency = adjacency, t0 = t0, tmax = tmax, tcovar = tcovar, compile_lna = F, messages = FALSE)

        # initialize stem object
        stem_object <- stem(dynamics = dynamics)

        # initialize model objects
        timestep <- (tmax - t0)/50

        stem_object$dynamics$tcovar <- build_tcovar_matrix(tcovar = stem_object$dynamics$tcovar, timestep = timestep, t0 = t0, tmax = tmax)
        stem_object$dynamics$tcovar_codes <- seq_len(ncol(stem_object$dynamics$tcovar) - 1)
        names(stem_object$dynamics$tcovar_codes) <- colnames(stem_object$dynamics$tcovar)[2:ncol(stem_object$dynamics$tcovar)]
        stem_object$dynamics$n_tcovar <- ncol(stem_object$dynamics$tcovar) - 1
        stem_object$dynamics$tcovar_changemat <- build_tcovar_changemat(stem_object$dynamics$tcovar)

        timecode <- which(names(stem_object$dynamics$tcovar_codes) == "TIME")

        rate_ptrs <- parse_rates(rates = stem_object$dynamics$rates, messages = FALSE)

        # objects with which to evaluate the rates
        rates_lumped   <- rep(0, 10) # lumped rate vector to be updated
        rates_unlumped <- rep(0, 10) # rate vector to be updated
        inds           <- rep(T, 10) # rates to update
        states         <- matrix(c(1:140), ncol = 14)
        eval_times     <- sort(runif(10, 0, 52))
        time_inds      <- findInterval(eval_times, stem_object$dynamics$tcovar[,1], rightmost.closed = TRUE)
        tcovars        <- stem_object$dynamics$tcovar[time_inds,]
        constants      <- stem_object$dynamics$constants
        parameters     <- stem_object$dynamics$parameters

        for(j in 1:nrow(states)) {
                state <- states[j,]
                tcov <- tcovars[j,]
                CALL_RATE_FCN(rates_lumped, inds, state, parameters, constants, tcov, rate_ptrs$lumped_ptr) # update lumped rate
                CALL_RATE_FCN(rates_unlumped, inds, state, parameters, constants, tcov, rate_ptrs$unlumped_ptr) # update unlumped rate

                testthat::expect_equivalent(rates_lumped[8], unname(parameters[2] * state[8]))
                testthat::expect_equivalent(rates_unlumped[1], unname(tcov[2] * state[5] + parameters[1]*sum(state[6:8]) + exp(parameters[4] + parameters[5] * tcov[3] + parameters[6] * sin(2 * constants[1] * tcov[3] / 52) + parameters[7] * cos(2 * constants[1] * tcov[3]/52))))
        }
})
