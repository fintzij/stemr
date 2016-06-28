library(stemr)

context("Calling rate functions")

# tests relate to CALL_RATE_FCN
test_that("Rates are computed properly for a simple system", {

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
        dynamics <- stem_dynamics(rates = rates, parameters = parameters,state_initializer = state_initializer, compartments=compartments, strata = strata, adjacency = adjacency, tcovar = tcovar)

        # initialize stem object
        stem_object <- stem(dynamics = dynamics)

        # initialize model objects
        timestep <- tmax - t0

        stem_object$dynamics$tcovar <- build_tcovar_matrix(tcovar = stem_object$dynamics$tcovar, timestep = timestep, t0 = t0, tmax = tmax)
        stem_object$dynamics$tcovar_codes <- seq_len(ncol(stem_object$dynamics$tcovar) - 1)
        names(stem_object$dynamics$tcovar_codes) <- colnames(stem_object$dynamics$tcovar)[2:ncol(stem_object$dynamics$tcovar)]

        # compile rate functions and get pointers
        rate_ptrs <- parse_rates(rates = stem_object$dynamics$rates,
                                 param_codes = stem_object$dynamics$param_codes,
                                 compartment_codes = stem_object$dynamics$comp_codes,
                                 const_codes = stem_object$dynamics$const_codes,
                                 tcovar_codes = stem_object$dynamics$tcovar_codes)

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

test_that("Rates are computed properly for a simple system with two time-varying covariates (time and one other)", {

        set.seed(12511)
        compartments <- c("S","I","R")
        rates <- list(rate("beta * I + iota * TIME", "S", "I", seasonality = seasonality(period = 52, s_params = c(SIN_52 = 2, COS_52 = 3))),
                      rate("mu * alpha", "I", "R"))
        state_initializer <- stem_initializer(c(S = 10, I = 1, R = 5), fixed = FALSE)
        parameters <- c(beta = 1, iota = 1, mu = 1/7, rho = 0.5)
        tcovar <- matrix(c(0,10,1, 0.5), ncol = 2); colnames(tcovar) <- c("t", "alpha")
        constants <- NULL
        strata <- NULL
        t0 <- 0; tmax <- 52

        # set dynamics
        dynamics <- stem_dynamics(rates = rates, parameters = parameters,state_initializer = state_initializer, compartments=compartments, strata = strata, adjacency = adjacency, tcovar = tcovar)

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

        rate_ptrs <- parse_rates(rates = stem_object$dynamics$rates,
                                 param_codes = stem_object$dynamics$param_codes,
                                 compartment_codes = stem_object$dynamics$comp_codes,
                                 const_codes = stem_object$dynamics$const_codes,
                                 tcovar_codes = stem_object$dynamics$tcovar_codes)

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
                                            unname(c((state[2] + tcov[3] + 2*sin(2*pi*tcov[3]/52) + 3*cos(2*pi*tcov[3]/52))*state[1], state[2]*tcov[2]/7)))
                testthat::expect_equivalent(rates_unlumped,
                                            unname(c(state[2] + tcov[3] + 2*sin(2*pi*tcov[3]/52) + 3*cos(2*pi*tcov[3]/52), tcov[2]/7)))
        }
})
