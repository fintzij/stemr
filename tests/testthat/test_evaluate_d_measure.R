library(stemr)

context("Evaluating the density of the measurement process")

test_that("Prevalence counts in a simple system are measured correctly", {
        skip_on_cran()

        set.seed(12511)
        compartments <- c("S","I","R")
        rates <- list(rate("beta * I", "S", "I"),
                      rate("mu", "I", "R"))
        state_initializer <- stem_initializer(c(S = 100, I = 10, R = 5), fixed = FALSE)
        parameters <- c(beta = 0.1, mu = 1/7, rho = 0.5)
        tcovar <- NULL
        constants <- NULL
        strata <- NULL
        t0 <- 0; tmax <- 52

        dynamics <- stem_dynamics(rates = rates, parameters = parameters, state_initializer = state_initializer, compartments=compartments, strata = strata, tcovar = tcovar, t0 = t0, compile_lna = F, tmax = 52, messages = FALSE)

# Poisson counts ----------------------------------------------------------

        emissions <- list(emission("I", "poisson", c("I * rho"), incidence = FALSE, obstimes = seq(0,20,by=2)))
        meas_process <- stem_measure(emissions = emissions, dynamics = dynamics, messages = FALSE)

        stem_object <- stem(dynamics = dynamics, measurement_process = meas_process)

        stemr_sim <- simulate_stem(stem_object, paths = TRUE, observations = T, tmax = 20, nsim = 1, paths_as_array = TRUE,
                                   datasets_as_array = TRUE, census_times = seq(0,20,by=2))

        emit_mat_stemr <- emit_mat_manual <- matrix(c(seq(0,20,by=2), seq(0, 20, by = 2)), ncol = 2)

        evaluate_d_measure(emitmat          = emit_mat_stemr,
                           obsmat           = stemr_sim$datasets[,,1],
                           statemat         = stemr_sim$paths[,,1],
                           measproc_indmat  = stem_object$measurement_process$measproc_indmat,
                           parameters       = stem_object$dynamics$parameters,
                           constants        = stem_object$dynamics$constants,
                           tcovar_censusmat = stem_object$measurement_process$tcovar_censmat,
                           d_meas_ptr       = stem_object$measurement_process$meas_pointers$d_measure_ptr)

        for(j in seq_len(nrow(emit_mat_manual))) {
                emit_mat_manual[j, 2] <- dpois(stemr_sim$datasets[j,2,1], stem_object$dynamics$parameters["rho"] * stemr_sim$paths[j,3,1], log = T)
        }

        expect_identical(emit_mat_stemr, emit_mat_manual)


# Binomial counts ---------------------------------------------------------

        emissions <- list(emission("I", "binomial", c("I", "rho"), incidence = FALSE, obstimes = seq(0,20,by=2)))
        meas_process <- stem_measure(emissions = emissions, dynamics = dynamics, messages = FALSE)

        stem_object <- stem(dynamics = dynamics, measurement_process = meas_process)

        stemr_sim <- simulate_stem(stem_object, paths = TRUE, observations = T, tmax = 20, nsim = 1, paths_as_array = TRUE,
                                   datasets_as_array = TRUE, census_times = seq(0,20,by=2))

        emit_mat_stemr <- emit_mat_manual <- matrix(c(seq(0,20,by=2), seq(0, 20, by = 2)), ncol = 2)

        evaluate_d_measure(emitmat          = emit_mat_stemr,
                           obsmat           = stemr_sim$datasets[,,1],
                           statemat         = stemr_sim$paths[,,1],
                           measproc_indmat  = stem_object$measurement_process$measproc_indmat,
                           parameters       = stem_object$dynamics$parameters,
                           constants        = stem_object$dynamics$constants,
                           tcovar_censusmat = stem_object$measurement_process$tcovar_censmat,
                           d_meas_ptr       = stem_object$measurement_process$meas_pointers$d_measure_ptr)

        for(j in seq_len(nrow(emit_mat_manual))) {
                emit_mat_manual[j, 2] <- dbinom(stemr_sim$datasets[j,2,1], stemr_sim$paths[j,3,1], stem_object$dynamics$parameters["rho"], log = T)
        }

        expect_identical(emit_mat_stemr, emit_mat_manual)

# Negative binomial counts ---------------------------------------------------------

        emissions <- list(emission("I", "negbinomial", c("I", "rho"), incidence = FALSE, obstimes = seq(0,20,by=2)))
        meas_process <- stem_measure(emissions = emissions, dynamics = dynamics, messages = FALSE)

        stem_object <- stem(dynamics = dynamics, measurement_process = meas_process)

        stemr_sim <- simulate_stem(stem_object, paths = TRUE, observations = T, tmax = 20, nsim = 1, paths_as_array = TRUE,
                                   datasets_as_array = TRUE, census_times = seq(0,20,by=2))

        emit_mat_stemr <- emit_mat_manual <- matrix(c(seq(0,20,by=2), seq(0, 20, by = 2)), ncol = 2)

        evaluate_d_measure(emitmat          = emit_mat_stemr,
                           obsmat           = stemr_sim$datasets[,,1],
                           statemat         = stemr_sim$paths[,,1],
                           measproc_indmat  = stem_object$measurement_process$measproc_indmat,
                           parameters       = stem_object$dynamics$parameters,
                           constants        = stem_object$dynamics$constants,
                           tcovar_censusmat = stem_object$measurement_process$tcovar_censmat,
                           d_meas_ptr       = stem_object$measurement_process$meas_pointers$d_measure_ptr)

        for(j in seq_len(nrow(emit_mat_manual))) {
                emit_mat_manual[j, 2] <- dnbinom(stemr_sim$datasets[j,2,1], stemr_sim$paths[j,3,1], stem_object$dynamics$parameters["rho"], log = T)
        }

        expect_identical(emit_mat_stemr, emit_mat_manual)

# Normally distributed densities counts ---------------------------------------------------------

        emissions <- list(emission("I", "gaussian", c("I", "rho"), incidence = FALSE, obstimes = seq(0,20,by=2)))
        meas_process <- stem_measure(emissions = emissions, dynamics = dynamics, messages = FALSE)

        stem_object <- stem(dynamics = dynamics, measurement_process = meas_process)

        stemr_sim <- simulate_stem(stem_object, paths = TRUE, observations = T, tmax = 20, nsim = 1, paths_as_array = TRUE,
                                   datasets_as_array = TRUE, census_times = seq(0,20,by=2))

        emit_mat_stemr <- emit_mat_manual <- matrix(c(seq(0,20,by=2), seq(0, 20, by = 2)), ncol = 2)

        evaluate_d_measure(emitmat          = emit_mat_stemr,
                           obsmat           = stemr_sim$datasets[,,1],
                           statemat         = stemr_sim$paths[,,1],
                           measproc_indmat  = stem_object$measurement_process$measproc_indmat,
                           parameters       = stem_object$dynamics$parameters,
                           constants        = stem_object$dynamics$constants,
                           tcovar_censusmat = stem_object$measurement_process$tcovar_censmat,
                           d_meas_ptr       = stem_object$measurement_process$meas_pointers$d_measure_ptr)

        for(j in seq_len(nrow(emit_mat_manual))) {
                emit_mat_manual[j, 2] <- dnorm(stemr_sim$datasets[j,2,1], stemr_sim$paths[j,3,1], stem_object$dynamics$parameters["rho"], log = T)
        }

        expect_identical(emit_mat_stemr, emit_mat_manual)

})

test_that("Prevalence and incidence counts in a stratified system with two measurement processes are measured correctly", {
        skip_on_cran()

        set.seed(12511)
        strata <- c("male", "female")
        compartments <- list(S = "ALL", I = "ALL", R = "ALL")
        rates <- list(rate("beta * I_SELF + gamma * comp_fcn(I_ADJ, sum)", "S", "I", "ALL", incidence = TRUE,
                           seasonality(period = 52, common_seasonality = TRUE, intercept = 5e-5, trend = 1e-4, s_params = c(0.5, 0.575e-4))),
                      rate("mu", "I", "R", "ALL"),
                      rate("phi", "R", "S", "ALL"))
        adjacency <- matrix(c(0,1,1,0), nrow = 2, ncol = 2); colnames(adjacency) <- rownames(adjacency) <- strata
        initial_state_vec <- c(S = 900, I = 10, R = 90)
        state_initializer <- list(stem_initializer(initial_state_vec, fixed = TRUE, strata = c("male", "female"), shared_params = TRUE))
        parameters <- c(beta = 0.000025, gamma = 0.0005, mu = 1/2, phi = 5, rho = 0.5)
        tcovar <- NULL
        constants <- NULL
        t0 <- 0; tmax <- 52

        dynamics <- stem_dynamics(rates = rates, parameters = parameters, state_initializer = state_initializer, compartments=compartments, strata = strata, tcovar = tcovar, t0 = 0, tmax = 52, adjacency = adjacency, compile_lna = F, messages = FALSE)

# Poisson counts ----------------------------------------------------------

        emissions <- list(emission("I_male", "poisson", c("I_male * rho"), incidence = TRUE, obstimes = seq(0,20,by=2)),
                          emission("I_female", "poisson", c("I_female * rho"), incidence = FALSE, obstimes = seq(1,21, by = 2)))
        meas_process <- stem_measure(emissions = emissions, dynamics = dynamics, messages = FALSE)

        stem_object <- stem(dynamics = dynamics, measurement_process = meas_process)

        stemr_sim <- simulate_stem(stem_object, paths = TRUE, observations = T, tmax = 20, nsim = 1, paths_as_array = TRUE,
                                   datasets_as_array = TRUE, census_times = 0:21)

        emit_mat_stemr <- emit_mat_manual <- matrix(0, nrow = nrow(stemr_sim$paths), ncol = 3)
        emit_mat_stemr[,1] <- emit_mat_manual[,1] <- 0:21

        evaluate_d_measure(emitmat          = emit_mat_stemr,
                           obsmat           = stemr_sim$datasets[,,1],
                           statemat         = stemr_sim$paths[,,1],
                           measproc_indmat  = stem_object$measurement_process$measproc_indmat,
                           parameters       = stem_object$dynamics$parameters,
                           constants        = stem_object$dynamics$constants,
                           tcovar_censusmat = stem_object$measurement_process$tcovar_censmat,
                           d_meas_ptr       = stem_object$measurement_process$meas_pointers$d_measure_ptr)

        for(j in seq_len(nrow(emit_mat_manual))) {
                if(j %% 2) {
                        emit_mat_manual[j, 2] <- dpois(stemr_sim$datasets[j,2,1], stem_object$dynamics$parameters["rho"] * stemr_sim$paths[j,8,1], log = T)
                } else {
                        emit_mat_manual[j, 3] <- dpois(stemr_sim$datasets[j,3,1], stem_object$dynamics$parameters["rho"] * stemr_sim$paths[j,5,1], log = T)
                }
        }

        expect_identical(emit_mat_stemr, emit_mat_manual)


# Binomial counts ---------------------------------------------------------

        emissions <- list(emission("I_male", "binomial", c("I_male", "rho"), incidence = TRUE, obstimes = seq(0,20,by=2)),
                          emission("I_female", "binomial", c("I_female", "rho"), incidence = FALSE, obstimes = seq(1,21, by = 2)))
        meas_process <- stem_measure(emissions = emissions, dynamics = dynamics, messages = FALSE)

        stem_object <- stem(dynamics = dynamics, measurement_process = meas_process)

        stemr_sim <- simulate_stem(stem_object, paths = TRUE, observations = T, tmax = 20, nsim = 1, paths_as_array = TRUE,
                                   datasets_as_array = TRUE, census_times = 0:21)

        emit_mat_stemr <- emit_mat_manual <- matrix(0, nrow = nrow(stemr_sim$paths), ncol = 3)
        emit_mat_stemr[,1] <- emit_mat_manual[,1] <- 0:21

        evaluate_d_measure(emitmat          = emit_mat_stemr,
                           obsmat           = stemr_sim$datasets[,,1],
                           statemat         = stemr_sim$paths[,,1],
                           measproc_indmat  = stem_object$measurement_process$measproc_indmat,
                           parameters       = stem_object$dynamics$parameters,
                           constants        = stem_object$dynamics$constants,
                           tcovar_censusmat = stem_object$measurement_process$tcovar_censmat,
                           d_meas_ptr       = stem_object$measurement_process$meas_pointers$d_measure_ptr)

        for(j in seq_len(nrow(emit_mat_manual))) {
                if(j %% 2) {
                        emit_mat_manual[j, 2] <- dbinom(stemr_sim$datasets[j,2,1], stemr_sim$paths[j,8,1], stem_object$dynamics$parameters["rho"], log = T)
                } else {
                        emit_mat_manual[j, 3] <- dbinom(stemr_sim$datasets[j,3,1], stemr_sim$paths[j,5,1], stem_object$dynamics$parameters["rho"], log = T)
                }
        }

        expect_identical(emit_mat_stemr, emit_mat_manual)

# Negative binomial counts ---------------------------------------------------------

        emissions <- list(emission("I_male", "negbinomial", c("I_male", "rho"), incidence = TRUE, obstimes = seq(0,20,by=2)),
                          emission("I_female", "negbinomial", c("I_female", "rho"), incidence = FALSE, obstimes = seq(1,21, by = 2)))
        meas_process <- stem_measure(emissions = emissions, dynamics = dynamics, messages = FALSE)

        stem_object <- stem(dynamics = dynamics, measurement_process = meas_process)

        stemr_sim <- simulate_stem(stem_object, paths = TRUE, observations = T, tmax = 20, nsim = 1, paths_as_array = TRUE,
                                   datasets_as_array = TRUE, census_times = 0:21)

        emit_mat_stemr <- emit_mat_manual <- matrix(0, nrow = nrow(stemr_sim$paths), ncol = 3)
        emit_mat_stemr[,1] <- emit_mat_manual[,1] <- 0:21

        evaluate_d_measure(emitmat          = emit_mat_stemr,
                           obsmat           = stemr_sim$datasets[,,1],
                           statemat         = stemr_sim$paths[,,1],
                           measproc_indmat  = stem_object$measurement_process$measproc_indmat,
                           parameters       = stem_object$dynamics$parameters,
                           constants        = stem_object$dynamics$constants,
                           tcovar_censusmat = stem_object$measurement_process$tcovar_censmat,
                           d_meas_ptr       = stem_object$measurement_process$meas_pointers$d_measure_ptr)

        for(j in seq_len(nrow(emit_mat_manual))) {
                if(j %% 2) {
                        emit_mat_manual[j, 2] <- dnbinom(stemr_sim$datasets[j,2,1], stemr_sim$paths[j,8,1], stem_object$dynamics$parameters["rho"], log = T)
                } else {
                        emit_mat_manual[j, 3] <- dnbinom(stemr_sim$datasets[j,3,1], stemr_sim$paths[j,5,1], stem_object$dynamics$parameters["rho"], log = T)
                }
        }

        expect_identical(emit_mat_stemr, emit_mat_manual)

# Normally distributed densities counts ---------------------------------------------------------

        emissions <- list(emission("I_male", "gaussian", c("I_male", "rho"), incidence = TRUE, obstimes = seq(0,20,by=2)),
                          emission("I_female", "gaussian", c("I_female", "rho"), incidence = FALSE, obstimes = seq(1,21, by = 2)))
        meas_process <- stem_measure(emissions = emissions, dynamics = dynamics, messages = FALSE)

        stem_object <- stem(dynamics = dynamics, measurement_process = meas_process)

        stemr_sim <- simulate_stem(stem_object, paths = TRUE, observations = T, tmax = 20, nsim = 1, paths_as_array = TRUE,
                                   datasets_as_array = TRUE, census_times = 0:21)

        emit_mat_stemr <- emit_mat_manual <- matrix(0, nrow = nrow(stemr_sim$paths), ncol = 3)
        emit_mat_stemr[,1] <- emit_mat_manual[,1] <- 0:21

        evaluate_d_measure(emitmat          = emit_mat_stemr,
                           obsmat           = stemr_sim$datasets[,,1],
                           statemat         = stemr_sim$paths[,,1],
                           measproc_indmat  = stem_object$measurement_process$measproc_indmat,
                           parameters       = stem_object$dynamics$parameters,
                           constants        = stem_object$dynamics$constants,
                           tcovar_censusmat = stem_object$measurement_process$tcovar_censmat,
                           d_meas_ptr       = stem_object$measurement_process$meas_pointers$d_measure_ptr)

        for(j in seq_len(nrow(emit_mat_manual))) {
                if(j %% 2) {
                        emit_mat_manual[j, 2] <- dnorm(stemr_sim$datasets[j,2,1], stemr_sim$paths[j,8,1], stem_object$dynamics$parameters["rho"], log = T)
                } else {
                        emit_mat_manual[j, 3] <- dnorm(stemr_sim$datasets[j,3,1], stemr_sim$paths[j,5,1], stem_object$dynamics$parameters["rho"], log = T)
                }
        }

        expect_identical(emit_mat_stemr, emit_mat_manual)

})