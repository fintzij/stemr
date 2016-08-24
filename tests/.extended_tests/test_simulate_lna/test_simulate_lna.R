# Extended test for LNA simulation functions. Tested against the 'issb' package and the Gillespie simulation functions.

# Simulation test to check that simulate_lna produces trajectories that
# approximate the MJP trajectories. Three models are simulated: simple SIR
# (simulation 1), SIRS with seasonality (simulation 2), SIRS with seasonality
# and two strata. Each of the simulations will produce a plot comparing the
# trajectories produced by simulate epimodel with the trajectories produced by
# GillespieSSA.

library(stemr)
library(ggplot2)
library(issb)

# Set up simulation objectx -----------------------------------------------

niter <- 500
census_times <- seq(0,52,by=0.1)

if(sim_num == 1){

        # set up the model in the issb package
        h = function(x, pars) {
                hazs = numeric(length(pars))
                hazs[1] = pars[1]*x[1]*x[2] # S->I
                hazs[2] = pars[2]*x[2]      # I->R
                return(hazs)
        }

        smat = matrix(c(-1,1,0,0,-1,1),nrow=3,ncol=2)
        rownames(smat) = c("S", "I", "R")

        initial = c(S = 1500, I = 10, R = 50)

        pars = c(0.00025,1/7)

        SIR_mod_issb <- create_model(smat, h, initial, pars)

        # Set up stemr objects
        compartments <- c("S","I","R")
        rates <- list(rate("beta * I", "S", "I"),
                      rate("mu", "I", "R"))
        state_initializer <- stem_initializer(c(S = 1500, I = 10, R = 50), fixed = TRUE)
        parameters <- c(beta = 0.00025, mu = 1/7, rho = 0.5)
        tcovar <- NULL
        constants <- NULL
        strata <- NULL
        t0 <- 0; tmax <- 52
        lna_scale <- "log"

        dynamics <- stem_dynamics(rates = rates, parameters = parameters,state_initializer = state_initializer, compartments=compartments, strata = strata, tcovar = tcovar, tmax = tmax, lna_scale = lna_scale, messages = TRUE, compile_rates = T, compile_lna = T, compile_ode = F)

        stem_object <- stem(dynamics = dynamics)

        # simulate paths
        sim1_gillespie_trajecs <- simulate_stem(stem_object = stem_object,
                                       nsim = niter,
                                       paths = TRUE,
                                       observations = FALSE,
                                       subject_paths = FALSE,
                                       method = "gillespie",
                                       t0 = t0,
                                       tmax = tmax,
                                       census_times = seq(0, 52, by = 0.1),
                                       paths_as_array = TRUE,
                                       messages = TRUE)$paths

        sim1_lna_norestart <- simulate_stem(stem_object = stem_object,
                                     nsim = niter,
                                     paths = TRUE,
                                     observations = FALSE,
                                     subject_paths = FALSE,
                                     method = "lna",
                                     lna_restart = FALSE,
                                     t0 = t0,
                                     tmax = tmax,
                                     census_times = seq(0, 52, by = 0.1),
                                     paths_as_array = TRUE,
                                     messages = TRUE)$paths

        # sim1_lna_restart <- simulate_stem(stem_object = stem_object,
        #                                     nsim = niter,
        #                                     paths = TRUE,
        #                                     observations = FALSE,
        #                                     subject_paths = FALSE,
        #                                     method = "lna",
        #                                     lna_restart = TRUE,
        #                                     t0 = t0,
        #                                     tmax = tmax,
        #                                     census_times = seq(0, 52, by = 0.1),
        #                                     paths_as_array = TRUE,
        #                                     messages = TRUE)$paths

        # sim1_ode <- stem_ode(stem_object = stem_object,
        #                      census_times = census_times,
        #                      init_states = c(S = 1500, I = 10, R = 50))

        sim1_ode_issb <- deterministic(SIR_mod_issb, 52, 0.1)

        path_quantiles <- data.frame(time = rep(census_times, 3),
                                     method = rep(c("gillespie", "lna_norestart",
                                                    # "lna_restart",
                                                    # "ode",
                                                    "ode_issb"), each = length(census_times)),
                                     mean = c(apply(sim1_gillespie_trajecs[,"I",,drop = FALSE], 1, mean),
                                              apply(sim1_lna_norestart[,"I",,drop = FALSE], 1, mean),
                                              # apply(sim1_lna_restart[,"I",,drop = FALSE], 1, mean),
                                              # sim1_ode[,"I"],
                                              sim1_ode_issb[,"I"]),
                                     lower = 0,
                                     upper = 0)
        path_quantiles$lower <- path_quantiles$mean - 1.96 * c(apply(sim1_gillespie_trajecs[,"I",,drop = FALSE], 1, sd),
                                                               apply(sim1_lna_norestart[,"I",,drop = FALSE], 1, sd),
                                                               # apply(sim1_lna_restart[,"I",,drop = FALSE], 1, sd),
                                                               rep(0, length(census_times))) / sqrt(niter)
        path_quantiles$upper <- path_quantiles$mean + 1.96 * c(apply(sim1_gillespie_trajecs[,"I",,drop = FALSE], 1, sd),
                                                               apply(sim1_lna_norestart[,"I",,drop = FALSE], 1, sd),
                                                               # apply(sim1_lna_restart[,"I",,drop = FALSE], 1, sd),
                                                               rep(0, length(census_times))) / sqrt(niter)
        path_quantiles[path_quantiles$method == "lna_norestart",3:5] = exp(path_quantiles[path_quantiles$method == "lna_norestart",3:5])
        pdf(paste0("test_simulate_lna",sim_num,"popsize_1560", ".pdf"))

        ggplot(path_quantiles, aes(x = time, y = mean, colour = method)) + geom_line(size = 1) + geom_ribbon(data = path_quantiles, aes(x = time, ymin = lower, ymax = upper, fill = method), alpha = 0.1, size = 0) + labs(y = "Average compartment counts", title = "S0 = 1500, I0 = 10, R0 = 50; beta = 0.00025, mu = 1/7") + theme_bw()

        dev.off()

} else if(sim_num == 2){

        # set up the model in the issb package
        h = function(x, pars) {
                hazs = numeric(length(pars))
                hazs[1] = pars[1]*x[1]*x[2] # S->I
                hazs[2] = pars[2]*x[2]      # I->R
                return(hazs)
        }

        smat = matrix(c(-1,1,0,0,-1,1),nrow=3,ncol=2)
        rownames(smat) = c("S", "I", "R")

        initial = c(S = 15000, I = 100, R = 500)

        pars = c(0.000025,1/7)

        SIR_mod_issb <- create_model(smat, h, initial, pars)

        # Set up stemr objects
        compartments <- c("S","I","R")
        rates <- list(rate("beta * I", "S", "I"),
                      rate("mu", "I", "R"))
        state_initializer <- stem_initializer(c(S = 15000, I = 100, R = 500), fixed = TRUE)
        parameters <- c(beta = 0.000025, mu = 1/7, rho = 0.5)
        tcovar <- NULL
        constants <- NULL
        strata <- NULL
        t0 <- 0; tmax <- 52

        dynamics <- stem_dynamics(rates = rates, parameters = parameters,state_initializer = state_initializer, compartments=compartments, strata = strata, tcovar = tcovar, tmax = tmax, messages = TRUE, compile_rates = T, compile_lna = T, compile_ode = T)

        stem_object <- stem(dynamics = dynamics)

        # simulate paths
        sim2_gillespie_trajecs <- simulate_stem(stem_object = stem_object,
                                                nsim = niter,
                                                paths = TRUE,
                                                observations = FALSE,
                                                subject_paths = FALSE,
                                                method = "gillespie",
                                                t0 = t0,
                                                tmax = tmax,
                                                census_times = seq(0, 52, by = 0.1),
                                                paths_as_array = TRUE,
                                                messages = TRUE)$paths

        sim2_lna_norestart <- simulate_stem(stem_object = stem_object,
                                            nsim = niter,
                                            paths = TRUE,
                                            observations = FALSE,
                                            subject_paths = FALSE,
                                            method = "lna",
                                            lna_restart = FALSE,
                                            t0 = t0,
                                            tmax = tmax,
                                            census_times = seq(0, 52, by = 0.1),
                                            paths_as_array = TRUE,
                                            messages = TRUE)$paths

        sim2_lna_restart <- simulate_stem(stem_object = stem_object,
                                          nsim = niter,
                                          paths = TRUE,
                                          observations = FALSE,
                                          subject_paths = FALSE,
                                          method = "lna",
                                          lna_restart = TRUE,
                                          t0 = t0,
                                          tmax = tmax,
                                          census_times = seq(0, 52, by = 0.1),
                                          paths_as_array = TRUE,
                                          messages = TRUE)$paths

        sim2_ode <- stem_ode(stem_object = stem_object,
                             census_times = census_times,
                             init_states = c(S = 15000, I = 100, R = 500))

        sim2_ode_issb <- deterministic(SIR_mod_issb, 52, 0.1)

        path_quantiles <- data.frame(time = rep(census_times, 5),
                                     method = rep(c("gillespie", "lna_norestart", "lna_restart", "ode", "ode_issb"), each = length(census_times)),
                                     mean = c(apply(sim2_gillespie_trajecs[,"I",,drop = FALSE], 1, mean),
                                              apply(sim2_lna_norestart[,"I",,drop = FALSE], 1, mean),
                                              apply(sim2_lna_restart[,"I",,drop = FALSE], 1, mean),
                                              sim2_ode[,"I"],
                                              sim2_ode_issb[,"I"]),
                                     lower = 0,
                                     upper = 0)
        path_quantiles$lower <- path_quantiles$mean - 1.96 * c(apply(sim2_gillespie_trajecs[,"I",,drop = FALSE], 1, sd),
                                                               apply(sim2_lna_norestart[,"I",,drop = FALSE], 1, sd),
                                                               apply(sim2_lna_restart[,"I",,drop = FALSE], 1, sd),
                                                               rep(0, length(census_times)*2)) / sqrt(niter)
        path_quantiles$upper <- path_quantiles$mean + 1.96 * c(apply(sim2_gillespie_trajecs[,"I",,drop = FALSE], 1, sd),
                                                               apply(sim2_lna_norestart[,"I",,drop = FALSE], 1, sd),
                                                               apply(sim2_lna_restart[,"I",,drop = FALSE], 1, sd),
                                                               rep(0, length(census_times)*2)) / sqrt(niter)

        pdf(paste0("test_simulate_lna",sim_num,"popsize_15600", ".pdf"))

        ggplot(path_quantiles, aes(x = time, y = mean, colour = method)) + geom_line(size = 1) + geom_ribbon(data = path_quantiles, aes(x = time, ymin = lower, ymax = upper, fill = method), alpha = 0.1, size = 0) + labs(y = "Average compartment counts", title = "S0 = 15000, I0 = 100, R0 = 500; beta = 0.000025, mu = 1/7") + theme_bw()

        dev.off()

} else if(sim_num == 3){

        # set up the model in the issb package
        h = function(x, pars) {
                hazs = numeric(length(pars))
                hazs[1] = pars[1]*x[1]*x[2] # S->I
                hazs[2] = pars[2]*x[2]      # I->R
                return(hazs)
        }

        smat = matrix(c(-1,1,0,1,0,-1,1,0),nrow=4,ncol=2)
        rownames(smat) = c("S", "I", "R", "I_INCIDENCE")

        initial = c(S = 1500, I = 10, R = 50, 10)

        pars = c(0.00025,1/7)

        SIR_mod_issb <- create_model(smat, h, initial, pars)

        # Set up stemr objects
        compartments <- c("S","I","R")
        rates <- list(rate("beta * I", "S", "I", incidence = TRUE),
                      rate("mu", "I", "R"))
        state_initializer <- stem_initializer(c(S = 1500, I = 10, R = 50), fixed = TRUE)
        parameters <- c(beta = 0.00025, mu = 1/7, rho = 0.5)
        tcovar <- NULL
        constants <- NULL
        strata <- NULL
        t0 <- 0; tmax <- 52

        dynamics <- stem_dynamics(rates = rates, parameters = parameters,state_initializer = state_initializer, compartments=compartments, strata = strata, tcovar = tcovar, tmax = tmax, messages = TRUE, compile_rates = T, compile_lna = T, compile_ode = T)

        stem_object <- stem(dynamics = dynamics)

        # simulate paths
        sim3_gillespie_trajecs <- simulate_stem(stem_object = stem_object,
                                                nsim = niter,
                                                paths = TRUE,
                                                observations = FALSE,
                                                subject_paths = FALSE,
                                                method = "gillespie",
                                                t0 = t0,
                                                tmax = tmax,
                                                timestep = 0.1,
                                                census_times = seq(0, 52, by = 0.1),
                                                paths_as_array = TRUE,
                                                messages = TRUE)$paths

        sim3_lna_norestart <- simulate_stem(stem_object = stem_object,
                                            nsim = 1,
                                            paths = TRUE,
                                            observations = FALSE,
                                            subject_paths = FALSE,
                                            method = "lna",
                                            lna_restart = FALSE,
                                            t0 = t0,
                                            tmax = tmax,
                                            census_times = seq(0, 52, by = 0.1),
                                            paths_as_array = TRUE,
                                            messages = TRUE)$paths

        sim3_lna_restart <- simulate_stem(stem_object = stem_object,
                                          nsim = niter,
                                          paths = TRUE,
                                          observations = FALSE,
                                          subject_paths = FALSE,
                                          method = "lna",
                                          lna_restart = TRUE,
                                          t0 = t0,
                                          tmax = tmax,
                                          census_times = seq(0, 52, by = 0.1),
                                          paths_as_array = TRUE,
                                          messages = TRUE)$paths

        sim3_ode <- stem_ode(stem_object = stem_object,
                             census_times = census_times,
                             init_states = c(S = 1500, I = 10, R = 50))

        sim3_ode_issb <- deterministic(SIR_mod_issb, 52, 0.1)

        path_quantiles <- data.frame(time = rep(census_times, 5),
                                     method = rep(c("gillespie", "lna_norestart", "lna_restart", "ode", "ode_issb"), each = length(census_times)),
                                     mean = c(apply(sim3_gillespie_trajecs[,"I",,drop = FALSE], 1, mean),
                                              apply(sim3_lna_norestart[,"I",,drop = FALSE], 1, mean),
                                              apply(sim3_lna_restart[,"I",,drop = FALSE], 1, mean),
                                              sim3_ode[,"I"],
                                              sim3_ode_issb[,"I"]),
                                     lower = 0,
                                     upper = 0)
        path_quantiles$lower <- path_quantiles$mean - 1.96 * c(apply(sim3_gillespie_trajecs[,"I",,drop = FALSE], 1, sd),
                                                               apply(sim3_lna_norestart[,"I",,drop = FALSE], 1, sd),
                                                               apply(sim3_lna_restart[,"I",,drop = FALSE], 1, sd),
                                                               rep(0, length(census_times)*2)) / sqrt(niter)
        path_quantiles$upper <- path_quantiles$mean + 1.96 * c(apply(sim3_gillespie_trajecs[,"I",,drop = FALSE], 1, sd),
                                                               apply(sim3_lna_norestart[,"I",,drop = FALSE], 1, sd),
                                                               apply(sim3_lna_restart[,"I",,drop = FALSE], 1, sd),
                                                               rep(0, length(census_times)*2)) / sqrt(niter)

        pdf(paste0("test_simulate_lna_3_popsize_1560", ".pdf"))

        ggplot(path_quantiles, aes(x = time, y = mean, colour = method)) + geom_line(size = 1) + geom_ribbon(data = path_quantiles, aes(x = time, ymin = lower, ymax = upper, fill = method), alpha = 0.1, size = 0) + labs(y = "Average compartment counts", title = "S0 = 1500, I0 = 10, R0 = 50; beta = 0.00025, mu = 1/7") + theme_bw()

        dev.off()

} else if(sim_num == 4){

        h = function(x, pars) {
                hazs = numeric(length(pars))
                hazs[1] = pars[1]*x[1]*x[2] # S->I
                hazs[2] = pars[2]*x[2]      # I->R
                return(hazs)
        }

        smat = matrix(c(-1,1,0,1,0,-1,1,0),nrow=4,ncol=2)
        rownames(smat) = c("S", "I", "R", "I_INCIDENCE")

        initial = c(S = 15000, I = 100, R = 500, 100)

        pars = c(0.000025,1/7)

        SIR_mod_issb <- create_model(smat, h, initial, pars)

        # Set up stemr objects
        compartments <- c("S","I","R")
        rates <- list(rate("beta * I", "S", "I", incidence = TRUE),
                      rate("mu", "I", "R"))
        state_initializer <- stem_initializer(c(S = 15000, I = 100, R = 500), fixed = TRUE)
        parameters <- c(beta = 0.000025, mu = 1/7, rho = 0.5)
        tcovar <- NULL
        constants <- NULL
        strata <- NULL
        t0 <- 0; tmax <- 52

        dynamics <- stem_dynamics(rates = rates, parameters = parameters,state_initializer = state_initializer, compartments=compartments, strata = strata, tcovar = tcovar, tmax = tmax, messages = TRUE, compile_rates = T, compile_lna = T, compile_ode = T)

        stem_object <- stem(dynamics = dynamics)

        # simulate paths
        sim4_gillespie_trajecs <- simulate_stem(stem_object = stem_object,
                                                nsim = niter,
                                                paths = TRUE,
                                                observations = FALSE,
                                                subject_paths = FALSE,
                                                method = "gillespie",
                                                t0 = t0,
                                                tmax = tmax,
                                                census_times = seq(0, 52, by = 0.1),
                                                paths_as_array = TRUE,
                                                messages = TRUE)$paths

        sim4_lna_norestart <- simulate_stem(stem_object = stem_object,
                                            nsim = niter,
                                            paths = TRUE,
                                            observations = FALSE,
                                            subject_paths = FALSE,
                                            method = "lna",
                                            lna_restart = FALSE,
                                            t0 = t0,
                                            tmax = tmax,
                                            census_times = seq(0, 52, by = 0.1),
                                            paths_as_array = TRUE,
                                            messages = TRUE)$paths

        sim4_lna_restart <- simulate_stem(stem_object = stem_object,
                                          nsim = niter,
                                          paths = TRUE,
                                          observations = FALSE,
                                          subject_paths = FALSE,
                                          method = "lna",
                                          lna_restart = TRUE,
                                          t0 = t0,
                                          tmax = tmax,
                                          census_times = seq(0, 52, by = 0.1),
                                          paths_as_array = TRUE,
                                          messages = TRUE)$paths

        sim4_ode <- stem_ode(stem_object = stem_object,
                             census_times = census_times,
                             init_states = c(S = 15000, I = 100, R = 500))

        sim4_ode_issb <- deterministic(SIR_mod_issb, 52, 0.1)

        path_quantiles <- data.frame(time = rep(census_times, 5),
                                     method = rep(c("gillespie", "lna_norestart", "lna_restart", "ode", "ode_issb"), each = length(census_times)),
                                     mean = c(apply(sim4_gillespie_trajecs[,"I",,drop = FALSE], 1, mean),
                                              apply(sim4_lna_norestart[,"I",,drop = FALSE], 1, mean),
                                              apply(sim4_lna_restart[,"I",,drop = FALSE], 1, mean),
                                              sim4_ode[,"I"],
                                              sim4_ode_issb[,"I"]),
                                     lower = 0,
                                     upper = 0)
        path_quantiles$lower <- path_quantiles$mean - 1.96 * c(apply(sim4_gillespie_trajecs[,"I",,drop = FALSE], 1, sd),
                                                               apply(sim4_lna_norestart[,"I",,drop = FALSE], 1, sd),
                                                               apply(sim4_lna_restart[,"I",,drop = FALSE], 1, sd),
                                                               rep(0, length(census_times)*2)) / sqrt(niter)
        path_quantiles$upper <- path_quantiles$mean + 1.96 * c(apply(sim4_gillespie_trajecs[,"I",,drop = FALSE], 1, sd),
                                                               apply(sim4_lna_norestart[,"I",,drop = FALSE], 1, sd),
                                                               apply(sim4_lna_restart[,"I",,drop = FALSE], 1, sd),
                                                               rep(0, length(census_times)*2)) / sqrt(niter)

        pdf(paste0("test_simulate_lna_4_popsize_15600", ".pdf"))

        ggplot(path_quantiles, aes(x = time, y = mean, colour = method)) + geom_line(size = 1) + geom_ribbon(data = path_quantiles, aes(x = time, ymin = lower, ymax = upper, fill = method), alpha = 0.1, size = 0) + labs(y = "Average compartment counts", title = "S0 = 15000, I0 = 100, R0 = 500; beta = 0.000025, mu = 1/7") + theme_bw()

        dev.off()

} else if(sim_num == 3) {

        s_params <- c(0.75, 0.575e-5)
        initial_state_vec <- c(S = 1500, I = 10, R = 50)
        compartments <- c("S","I","R")
        rates <- list(rate("beta * I", "S", "I", seasonality = seasonality(period = 52, intercept = 5e-5, trend = 1e-4, s_params = s_params), incidence = T),
                      rate("mu", "I", "R"),
                      rate("phi", "R", "S"))
        state_initializer <- stem_initializer(initial_state_vec, fixed = TRUE)
        parameters <- c(beta = 0.00025, mu = 1/2, phi = 3)
        tcovar <- NULL
        constants <- NULL
        strata <- NULL
        t0 <- 0; tmax <- 52

        dynamics <- stem_dynamics(rates = rates, parameters = parameters,state_initializer = state_initializer, compartments=compartments, strata = strata, tcovar = tcovar, t0 = 0, tmax = 52, messages = FALSE)

        stem_object <- stem(dynamics = dynamics)

        # simulate paths
        sim2_gillespie_trajecs <- simulate_stem(stem_object = stem_object,
                                       nsim = niter,
                                       paths = TRUE,
                                       observations = FALSE,
                                       subject_paths = FALSE,
                                       method = "gillespie",
                                       t0 = t0,
                                       tmax = tmax,
                                       census_times = census_times,
                                       paths_as_array = FALSE,
                                       messages = FALSE)$paths

        sim2_lna_norestart <- simulate_stem(stem_object = stem_object,
                                           nsim = niter,
                                           paths = TRUE,
                                           observations = FALSE,
                                           subject_paths = FALSE,
                                           method = "lna",
                                           lna_restart = FALSE,
                                           t0 = t0,
                                           tmax = tmax,
                                           census_times = census_times,
                                           paths_as_array = FALSE,
                                           messages = FALSE)$paths


} else if(sim_num == 3) {

        s_params <- c(0.5, 0.575e-4)
        strata <- c("male", "female")
        compartments <- list(S = "ALL", I = "ALL", R = "ALL")
        rates <- list(rate("beta * I_SELF + gamma * comp_fcn(I_ADJ, sum)", "S", "I", "ALL",
                           seasonality(period = 52, common_seasonality = TRUE, intercept = 5e-5, trend = 1e-4, s_params = s_params)),
                      rate("mu", "I", "R", "ALL"),
                      rate("phi", "R", "S", "ALL"))
        adjacency <- matrix(c(0,1,1,0), nrow = 2, ncol = 2); colnames(adjacency) <- rownames(adjacency) <- strata
        initial_state_vec <- c(S = 900, I = 10, R = 90)
        state_initializer <- list(stem_initializer(initial_state_vec, fixed = TRUE, strata = c("male", "female"), shared_params = TRUE))
        parameters <- c(beta = 0.000025, gamma = 0.0005, mu = 1/2, phi = 5)
        tcovar <- NULL
        constants <- NULL
        t0 <- 0; tmax <- 52

        dynamics <- stem_dynamics(rates = rates, parameters = parameters,state_initializer = state_initializer, compartments=compartments, strata = strata, tcovar = tcovar, t0 = 0, tmax = 52, adjacency = adjacency, messages = T)

        stem_object <- stem(dynamics = dynamics)

        # Set up GillespieSSA objects
        gillespie_tcovar <- data.frame(time = seq(0, 52, by = 52/50),
                                       timevar = 5e-5  + 1e-4 * seq(0, 52, by = 52/50),
                                       SIN52 = s_params[1]*sin(2*pi*seq(0, 52, by = 52/50)/52),
                                       COS52 = s_params[2]*cos(2*pi*seq(0, 52, by = 52/50)/52))

        a = c("(beta * I_male + gamma * I_female + exp(timevar + SIN52 + COS52)) * S_male",
              "(beta * I_female + gamma * I_male + exp(timevar + SIN52 + COS52)) * S_female",
              "mu * I_male", "mu * I_female",
              "phi * R_male", "phi * R_female")
        nu <- matrix(c(-1, 0, 1, 0, 0, 0,
                       0, -1, 0, 1, 0, 0,
                       0, 0, -1, 0, 1, 0,
                       0, 0, 0, -1, 0, 1,
                       1, 0, 0, 0, -1, 0,
                       0, 1, 0, 0, 0, -1), nrow = 6)

        # simulate paths
        stemr_trajecs <- simulate_stem(stem_object = stem_object,
                                       nsim = niter,
                                       paths = TRUE,
                                       observations = FALSE,
                                       subject_paths = FALSE,
                                       method = "gillespie",
                                       t0 = t0,
                                       tmax = tmax,
                                       census_times = census_times,
                                       paths_as_array = FALSE,
                                       messages = FALSE)$paths

        # nclust = 15
        # cl <- makeForkCluster(nclust)
        cl <- makePSOCKcluster(2)
        registerDoParallel(cl)

        gSSA_trajecs <- foreach(j = 1:niter, .packages = c("stemr", "GillespieSSA")) %dopar% {

                gillespie_initstate <- c(S_male = 900, S_female = 900, I_male = 10, I_female = 10, R_male = 90, R_female = 90)
                gSSA_sim <- matrix(0, nrow = 1e6, ncol = 7)
                gSSA_sim[1,] <- c(0, gillespie_initstate)
                colnames(gSSA_sim) <- c("time", "S_male", "S_female", "I_male", "I_female", "R_male", "R_female")

                subsim_ind <- 2

                for(k in 1:(nrow(gillespie_tcovar) - 1)) {
                        t0   <- gillespie_tcovar[k, 1]
                        tmax <- gillespie_tcovar[k+1, 1]

                        params <- c(beta = 0.000025,
                                    gamma = 0.0005,
                                    mu = 1/2,
                                    phi = 5,
                                    timevar = gillespie_tcovar[k, 2],
                                    SIN52 = gillespie_tcovar[k, 3],
                                    COS52 = gillespie_tcovar[k, 4])

                        initstate <- gSSA_sim[subsim_ind-1,-1]

                        subsim <- ssa(initstate, a, nu, params, tf = tmax-t0)$data
                        subsim <- subsim[subsim[,1] < (tmax-t0),, drop = FALSE]
                        subsim[,1] <- subsim[,1] + t0

                        gSSA_sim[subsim_ind + (seq_len(nrow(subsim)-1)-1),] <- subsim[-1, , drop = FALSE]

                        subsim_ind <- subsim_ind + nrow(subsim) - 1
                }

                gSSA_sim <- gSSA_sim[1:(subsim_ind-1),]

                gSSA_traj <- build_census_path(gSSA_sim, census_times, 1:6)

                return(gSSA_traj)
        }
}