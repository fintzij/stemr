# Extended test for Gillespie forward simulation functions. Tested against the 'GillespieSSA' package.

# Simulation test to check that simulate_stem produces the same results as the
# GillespieSSA package. Three models are simulated: simple SIR (simulation 1),
# SIRS with seasonality (simulation 2), SIRS with seasonality and two strata.
# Each of the simulations will produce a plot comparing the trajectories
# produced by simulate epimodel with the trajectories produced by GillespieSSA.

library(stemr)
library(GillespieSSA)
library(doParallel)
library(doMC)
library(itertools)
library(foreach)
library(iterators)
library(ggplot2)

args <- commandArgs(TRUE)
print(args)
sim_num <- as.numeric(args[1])

# Set up simulation objects ----------------------------------------------------------------

niter <- 100; census_times <- seq(0,52,by=1)

# test gillespie simulation for large population
compartments <- c("S","I","R")
rates <- list(rate("beta * I", "S", "I", incidence = TRUE),
              rate("mu", "I", "R"))
init_state = c(S = 500, I = 10, R = 50)
state_initializer <- stem_initializer(c(S = 50000, I = 10, R = 50), fixed = T)
parameters <- c(beta = 0.001, mu = 1/7, rho = 0.5)
tcovar <- NULL
constants <- NULL
strata <- NULL
t0 <- 0; tmax <- 52
timestep <- NULL
adjacency <- NULL
messages <- T
census_times = 0:tmax

# compile dynamics
dynamics <- stem_dynamics(rates = rates, parameters = parameters, tmax = tmax, state_initializer = state_initializer, compartments=compartments, strata = strata, tcovar = tcovar, messages = TRUE, compile_ode = F, compile_rates = T)

stem_object <- stem(dynamics = dynamics)

stemr_trajecs <-
        simulate_stem(
                stem_object = stem_object,
                method = "gillespie",
                paths = TRUE,
                observations = F,
                tmax = tmax,
                nsim = niter,
                census_times = 0:tmax
        )$paths

# Set up GillespieSSA objects
gillespie_initstate <- c(infections = 1, recoveries = 0)
gillespie_params <- c(beta = 0.001, mu = 1/7, initsusc = 500, initinfecs = 10)
a = c("beta*(initsusc - infections)*(initinfecs + infections - recoveries)", "mu*(initinfecs + infections - recoveries)")
nu <- diag(1,2)

gSSA_sim <- ssa(gillespie_initstate, a, nu, gillespie_params, tf = tmax, censusInterval = 1)$data
gSSA_traj <- build_census_path(gSSA_sim, census_times, 1:2)


if(sim_num == 1){

        # Set up stemr objects
        compartments <- c("S","I","R")
        rates <- list(rate("beta * I", "S", "I"),
                      rate("mu", "I", "R"))
        state_initializer <- stem_initializer(c(S = 150, I = 10, R = 5), fixed = TRUE)
        parameters <- c(beta = 0.0025, mu = 1/7, rho = 0.5)
        tcovar <- NULL
        constants <- NULL
        strata <- NULL
        t0 <- 0; tmax <- 52

        dynamics <- stem_dynamics(rates = rates, parameters = parameters,state_initializer = state_initializer, compartments=compartments, strata = strata, tcovar = tcovar, messages = FALSE)

        stem_object <- stem(dynamics = dynamics)

        # Set up GillespieSSA objects
        gillespie_initstate <- c(S = 150, I = 10, R = 5)
        a = c("beta*S*I", "mu*I")
        nu <- matrix(c(-1,+1,0, 0, -1, +1),nrow=3)

        # simulate paths
        stemr_trajecs <- simulate_stem(stem_object = stem_object,
                                       nsim = niter,
                                       paths = TRUE,
                                       observations = FALSE,
                                       subject_paths = FALSE,
                                       method = "gillespie",
                                       t0 = t0,
                                       tmax = tmax,
                                       census_times = seq(0, 52, by = 0.1),
                                       paths_as_array = FALSE,
                                       messages = FALSE)$paths

        nclust = 15
        cl <- makeForkCluster(nclust)
        # cl <- makePSOCKcluster(2)
        registerDoParallel(cl)

        gSSA_trajecs <- foreach(j = 1:niter, .packages = c("stemr", "GillespieSSA")) %dopar% {

                gSSA_sim <- ssa(gillespie_initstate, a, nu, parameters, tf = tmax)$data

                gSSA_traj <- build_census_path(gSSA_sim, census_times, 1:3)

                return(gSSA_traj)
        }

} else if(sim_num == 2) {

        s_params <- c(0.75, 0.575e-5)
        initial_state_vec <- c(S = 150, I = 10, R = 5)
        compartments <- c("S","I","R")
        rates <- list(rate("beta * I", "S", "I", seasonality = seasonality(period = 52, intercept = 5e-5, trend = 1e-4, s_params = s_params)),
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


        # Set up GillespieSSA objects
        gillespie_tcovar <- data.frame(time = seq(0, 52, by = 52/50),
                                       timevar = 5e-5 + 1e-4 * seq(0, 52, by = 52/50),
                                       SIN52 = s_params[1]*sin(2*pi*seq(0, 52, by = 52/50)/52),
                                       COS52 = s_params[2]*cos(2*pi*seq(0, 52, by = 52/50)/52))

        a = c("(beta * I + exp(timevar + SIN52 + COS52)) * S", "mu * I", "phi * R")
        nu <- matrix(c(-1,+1,0, 0, -1, +1, 1, 0, -1),nrow=3)

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

        nclust = 15
        cl <- makeForkCluster(nclust)
        # cl <- makePSOCKcluster(2)
        registerDoParallel(cl)

        gSSA_trajecs <- foreach(j = 1:niter, .packages = c("stemr", "GillespieSSA")) %dopar% {

                gillespie_initstate <- initial_state_vec
                gSSA_sim <- matrix(0, nrow = 1e5, ncol = 4)
                gSSA_sim[1,] <- c(0, gillespie_initstate)
                colnames(gSSA_sim) <- c("time", "S", "I", "R")

                subsim_ind <- 2

                for(k in 1:(nrow(gillespie_tcovar) - 1)) {
                        t0   <- gillespie_tcovar[k, 1]
                        tmax <- gillespie_tcovar[k+1, 1]

                        params <- c(beta = 0.00025,
                                    mu = 1/2,
                                    phi = 3,
                                    timevar = gillespie_tcovar[k, 2],
                                    SIN52 = gillespie_tcovar[k, 3],
                                    COS52 = gillespie_tcovar[k, 4])

                        print(params["timevar"])

                        initstate <- gSSA_sim[subsim_ind-1,-1]

                        subsim <- ssa(initstate, a, nu, params, tf = tmax-t0)$data
                        subsim <- subsim[subsim[,1] < (tmax-t0),, drop = FALSE]
                        subsim[,1] <- subsim[,1] + t0

                        gSSA_sim[subsim_ind + (seq_len(nrow(subsim)-1)-1),] <- subsim[-1, , drop = FALSE]

                        subsim_ind <- subsim_ind + nrow(subsim) - 1
                }

                gSSA_sim <- gSSA_sim[1:(subsim_ind-1),]

                gSSA_traj <- build_census_path(gSSA_sim, census_times, 1:3)

                return(gSSA_traj)
        }


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


#  Generate plots of the results ------------------------------------------

if(sim_num != 3) {
        # convert lists of censused paths to arrays
        n_censustimes <- length(census_times)
        stemr_trajecs <- array(unlist(stemr_trajecs), dim = c(n_censustimes, 4, niter))
        gSSA_trajecs  <- array(unlist(gSSA_trajecs), dim = c(n_censustimes, 4, niter))

        colnames(stemr_trajecs) <- colnames(gSSA_trajecs) <- c("time", "S", "I", "R")

        results <- data.frame(method = rep(c("stemr", "GillespieSSA"), each = n_censustimes*3),
                              compartment = rep(c("S","I","R"), each = n_censustimes),
                              time = census_times,
                              estimate = 0,
                              lower = 0,
                              upper = 0)

        for(s in c("S", "I", "R")) {
                results$estimate[results$compartment == s & results$method == "stemr"] <- stemr_ests <- apply(stemr_trajecs[,s,,drop = FALSE], 1, mean)
                results$estimate[results$compartment == s & results$method == "GillespieSSA"] <- gSSA_ests <- apply(gSSA_trajecs[,s,,drop = FALSE], 1, mean)

                results$lower[results$compartment == s & results$method == "stemr"] <-
                        pmax(0, stemr_ests - 1.96 * apply(stemr_trajecs[,s,,drop = FALSE], 1, sd)/sqrt(niter))
                results$lower[results$compartment == s & results$method == "GillespieSSA"] <-
                        pmax(0, gSSA_ests - 1.96 * apply(gSSA_trajecs[,s,,drop = FALSE], 1, sd)/sqrt(niter))

                results$upper[results$compartment == s & results$method == "stemr"] <- stemr_ests + 1.96 * apply(stemr_trajecs[,s,,drop = FALSE], 1, sd)/sqrt(niter)
                results$upper[results$compartment == s & results$method == "GillespieSSA"] <- gSSA_ests + 1.96 * apply(gSSA_trajecs[,s,,drop = FALSE], 1, sd)/sqrt(niter)
        }

        pdf(paste("test_simulate_gillespie",sim_num,".pdf", sep = ""))

        ggplot(results, aes(x = time, y = estimate, colour = method, linetype = compartment)) + geom_line(size = 1) + geom_ribbon(data = results, aes(x = time, ymin = lower, ymax = upper, fill = method), alpha = 0.1, size = 0) + labs(y = "Average compartment counts") + theme_bw()

        dev.off()

} else {
        # convert lists of censused paths to arrays
        n_censustimes <- length(census_times)
        stemr_trajecs <- array(unlist(stemr_trajecs), dim = c(n_censustimes, 7, niter))
        gSSA_trajecs  <- array(unlist(gSSA_trajecs), dim = c(n_censustimes, 7, niter))

        colnames(stemr_trajecs) <- colnames(gSSA_trajecs) <- c("time", "S_male", "S_female","I_male", "I_female","R_male", "R_female")

        results <- data.frame(method = rep(c("stemr", "GillespieSSA"), each = n_censustimes*6),
                              compartment = rep(c("S_male", "S_female","I_male", "I_female","R_male", "R_female"), each = n_censustimes),
                              time = census_times,
                              estimate = 0,
                              lower = 0,
                              upper = 0)

        for(s in c("S_male", "S_female","I_male", "I_female","R_male", "R_female")) {
                results$estimate[results$compartment == s & results$method == "stemr"] <- stemr_ests <- apply(stemr_trajecs[,s,,drop = FALSE], 1, mean)
                results$estimate[results$compartment == s & results$method == "GillespieSSA"] <- gSSA_ests <- apply(gSSA_trajecs[,s,,drop = FALSE], 1, mean)

                results$lower[results$compartment == s & results$method == "stemr"] <-
                        pmax(0, stemr_ests - 1.96 * apply(stemr_trajecs[,s,,drop = FALSE], 1, sd)/sqrt(niter))
                results$lower[results$compartment == s & results$method == "GillespieSSA"] <-
                        pmax(0, gSSA_ests - 1.96 * apply(gSSA_trajecs[,s,,drop = FALSE], 1, sd)/sqrt(niter))

                results$upper[results$compartment == s & results$method == "stemr"] <- stemr_ests + 1.96 * apply(stemr_trajecs[,s,,drop = FALSE], 1, sd)/sqrt(niter)
                results$upper[results$compartment == s & results$method == "GillespieSSA"] <- gSSA_ests + 1.96 * apply(gSSA_trajecs[,s,,drop = FALSE], 1, sd)/sqrt(niter)
        }


        pdf(paste("test_simulate_gillespie",sim_num,".pdf", sep = ""))

        ggplot(results, aes(x = time, y = estimate, colour = method, linetype = compartment)) + geom_line(size = 1) + geom_ribbon(data = results, aes(x = time, ymin = lower, ymax = upper, fill = method), alpha = 0.1, size = 0) + labs(y = "Average compartment counts") + theme_bw()

         dev.off()
}

