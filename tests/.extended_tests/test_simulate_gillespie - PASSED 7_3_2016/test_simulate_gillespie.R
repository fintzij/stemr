# Extended test for Gillespie forward simulation functions. Tested against the 'GillespieSSA' package.

# Simulation test to check that simulate_epimodel produces the same results as
# the GillespieSSA package. Two models are simulated: simple SIR (simulation 1),
# and SIRS with seasonality (simulation 2). Each of the simulations will produce
# a plot comparing the trajectories produced by simulate epimodel with the
# trajectories produced by GillespieSSA.

library(stemr)
library(GillespieSSA)
library(doParallel)
library(doMC)
library(itertools)
library(foreach)
library(iterators)
library(ggplot2)
library(batch)

parseCommandArgs()

# Set up simulation objects ----------------------------------------------------------------

niter <- 100000; census_times <- seq(0,52,by=0.1)

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
                                       as_array = FALSE,
                                       messages = FALSE)

        # nclust = 15
        # cl <- makeForkCluster(nclust)
        cl <- makePSOCKcluster(2)
        registerDoParallel(cl)

        gSSA_trajecs <- foreach(j = 1:niter, .packages = c("stemr", "GillespieSSA")) %dopar% {

                gSSA_sim <- ssa(gillespie_initstate, a, nu, parameters, tf = tmax)$data

                gSSA_traj <- get_census_path(gSSA_sim, census_times, 1:3)

                return(gSSA_traj)
        }

} else if(sim_num == 2) {

        s_params <- c(0.75e-4, 0.575e-5)
        # s_params <- c(0,0)
        initial_state_vec <- c(S = 150, I = 10, R = 5)
        compartments <- c("S","I","R")
        rates <- list(rate("beta * I + exp(delta + iota * TIME)", "S", "I", seasonality = seasonality(period = 52, s_params)),
                      rate("mu", "I", "R"),
                      rate("phi", "R", "S"))
        state_initializer <- stem_initializer(initial_state_vec, fixed = TRUE)
        parameters <- c(beta = 0.0025, delta = 0.001, iota = 0.002, mu = 1/7, phi = 0.075)
        tcovar <- NULL
        constants <- NULL
        strata <- NULL
        t0 <- 0; tmax <- 52

        dynamics <- stem_dynamics(rates = rates, parameters = parameters,state_initializer = state_initializer, compartments=compartments, strata = strata, tcovar = tcovar, t0 = 0, tmax = 52, messages = FALSE)

        stem_object <- stem(dynamics = dynamics)


        # Set up GillespieSSA objects
        gillespie_tcovar <- data.frame(time = seq(0, 52, by = 52/50),
                                       timevar = exp(0.001 + seq(0, 52, by = 52/50) * parameters["iota"]),
                                       SIN52 = s_params[1]*sin(2*pi*seq(0, 52, by = 52/50)/52),
                                       COS52 = s_params[2]*cos(2*pi*seq(0, 52, by = 52/50)/52))

        a = c("(beta*I + timevar + exp(SIN52 + COS52))*S", "mu*I", "phi * R")
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
                                       as_array = FALSE,
                                       messages = FALSE)

        # nclust = 15
        # cl <- makeForkCluster(nclust)
        cl <- makePSOCKcluster(2)
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

                        params <- c(beta = 0.0025,
                                    mu = 1/7,
                                    phi = 0.075,
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

                gSSA_traj <- get_census_path(gSSA_sim, census_times, 1:3)

                return(gSSA_traj)
        }
}


#  Generate plots of the results ------------------------------------------

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

