library(stemr)
library(extraDistr)
library(foreach)
library(doRNG)
library(doParallel)

# no strata the stemr object --------------------------------------------------
set.seed(12511)
strata <- NULL
compartments <- c("S", "I", "R")
rates <- list(rate("beta_t * I", "S", "I", incidence = T),
              rate("mu", "I", "R", incidence = T),
              rate("omega", "R", "S", incidence = T))
state_initializer <-
   list(stem_initializer(c(S = 1e5-15, I = 10, R = 5),
                         fixed = T, 
                         prior = c(1e5-15, 10, 5)))
adjacency <- NULL
tcovar <- cbind(time = 0:30,
                beta_t = (2.1 + 0.35 * sin(seq(0,30) / 3.5) + rnorm(31, 0, 0.1)) / 1e5)
parameters = c(mu = 1, omega = 1/4, rho = 0.25, phi = 50)
constants <- c(t0 = 0)
t0 <- 0; tmax <- 30;

dynamics <-
      stem_dynamics(
            rates = rates,
            tmax = tmax,
            parameters = parameters,
            state_initializer = state_initializer,
            compartments = compartments,
            constants = constants,
            strata = strata,
            adjacency = adjacency,
            tcovar = tcovar,
            messages = T,
            compile_ode = T,
            compile_rates = T,
            compile_lna = T,
            rtol = 1e-6,
            atol = 1e-6,
            step_size = 1e-6
      )

emissions <-
   list(emission("S2I", "negbinomial", c("phi", "S2I * rho"), 
                 incidence = TRUE,
                 obstimes = seq(1, tmax, by = 1)))

measurement_process <- stem_measure(emissions = emissions, dynamics = dynamics, messages = T)

stem_object <- make_stem(dynamics = dynamics, measurement_process = measurement_process)

stem_data <- simulate_stem(stem_object  = stem_object,
                           method       = "gillespie",
                           paths        = TRUE,
                           observations = T,
                           nsim         = 1,
                           census_times = unique(c(0:tmax)))

# grab the dataset
true_path <- stem_data$paths[[1]]
dat <- stem_data$datasets[[1]]

# Set up model to be fit --------------------------------------------------

set.seed(12511)
strata <- NULL
compartments <- c("S", "I", "R")
rates <- list(rate("beta_t * I", "S", "I", incidence = T),
              rate("mu", "I", "R", incidence = T),
              rate("omega", "R", "S", incidence = TRUE))
state_initializer <-
   list(stem_initializer(c(S = 1e5 - 15, I = 10, R = 5),
                         fixed = F,
                         prior = c(1e5 - 15, 10, 5)))
adjacency <- NULL
tcovar <- NULL
parameters = c(
      log_R0_init = log(2),
      mu = 1, 
      omega = 1/4,
      rho = 0.25,
      phi = 50,
      sigma = 0.05)

constants <- c(t0 = 0)
popsize = 1e5; log_popsize = log(popsize)
t0 <- 0; tmax <- 30;

# RW in terms of log(R0), parameterized by differences
foi_rw1 <- function(parameters, draws, log_pop = log_popsize) {
      
      log_R0_t <- numeric(length = 1 + length(draws))
      log_R0_t[1] <- parameters["log_R0_init"] 
      
      for(t in 2:(length(log_R0_t))) {
            log_R0_t[t] <- log_R0_t[t-1] + draws[t-1] * parameters["sigma"]
      }
      
      return(exp(log_R0_t - log_pop + log(parameters["mu"])))
}

tparam <- 
   list(tpar(tparam_name = "beta_t",
             times = 0:(tmax-1), 
             n_draws = tmax,
             draws2par = foi_rw1))

dynamics <-
      stem_dynamics(
            rates = rates,
            tmax = tmax,
            timestep = NULL,
            parameters = parameters,
            state_initializer = state_initializer,
            compartments = compartments,
            constants = constants,
            strata = strata,
            adjacency = adjacency,
            tcovar = NULL,
            tparam = tparam,
            messages = T,
            compile_ode = T,
            compile_rates = T,
            compile_lna = T,
            rtol = 1e-6,
            atol = 1e-6,
            step_size = 1e-6
      )

emissions <- list(emission("S2I", "negbinomial", c("phi","S2I * rho"), incidence = TRUE, obstimes = seq(1, tmax, by = 1)))

measurement_process <- stem_measure(emissions = emissions, dynamics = dynamics, data = dat, messages = T)

stem_object <- make_stem(dynamics = dynamics, measurement_process = measurement_process)

# initialize the inference
to_estimation_scale <- function(params_nat) {
      return(c(
            params_nat[1],              # log_R0_0      -> log_R0_0
            -log(params_nat[2]),        # mu            -> log(mu)
            -log(params_nat[3]),        # omega         -> log(omega)
            qlogis(params_nat[4]),      # rho           -> logit(rho)
            - 0.5 * log(params_nat[5]), # phi           -> log(1/sqrt(phi))
            log(params_nat[6])
      ))
}

from_estimation_scale <- function(params_est) {
      return(c(
            params_est[1],                         # log_R0_0       -> log_R0_0
            exp(-params_est[2]),                   # log(mu)        -> mu
            exp(-params_est[3]),                   # log(omega)     -> omega
            plogis(params_est[4]),                 # logit(rho)     -> rho
            exp(-2 * params_est[5]),               # log(1/sqrt(phi)) -> phi
            exp(params_est[6])
      )) 
}

# Priors
# log_R0_0                  ~ norm(log(2), 0.5): R0_0 = 0.88, 1.05, 1.43, 2, 2.8, 3.8, 4.55
# log(1/mu)                 ~ norm(0.05, 0.25): mean infectious period = 0.76, 0.89, 1.05, 1.24, 1.45
# log(1/omega)              ~ norm(2, 0.25): mean immune period = 5.4, 6.2, 7.4, 8.7, 10.2
# logit(rho_0)              ~ norm(mean = logit(0.225), sd = 0.5): rho = 0.13, 0.17, 0.225, 0.29, 0.36
# 1/sqrt(phi)               ~ exp(3.5): phi = 1.1e3, 150, 25.5, 6.37, 2.3 
# sigma                     ~ exp(20): sigma = 0.0026, 0.005, 0.014, 0.035, 0.069, 0.115, 0.15

logprior = function(params_est) {
   dnorm(params_est[1], log(2), 0.5, log = T) + 
      dnorm(params_est[2], 0.05, 0.5, log = TRUE) +                    # log(1/mu) 
      dnorm(params_est[3], 1.35, 0.5, log = TRUE) +                    # log(1/omega)
      dnorm(params_est[4], mean = logit(0.3), sd = 0.5, log = TRUE) +  # logit(rho)
      dexp(exp(params_est[5]), 5, log = TRUE) + params_est[5] + 
      dexp(exp(params_est[6]), 20, log = TRUE) + params_est[6]
}

priors <- 
      list(logprior = logprior,
           to_estimation_scale = to_estimation_scale,
           from_estimation_scale = from_estimation_scale)

par_initializer = function() {
   priors$from_estimation_scale(priors$to_estimation_scale(parameters) + 
                                   rnorm(6, 0, 0.1))
}

# specify the kernel
mcmc_kern <-
   mcmc_kernel(
      parameter_blocks = 
         list(parblock(
            pars_nat = c("log_R0_init", "mu", "omega", "rho", "phi", "sigma"),
            pars_est = c("log_R0_init_est", "log_mu", "log_omega", "logit_rho", "log_phi", "log_sigma"),
            priors = priors,
            # alg = "mvnss",
            alg = "mvnmh",
            sigma = diag(0.01, 6),
            initializer = par_initializer,
            control = 
               # mvnss_control(stop_adaptation = 1e2))),
               mvnmh_control(stop_adaptation = 1e2))),
      lna_ess_control = lna_control(bracket_update_iter = 50,
                                    joint_initdist_update = FALSE),
      tparam_ess_control = tpar_control(bracket_update_iter = 50))

res <-
   fit_stem(stem_object = stem_object,
            method = "lna",
            mcmc_kern = mcmc_kern,
            iterations = 3.5e2)
