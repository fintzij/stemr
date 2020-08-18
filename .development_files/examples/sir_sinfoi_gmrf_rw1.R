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
state_initializer <- list(stem_initializer(c(S = 1e5-15, I = 10, R = 5), fixed = T, prior = c(1e5-15, 10, 5)))
adjacency <- NULL
tcovar <- cbind(time = 0:30,
                beta_t = (2 + 0.35 * sin(seq(0,30) / 3.5) + rnorm(31, 0, 0.1)) / 1e5)
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
            compile_ode = F,
            compile_rates = T,
            compile_lna = T,
            rtol = 1e-6,
            atol = 1e-6,
            step_size = 1e-6
      )

emissions <- list(emission("S2I", "negbinomial", c("phi", "S2I * rho"), incidence = TRUE, obstimes = seq(1, tmax, by =1)))

measurement_process <- stem_measure(emissions = emissions, dynamics = dynamics, messages = T)

stem_object <- stem(dynamics = dynamics, measurement_process = measurement_process)

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
state_initializer <- list(stem_initializer(c(S = 1e5-15, I = 10, R = 5), fixed = F, prior = c(1e5-15, 10, 5)))
adjacency <- NULL
tcovar <- NULL
parameters = c(
      log_R0_0 = log(2),
      mu = 1, 
      omega = 1/4,
      rho = 0.25,
      phi = 50)

constants <- c(t0 = 0)
popsize = 1e5; log_popsize = log(popsize)
t0 <- 0; tmax <- 30;

# RW in terms of log(R0), parameterized by differences
foi_rw1 <- function(parameters, draws) {
      
      sig <- exp(-2.75 + 0.25 * draws[length(draws)])   
      log_R0_t <- numeric(length = length(draws)-1)
      log_R0_t[1] <- parameters[1] 
      
      for(t in 2:(length(draws)-1)) {
            log_R0_t[t] <- log_R0_t[t-1] + draws[t-1] * sig
      }
      
      return(c(exp(log_R0_t - log_popsize + log(parameters[2])), sig))
}

tparam <- list(tpar(tparam_name = "beta_t", times = 0:30, draws2par = foi_rw1))

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

stem_object <- stem(dynamics = dynamics, measurement_process = measurement_process)

# initialize the inference
to_estimation_scale <- function(params_nat) {
      return(c(
            log(expm1(params_nat[1])),  # log_R0_0      -> log_R0_0
            -log(params_nat[2]),        # mu            -> log(mu)
            -log(params_nat[3]),        # omega         -> log(omega)
            logit(params_nat[4]),       # rho           -> logit(rho)
            - 0.5 * log(params_nat[5])  # phi           -> log(1/sqrt(phi))
      ))
}

from_estimation_scale <- function(params_est) {
      return(c(
            log1p(exp(params_est[1])),             # log_R0_0       -> log_R0_0
            exp(-params_est[2]),                   # log(mu)        -> mu
            exp(-params_est[3]),                   # log(omega)     -> omega
            expit(params_est[4]),                  # logit(rho)     -> rho
            exp(-2 * params_est[5])                  # log(1/sqrt(phi)) -> phi
      )) 
}

# Priors
# log_R0_0                  ~ norm(0.6, 0.3): R0_0 = 1.54, 1.79, 2.12, 2.51, 2.92
# log(1/mu)                 ~ norm(0.05, 0.25): mean infectious period = 0.76, 0.89, 1.05, 1.24, 1.45
# log(1/omega)              ~ norm(2, 0.25): mean immune period = 5.4, 6.2, 7.4, 8.7, 10.2
# logit(rho_0)              ~ norm(mean = logit(0.225), sd = 0.5): rho = 0.13, 0.17, 0.225, 0.29, 0.36
# 1/sqrt(phi)               ~ exp(3.5): phi = 1.1e3, 150, 25.5, 6.37, 2.3 

priors <- 
      list(
            prior_density =
                  function(params_nat, params_est) {
                        dnorm(params_est[1], 0, 0.5, log = TRUE) +                             # log_R0_0
                              dnorm(params_est[2], 0.05, 0.5, log = TRUE) +                    # log(1/mu) 
                              dnorm(params_est[3], 1.35, 0.5, log = TRUE) +                    # log(1/omega)
                              dnorm(params_est[4], mean = logit(0.3), sd = 0.5, log = TRUE) +  # logit(rho)
                              dexp(exp(params_est[5]), 5, log = TRUE) + params_est[5]
                  },
            to_estimation_scale = to_estimation_scale,
            from_estimation_scale = from_estimation_scale)

covmat_names <- c(
      "log_R0_0_m1",
      "log_mu_inv",
      "log_omega",
      "logit_rho",
      "log_sqrt_phi_inv"
)
covmat <- diag(0.01, length(parameters))
rownames(covmat) <- colnames(covmat) <- covmat_names

mcmc_kernel <-
      kernel(
            method = "mvnss",
            sigma = covmat,
            scale_constant = 0.5,
            scale_cooling = 0.99,
            stop_adaptation = 5e4,
            step_size = 0.01,
            harss_warmup = 0,
            nugget = 0.001,
            mvnss_setting_list = 
                  mvnss_settings(n_mvnss_updates = 1,
                                 initial_bracket_width = 0.5,
                                 bracket_limits = c(0.125, 8),
                                 nugget_cooling = 0.99, 
                                 nugget_step_size = 0.01),
            messages = FALSE
      )

stem_object$dynamics$parameters <- function() {
      setNames(from_estimation_scale(to_estimation_scale(parameters) + rnorm(length(parameters), 0, 0.001)),
               names(parameters))
}

registerDoParallel(5)

results <- foreach(chain = 1:5,
                   .packages=c("stemr"),
                   .options.RNG = 52787,
                   .export = ls(all.names = T)) %dorng% {
                         
                         chain_res <- stem_inference(stem_object = stem_object,
                                                     method = "lna",
                                                     iterations = 1e5,
                                                     thin_params = 10,
                                                     thin_latent_proc = 10,
                                                     initialization_attempts = 500,
                                                     priors = priors,
                                                     mcmc_kernel = mcmc_kernel,
                                                     print_progress = 1e3,
                                                     t0_kernel = NULL,
                                                     messages = FALSE,
                                                     ess_args = ess_settings(n_ess_updates = 5,
                                                                             n_tparam_updates = 5,
                                                                             ess_warmup = 500, 
                                                                             initdist_bracket_width = 2*pi,
                                                                             initdist_bracket_update = 5e3,
                                                                             tparam_bracket_width = 2*pi,
                                                                             tparam_bracket_update = 5e3,
                                                                             lna_bracket_width = 2*pi,
                                                                             lna_bracket_update = 5e3))
                         return(chain_res)
                   }

results$true_path <- true_path
save(results, file = "sir_sinfoi_rw1_lna.Rdata")
