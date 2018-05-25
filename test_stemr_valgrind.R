library(stemr)
library(extraDistr)
# library(foreach)
# library(doRNG)
# library(doParallel)

# no strata the stemr object --------------------------------------------------
set.seed(12511)
strata <- NULL
compartments <- c("S", "I", "R")
rates <- list(rate("beta * I", "S", "I", incidence = T),
              rate("mu", "I", "R"))
state_initializer <- list(stem_initializer(c(S = 50000, I = 10, R = 5), fixed = T, prior = c(50000, 10, 5)))
adjacency <- NULL
tcovar <- NULL
parameters <- c(beta = 0.00001, mu = 1/7, rho = 0.5, phi = 10)
constants <- c(t0 = 0)
t0 <- 0; tmax <- 52; 

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
            atol = 1e-6
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

# set the seed and recompile the stemr object
init_state <- extraDistr::rdirichlet(1, c(50000, 10, 5)) * 50015
state_initializer <- list(stem_initializer(c(S = init_state[1], I = init_state[2], R = init_state[3]),
                                           fixed = FALSE, prior = c(50000, 10, 5)))
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
            step_size = 1e-6,
            rtol = 1e-6,
            atol = 1e-6
      )

measurement_process <- stem_measure(emissions = emissions, dynamics = dynamics, data = dat)
stem_object <- stem(dynamics = dynamics, measurement_process = measurement_process)

# initialize the inference
### Parameterization in terms of log(R0) and log(mu)
## Priors for log(R0), log(mu), rho, log(phi)
# Parameters (natural scale): beta, mu, rho, phi
# Parameters (estimation scale): log(beta * N / mu), log(mu), logit(rho), log(phi)

priors <- list(prior_density =
                     function(params_nat, params_est) {
                           sum(dnorm(params_est[1], mean = 1.2, sd = 1, log = TRUE),      # log(R0)
                               dnorm(params_est[2], mean = -2, sd = 1, log = TRUE),       # log(mu)
                               dnorm(params_est[3], mean = -0.1, sd = 1, log = TRUE),      # logit(rho)
                               dnorm(params_est[4], mean = 1.5, sd=1, log = TRUE))        # log(phi)
                     },
               to_estimation_scale = function(params) {
                     c(log(50015) + log(params[1]) - log(params[2]), # (beta,mu,N) -> log(R0)
                       log(params[2]),                     # mu -> log(mu)
                       logit(params[3]),                   # rho -> logit(rho)
                       log(params[4]))                     # phi -> log(phi)
               },
               from_estimation_scale = function(params) {
                     c(exp(params[1] + params[2] - log(50015)), # (log(R0), log(mu), N) -> beta = exp(log(R0) + log(mu) - log(N))
                       exp(params[2]), # log(mu) -> mu
                       expit(params[3]), # logit(rho) -> rho
                       exp(params[4])) # log(phi) -> phi
               })

t0_kernel <- NULL
covmat <- diag(1e-2, 4);
rownames(covmat) <- colnames(covmat) <- c("log_R0", "log_mu", "logit_rho", "log_phi")
mcmc_kernel <-
      kernel(
            method = "mvnss",
            sigma = covmat,
            scale_cooling = 0.50001,
            stop_adaptation = 5e3,
            harss_warmup = 1e3,
            nugget = 0.00001,
            mvnss_setting_list = mvnss_settings(n_mvnss_updates = 1, cov_update_interval = 10),
            harss_setting_list = harss_settings(n_harss_updates = 1),
            afss_setting_list = afss_settings(factor_update_interval = 100,
                                              prob_update_interval = 100,
                                              first_factor_update = 1,
                                              first_prob_update = 1,
                                              n_afss_updates = 2,
                                              target_prop_totsd = 0.95),
            messages = FALSE
      )

stem_object$dynamics$parameters <- function() {
      priors$from_estimation_scale(priors$to_estimation_scale(parameters) + rnorm(4, 0, 0.01))
}

chain_res <- stem_inference(stem_object = stem_object,
                            method = "ode",
                            iterations = 1e3,
                            thin_params = 10,
                            thin_latent_proc = 10,
                            initialization_attempts = 500,
                            priors = priors,
                            mcmc_kernel = mcmc_kernel,
                            t0_kernel = t0_kernel,
                            messages = FALSE,
                            ess_args = ess_settings(ess_warmup = 1e3))