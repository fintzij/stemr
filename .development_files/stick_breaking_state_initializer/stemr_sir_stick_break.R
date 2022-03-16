library(tidyverse)
library(tidybayes)
require(cowplot)
library(stemr)
set.seed(12511)
popsize = 1e4 # population size

true_pars =
  c(R0     = 1.5,  # basic reproduction number
    mu_inv = 2,    # infectious period duration = 2 days
    rho    = 0.5,  # case detection rate
    phi    = 10)   # negative binomial overdispersion

# initialize model compartments and rates
strata <- "" # no strata
compartments <- c("S", "I", "R")

# rates initialized as a list of rate lists
rates <-
  list(rate(rate = "beta * I", # individual level rate (unlumped)
            from = "S",        # source compartment
            to   = "I",        # destination compartment
            incidence = T),    # compute incidence of S2I transitions, required for simulating incidence data
       rate(rate = "mu",       # individual level rate
            from = "I",        # source compartment
            to   = "R",        # destination compartment
            incidence = TRUE)) # compute incidence of I2R transitions (not required for simulating data)

# list used for simulation/inference for the initial state, initial counts fixed.
# state initializer a list of stem_initializer lists.
state_initializer <-
  list(stem_initializer(
    init_states = c(S = popsize-10, I = 10, R = 0), # must match compartment names
    prior = c(4.5, 36, 0.5, 0.5),
    dist = "sbln",
    fixed = F)) # initial state fixed for simulation, we'll change this later

# state_initializer <-
#   list(stem_initializer(
#     init_states = c(S = popsize-10, I = 10, R = 0), # must match compartment names
#     prior = c(popsize-10, -10, 0),
#     dist = "multinomial",
#     fixed = F)) # initial state fixed for simulation, we'll change this later

# state_initializer <-
#   list(
#     stem_initializer(
#       init_states = c(S_guin = 19970, E_guin = 15, I_guin = 10, R_guin = 5),
#       fixed = TRUE,
#       strata = "guin",
#       prior = c(11799970, 15, 10, 5),
#       dist = "multinom"
#       ),
#     stem_initializer(
#       init_states = c(S_lib = 34970, E_lib = 15, I_lib = 10, R_lib = 5),
#       fixed = TRUE,
#       strata = "lib",
#       prior = c(4399970, 15, 10, 5),
#       dist = "multinom"
#     ),
#     stem_initializer(
#       init_states = c(S_sln = 24970, E_sln = 15, I_sln = 10, R_sln = 5),
#       fixed = TRUE,
#       strata = "sln",
#       prior = c(7099970, 15, 10, 5),
#       dist = "multinom"
#   )
# )

# set the parameter values - must be a named vector
parameters =
  c(true_pars["R0"] / popsize / true_pars["mu_inv"], # R0 = beta * P / mu
    1/true_pars["mu_inv"],
    true_pars["rho"],
    true_pars["phi"])
names(parameters) <- c("beta", "mu", "rho", "phi")

# declare the initial time to be constant
constants <- c(t0 = 0)
t0 <- 0; tmax <- 40

# compile the model
dynamics <-
  stem_dynamics(
    rates = rates,
    tmax = tmax,
    parameters = parameters,
    state_initializer = state_initializer,
    compartments = compartments,
    constants = constants,
    compile_ode = T,   # compile ODE functions
    compile_rates = T, # compile MJP functions for Gillespie simulation
    compile_lna = T,   # compile LNA functions
    messages = F       # don't print messages
  )

# list of emission distribution lists (analogous to rate specification)
emissions <-
  list(emission(meas_var = "S2I", # transition or compartment being measured (S->I transitions)
                distribution    = "negbinomial",         # emission distribution
                emission_params = c("phi", "S2I * rho"), # distribution pars, here overdispersion and mean
                incidence       = TRUE,                  # is the data incidence
                obstimes        = seq(1, tmax, by =1)))  # vector of observation times

# compile the measurement process
measurement_process <-
  stem_measure(emissions = emissions,
               dynamics  = dynamics,
               messages  = F)

# put it all together into a stochastic epidemic model object
stem_object <-
  make_stem(
    dynamics = dynamics,
    measurement_process = measurement_process)



# Need check storage of initial compartments ------------------------------
set.seed(100)
sim_mjp_200 <- simulate_stem(stem_object = stem_object, method = "gillespie", nsim = 20)
set.seed(100)
sim_ode_200 <- simulate_stem(stem_object = stem_object, method = "ode", nsim = 20)
set.seed(100)
sim_lna_200 <- simulate_stem(stem_object = stem_object, method = "lna", nsim = 20)


to_estimation_scale = function(params_nat) {
  c(log(params_nat[1] * popsize / params_nat[2] - 1), # (beta,mu,N) -> log(R0-1)
    log(params_nat[2]),                     # mu -> log(mu)
    logit(params_nat[3]),                   # rho -> logit(rho)
    log(params_nat[4]))                     # phi -> log(phi)
}

# function to take params_est and return params_nat
from_estimation_scale = function(params_est) {
  c(exp(log(exp(params_est[1])+1) + params_est[2] - log(popsize)), # (log(R0), log(mu), N) -> beta = exp(log(R0) + log(mu) - log(N))
    exp(params_est[2]), # log(mu) -> mu
    expit(params_est[3]), # logit(rho) -> rho
    exp(params_est[4])) # log(phi) -> phi
}

# calculate the log prior density. note the jacobian for phi
logprior =
  function(params_est) {
    sum(dnorm(params_est[1], 0, 1, log = TRUE),
        dnorm(params_est[2], -0.7, 0.35, log = TRUE),
        dnorm(params_est[3], 0, 1, log = TRUE),
        dexp(exp(params_est[4]), 0.1, log = TRUE) + params_est[4])
  }

# return all three functions in a list
priors <- list(logprior = logprior,
               to_estimation_scale = to_estimation_scale,
               from_estimation_scale = from_estimation_scale)


#'
#' We now specify the MCMC transition kernel. In this simple example, we'll update
#' the model hyperparameters using a multivariate Metropolis algorithm. We'll tune
#' the algorithm using a global adaptive scheme (algorithm 4 in Andrieu and Thoms). We'll also initialize the parameters at random values, which is done by replacing the vector of parameters in the stem object with a function that returns a named vector of parameters.
#'
## ----mcmc_kern, echo = TRUE----------------------------------------------
# specify the initial proposal covariance matrix with row and column names
# corresponding to parameters on their estimation scales

par_initializer = function() {
  priors$from_estimation_scale(priors$to_estimation_scale(parameters) +
                                 rnorm(4, 0, 0.1))
}

measurement_process <-
  stem_measure(emissions = emissions,
               dynamics = dynamics,
               data = sim_mjp_200$datasets[[1]])
stem_object <- make_stem(dynamics = dynamics,
                         measurement_process = measurement_process)

# specify the kernel
mcmc_kern <-
  mcmc_kernel(
    parameter_blocks =
      list(parblock(
        pars_nat = c("beta", "mu", "rho", "phi"),
        pars_est = c("log_R0", "log_mu", "logit_rho", "log_phi"),
        priors = priors,
        # alg = "mvnss",
        alg = "mvnmh",
        sigma = diag(0.01, 4),
        initializer = par_initializer,
        control =
          # mvnss_control(stop_adaptation = 2.5e4))),
          mvnmh_control(stop_adaptation = 2.5e2))),
    lna_ess_control = lna_control(bracket_update_iter = 50,
                                  joint_initdist_update = TRUE))


# debug(initdist_update)
res <-
  fit_stem(stem_object = stem_object,
           method = "ode",
           mcmc_kern = mcmc_kern,
           iterations = 5e4,
           # iterations = 5e2,
           print_progress = 1000)

write_rds(res, "res.rds")
res <- read_rds("res.rds")
res$results$posterior$initdist_samples

source("/Users/damon/Documents/uci_covid_modeling2/code/stemr_functions.R")


# Posterior Predictive looks good -----------------------------------------

simulation_parameters_list <- split_along_dim(res$results$posterior$parameter_samples_nat, 1)

init_dist_list <- split_along_dim(res$results$posterior$initdist_samples, 1)

sim_results_tbl <- map2_dfr(simulation_parameters_list, init_dist_list,
                            function(simulation_parameters, init_dist) {
                              stem_object$dynamics$fixed_inits <- T
                              stem_object$dynamics$initdist_params <-
                                stem_object$dynamics$initdist_priors <-
                                stem_object$dynamics$initializer[[1]]$init_states <-
                                stem_object$dynamics$initializer[[1]]$prior <-
                                init_dist
                              stem_object$dynamics$parameters <- simulation_parameters
                              tibble(sim_result = list(
                                unlist(
                                  simulate_stem(stem_object = stem_object,
                                                method = "ode"),
                                  recursive = F)))
                            }) %>%
  unnest_wider(sim_result)




pp <- sim_results_tbl %>%
  select(datasets) %>%
  # mutate(.iteration = 1:n()) %>%
  mutate(datasets = map(datasets, as_tibble)) %>%
  unnest(datasets) %>%
  group_by(time) %>%
  median_qi(.width = c(0.5, 0.8, 0.95))

ggplot() +
  geom_lineribbon(data = pp, aes(time, S2I, ymin = .lower, ymax = .upper)) +
  geom_point(data = as_tibble(measurement_process$data), aes(time, S2I))

