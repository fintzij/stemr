# Overview
# 
# Code for fitting an SEIR model to weekly and monthly incidence data simulated
# from an SEIR model with COVID-like dynamics. 
#
# Settings: 
# - R0 = 2.5
# - Mean latent period duration = 1 week
# - Mean infectious period duration = 10 days
# - P  = 5000
# - X(t_0) = (S_0 = 9990, E_0 = 5, I_0 = 5, R_0 = 0)
# - Mean case detection rate = 0.1

library(stemr)
library(coda)
library(foreach)
library(doParallel)
library(doRNG)
library(future)

# load simulation arguments
args <- commandArgs(TRUE)
print(args)
replication  <- as.numeric(args[1])
set.seed(12511 + replication)

# create stem model
popsize = 1e4 # population size

true_pars =
  c(R0         = 2.5, # basic reproduction number
    latent_dur = 5/7, # latent period duration = 5 days
    infec_dur  = 10/7, # infectious period duration = 10 days
    rho        = 0.1, # case detection rate
    phi        = 36)  # negative binomial overdispersion

# initialize model compartments and rates
strata <- NULL # no strata
compartments <- c("S", "E", "I", "R")

# rates initialized as a list of rate lists
rates <-
  list(rate(rate = "beta * I", # individual level rate (unlumped)
            from = "S",        # source compartment
            to   = "E",        # destination compartment
            incidence = FALSE),# don't compute incidence of S2E transitions, not required for simulating incidence data
       rate(rate = "omega",    # individual level rate
            from = "E",        # source compartment
            to   = "I",        # destination compartment
            incidence = TRUE), # compute incidence of IE2I transitions (required for simulating data)
       rate(rate = "mu",       # individual level rate
            from = "I",        # source compartment
            to   = "R",        # destination compartment
            incidence = FALSE)) # don't compute incidence of I2R transitions (not required for simulating data)

# list used for simulation/inference for the initial state, initial counts fixed.
# state initializer a list of stem_initializer lists.
state_initializer <-
  list(stem_initializer(
          init_states = c(S = popsize-10, E = 5, I = 5, R = 0), # must match compartment names
          fixed = T)) # initial state fixed for simulation, we'll change this later

# set the parameter values - must be a named vector
parameters =
  c(true_pars["R0"] / popsize / true_pars["infec_dur"], # R0 = beta * P / mu
    1/true_pars["latent_dur"],
    1/true_pars["infec_dur"],
    true_pars["rho"],
    true_pars["phi"])
names(parameters) <- c("beta", "omega", "mu", "rho", "phi")

# declare the initial time to be constant
constants <- c(t0 = 0)
t0 <- 0; tmax <- 52

# compile the model
dynamics <-
      stem_dynamics(
            rates = rates,
            tmax = tmax,
            parameters = parameters,
            state_initializer = state_initializer,
            compartments = compartments,
            constants = constants,
            compile_ode = F,   # compile ODE functions
            compile_rates = T, # compile MJP functions for Gillespie simulation
            compile_lna = T,   # compile LNA functions
            messages = F       # don't print messages
      )

# list of emission distribution lists (analogous to rate specification)
emissions <-
  list(emission(meas_var = "E2I", # transition or compartment being measured (S->I transitions)
                distribution    = "negbinomial",         # emission distribution
                emission_params = c("phi", "E2I * rho"), # distribution pars, here overdispersion and mean
                incidence       = TRUE,                  # is the data incidence
                obstimes        = seq(1, tmax, by = 1)))  # vector of observation times

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

 
# Simulating an outbreak and data 
sim_mjp <- simulate_stem(stem_object = stem_object, method = "gillespie")

while(max(sim_mjp$datasets[[1]][,"E2I"]) < 15) {
    sim_mjp <- simulate_stem(stem_object = stem_object, method = "gillespie")
}

# grab the true path and the dataset
dat_sim <- sim_mjp$datasets[[1]]
true_path <- sim_mjp$paths[[1]]

# chop off the trailing data after the outbreak has ended (four consecutive zeros)
rle_dat <- rle(dat_sim[-c(1:8),2] < 5) # run length encoding representation
runs    <- which(rle_dat$values & rle_dat$lengths >= 4) # find runs of zeros of length >= 4

# if a run of zeros was found, truncate data
if(length(runs) != 0) {
  
  inds <- 
    if(runs[1] == 1) {
      9:12
    } else {
      c(seq(9, length.out = sum(rle_dat$length[seq_len(runs[1] - 1)])),
        seq(9 + sum(rle_dat$length[seq_len(runs[1] - 1)]), length.out = 4))
    }
  
  # truncate at one year
  if(any(inds > 52)) inds <- 9:52
  
  dat <- rbind(dat_sim[1:8,],
               dat_sim[inds,])
  
} else {
  dat <- dat_sim
}

# recompile the model for inference
t0 <- 0; tmax <- max(dat[,"time"])

# compile the model
dynamics <-
  stem_dynamics(
    rates = rates,
    tmax = tmax,
    parameters = parameters,
    state_initializer = state_initializer,
    compartments = compartments,
    constants = constants,
    compile_ode = F,   # compile ODE functions
    compile_rates = T, # compile MJP functions for Gillespie simulation
    compile_lna = T,   # compile LNA functions
    messages = F       # don't print messages
  )

# recompile emission lists with new tmax
emissions <-
  list(emission(meas_var = "E2I", # transition or compartment being measured (S->I transitions)
                distribution    = "negbinomial",         # emission distribution
                emission_params = c("phi", "E2I * rho"), # distribution pars, here overdispersion and mean
                incidence       = TRUE,                  # is the data incidence
                obstimes        = dat[,1]))              # vector of observation times

measurement_process <-
  stem_measure(emissions = emissions,
               dynamics = dynamics,
               data = dat)

stem_object <- make_stem(dynamics = dynamics, 
                         measurement_process = measurement_process)

## Priors for log(R0), log(omega/mu), log(1/mu), logit(rho), phi
# Parameters (natural scale): beta, mu, rho, phi
# Parameters (estimation scale): log(beta * N / mu), log(omega/mu), log(mu), logit(rho), log(phi)

# function to take params_nat and return params_est
to_estimation_scale = function(params_nat) {
  c(log(params_nat[1] * popsize / params_nat[3]), # (beta,mu,N) -> log(R0)
    -log(params_nat[2]),                          # omega, mu -> log(omega/mu)
    -log(params_nat[3]),                          # mu -> log(mu)
    qlogis(params_nat[4]),                        # rho -> logit(rho)
    -0.5 * log(params_nat[5]))                    # phi -> log(1/sqrt(phi))
}

# function to take params_est and return params_nat
from_estimation_scale = function(params_est) {
  c(exp(params_est[1] - params_est[3] - log(popsize)), # (log(R0), log(mu), N) -> beta = exp(log(R0) + log(mu) - log(N))
    exp(-params_est[2]),                               # log(omega/mu) -> omega, mu
    exp(-params_est[3]),                               # log(mu) -> mu
    plogis(params_est[4]),                             # logit(rho) -> rho
    exp(-2 * params_est[5]))                           # -0.5 * log(phi) -> phi
}

# calculate the log prior density. note the jacobian for rho and phi
logprior =
      function(params_est) {
            dnorm(params_est[1], log(2.4), 0.5, log = TRUE) + 
                dnorm(params_est[2], log(6/7), 0.82, log = TRUE) + 
                dnorm(params_est[3], log(9/7), 0.82, log = TRUE) + 
                params_est[4] - 2 * log1p(exp(params_est[4])) + 
                dexp(exp(params_est[5]), 3, log = TRUE) + params_est[5]
      }

# return all three functions in a list
priors <- list(logprior = logprior,
               to_estimation_scale = to_estimation_scale,
               from_estimation_scale = from_estimation_scale)

# specify the MCMC transition kernel
par_initializer = function() {
  priors$from_estimation_scale(priors$to_estimation_scale(parameters) + 
                                 rnorm(length(parameters), 0, 0.1))
}
  
# specify the kernel
mcmc_kern <-
        mcmc_kernel(
          parameter_blocks = 
              list(parblock(
                  pars_nat = c("beta", "omega", "mu", "rho", "phi"),
                  pars_est = c("log_R0", "log_latent_dur", "log_infec_dur", "logit_rho", "log_sqrt_phi_inv"),
                  priors = priors,
                  alg = "mvnmh",
                  sigma = diag(0.1, length(parameters)),
                  initializer = par_initializer,
                  control = 
                      mvnmh_control(stop_adaptation = 2.5e4))),
          lna_ess_control = lna_control(bracket_update_iter = 1e4))

# run the MCMC algorithm to fit the model - 5 chains in parallel, 5e4 posterior samples each
registerDoParallel(cores = future::availableCores())
set.seed(52787 + replication)

fits <- 
    foreach(i = seq_len(5),
            .export = ls(),
            .packages = c("stemr")) %dorng% { 
                
                fit_stem(stem_object = stem_object,
                         method = "lna",
                         mcmc_kern = mcmc_kern,
                         iterations = 7.5e4, 
                         thinning_interval = 50)
                }

# extract results
posts <- 
  do.call(rbind, 
          lapply(fits, 
                 function(x) x$results$posterior$parameter_samples_est))
post_quants <- apply(posts, 2, quantile, c(0.025, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.975)) 

posts_mcmc <- 
  as.mcmc.list(lapply(fits, 
                      function(x) 
                        as.mcmc(cbind(
                          logpost = 
                            x$results$posterior$data_log_lik + 
                            x$results$posterior$lna_log_lik + 
                            x$results$posterior$params_log_prior,
                          x$results$posterior$parameter_samples_est))))
psrf              <- gelman.diag(posts_mcmc)
effective_samples <- do.call(rbind, lapply(posts_mcmc, effectiveSize))  

results <- 
  list(true_pars = true_pars,
       dat = dat,
       post_quants = post_quants,
       effective_samples = effective_samples,
       psrf = psrf,
       times = sapply(fits, function(x) x$results$runtime))

saveRDS(results, file = paste0("seir_covid_",replication,".Rds"))
