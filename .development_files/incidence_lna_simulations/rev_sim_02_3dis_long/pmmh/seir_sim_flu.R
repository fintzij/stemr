# Overview
# 
# Code for fitting an SEIR model to weekly and monthly incidence data simulated
# from an SEIR model with influenza-like dynamics. 
#
# Settings: 
# - R0 = 1.3
# - Mean latent period duration = 1.6 days
# - Mean infectious period duration = 2 days
# - P  = 100,000
# - X(t_0) = (S_0 = 999000, E_0 = 50, I_0 = 50, R_0 = 0)
# - Mean case detection rate = 0.01

library(pomp)
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
popsize = 5e4 # population size

true_pars =
  c(R0         = 1.3, # basic reproduction number
    latent_dur = 2/7, # latent period duration = 2 days
    infec_dur  = 2.5/7, # infectious period duration = 2.5 days
    rho        = 0.02, # case detection rate
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
       rate(rate = "omega",       # individual level rate
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
          init_states = c(S = popsize-50, E = 25, I = 25, R = 0), # must match compartment names
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
  list(emission(meas_var = "cases", # transition or compartment being measured (S->I transitions)
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

while(max(sim_mjp$datasets[[1]][,"cases"]) < 15) {
    sim_mjp <- simulate_stem(stem_object = stem_object, method = "gillespie")
}

# grab the true path and the dataset
dat_sim <- sim_mjp$datasets[[1]]
true_path <- sim_mjp$paths[[1]]

# chop off the trailing data after the outbreak has ended (eight consecutive zeros)
rle_dat <- rle(dat_sim[-c(1:8),2] == 0) # run length encoding representation
runs    <- which(rle_dat$values & rle_dat$lengths >= 8) # find runs of zeros of length >= 4

# if a run of zeros was found, truncate data
if(length(runs) != 0) {
  
  inds <- 
    if(runs[1] == 1) {
      9:16
    } else {
      c(seq(9, length.out = sum(rle_dat$length[seq_len(runs[1] - 1)])),
        seq(9 + sum(rle_dat$length[seq_len(runs[1] - 1)]), length.out = 8))
    }
  
  # truncate at one year
  if(any(inds > 52)) inds <- 9:52
  
  dat <- rbind(dat_sim[1:8,],
               dat_sim[inds,])
  
} else {
  dat <- dat_sim
}

# Set up pomp objects -----------------------------------------------------

# Measurement process objects
rmeas <- "
  cases=rnbinom_mu(exp(-2*log_sqrt_phi_inv),E2I/(1+exp(-logit_rho)));    // simulates the data
"

dmeas <- "
  lik=dnbinom_mu(cases,exp(-2*log_sqrt_phi_inv),E2I/(1+exp(-logit_rho)),give_log); // emission density
"

# Define the priors
seir.dprior = 
  "double loglik;
   loglik = 
      dnorm(log_R0, log(1.5), 0.5, 1) + 
      dnorm(log_latent_dur, log(2.25) - log(7), 0.82, 1) + 
      dnorm(log_infec_dur, log(2.25) -log(7), 0.82, 1) + 
      logit_rho - 2 * log1p(exp(logit_rho)) + 
      dexp(exp(log_sqrt_phi_inv), exp(-log(3)), 1) + log_sqrt_phi_inv;
  lik = give_log == 1? loglik : exp(loglik);
"

# define the stepper
seir.step<-"
  double rate[3];
  double dN[3];
  rate[0]=exp(log_R0 - log_infec_dur - log(*get_userdata_double(\"popsize\")))*I; // Infection rate
  rate[1]=exp(-log_latent_dur);                // rate of E-I transitions
  rate[2]=exp(-log_infec_dur);                 // recovery rate
  reulermultinom(1,S,&rate[0],dt,&dN[0]);      // generate the number of newly exposed people
  reulermultinom(1,E,&rate[1],dt,&dN[1]);      // generate the number of newly infected people
  reulermultinom(1,I,&rate[2],dt,&dN[2]);      // generate the number of newly recovered people
  S+=-dN[0];                                   // update the number of Susceptibles
  E+=dN[0]-dN[1];                              // update the number of Infections
  I+=dN[1]-dN[2];                              // update the number of Infections
  R+=dN[2];                                    // update the number of Recoveries
  E2I += dN[1];
"

# instatiate the euler stepper function
SEIR_sim <- euler(step.fun = Csnippet(seir.step), delta.t = 1/24/7)

# instatiate the pomp object
parnames = c("log_R0", "log_latent_dur", "log_infec_dur", "logit_rho", "log_sqrt_phi_inv")
pars = c(log(true_pars["R0"]),
         log(true_pars["latent_dur"]),
         log(true_pars["infec_dur"]),
         qlogis(true_pars["rho"]),
         -0.5 * log(true_pars["phi"]))
names(pars) = parnames

seir_mod <- pomp(
  data = as.data.frame(dat),  #"cases" is the dataset, "time" is the observation time
  times = "time",
  t0 = 0,                      # initial time point
  popsize = popsize,
  dmeasure = Csnippet(dmeas),  # evaluates the density of the measurement process
  rmeasure = Csnippet(rmeas),  # simulates from the measurement process
  rprocess = SEIR_sim,          # simulates from the latent process
  statenames = c("S", "E", "I", "R", "E2I"),  #state space variable name
  params = pars,
  paramnames = parnames, #parameters name
  accumvars = c("E2I"),
  dprior = Csnippet(seir.dprior),
  rinit = function(params, t0, ...) {
    return(c(S = popsize-50, E = 25, I = 25, R = 0, E2I = 0))
  }
)

# Metropolis kernel proposal covariance matrix
cov_init <- diag(1e-2, length(true_pars)); 
colnames(cov_init) <- rownames(cov_init) <- 
  c("log_R0", "log_latent_dur", "log_infec_dur", "logit_rho", "log_sqrt_phi_inv")

# run the MCMC
registerDoParallel(cores = future::availableCores())
set.seed(52787 + replication)

fits <- 
  foreach(i = seq_len(5),
          .export = ls(),
          .packages = c("pomp")) %dorng% {
            
            init <- pars + rnorm(length(parameters), 0, 0.1)
            
            adapt_res_1 <- 
              pmcmc(seir_mod,
                    Nmcmc = 1e3,
                    Np = 500,
                    params = init,
                    proposal =
                      mvn.rw.adaptive(
                        rw.sd = diag(cov_init),
                        scale.start = 100,
                        shape.start = 100,
                        scale.cooling = 0.999))
            
            while(adapt_res_1@accepts < 100 | adapt_res_1@accepts > 800) {
              
              init <- pars + rnorm(length(parameters), 0, 0.1)
              
              adapt_res_1 <- 
                pmcmc(seir_mod,
                      Nmcmc = 1e3,
                      Np = 500,
                      params = init,
                      proposal =
                        mvn.rw.adaptive(
                          rw.sd = diag(cov_init),
                          scale.start = 100,
                          shape.start = 100,
                          scale.cooling = 0.999))
            }
            
            start.time <- Sys.time()
            adapt_res_2 <- 
              pmcmc(adapt_res_1,
                    Nmcmc = 1e4,
                    Np = 500,
                    proposal =
                      mvn.rw.adaptive(
                        rw.var = covmat(adapt_res_1) + diag(0.01, length(true_pars)),
                        scale.start = 100,
                        shape.start = 100,
                        scale.cooling = 0.999))
            
            while(adapt_res_2@accepts < 500 | adapt_res_2@accepts > 8e3) {
              
              init <- pars + rnorm(length(parameters), 0, 0.1)
              
              adapt_res_2 <- 
                pmcmc(adapt_res_1,
                      Nmcmc = 1e4,
                      Np = 500,
                      params = init,
                      proposal =
                        mvn.rw.adaptive(
                          rw.var = covmat(adapt_res_1) + diag(0.01, length(true_pars)),
                          scale.start = 100,
                          shape.start = 100,
                          scale.cooling = 0.999))
            }
            
            # stop the adaptation and do some more MCMC
            fit <- 
              pmcmc(adapt_res_2, 
                    Nmcmc = 5e4, 
                    Np = 500,
                    proposal = mvn.rw(rw.var = covmat(adapt_res_2)))
            
            end.time <- Sys.time()
            
            return(list(fit = fit, runtime = difftime(end.time, start.time, units = "hours")))
          }

# extract results
posts <- 
  do.call(rbind, 
          lapply(fits, 
                 function(x) traces(x$fit, pars = parnames)))
post_quants <- apply(posts, 2, quantile, c(0.025, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.975)) 

posts_mcmc <- 
  as.mcmc.list(lapply(fits, 
                      function(x) 
                        as.mcmc(cbind(
                          logpost = 
                            rowSums(traces(x$fit, pars = c("loglik", "log.prior"))),
                          traces(x$fit, pars = parnames)))))
psrf              <- gelman.diag(posts_mcmc)
effective_samples <- do.call(rbind, lapply(posts_mcmc, effectiveSize))  

results <- 
  list(true_pars = true_pars,
       dat = dat,
       post_quants = post_quants,
       effective_samples = effective_samples,
       psrf = psrf,
       times = sapply(fits, function(x) x$runtime))

saveRDS(results, file = paste0("seir_flu_",replication,".Rds"))
