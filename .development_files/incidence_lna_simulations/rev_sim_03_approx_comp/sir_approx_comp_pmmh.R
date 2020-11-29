# Overview
# 
# Code for fitting an SEIR model to weekly incidence data simulated
# from an SEIR model.
#
# Settings: 
# - R0 = 2
# - Mean latent period duration = 1 week
# - Mean infectious period duration = 1 week
# - P = 2000, 10000, or 50000
# - X(t_0) = (S_0 = 0.999 * P, E_0 = 0.0005 * P, I_0 = 0.0005 * P, R_0 = 0)
# - Mean case detection rate = 0.5

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
# popsize <- as.numeric(args[1])
# replication  <- as.numeric(args[2])
replication = 994; popsize = 2e3
set.seed(12511 + replication)

# create stem model
true_pars =
  c(R0         = 1 + rlnorm(1, 0, 0.56), # basic reproduction number
    infec_dur  = rlnorm(1, 0, 0.354), # infectious period duration = 1 week
    rho        = plogis(rnorm(1, 0, 1)), # case detection rate
    sqrt_phi_inv = rgamma(1, 2, 4))  # negative binomial overdispersion

# initialize model compartments and rates
strata <- NULL # no strata
compartments <- c("S", "I", "R")

# rates initialized as a list of rate lists
rates <-
  list(rate(rate = "beta * I", # individual level rate (unlumped)
            from = "S",        # source compartment
            to   = "I",        # destination compartment
            incidence = TRUE), # compute incidence of S2E transitions, required for simulating incidence data
       rate(rate = "mu",       # individual level rate
            from = "I",        # source compartment
            to   = "R",        # destination compartment
            incidence = FALSE)) # don't compute incidence of I2R transitions (not required for simulating data)

# list used for simulation/inference for the initial state, initial counts fixed.
# state initializer a list of stem_initializer lists.
state_initializer <-
  list(stem_initializer(
    init_states = c(S = 0.999 * popsize, I = 0.001 * popsize, R = 0), # must match compartment names
    fixed = T)) # initial state fixed for simulation, we'll change this later

# set the parameter values - must be a named vector
parameters =
  c(true_pars["R0"] / popsize / true_pars["infec_dur"], # R0 = beta * P / mu
    1/true_pars["infec_dur"],
    true_pars["rho"],
    1/true_pars["sqrt_phi_inv"]^2)
names(parameters) <- c("beta", "mu", "rho", "phi")

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
    compile_ode = TRUE,   # compile ODE functions
    compile_rates = TRUE, # compile MJP functions for Gillespie simulation
    compile_lna = TRUE,   # compile LNA functions
    messages = F       # don't print messages
  )

# list of emission distribution lists (analogous to rate specification)
emissions <-
  list(emission(meas_var = "cases", # cases = observed transitions (E->I transitions)
                distribution    = "negbinomial",         # emission distribution
                emission_params = c("phi", "S2I * rho"), # distribution pars, here overdispersion and mean
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

# Set up pomp objects -----------------------------------------------------

# Measurement process objects
rmeas <- "
  cases=rnbinom_mu(exp(-2*log_sqrt_phi_inv),exp(logit_rho)*S2I/(1+exp(logit_rho)));    // simulates the data
"

dmeas<-"
  lik=dnbinom_mu(cases,exp(-2*log_sqrt_phi_inv),exp(logit_rho)*S2I/(1+exp(logit_rho)),give_log); // emission density
"

# Define the priors
sir.dprior = 
  "double loglik;
   loglik = 
      dnorm(log_R0_m1, 0, 0.56, 1) + 
      dnorm(log_infec_dur, 0, 0.354, 1) + 
      dnorm(logit_rho, 0, 1, 1) + 
      dgamma(exp(log_sqrt_phi_inv), 2, exp(-log(4)), 1) + log_sqrt_phi_inv;
  lik = give_log == 1? loglik : exp(loglik);
"

# define the stepper
sir.step<-"
  double rate[2];
  double dN[2];
  rate[0]=exp(log(exp(log_R0_m1)+1) - log_infec_dur - log(*get_userdata_double(\"popsize\")))*I; // Infection rate
  rate[1]=exp(-log_infec_dur);                 // recovery rate
  reulermultinom(1,S,&rate[0],dt,&dN[0]);      // generate the number of newly infected people
  reulermultinom(1,I,&rate[1],dt,&dN[1]);      // generate the number of newly recovered people
  S+=-dN[0];                                   // update the number of Susceptibles
  I+=dN[0]-dN[1];                              // update the number of Infections
  R+=dN[1];                                    // update the number of Recoveries
  S2I += dN[0];
"

# instatiate the euler stepper function
SIR_sim <- euler(step.fun = Csnippet(sir.step), delta.t = 1/7/24)

# instatiate the pomp object
parnames = c("log_R0_m1", "log_infec_dur", "logit_rho", "log_sqrt_phi_inv")
pars = c(log(true_pars["R0"] - 1),
         log(true_pars["infec_dur"]),
         qlogis(true_pars["rho"]),
         log(true_pars["sqrt_phi_inv"]))
names(pars) = parnames

sir_mod <- pomp(
  data = as.data.frame(dat),  #"cases" is the dataset, "time" is the observation time
  times = "time",
  t0 = 0,                      # initial time point
  popsize = popsize,
  dmeasure = Csnippet(dmeas),  # evaluates the density of the measurement process
  rmeasure = Csnippet(rmeas),  # simulates from the measurement process
  rprocess = SIR_sim,          # simulates from the latent process
  statenames = c("S", "I", "R", "S2I"),  #state space variable name
  params = pars,
  paramnames = parnames, #parameters name
  accumvars = c("S2I"),
  dprior = Csnippet(sir.dprior),
  rinit = function(params, t0, ...) {
    return(c(S = 0.999 * popsize, I = 0.001 * popsize, R = 0, S2I = 0))
  }
)

# Metropolis kernel proposal covariance matrix
cov_init <- diag(1e-2, length(parameters)); 
colnames(cov_init) <- rownames(cov_init) <- parnames

# run the MCMC
registerDoParallel(cores = future::availableCores())
set.seed(52787 + replication)

fits <- 
  foreach(i = seq_len(5),
          .export = ls(),
          .packages = c("pomp")) %dorng% {
            
            init <- pars + rnorm(length(true_pars), 0, 0.1)
            
            start.time <- Sys.time()
            adapt_res_1 <- 
              pmcmc(sir_mod,
                    Nmcmc = 1e3,
                    Np = 250,
                    params = init,
                    proposal =
                      mvn.rw.adaptive(
                        rw.sd = diag(cov_init),
                        scale.start = 100,
                        shape.start = 100,
                        scale.cooling = 0.999))
            
            while(adapt_res_1@accepts < 1e2) {
              
              init <- pars + rnorm(length(true_pars), 0, 0.1)
              
              adapt_res_1 <- 
                pmcmc(sir_mod,
                      Nmcmc = 1e4,
                      Np = 250,
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
                    Nmcmc = 2.5e4,
                    Np = 250,
                    proposal =
                      mvn.rw.adaptive(
                        rw.var = covmat(adapt_res_1) + diag(0.01, length(true_pars)),
                        scale.start = 100,
                        shape.start = 100,
                        scale.cooling = 0.999))
            
            while(adapt_res_2@accepts < 500) {
              
              init <- pars + rnorm(length(true_pars), 0, 0.1)
              
              adapt_res_2 <- 
                pmcmc(adapt_res_1,
                      Nmcmc = 1e4,
                      Np = 250,
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
                    Np = 250,
                    proposal = mvn.rw(rw.var = covmat(adapt_res_2))) 
            
            end.time <- Sys.time()
            
            return(list(fit = fit, adapt_res_1 = adapt_res_1, adapt_res_2 = adapt_res_2, runtime = difftime(end.time, start.time, units = "hours")))
          }

# # extract results
# posts <- 
#   do.call(rbind, 
#           lapply(fits, 
#                  function(x) traces(x$fit, pars = parnames)))
# post_quants <- apply(posts, 2, quantile, c(0.025, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.975)) 
# 
# posts_mcmc <- 
#   as.mcmc.list(lapply(fits, 
#                       function(x) 
#                         as.mcmc(cbind(
#                                    logpost = 
#                                      rowSums(traces(x$fit, pars = c("loglik", "log.prior"))),
#                                    traces(x$fit, pars = parnames)))))
# psrf              <- gelman.diag(posts_mcmc)
# effective_samples <- do.call(rbind, lapply(posts_mcmc, effectiveSize))  
# 
# results <- 
#   list(true_pars = true_pars,
#        dat = dat,
#        post_quants = post_quants,
#        effective_samples = effective_samples,
#        psrf = psrf,
#        times = sapply(fits, function(x) x$runtime))
# 
# saveRDS(results, file = paste0("sir_pmmh_popsize_",popsize,"_",replication,".Rds"))
