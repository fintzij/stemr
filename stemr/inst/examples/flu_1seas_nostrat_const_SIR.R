library(stemr)
library(extraDistr)
library(foreach)
library(doParallel)
library(doRNG)

set.seed(12511)

# Load the data -----------------------------------------------------------
popsize <- 5351427
log_popsize <- log(popsize)

dat <- data.frame(week = 1:113,
                  cases_young = 
                        c(0,0,0,0,0,0,0,0,0,2,5,6,15,20,18,15,10,5,4,5,2,20,16,11,8,2,28,66,
                          159,414,1168,1524,1034,437,120,19,10,2,1,0,2,1,0,0,1,0,0,0,0,0,0,0,
                          0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,2,0,2,
                          8,17,66,38,36,30,56,81,68,59,46,43,29,35,31,21,15,6,8,5,2,2,1,1,0,0,0),
                  cases_adult = 
                        c(0,0,0,0,2,0,5,3,0,4,8,5,26,17,23,13,6,3,0,7,5,3,7,5,7,5,22,49,160,251,
                          491,670,480,231,80,40,15,13,5,5,6,0,2,1,0,0,0,2,0,0,0,0,0,0,0,0,0,0,0,
                          0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,1,1,1,0,4,16,48,60,97,84,
                          137,117,151,152,118,110,106,92,46,36,33,15,13,14,6,5,3,0,0,0,0))

dat <- data.frame(week = dat[,1], cases = rowSums(dat[,2:3]))

vaccinations <- 
      data.frame(week = 1:113,
                 doses_young = 
                       c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9,22,8,124,17,35,43,190,4558,
                         47901,111321,174252,198130,168458,112546,30541,2198,1735,1851,3715,3609,
                         3157,3111,2614,1852,1826,1514,1115,933,1137,525,423,575,555,477,385,267,
                         492,834,340,221,199,138,163,89,92,65,82,86,67,88,17,0,0,0,0,0,0,0,0,0,0,
                         0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
                 doses_adult = 
                       c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,43,536,26,37,233,67,251,12732,
                         90590,159005,172237,61612,22662,20411,40267,99941,21796,22854,50575,
                         134735,163487,164602,175421,152297,104068,64992,35803,25016,18741,19059,
                         8328,4985,5984,5042,5235,3014,2521,1989,3218,1431,488,540,308,237,185,
                         330,421,138,228,268,349,113,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                         0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0))

vaccinations <- data.frame(week = vaccinations[,1], doses = rowSums(vaccinations[,2:3]))

# get rid of the second season
dat <- dat[1:52,]
vaccinations <- vaccinations[1:52,]

vaccprop_s1 <- sum(vaccinations[1:52,2])/popsize


# build the stemr object --------------------------------------------------
# model compartments
strata <- c("u", "v")
compartments <- list(S = "ALL", I = "ALL", R = "ALL", D = "ALL")

# model rates and forcings (deterministic compartment flows)
rates <- list(rate("(alpha_t + beta_t * (I_u + I_v)) * S_u", from = "S", to = "I", strata = "u", incidence = T, lumped = TRUE),
              rate("VE_susc * (alpha_t + beta_t * (I_u + I_v)) * S_v", from = "S", to = "I", strata = "v", incidence = T, lumped = TRUE),
              rate("mu", "I_u", "R_u", "u", incidence = TRUE),
              rate("mu", "I_v", "R_v", "v", incidence = TRUE))

forcings <- list(forcing("doses", c("S_u", "I_u", "R_u","D_u"), c("S_v", "I_v", "R_v", "D_v")))
constants <- c(t0 = 0, popsize = popsize)
tcovar    <- vaccinations
t0 <- 0; tmax <- nrow(dat);

## set the initial parameters
inits_u <- c(S = popsize * 0.25, I = 0, R = 0, D = popsize * 0.75)
inits_v <- c(S = 0, I = 0, R = 0, D = 0)

state_initializer <- list(stem_initializer(inits_u, fixed = F, strata = "u", prior = c(25, 0, 0, 75), dist = "dirmultinom"),
                          stem_initializer(inits_v, fixed = T, strata = "v"))

# adjacency matrix
adjacency <- matrix(1, nrow = length(strata), ncol = length(strata)); diag(adjacency) <- 0
colnames(adjacency) = rownames(adjacency) = strata

# parameters and other objects
parameters = c(mu       = 3,
               VE_susc  = 0.25,
               rho      = 0.0125,
               phi      = 1)

foi_rw2 <- function(parameters, draws) {
      
      log_R0_t <- numeric(2)
      log_R0_t[1] <- log(1.15) + 0.06 * draws[1]
      log_R0_t[2] <- log(1.4)  + 0.145 * draws[2]
      
      return(c(exp(log_R0_t + log(parameters["mu"]) - log(parameters["S_u_0"]))))
}

alpha_rw1 <- function(parameters, draws) {
      
      log_alpha_t <- numeric(2)
      log_alpha_t[1] <- 1.6 + 0.25 * draws[1] 
      log_alpha_t[2] <- 2.7 + 0.42 * draws[2]
      
      # log_alpha_t_x_N ~ RW1(sig)
      return(c(exp(log_alpha_t - log(parameters["S_u_0"]))))  
}

tparam <- list(tpar(tparam_name = "beta_t", times = c(0,20), draws2par = foi_rw2),
               tpar(tparam_name = "alpha_t", times = c(0,20), draws2par = alpha_rw1))

dynamics <-
      stem_dynamics(
            rates = rates,
            tmax = tmax,
            parameters = parameters,
            state_initializer = state_initializer,
            compartments = compartments,
            strata = strata,
            adjacency = adjacency,
            constants = constants,
            tcovar = tcovar,
            tparam = tparam,
            forcings = forcings,
            messages = F,
            compile_ode = T,
            compile_rates = T,
            compile_lna = T,
            rtol = 1e-6,
            atol = 1e-6,
            step_size = 1e-6,
            stepper = "rk54_a"
      )

emissions <- list(emission("cases", "negbinomial", c("phi", "(S_u2I_u + S_v2I_v) * rho"),
                           incidence = TRUE, obstimes = seq(1,tmax, by =1)))

measurement_process <- stem_measure(emissions = emissions, dynamics = dynamics, data = dat, messages = T)

stem_object <- stem(dynamics = dynamics, measurement_process = measurement_process)

#### initialize the inference
## Parameters
# alpha:    Baseline FOI from outside the population
# beta_0:   Initial reproduction number
# mu:       Recovery rate for unvaccinated individuals
# VE_susc:  Vaccine efficacy for susceptibility (multiplicative change in force of infection for vaccinated individuals).
# rho:      Negative binomial case detection probability.
# phi:      Negative binomial overdispersion.

## Converting to and from the estimation scale
log_popsize <- log(popsize)
to_estimation_scale <- function(params_nat) {
      return(c(
            log(expm1(log(7)-log(params_nat[1]))), #-------------------------------# mu       -> log(1/mu)
            logit(params_nat[2]), #------------------------------------------------# VE_susc  -> log(VE_susc)
            logit(params_nat[3]), #------------------------------------------------# rho      -> logit(rho)
            -0.5 * log(params_nat[4]) #--------------------------------------------# phi      -> log(1/sqrt(phi))
      ))
}

from_estimation_scale <- function(params_est) {
      return(c(
            exp(log(7) - log1p(exp(params_est[1]))), #-----------------------------# log(1/mu)          -> mu
            expit(params_est[2]), #------------------------------------------------# logit(VE_susc)     -> VE_susc
            expit(params_est[3]), #------------------------------------------------# logit(rho)         -> rho
            exp(-2 * params_est[4]) #----------------------------------------------# log(phi)           -> phi
      ))
}

## Priors 
priors <- list(prior_density =
                     function(params_nat, params_est) {
                                 dnorm(params_est[1], 0.41, 0.175, log = TRUE) + 
                                 dnorm(params_est[2], -1.1, 0.5, log = TRUE) + 
                                 dnorm(params_est[3], logit(0.0125), 0.42, log = TRUE) + 
                                 dexp(exp(params_est[4]), 1, log = TRUE) + params_est[4]
                     },
               to_estimation_scale   = to_estimation_scale,
               from_estimation_scale = from_estimation_scale)

covmat <- diag(0.01, length(parameters));
rownames(covmat) <-
      colnames(covmat) <- 
      c(
            "log_mu_inv_m1",
            "logit_VE_susc",
            "logit_rho",
            "log_sqrt_phi_inv"
      )

mcmc_kernel <-
      kernel(
            method = "mvnss",
            sigma = covmat,
            scale_constant = 0.5,
            scale_cooling = 0.99,
            stop_adaptation = 2.5e4,
            step_size = 0.025,
            harss_warmup = 0,
            nugget = 1e-5,
            mvnss_setting_list = 
                  mvnss_settings(n_mvnss_updates = 1,
                                 initial_bracket_width = 0.5,
                                 bracket_limits = c(0.125, 6),
                                 nugget_cooling = 0.99, 
                                 nugget_step_size = 0.01),
            messages = FALSE
      )

stem_object$dynamics$parameters <- function() {
      setNames(from_estimation_scale(to_estimation_scale(parameters) + rnorm(length(parameters), 0, 0.01)),
               names(parameters))
}

registerDoParallel(5)

results <- foreach(chain = 1:5,
                   .packages=c("stemr"),
                   .options.RNG = 52787,
                   .export = ls(all.names = T)) %dorng% {
                         
                         chain_res <- stem_inference(stem_object = stem_object,
                                                     method = "lna",
                                                     iterations = 7.5e4,
                                                     thin_params = 25,
                                                     thin_latent_proc = 25,
                                                     initialization_attempts = 500,
                                                     priors = priors,
                                                     mcmc_kernel = mcmc_kernel,
                                                     print_progress = 1e3,
                                                     t0_kernel = NULL,
                                                     messages = FALSE,
                                                     ess_args = ess_settings(n_ess_updates = 1,
                                                                             n_tparam_updates = 1,
                                                                             ess_warmup = 500, 
                                                                             initdist_bracket_width = 2*pi,
                                                                             initdist_bracket_update = 5e3,
                                                                             joint_initdist_update = FALSE,
                                                                             tparam_bracket_width = 2*pi,
                                                                             tparam_bracket_update = 5e3,
                                                                             joint_tparam_update = FALSE,
                                                                             lna_bracket_width = 2*pi,
                                                                             lna_bracket_update = 5e3))
                         return(chain_res)
                   }

save(results, file = paste0("flu_1seas_nostrat_const_SIR_lna.Rdata"))
