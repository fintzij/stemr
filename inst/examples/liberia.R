library(stemr)
library(extraDistr)
library(foreach)
library(doParallel)
library(doRNG)


# Data and population size ------------------------------------------------
stopImplicitCluster()
rm(list=ls())
popsize <- 4.4e6
log_popsize <- log(popsize)

dat <- 
      cbind(
            time = 
                  1:64,
            cases = 
                  c(0,0,1,1,3,1,1,0,1,4,0,1,0,3,18,28,27,26,38,36,86,106,208,289,
                    394,362,431,454,480,404,291,272,192,118,122,130,100,82,74,66,
                    34,15,25,9,4,14,7,12,4,11,4,1,1,0,2,1,0,0,0,0,0,0,0,0)
      )

# no strata the stemr object --------------------------------------------------
set.seed(12511)
strata <- NULL
compartments <- c("S", "E", "I", "R")
rates <- list(rate("(alpha + beta * I) * (S - effpop)", "S", "E", lumped = TRUE, incidence = T),
              rate("omega", "E", "I", incidence = T),
              rate("mu", "I", "R", incidence = T))
state_initializer <- list(stem_initializer(c(S = popsize-30, E = 15, I = 10, R = 5), fixed = F, prior = c(popsize-30, 15, 10, 5)))
adjacency <- NULL
tcovar <- NULL
parameters = c(alpha = 0.05 / (popsize - 4.35e6), beta = 2.25 / (popsize - 4.35e6), omega = 1, mu = 1, rho = 0.75, phi = 5, effpop = 4.35e6)
constants <- c(t0 = 0)
t0 <- 0; tmax <- nrow(dat);

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
            tcovar = tcovar,
            messages = T,
            compile_ode = T,
            compile_rates = T,
            compile_lna = T,
            rtol = 1e-6,
            atol = 1e-6,
            step_size = 1e-6
      )

emissions <- list(emission("cases", "negbinomial", c("phi","E2I * rho"), incidence = TRUE, obstimes = seq(1, tmax, by =1)))

measurement_process <- stem_measure(data = dat, emissions = emissions, dynamics = dynamics, messages = T)

stem_object <- stem(dynamics = dynamics, measurement_process = measurement_process)

#### initialize the inference
# {alpha, beta, omega, mu, rho, phi, Neff} 
# -> {log(Rext), log(Reff-1) + log(rho * Neff), log(omega / mu), log(1/mu), logit(rho), log(phi), log(rho * Neff)}
to_estimation_scale <- 
      function (params_nat) {
            
            l_effpop     <- log(popsize - params_nat[7])
            l_Neff_x_rho <- l_effpop + log(params_nat[5]) 
            l_infecdur   <- -log(params_nat[4])
            l_Reff_m1    <- log(exp(log(params_nat[2]) + l_effpop + l_infecdur)-1)

            return(
                  c(
                        log(params_nat[1]) + log(1000),
                        l_Reff_m1 + l_Neff_x_rho,
                        log(params_nat[3]) + l_infecdur, 
                        l_infecdur,
                        logit(params_nat[5]),
                        log(params_nat[6]),
                        l_Neff_x_rho
                  )
            )
      }

from_estimation_scale <- 
      function (params_est) {
            
            rho <- expit(params_est[5])
            l_effpop <- params_est[7] - log(rho)
            
            return(c(
                  exp(params_est[1] - log(1000)),
                  exp(log(exp(params_est[2] - params_est[7])+1) - l_effpop - params_est[4]),
                  exp(params_est[3] - params_est[4]),
                  exp(-params_est[4]),
                  rho,
                  exp(params_est[6]),
                  popsize - exp(l_effpop)
            ))
      }

## Priors
priors <- list(prior_density =
                  function(params_nat, params_est) {
                     
                     l_effpop <- params_est[7] - log(expit(params_est[5]))
                     l_Reff_m_1 <- params_est[2] - params_est[7]
                     
                     sum(
                        dexp(exp(params_est[1]), 40, log = TRUE) + params_est[1], 
                        dnorm(l_Reff_m_1, log(0.5), 1.08, log = TRUE),
                        dnorm(params_est[3], 0, 0.3, log = TRUE),
                        dnorm(params_est[4], 0, 0.3, log = TRUE),
                        dnorm(params_est[5], 0.85, 0.75, log = TRUE),
                        dexp(exp(params_est[6]), 0.69, log = TRUE) + params_est[6],
                        dnorm(l_effpop, 9.9, 0.62, log = TRUE)
                     )
                  },
               to_estimation_scale   = to_estimation_scale,
               from_estimation_scale = from_estimation_scale)

covmat_names <- c(
   "log_Reff_ext",
   "log_Reff_m_1_o",
   "log_omega_d_mu",
   "log_carriage_dur",
   "logit_rho",
   "log_phi",
   "log_effpop_o"
)
covmat <- diag(0.01, length(parameters))
rownames(covmat) <- colnames(covmat) <- covmat_names

mcmc_kernel <-
   kernel(
       method = "mvnss",
       sigma = covmat,
       scale_constant = 0.5,
       scale_cooling = 0.7,
       stop_adaptation = 5e4,
       step_size = 0.5,
       nugget = 1e-5, 
       harss_warmup = 0,
       mvnss_setting_list = 
           mvnss_settings(n_mvnss_updates = 1, 
                          initial_bracket_width = 0.5,
                          bracket_limits = c(0.001, 5),
                          nugget_cooling = 0.99, 
                          nugget_step_size = 0.001),
      messages = FALSE
   )

stem_object$dynamics$parameters <- function() {
   setNames(from_estimation_scale(to_estimation_scale(parameters) + rnorm(length(parameters), 0, 0.1)),
            names(parameters))
}

# registerDoParallel(3)
# 
# results <- foreach(chain = 1:5,
#                    .packages=c("stemr"),
#                    .options.RNG = 52787,
#                    .export = ls(all.names = T)) %dorng% {
#                       
#                       chain_res <- stem_inference(stem_object = stem_object,
#                                                   method = "ode",
#                                                   iterations = 1.5e5,
#                                                   thin_params = 50,
#                                                   thin_latent_proc = 50,
#                                                   initialization_attempts = 500,
#                                                   priors = priors,
#                                                   mcmc_kernel = mcmc_kernel,
#                                                   t0_kernel = t0_kernel,
#                                                   ess_args = ess_settings(n_ess_updates = 1,
#                                                                           ess_warmup = 100, 
#                                                                           initdist_bracket_width = 2*pi,
#                                                                           initdist_bracket_update = 5e3,
#                                                                           lna_bracket_width = 2*pi,
#                                                                           lna_bracket_update = 5e3,
#                                                                           joint_strata_update = FALSE),
#                                                   print_progress = 5e3,
#                                                   messages = F)
#                       return(chain_res)
#                    }
# 
# save(results, file = paste0("liberia_ODE.Rdata"))
# 
# # grab the initial covariance matrix
# covs <- sapply(results, function(x) cov(x$stem_settings$mcmc_kernel$sigma))
# covmat <- matrix(rowMeans(covs), length(parameters), length(parameters), dimnames = list(covmat_names, covmat_names))
# 
# mcmc_kernel <-
#    kernel(
#        method = "mvnss",
#        sigma = covmat,
#        scale_constant = 0.5,
#        scale_cooling = 0.7,
#        stop_adaptation = 5e4,
#        step_size = 0.5,
#        nugget = 1e-5, 
#        harss_warmup = 0,
#        mvnss_setting_list = 
#            mvnss_settings(n_mvnss_updates = 1, 
#                           initial_bracket_width = 0.5,
#                           bracket_limits = c(0.001, 5),
#                           nugget_cooling = 0.99, 
#                           nugget_step_size = 0.001),
#       messages = FALSE
#    )
# 
# # grab the initial parameters and compartment volumes
# init_pars <- lapply(results, function(x) x$dynamics$parameters)
# init_states <- lapply(results, function(x) x$dynamics$initdist_params)
# 
# rm(results)

registerDoParallel(3)

results <- foreach(chain = 1:5,
                   .packages="stemr",
                   .options.RNG = 52787,
                   .export = ls(all.names = T)) %dorng% {
                      
                      # stem_object$dynamics$parameters <- init_pars[[chain]]
                      # stem_object$dynamics$initdist_params <- init_states[[chain]]
                      
                      chain_res <- stem_inference(stem_object = stem_object,
                                                  method = "lna",
                                                  iterations = 1.5e5,
                                                  thin_params = 50,
                                                  thin_latent_proc = 50,
                                                  initialization_attempts = 500,
                                                  priors = priors,
                                                  mcmc_kernel = mcmc_kernel,
                                                  t0_kernel = t0_kernel,
                                                  ess_args = ess_settings(n_ess_updates = 1,
                                                                          ess_warmup = 100, 
                                                                          initdist_bracket_width = 2*pi,
                                                                          initdist_bracket_update = 5e3,
                                                                          lna_bracket_width = 2*pi,
                                                                          lna_bracket_update = 5e3,
                                                                          joint_strata_update = FALSE),    
                                                  messages = F,
                                                  print_progress = 5e3)
                      return(chain_res)
                   }

save(results, file = paste0("liberia_LNA.Rdata"))
