library(stemr)
library(extraDistr)
library(foreach)
library(doParallel)
library(doRNG)

# Simulate data -----------------------------------------------------------

popsize_guin <- 11.8e6
popsize_lib <- 4.4e6
popsize_sln <- 7.1e6

# effective population sizes
ep_guin <- 2e4
ep_lib <- 3.5e4
ep_sln <- 2.5e4

# dynamics
set.seed(52787)
strata <- c("guin", "lib", "sln")
compartments <- list(S = "ALL", E = "ALL", I = "ALL", R = "ALL")

rates <- list(rate("beta_guin * (I_guin + transmission_lib * alpha_lib2guin * I_lib + transmission_sln * alpha_sln2guin * I_sln) * S_guin", from = "S", to = "E", strata = "guin", lumped = TRUE, incidence = T),
              rate("transmission_lib * beta_lib * (I_lib + alpha_guin2lib * I_guin + transmission_sln * alpha_sln2lib * I_sln) * S_lib", from = "S", to = "E", strata = "lib", lumped = TRUE, incidence = T),
              rate("transmission_sln * beta_sln * (I_sln + alpha_guin2sln * I_guin + transmission_lib * alpha_lib2sln * I_lib) * S_sln", from = "S", to = "E", strata = "sln", lumped = TRUE, incidence = T),
              rate("omega_guin", from = "E", to = "I", strata = "guin", incidence = T),
              rate("transmission_lib * omega_lib", from = "E", to = "I", strata = "lib", incidence = T),
              rate("transmission_sln * omega_sln", from = "E", to = "I", strata = "sln", incidence = T),
              rate("mu_guin", from = "I", to = "R", strata = "guin", incidence = T),
              rate("transmission_lib * mu_lib", from = "I", to = "R", strata = "lib", incidence = T),
              rate("transmission_sln * mu_sln", from = "I", to = "R", strata = "sln", incidence = T))

state_initializer <- list(stem_initializer(c(S_guin = ep_guin - 30, E_guin = 15, I_guin = 10, R_guin = 5),
                                           fixed = TRUE, 
                                           strata = "guin",
                                           prior = c(popsize_guin-30, 15, 10, 5)),
                          stem_initializer(c(S_lib = ep_lib - 30, E_lib = 15, I_lib = 10, R_lib = 5), 
                                           fixed = TRUE,
                                           strata = "lib",
                                           prior = c(popsize_lib - 30, 15, 10, 5)),
                          stem_initializer(c(S_sln = ep_sln - 30, E_sln = 15, I_sln = 10, R_sln = 5),
                                           fixed = TRUE,
                                           strata = "sln",
                                           prior = c(popsize_sln - 30, 15, 10, 5)))
constants <- c(t0 = 0)
t0 <- 0; tmax <- 70;
adjacency <- NULL
tcovar <- cbind(time = 0:tmax,
                transmission_lib = c(rep(0, 10), rep(1, (tmax + 1) - 10)),
                transmission_sln = c(rep(0, 19), rep(1, (tmax + 1) - 19)))

mu_guin <- 0.9
mu_lib <- 1.1
mu_sln <- 1

parameters = c(beta_guin = 1.2 / ep_guin * mu_guin, 
               beta_lib = 1.35 / ep_lib * mu_lib,
               beta_sln = 1.45 / ep_sln * mu_sln,
               alpha_guin2lib = 0.02 / ep_lib * mu_guin,
               alpha_guin2sln = 0.02 / ep_sln * mu_guin,
               alpha_lib2guin = 0.02 / ep_guin * mu_lib,
               alpha_lib2sln  = 0.02 / ep_sln * mu_lib,
               alpha_sln2guin = 0.02 / ep_guin * mu_sln,
               alpha_sln2lib  = 0.02 / ep_lib * mu_sln,
               omega_guin = 1.2,
               omega_lib = 1,
               omega_sln = 0.8,
               mu_guin = mu_guin,
               mu_lib = mu_lib,
               mu_sln = mu_sln,
               rho_guin = 100/150,
               rho_lib = 100/175,
               rho_sln = 100/125,
               phi_guin = 100,
               phi_lib = 100,
               phi_sln = 100)
true_pars <- parameters

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
            rtol = 1e-5,
            atol = 1e-5,
            step_size = 1e-7
      )

emissions <- list(emission("guin_cases", "negbinomial", c("phi_guin", "E_guin2I_guin * rho_guin"), incidence = TRUE, obstimes = seq(1, tmax, by =1)),
                  emission("lib_cases", "negbinomial", c("phi_lib", "E_lib2I_lib * rho_lib"), incidence = TRUE, obstimes = seq(10, tmax, by =1)),
                  emission("sln_cases", "negbinomial", c("phi_sln", "E_sln2I_sln * rho_sln"), incidence = TRUE, obstimes = seq(19, tmax, by =1)))

measurement_process <- stem_measure(emissions = emissions, dynamics = dynamics, messages = T)

stem_object <- stem(dynamics = dynamics, measurement_process = measurement_process)

sim <- simulate_stem(stem_object)

true_path   <- sim$paths[[1]]
dat         <- sim$datasets[[1]]
total_cases <- colSums(dat[,-1])

# Data and population size ------------------------------------------------

log_popsize_guin <- log(popsize_guin)
log_popsize_lib <- log(popsize_lib)
log_popsize_sln <- log(popsize_sln)

# stemr object separate dynamics for all countries --------------------------------------------------
set.seed(12511)
strata <- c("guin", "lib", "sln")
compartments <- list(S = "ALL", E = "ALL", I = "ALL", R = "ALL")
rates <- list(rate("beta_guin * (I_guin + transmission_lib * alpha_lib2guin * I_lib + transmission_sln * alpha_sln2guin * I_sln) * (S_guin - effpop_guin)", from = "S", to = "E", strata = "guin", lumped = TRUE, incidence = T),
              rate("transmission_lib * beta_lib * (I_lib + alpha_guin2lib * I_guin + transmission_sln * alpha_sln2lib * I_sln) * (S_lib - effpop_lib)", from = "S", to = "E", strata = "lib", lumped = TRUE, incidence = T),
              rate("transmission_sln * beta_sln * (I_sln + alpha_guin2sln * I_guin + transmission_lib * alpha_lib2sln * I_lib) * (S_sln - effpop_sln)", from = "S", to = "E", strata = "sln", lumped = TRUE, incidence = T),
              rate("omega_guin", from = "E", to = "I", strata = "guin", incidence = T),
              rate("transmission_lib * omega_lib", from = "E", to = "I", strata = "lib", incidence = T),
              rate("transmission_sln * omega_sln", from = "E", to = "I", strata = "sln", incidence = T),
              rate("mu_guin", from = "I", to = "R", strata = "guin", incidence = T),
              rate("transmission_lib * mu_lib", from = "I", to = "R", strata = "lib", incidence = T),
              rate("transmission_sln * mu_sln", from = "I", to = "R", strata = "sln", incidence = T))

state_initializer <- list(stem_initializer(c(S_guin = popsize_guin - 30, E_guin = 15, I_guin = 10, R_guin = 5),
                                           fixed = FALSE, 
                                           strata = "guin",
                                           prior = c(popsize_guin-30, 15, 10, 5)),
                          stem_initializer(c(S_lib = popsize_lib - 30, E_lib = 15, I_lib = 10, R_lib = 5), 
                                           fixed = FALSE,
                                           strata = "lib",
                                           prior = c(popsize_lib - 30, 15, 10, 5)),
                          stem_initializer(c(S_sln = popsize_sln - 30, E_sln = 15, I_sln = 10, R_sln = 5),
                                           fixed = FALSE,
                                           strata = "sln",
                                           prior = c(popsize_sln - 30, 15, 10, 5)))

adjacency <- NULL
tcovar <- cbind(time = 0:tmax,
                transmission_lib = c(rep(0, 10), rep(1, (tmax+1) - 10)),
                transmission_sln = c(rep(0, 19), rep(1, (tmax+1) - 19)))

parameters = c(beta_guin = 1.25 / ep_guin * mu_guin, 
               beta_lib = 1.35 / ep_lib * mu_lib,
               beta_sln = 1.45 / ep_sln * mu_sln,
               alpha_guin2lib = 0.02 / ep_lib * mu_guin,
               alpha_guin2sln = 0.02 / ep_sln * mu_guin,
               alpha_lib2guin = 0.02 / ep_guin * mu_lib,
               alpha_lib2sln  = 0.02 / ep_sln * mu_lib,
               alpha_sln2guin = 0.02 / ep_guin * mu_sln,
               alpha_sln2lib  = 0.02 / ep_lib * mu_sln,
               effpop_guin = popsize_guin - ep_guin,
               effpop_lib = popsize_lib - ep_lib,
               effpop_sln = popsize_sln - ep_sln,
               omega_guin = 1,
               omega_lib = 1,
               omega_sln = 1,
               mu_guin = mu_guin,
               mu_lib = mu_lib,
               mu_sln = mu_sln,
               rho_guin = 100/150,
               rho_lib = 100/175,
               rho_sln = 100/125,
               phi_guin = 50,
               phi_lib = 50,
               phi_sln = 50)

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
            rtol = 1e-5,
            atol = 1e-5,
            step_size = 1e-7
      )

emissions <- list(emission("guin_cases", "negbinomial", c("phi_guin", "E_guin2I_guin * rho_guin"), incidence = TRUE, obstimes = seq(1, tmax, by =1)),
                  emission("lib_cases", "negbinomial", c("phi_lib", "E_lib2I_lib * rho_lib"), incidence = TRUE, obstimes = seq(10, tmax, by =1)),
                  emission("sln_cases", "negbinomial", c("phi_sln", "E_sln2I_sln * rho_sln"), incidence = TRUE, obstimes = seq(19, tmax, by =1)))

# recompile stemr object
measurement_process <- stem_measure(data = dat, emissions = emissions, dynamics = dynamics, messages = T)

stem_object <- stem(dynamics = dynamics, measurement_process = measurement_process)

#### initialize the inference
popsizes <- c(popsize_guin, popsize_lib, popsize_sln)

# to and from the MCMC estimation scale
to_estimation_scale <- function(params_nat) {
      
      l_effpops      <- log(popsizes - params_nat[10:12])
      l_Neffs_x_rhos <- l_effpops + log(params_nat[19:21])
      l_infecdurs    <- -log(params_nat[16:18])
      l_Reffs_m1     <- log(expm1(log(params_nat[1:3]) + l_effpops + l_infecdurs))
      
      return(c(
            l_Reffs_m1 + l_Neffs_x_rhos,
            log(params_nat[4:9]) + l_effpops[c(2,3,1,3,1,2)] + l_infecdurs[c(1,1,2,2,3,3)],
            l_Neffs_x_rhos,
            log(params_nat[13:15]) + l_infecdurs,
            l_infecdurs,
            logit(params_nat[19:21]),
            -0.5 * log(params_nat[22:24])
      ))
}

from_estimation_scale <- function(params_est) {
      
      rhos       <- expit(params_est[19:21])
      l_effpops  <- params_est[10:12] - log(rhos)
      l_Reffs_m1 <- params_est[1:3] - params_est[10:12]
      
      return(c(
            exp(log1p(exp(l_Reffs_m1)) - l_effpops - params_est[16:18]),
            exp(params_est[4:9] - l_effpops[c(2,3,1,3,1,2)] - params_est[c(16, 16, 17, 17, 18, 18)]),                               
            popsizes - exp(l_effpops), 
            exp(params_est[13:15] - params_est[16:18]),
            exp(-params_est[16:18]),                                                   
            expit(params_est[19:21]),                                     
            exp(-2 * params_est[22:24])                                        
      ))
}

priors <- list(prior_density =
                     function(params_nat, params_est) {
                           
                           rhos       <- expit(params_est[19:21])
                           l_effpops  <- params_est[10:12] - log(rhos)
                           
                           sum(
                                 dnorm(params_est[1:3] - params_est[10:12], log(0.5), 1.08, log = TRUE), 
                                 dexp(exp(params_est[4:9]), 40, log = TRUE) + params_est[4:9],
                                 dnorm(l_effpops, c(9.8, 10.5, 10.6), c(0.62, 0.62, 0.62), log = TRUE),
                                 dnorm(params_est[13:15], 0, 0.3, log = TRUE),
                                 dnorm(params_est[16:18], 0, 0.3, log = TRUE), 
                                 dnorm(params_est[19:21], 0.85, 0.75, log = TRUE),
                                 dexp(exp(params_est[22:24]), log = TRUE) + params_est[22:24])
                     },
               to_estimation_scale   = to_estimation_scale,
               from_estimation_scale = from_estimation_scale)

covmat_names <- c(
      "log_Reff_guin_o",
      "log_Reff_lib_o",
      "log_Reff_sln_o",
      "log_Rext_g2l",
      "log_Rext_g2s",
      "log_Rext_l2g",
      "log_Rext_l2s",
      "log_Rext_s2g",
      "log_Rext_s2l",
      "log_Neff_x_rho_guin",
      "log_Neff_x_rho_lib",
      "log_Neff_x_rho_sln",
      "log_omega_d_mu_guin",
      "log_omega_d_mu_lib",
      "log_omega_d_mu_sln",
      "log_infecdur_guin",
      "log_infecdur_lib",
      "log_infecdur_sln",
      "logit_rho_guin",
      "logit_rho_lib",
      "logit_rho_sln",
      "log_sqrt_phi_inv_guin",
      "log_sqrt_phi_inv_lib",
      "log_sqrt_phi_inv_sln"
)

covmat <- diag(c(rep(0.05,3),  # l_Reffs
                 rep(0.1, 6),  # Rext
                 rep(1, 3),    # log_effpop_x_rho
                 rep(0.05, 3), # log_omega_d_mu
                 rep(0.05, 3), # log_carriage_dur
                 rep(0.5, 3),  # logit_rho
                 rep(0.05, 3))) # log_phi
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
            parameter_blocks = 
                  list(parblock(c("log_Reff_guin_o",
                                  "log_Rext_l2g", 
                                  "log_Rext_s2g",
                                  "log_Neff_x_rho_guin",
                                  "log_omega_d_mu_guin",
                                  "log_infecdur_guin",
                                  "logit_rho_guin",
                                  "log_sqrt_phi_inv_guin")),
                       parblock(c("log_Reff_lib_o", 
                                  "log_Rext_g2l", 
                                  "log_Rext_s2l",
                                  "log_Neff_x_rho_lib",
                                  "log_omega_d_mu_lib",
                                  "log_infecdur_lib",
                                  "logit_rho_lib",
                                  "log_sqrt_phi_inv_lib")),
                       parblock(c("log_Reff_sln_o", 
                                  "log_Rext_g2s", 
                                  "log_Rext_l2s",
                                  "log_Neff_x_rho_sln",
                                  "log_omega_d_mu_sln",
                                  "log_infecdur_sln",
                                  "logit_rho_sln",
                                  "log_sqrt_phi_inv_sln"))),
            joint_block_update = FALSE,
            messages = FALSE
      )

stem_object$dynamics$parameters <- function() {
      setNames(from_estimation_scale(to_estimation_scale(parameters) + rnorm(length(parameters), 0, 0.1)),
               names(parameters))
}

registerDoParallel(5)

results <- foreach(chain = 1:5,
                   .packages=c("stemr"),
                   .options.RNG = 52787,
                   .export = ls(all.names = T)) %dorng% {
                         
                         chain_res <- stem_inference(stem_object = stem_object,
                                                     method = "ode",
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
                                                                             joint_strata_update = FALSE,
                                                                             joint_initdist_update = FALSE),
                                                     print_progress = 1e3,
                                                     messages = F)
                         return(chain_res)
                   }

save(results, file = paste0("ebola_synth_ODE.Rdata"))

# grab the initial covariance matrix
covs <- sapply(results, function(x) x$stem_settings$mcmc_kernel$sigma)
covmat <- matrix(rowMeans(covs), length(parameters), length(parameters), dimnames = list(covmat_names, covmat_names))

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
            parameter_blocks = 
                  list(parblock(c("log_Reff_guin_o",
                                  "log_Rext_l2g", 
                                  "log_Rext_s2g",
                                  "log_Neff_x_rho_guin",
                                  "log_omega_d_mu_guin",
                                  "log_infecdur_guin",
                                  "logit_rho_guin",
                                  "log_sqrt_phi_inv_guin")),
                       parblock(c("log_Reff_lib_o", 
                                  "log_Rext_g2l", 
                                  "log_Rext_s2l",
                                  "log_Neff_x_rho_lib",
                                  "log_omega_d_mu_lib",
                                  "log_infecdur_lib",
                                  "logit_rho_lib",
                                  "log_sqrt_phi_inv_lib")),
                       parblock(c("log_Reff_sln_o", 
                                  "log_Rext_g2s", 
                                  "log_Rext_l2s",
                                  "log_Neff_x_rho_sln",
                                  "log_omega_d_mu_sln",
                                  "log_infecdur_sln",
                                  "logit_rho_sln",
                                  "log_sqrt_phi_inv_sln"))),
            joint_block_update = FALSE,
            messages = FALSE
      )

# grab the initial parameters and compartment volumes
init_pars <- lapply(results, 
                    function(x) setNames(x$dynamics$parameters, names(parameters)))
init_states <- lapply(results, 
                      function(x) setNames(x$dynamics$initdist_params, names(stem_object$dynamics$initdist_params)))

rm(results)

registerDoParallel(5)

results <- foreach(chain = 1:5,
                   .packages="stemr",
                   .options.RNG = 52787,
                   .export = ls(all.names = T)) %dorng% {
                         
                         stem_object$dynamics$parameters <- init_pars[[chain]]
                         stem_object$dynamics$initdist_params <- init_states[[chain]]
                         
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
                                                                             initdist_bracket_width = pi/6,
                                                                             initdist_bracket_update = 5e3,
                                                                             lna_bracket_width = pi/6,
                                                                             lna_bracket_update = 5e3,
                                                                             joint_strata_update = FALSE,
                                                                             joint_initdist_update = FALSE),
                                                     print_progress = 5e3,
                                                     messages = F)
                         return(chain_res)
                   }

results$true_path <- true_path
results$true_pars <- true_pars

save(results, file = paste0("ebola_synth_LNA.Rdata"))
