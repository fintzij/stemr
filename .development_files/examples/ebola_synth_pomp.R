
# Load libraries ----------------------------------------------------------

library(pomp)
library(stemr)
library(foreach)
library(doParallel)
library(doRNG)

# Load data ------------------------------------------------

popsize_guin <- 11.8e6
popsize_lib <- 4.4e6
popsize_sln <- 7.1e6

load("ebola_synth_dat.Rdata")
dat <- as.data.frame(dat)

# set zero measurements to NA
dat[1:10,3] = NA
dat[1:19,4] = NA

# time-varying covariates for transmission
tcovar_pomp <-
      data.frame(time = seq(0, max(dat[,1]), by = 1/7))
tcovar_pomp$transmission_lib <- ifelse(tcovar_pomp$time < 10, 0.0, 1.0)
tcovar_pomp$transmission_sln <- ifelse(tcovar_pomp$time < 19, 0.0, 1.0)

# Set up pomp objects -----------------------------------------------------
# 
# # Measurement process objects
rmeas <- "
guin_cases = rnbinom_mu(exp(-2 * log_sqrt_phi_inv_guin), E2I_guin / (1 + exp(-logit_rho_guin)));
lib_cases  = rnbinom_mu(exp(-2 * log_sqrt_phi_inv_lib), E2I_lib / (1 + exp(-logit_rho_lib)));
sln_cases  = rnbinom_mu(exp(-2 * log_sqrt_phi_inv_sln), E2I_sln / (1 + exp(-logit_rho_sln)));
"

dmeas <- "
      double loglik_guin;
      double loglik_lib;
      double loglik_sln;

      loglik_guin = dnbinom_mu(guin_cases, exp(-2 * log_sqrt_phi_inv_guin), E2I_guin / (1 + exp(-logit_rho_guin)), 1);

      if(ISNA(lib_cases)) {
            loglik_lib = 0.0;
      } else {
            loglik_lib = dnbinom_mu(lib_cases, exp(-2 * log_sqrt_phi_inv_lib), E2I_lib / (1 + exp(-logit_rho_lib)), 1);
      }

      if(ISNA(sln_cases)) {
            loglik_sln = 0.0;
      } else {
            loglik_sln = dnbinom_mu(sln_cases, exp(-2 * log_sqrt_phi_inv_sln), E2I_sln / (1 + exp(-logit_rho_sln)), 1);
      }

      lik = loglik_guin + loglik_lib + loglik_sln;
      lik = (give_log) ? lik : exp(lik);
"

# define the stepper
ebola.step<-"
      // per contact infection
      double beta_guin;
      double beta_lib;
      double beta_sln;

      // cross country transmission
      double alpha_guin2lib;
      double alpha_guin2sln;
      double alpha_lib2guin;
      double alpha_lib2sln;
      double alpha_sln2guin;
      double alpha_sln2lib;

      // log effective population sizes
      double log_effpop_guin;
      double log_effpop_lib;
      double log_effpop_sln;

      // effective population sizes
      double effpop_guin;
      double effpop_lib;
      double effpop_sln;
      
      // effective numbers of susceptibles
      double S_eff_guin;
      double S_eff_lib;
      double S_eff_sln;

      // latent to infectious
      double omega_guin;
      double omega_lib;
      double omega_sln;

      // infectious to recovered
      double mu_guin;
      double mu_lib;
      double mu_sln;

      // detection rates
      double rho_guin;
      double rho_lib;
      double rho_sln;

      // rates and state transitions
      double rates[9];
      double dN[9];

      // get the parameters on their natural scales
      rho_guin = 1 / (1 + exp(-logit_rho_guin));
      rho_lib  = 1 / (1 + exp(-logit_rho_lib));
      rho_sln  = 1 / (1 + exp(-logit_rho_sln));

      log_effpop_guin = log_Neff_x_rho_guin - log(rho_guin);
      log_effpop_lib  = log_Neff_x_rho_lib - log(rho_lib);
      log_effpop_sln  = log_Neff_x_rho_sln - log(rho_sln);

      effpop_guin = 11800000 - exp(log_effpop_guin);
      effpop_lib  = 4400000 - exp(log_effpop_lib);
      effpop_sln  = 7100000 - exp(log_effpop_sln);

      beta_guin = exp(log(1 + exp(log_Reff_m1_o_guin - log_Neff_x_rho_guin)) - log_effpop_guin - log_mu_inv_guin);
      beta_lib  = exp(log(1 + exp(log_Reff_m1_o_lib - log_Neff_x_rho_lib)) - log_effpop_lib - log_mu_inv_lib);
      beta_sln  = exp(log(1 + exp(log_Reff_m1_o_sln - log_Neff_x_rho_sln)) - log_effpop_sln - log_mu_inv_sln);

      alpha_guin2lib = exp(log_Rext_g2l - log_effpop_lib  - log_mu_inv_guin);
      alpha_guin2sln = exp(log_Rext_g2s - log_effpop_sln  - log_mu_inv_guin);
      alpha_lib2guin = exp(log_Rext_l2g - log_effpop_guin - log_mu_inv_lib);
      alpha_lib2sln  = exp(log_Rext_l2s - log_effpop_sln  - log_mu_inv_lib);
      alpha_sln2guin = exp(log_Rext_s2g - log_effpop_guin - log_mu_inv_sln);
      alpha_sln2lib  = exp(log_Rext_s2l - log_effpop_lib  - log_mu_inv_sln);

      omega_guin = exp(log_omega_d_mu_guin - log_mu_inv_guin);
      omega_lib  = exp(log_omega_d_mu_lib - log_mu_inv_lib);
      omega_sln  = exp(log_omega_d_mu_sln - log_mu_inv_sln);

      mu_guin = exp(-log_mu_inv_guin);
      mu_lib  = exp(-log_mu_inv_lib);
      mu_sln  = exp(-log_mu_inv_sln);
      
      // compute effective number of susceptibles
      S_eff_guin = S_guin - effpop_guin;
      S_eff_lib = S_lib - effpop_lib;
      S_eff_sln = S_sln - effpop_sln;
      
      if(S_eff_guin < 0.0) S_eff_guin = 0.0;
      if(S_eff_lib < 0.0) S_eff_lib = 0.0;
      if(S_eff_sln < 0.0) S_eff_sln = 0.0;

      // compute rates
      rates[0] = beta_guin * (I_guin + transmission_lib * alpha_lib2guin * I_lib + transmission_sln * alpha_sln2guin * I_sln) / S_guin * S_eff_guin;
      rates[1] = transmission_lib * beta_lib * (I_lib + alpha_guin2lib * I_guin + transmission_sln * alpha_sln2lib * I_sln) / S_lib * S_eff_lib;
      rates[2] = transmission_sln * beta_sln * (I_sln + alpha_guin2sln * I_guin + transmission_lib * alpha_lib2sln * I_lib) / S_sln * S_eff_sln;
      rates[3] = omega_guin;
      rates[4] = transmission_lib * omega_lib;
      rates[5] = transmission_sln * omega_sln;
      rates[6] = mu_guin;
      rates[7] = transmission_lib * mu_lib;
      rates[8] = transmission_sln * mu_sln;

      // generate the incidence increments
      reulermultinom(1, S_guin, &rates[0], dt, &dN[0]);
      reulermultinom(1, S_lib, &rates[1], dt, &dN[1]);
      reulermultinom(1, S_sln, &rates[2], dt, &dN[2]);
      reulermultinom(1, E_guin, &rates[3], dt, &dN[3]);
      reulermultinom(1, E_lib,  &rates[4], dt, &dN[4]);
      reulermultinom(1, E_sln,  &rates[5], dt, &dN[5]);
      reulermultinom(1, I_guin, &rates[6], dt, &dN[6]);
      reulermultinom(1, I_lib,  &rates[7], dt, &dN[7]);
      reulermultinom(1, I_sln,  &rates[8], dt, &dN[8]);

      // increment the compartment counts
      S_guin += -dN[0];
      S_lib  += -dN[1];
      S_sln  += -dN[2];

      E_guin += dN[0] - dN[3];
      E_lib  += dN[1] - dN[4];
      E_sln  += dN[2] - dN[5];

      I_guin += dN[3] - dN[6];
      I_lib  += dN[4] - dN[7];
      I_sln  += dN[5] - dN[8];

      R_guin += dN[6];
      R_lib  += dN[7];
      R_sln  += dN[8];

      // increment the incidence
      E2I_guin += dN[3];
      E2I_lib  += dN[4];
      E2I_sln  += dN[5];
"

# instatiate the euler stepper function
ebola_sim <- euler.sim(step.fun = Csnippet(ebola.step), delta.t = 1/7)

# Define the priors
prior_density = function(params,..., log) {
      lik <- sum(
            # endogenous effective reproduction numbers - 1
            dnorm(params[c("log_Reff_m1_o_guin", "log_Reff_m1_o_lib", "log_Reff_m1_o_sln")] -
                        params[c("log_Neff_x_rho_guin", "log_Neff_x_rho_lib", "log_Neff_x_rho_sln")], log(0.5), 1.08, log = TRUE) +

            # exogenous effective reproduction numbers
            dexp(exp(params[c("log_Rext_g2l", "log_Rext_g2s", "log_Rext_l2g", "log_Rext_l2s", "log_Rext_s2g", "log_Rext_s2l")]), 40, log = TRUE) +
            params[c("log_Rext_g2l", "log_Rext_g2s", "log_Rext_l2g", "log_Rext_l2s", "log_Rext_s2g", "log_Rext_s2l")] +

            # effective population sizes
            dnorm(params[c("log_Neff_x_rho_guin", "log_Neff_x_rho_lib", "log_Neff_x_rho_sln")] -
                        log(expit(params[c("logit_rho_guin", "logit_rho_lib", "logit_rho_sln")])), c(9.8, 10.5, 10.6), c(0.62, 0.62, 0.62), log = TRUE) +

            # ratios of latent and infectious period durations
            dnorm(params[c("log_omega_d_mu_guin", "log_omega_d_mu_lib", "log_omega_d_mu_sln",
                           "log_mu_inv_guin", "log_mu_inv_lib", "log_mu_inv_sln")], 0, 0.3, log = TRUE) +

            # mean case detection rates
            dnorm(params[c("logit_rho_guin", "logit_rho_lib", "logit_rho_sln")], 0.85, 0.75, log = TRUE) +

            # negative binomial overdispersions
            dexp(exp(params[c("log_sqrt_phi_inv_guin", "log_sqrt_phi_inv_lib", "log_sqrt_phi_inv_sln")]), log = TRUE))

      if(!log) lik <- exp(lik)
      return(lik)
}

ebola_mod <- pomp(
      data = dat,
      times = "time",
      covar = tcovar_pomp,
      tcovar = "time",
      t0 = 0,                      # initial time point
      dmeasure = Csnippet(dmeas),  # evaluates the density of the measurement process
      rmeasure = Csnippet(rmeas),  # simulates from the measurement process
      rprocess = ebola_sim,        # simulates from the latent process
      obsnames = c("guin_cases", "lib_cases", "sln_cases"),
      statenames = c("S_guin", "E_guin", "I_guin", "R_guin",
                     "S_lib", "E_lib", "I_lib", "R_lib",
                     "S_sln", "E_sln", "I_sln", "R_sln",
                     "E2I_guin", "E2I_lib", "E2I_sln"),
      paramnames = c("log_Reff_m1_o_guin", "log_Reff_m1_o_lib", "log_Reff_m1_o_sln",
                     "log_Rext_g2l", "log_Rext_g2s", "log_Rext_l2g", "log_Rext_l2s", "log_Rext_s2g", "log_Rext_s2l",
                     "log_Neff_x_rho_guin", "log_Neff_x_rho_lib", "log_Neff_x_rho_sln",
                     "log_omega_d_mu_guin", "log_omega_d_mu_lib", "log_omega_d_mu_sln",
                     "log_mu_inv_guin", "log_mu_inv_lib", "log_mu_inv_sln",
                     "logit_rho_guin", "logit_rho_lib", "logit_rho_sln",
                     "log_sqrt_phi_inv_guin", "log_sqrt_phi_inv_lib", "log_sqrt_phi_inv_sln"),
      covarnames = c("transmission_lib", "transmission_sln"),
      zeronames = c("E2I_guin", "E2I_lib", "E2I_sln"),
      initializer = function (params, t0, ...)
      {
        init_counts <-
          setNames(
            c(rmultinom(1, popsize_guin, c(popsize_guin - 30, 15, 10, 5)),
              rmultinom(1, popsize_lib, c(popsize_lib - 30, 15, 10, 5)),
              rmultinom(1, popsize_sln, c(popsize_sln - 30, 15, 10, 5)),
              rep(0, 3)),
            c("S_guin","E_guin","I_guin","R_guin",
              "S_lib","E_lib","I_lib","R_lib",
              "S_sln","E_sln","I_sln","R_sln",
              "E2I_guin","E2I_lib","E2I_sln")
          )
        while (all(init_counts[c("E_guin", "I_guin")] == 0)) {
          init_counts <-
            setNames(
              c(rmultinom(1, popsize_guin, c(popsize_guin - 30, 15, 10, 5)),
                rmultinom(1, popsize_lib, c(popsize_lib - 30, 15, 10, 5)),
                rmultinom(1, popsize_sln, c(popsize_sln - 30, 15, 10, 5)),
                rep(0, 3)),
              c("S_guin","E_guin","I_guin","R_guin",
                "S_lib","E_lib","I_lib","R_lib",
                "S_sln","E_sln","I_sln","R_sln",
                "E2I_guin","E2I_lib","E2I_sln")
            )
        }
        return(init_counts)
      },
      params = c(log_Reff_m1_o_guin = 8.1,
                 log_Reff_m1_o_lib = 8.9,
                 log_Reff_m1_o_sln = 9.1,
                 log_Rext_g2l = -3.9,
                 log_Rext_g2s = -3.9,
                 log_Rext_l2g = -3.9,
                 log_Rext_l2s = -3.9,
                 log_Rext_s2g = -3.9,
                 log_Rext_s2l = -3.9,
                 log_Neff_x_rho_guin = 9.5,
                 log_Neff_x_rho_lib = 9.9,
                 log_Neff_x_rho_sln = 9.9,
                 log_omega_d_mu_guin = 0.1,
                 log_omega_d_mu_lib = -0.1,
                 log_omega_d_mu_sln = 0,
                 log_mu_inv_guin = 0.1,
                 log_mu_inv_lib = -0.1,
                 log_mu_inv_sln = 0,
                 logit_rho_guin = 0.7,
                 logit_rho_lib = 0.3,
                 logit_rho_sln = 1.4,
                 log_sqrt_phi_inv_guin = -2,
                 log_sqrt_phi_inv_lib = -2,
                 log_sqrt_phi_inv_sln = -2),
      dprior = prior_density
)

# initial values
params_init <- c(log_Reff_m1_o_guin = 8.1,
                 log_Reff_m1_o_lib = 8.9,
                 log_Reff_m1_o_sln = 9.1,
                 log_Rext_g2l = -3.9,
                 log_Rext_g2s = -3.9,
                 log_Rext_l2g = -3.9,
                 log_Rext_l2s = -3.9,
                 log_Rext_s2g = -3.9,
                 log_Rext_s2l = -3.9,
                 log_Neff_x_rho_guin = 9.5,
                 log_Neff_x_rho_lib = 9.9,
                 log_Neff_x_rho_sln = 9.9,
                 log_omega_d_mu_guin = 0.1,
                 log_omega_d_mu_lib = -0.1,
                 log_omega_d_mu_sln = 0,
                 log_mu_inv_guin = 0.1,
                 log_mu_inv_lib = -0.1,
                 log_mu_inv_sln = 0,
                 logit_rho_guin = 0.7,
                 logit_rho_lib = 0.3,
                 logit_rho_sln = 1.4,
                 log_sqrt_phi_inv_guin = -2,
                 log_sqrt_phi_inv_lib = -2,
                 log_sqrt_phi_inv_sln = -2)

starts <- function() params_init + rnorm(length(params_init), 0, 0.1)

# run the MCMC
covmat <- diag(0.01, length(params_init))
rownames(covmat) <- colnames(covmat) <- names(params_init)

registerDoParallel(5)

pomp_results <- foreach(chain = 1:5,
                        .packages=c("pomp"),
                        .options.RNG = 52787,
                        .export = ls(all.names = T)) %dorng% { 
                              
                          start_time = Sys.time()
                          res_adapt <-
                                    pmcmc(ebola_mod,
                                          Nmcmc  = 5e4,
                                          Np     = 1e3,
                                          start  = starts(),
                                          proposal =
                                                mvn.rw.adaptive(rw.var = covmat,
                                                                scale.cooling = 0.99975,
                                                                scale.start = 100,
                                                                shape.start = 100
                                                ))
                              
                              # freeze adaptation and do final run
                              res <-
                                    continue(res_adapt,
                                             Nmcmc = 3e5,
                                             proposal =
                                                   mvn.rw(covmat(res_adapt)))
         
                              end_time = Sys.time()
                              runtime = difftime(end_time, start_time)
                              
                              # return results
                              return(list(res = res, runtime = runtime))
                        }

save(pomp_results, file = "ebola_results_pomp_long.Rdata")
