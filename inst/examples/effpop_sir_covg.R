library(stemr)
library(extraDistr)
library(foreach)
library(doRNG)
library(doParallel)
library(coda)

args <- commandArgs(TRUE)
print(args)
replication <- as.numeric(args[1])

popsize <- 1e5
S0 <- 1e5-10
I0 <- 10
R0 <- 0

set.seed(100817 + replication)

true_pars <-
      c(R0 = 1 + exp(rnorm(1, 0, 0.5)),
        mu = exp(rnorm(1, -0.7, 0.35)),
        rho = expit(rnorm(1, 0, 1.4)),
        phi = rexp(1, 0.1),
        effpop = runif(1, 5e3,5e4))

# no strata the stemr object --------------------------------------------------

strata <- NULL
compartments <- c("S", "I", "R")
rates <- list(rate("beta * I * (S - effpop)", from = "S", to = "I", incidence = T, lumped = TRUE),
              rate("mu", "I", "R"))
state_initializer <- list(stem_initializer(c(S = S0, I = I0, R = R0), fixed = T))
adjacency <- NULL
tcovar <- NULL
parameters = c(true_pars["R0"] / true_pars["effpop"] * true_pars["mu"], true_pars["mu"], true_pars["rho"], true_pars["phi"], popsize - true_pars["effpop"])
names(parameters) <- c("beta", "mu", "rho", "phi", "effpop")
constants <- c(t0 = 0)
t0 <- 0; tmax <- 100

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
            rtol = 1e-6,
            atol = 1e-6
      )

emissions <- list(emission("S2I", "negbinomial", c("phi", "S2I * rho"), incidence = TRUE, obstimes = seq(1, tmax, by =1)))

measurement_process <- stem_measure(emissions = emissions, dynamics = dynamics, messages = T)

stem_object <- stem(dynamics = dynamics, measurement_process = measurement_process)

dat <- matrix(0.0, nrow = tmax, ncol = 2)

while(max(dat[,2]) < 15) {
      
        stem_data <- simulate_stem(stem_object  = stem_object,
                                   method       = "gillespie",
                                   paths        = TRUE,
                                   observations = T,
                                   nsim         = 1,
                                   census_times = unique(c(0:tmax)))

        # grab the dataset
        true_path <- stem_data$paths[[1]]
        dat <- stem_data$datasets[[1]]
}

g <- which(true_path[-1,5] < 25)
if(length(g) != 0) {
      tmax <- g[which(g >= 15)[1]]
      if(tmax <= 15 & !which.max(true_path[,5])>tmax) {
            tmax <- 15
      } else if(tmax > 50 | (tmax <= 15 & which.max(true_path[,5])>tmax)) {
            tmax <- 50
      }
      dat <- dat[1:tmax,]
      true_path <- true_path[1:tmax,]  
}

emissions <- list(emission("S2I", "negbinomial", c("phi", "S2I * rho"), incidence = TRUE, obstimes = seq(1, tmax, by =1)))
measurement_process <- stem_measure(emissions = emissions, dynamics = dynamics, data = dat)
stem_object <- stem(dynamics = dynamics, measurement_process = measurement_process)

# initialize the inference
### Parameterization in terms of log(R0) and log(mu)
## Priors for log(R0), log(mu), rho, log(phi)
# Parameters (natural scale): beta, mu, rho, phi
# Parameters (estimation scale): log(beta * N / mu), log(mu), logit(rho), log(phi)
to_estimation_scale = function(params_nat) {
      
      c(log(params_nat[1] * (popsize - params_nat[5]) / params_nat[2] - 1), # (beta,mu,Neff) -> log(R0-1)
        log(params_nat[2]),                               # mu -> log(mu)
        logit(params_nat[3]),                             # rho -> logit(rho)
        log(params_nat[4]),                               # phi -> log(phi)
        log(popsize - params_nat[5]) + log(params_nat[3]))                     
}

from_estimation_scale = function(params_est) {
      
      rho <- expit(params_est[3])
      l_effpop <- params_est[5] - log(rho)
      
      c(exp(log(exp(params_est[1])+1) + params_est[2] - l_effpop), # (log(R0), log(mu), N) -> beta = exp(log(R0) + log(mu) - log(N))
        exp(params_est[2]),                # log(mu) -> mu
        rho,                               # logit(rho) -> rho
        exp(params_est[4]),                # log(phi) -> phi
        popsize - exp(l_effpop))           # log(effpop) -> effpop 
}

prior_density = 
      function(params_nat, params_est) {
            
            l_effpop <- params_est[5] - log(expit(params_est[3]))
            
            sum(dnorm(params_est[1], 0, 0.5, log = TRUE),
                dnorm(params_est[2], -0.7, 0.35, log = TRUE),
                dnorm(params_est[3], 0, 1.4, log = TRUE),
                dexp(exp(params_est[4]), 0.1, log = TRUE) + params_est[4],
                dunif(exp(l_effpop), 5e3, 5e4, log = T) + l_effpop) 
      }

priors <- list(prior_density = prior_density,
               to_estimation_scale = to_estimation_scale,
               from_estimation_scale = from_estimation_scale)

covmat <- diag(0.01, 5)
rownames(covmat) <- colnames(covmat) <- c("log_R0_m1", "log_mu", "logit_rho", "log_phi", "log_Neff_x_rho")
mcmc_kernel <-
        kernel(
                method = "mvn_g_adaptive",
                stop_adaptation = 1e4,
                sigma = covmat,
                scale_constant = 0.5,
                scale_cooling = 0.99,
                step_size = 0.1,
                nugget = 1e-5,
                messages = FALSE
        )

stem_object$dynamics$parameters <- function() {
      priors$from_estimation_scale(priors$to_estimation_scale(parameters) + rnorm(5, 0, 0.01))
}

# register the cluster and set the seed
registerDoParallel(5)

# # Estimate an empirical covariance matrix
results <- foreach(chain = 1:5,
                   .packages="stemr",
                   .options.RNG = 52787,
                   .export = ls(all.names = T)) %dorng% {
                         
                         chain_res <- stem_inference(stem_object = stem_object,
                                                     method = "lna",
                                                     iterations = 3.5e4,
                                                     thin_params = 10,
                                                     thin_latent_proc = 10,
                                                     initialization_attempts = 500,
                                                     priors = priors,
                                                     mcmc_kernel = mcmc_kernel,
                                                     t0_kernel = t0_kernel,
                                                     messages = FALSE)
                         return(chain_res)
                   }

# collect the results
posterior_samples   <- cbind(chain = rep(1:5, each = 2.5e3), iter = seq(1,25000,by=10), do.call(rbind, lapply(results, function(x) x$results$MCMC_results[-c(1:1001),])))
posterior_quantiles <- apply(posterior_samples[,c("log_R0_m1", "log_mu", "logit_rho", "log_phi", "log_Neff_x_rho")], 2, quantile, c(0.025, 0.1, 0.25, 0.5, 0.75, 0.9, 0.975))
posterior_quantiles <- cbind(posterior_quantiles, log_Neff = quantile(log(1e5 - posterior_samples$effpop), c(0.025, 0.1, 0.25, 0.5, 0.75, 0.9, 0.975)))

latent_posts <- lapply(results, function(x) x$results$lna_paths[,-1,-1])
latent_S2I   <- do.call(cbind, lapply(latent_posts, function(x) x[,1,]))
latent_I2R   <- do.call(cbind, lapply(latent_posts, function(x) x[,2,]))

latent_quantiles <- list(S2I = apply(latent_S2I, 1, quantile, c(0.025, 0.1, 0.25, 0.5, 0.75, 0.9, 0.975)),
                         I2R = apply(latent_I2R, 1, quantile, c(0.025, 0.1, 0.25, 0.5, 0.75, 0.9, 0.975)))

posterior_res_mcmc <- as.mcmc.list(lapply(results, 
                                          function(x) mcmc(
                                                cbind(logpost = rowSums(x$results$MCMC_results[-c(1:1001),1:3]),
                                                      x$results$MCMC_results[-c(1:1001),c("log_R0_m1", "log_mu", "logit_rho", "log_phi", "log_Neff_x_rho")]))))
psrf              <- gelman.diag(posterior_res_mcmc)
effective_samples <- do.call(rbind, lapply(posterior_res_mcmc, effectiveSize))

posterior_results <- list(true_pars = true_pars,
                          true_path = true_path,
                          dat = dat,
                          posterior_quantiles = posterior_quantiles, 
                          latent_quantiles = latent_quantiles, 
                          effective_samples = effective_samples,
                          psrf = psrf,
                          times = sapply(results, function(x) x$results$time))

save(posterior_results, file = paste0("sir_effpop_",replication,".Rdata"))
