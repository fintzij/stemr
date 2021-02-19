require(ggplot2)
require(cowplot)
set.seed(12511)
library(stemr)
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
    prior = c(-5.05439386925471, -3.96216160104193, 0.508799264231241, 0.447697131924911, 10000),
    dist = "rsbln",
    fixed = F)) # initial state fixed for simulation, we'll change this later

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











## ----sim_paths, echo = TRUE----------------------------------------------
sim_mjp <- simulate_stem(stem_object = stem_object, method = "gillespie", full_paths = T)
sim_lna <- simulate_stem(stem_object = stem_object, method = "lna", lna_method = "approx")
sim_ode <- simulate_stem(stem_object = stem_object, method = "ode")


## ----plot_sims, echo = FALSE, fig.width=8, fig.height=4------------------
sim_paths =
  expand.grid(time = 0:tmax,
              Method = c("Gillespie", "LNA", "ODE"),
              Compartment = c("S","I","R","S2I","I2R"),
              Type = c("Prevalence","Incidence"))
sim_paths =
  sim_paths[!((
    sim_paths$Compartment %in% c("S", "I", "R") &
      sim_paths$Type == "Incidence"
  ) |
    (
      sim_paths$Compartment %in% c("S2I", "I2R") &
        sim_paths$Type == "Prevalence"
    )
  ), ]
sim_paths$Compartment = factor(sim_paths$Compartment, levels = c("S", "I", "R", "S2I", "I2R"))
sim_paths = sim_paths[with(sim_paths, order(Method, Compartment, Type, time)), ]
sim_paths$Count =
  c(
    sim_mjp$paths[[1]][, -1],
    sim_lna$natural_paths[[1]][, -1],
    sim_lna$paths[[1]][, -1],
    sim_ode$natural_paths[[1]][, -1],
    sim_ode$paths[[1]][, -1]
  )

mjp_prev =
  data.frame(time = sim_mjp$full_paths[[1]][,1],
             Compartment = rep(c("S","I","R"), each = nrow(sim_mjp$full_paths[[1]])),
             Count = c(sim_mjp$full_paths[[1]][,3:5]))

mjp_counts =
  ggplot(mjp_prev, aes(x = time, y = Count,
                       colour = Compartment,
                       group = Compartment)) +
  geom_step() +
  theme_minimal() +
  scale_color_brewer(type = "qual", palette = 6) +
  scale_y_continuous(trans = "sqrt",
                     breaks = c(0,50, 250, 1000, 2.5e3, 5e3,7.5e3,1e4),
                     expand = c(0,0)) +
  labs(title = "Compartment counts", subtitle = "MJP")

mjp_incid =
  ggplot(subset(sim_paths, Method == "Gillespie" & Type == "Incidence"),
         aes(x = time, y = Count,
             colour = Compartment,
             group = Compartment)) +
  geom_point() +
  theme_minimal() +
  scale_color_brewer("Transition", type = "qual", palette = 2) +
  labs(title = "Incident transition events",subtitle = "MJP")

lna_prev =
  ggplot(subset(sim_paths, Method == "LNA" & Type == "Prevalence"),
         aes(x = time, y = Count,
             colour = Compartment,
             group = Compartment)) +
  geom_line(linetype = 2) +
  theme_minimal() +
  scale_color_brewer("Compartment", type = "qual", palette = 6) +
  scale_y_continuous(trans = "sqrt",
                     breaks = c(0,50, 250, 1000, 2.5e3, 5e3,7.5e3,1e4),
                     expand = c(0,0)) +
  labs(title = "", subtitle = "LNA")

lna_incid =
  ggplot(subset(sim_paths, Method == "LNA" & Type == "Incidence"),
         aes(x = time, y = Count,
             colour = Compartment,
             group = Compartment)) +
  geom_point(shape = 2) +
  theme_minimal() +
  scale_color_brewer("Transition", type = "qual", palette = 2) +
  labs(title = "",subtitle = "LNA")

ode_prev =
  ggplot(subset(sim_paths, Method == "ODE" & Type == "Prevalence"),
         aes(x = time, y = Count,
             colour = Compartment,
             group = Compartment)) +
  geom_line(linetype = 3) +
  theme_minimal() +
  scale_color_brewer("Compartment", type = "qual", palette = 6) +
  scale_y_continuous(trans = "sqrt",
                     breaks = c(0,50, 250, 1000, 2.5e3, 5e3,7.5e3,1e4),
                     expand = c(0,0)) +
  labs(title = "", subtitle = "ODE")

ode_incid =
  ggplot(subset(sim_paths, Method == "ODE" & Type == "Incidence"),
         aes(x = time, y = Count,
             colour = Compartment,
             group = Compartment)) +
  geom_point(shape = 3) +
  theme_minimal() +
  scale_color_brewer("Transition", type = "qual", palette = 2) +
  labs(title = "", subtitle = "ODE")

cowplot::plot_grid(
  mjp_counts, lna_prev, ode_prev, mjp_incid, lna_incid, ode_incid
)







## ---- plot_dat-----------------------------------------------------------
ggplot(data = as.data.frame(sim_mjp$datasets[[1]]),
       aes(x=time, y = S2I)) +
  geom_point() +
  theme_minimal() +
  labs(x = "Week", y = "Count", title = "Observed Incidence")









## ----recompile_meas, echo = TRUE-----------------------------------------
state_initializer <-
  list(
    stem_initializer(
      init_states = c(S = popsize-10, I = 10, R = 0), # must match compartment names
      fixed = TRUE,
      prior = c(popsize, 10, 0)/10)) # we now do inference on the initial compartment counts

dynamics <-
  stem_dynamics(
    rates = rates,
    tmax = tmax,
    parameters = parameters,
    state_initializer = state_initializer,
    compartments = compartments,
    constants = constants,
    compile_ode = T,   # compile ODE functions
    compile_rates = F, # compile MJP functions for Gillespie simulation
    compile_lna = T,   # compile LNA functions
    messages = F       # don't print messages
  )

measurement_process <-
  stem_measure(emissions = emissions,
               dynamics = dynamics,
               data = sim_mjp$datasets[[1]])
stem_object <- make_stem(dynamics = dynamics,
                         measurement_process = measurement_process)




## ----priors_est_scale_fcns, echo = TRUE----------------------------------

### Parameterization in terms of log(R0) and log(mu)
## Priors for log(R0), log(mu), logit(rho), phi
# Parameters (natural scale): beta, mu, rho, phi
# Parameters (estimation scale): log(beta * N / mu), log(mu), logit(rho), log(phi)

# function to take params_nat and return params_est
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







## ----mcmc_kern, echo = TRUE----------------------------------------------
# specify the initial proposal covariance matrix with row and column names
# corresponding to parameters on their estimation scales

par_initializer = function() {
  priors$from_estimation_scale(priors$to_estimation_scale(parameters) +
                                 rnorm(4, 0, 0.1))
}

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




## ----fit_mod, echo = TRUE------------------------------------------------
res <-
  fit_stem(stem_object = stem_object,
           method = "ode",
           mcmc_kern = mcmc_kern,
           iterations = 5e2,
           print_progress = 0)







## ----access_res, echo = TRUE---------------------------------------------
runtime = res$results$time
mcmc_samples = res$results$MCMC_results
ode_paths = res$results$ode_paths # or res$results$lna_paths if using the LNA for inference
adapt_par = res$results$adaptation_record$adaptation_scale_record






## ----pomp_fit, echo = TRUE, warning=FALSE, eval = FALSE------------------
## require(pomp)
## S0 <- 1e4 - 10
## I0 <- 10
## R0 <- 0
## popsize = 1e4
## cases <- sim_mjp$datasets[[1]][,2]
##
## # Set up pomp objects -----------------------------------------------------
##
## # Measurement process objects
## rmeas <- "
##   cases=rnbinom_mu(exp(log_phi),exp(logit_rho)*S2I/(1+exp(logit_rho)));    // simulates the data
## "
##
## dmeas<-"
##   lik=dnbinom_mu(cases,exp(log_phi),exp(logit_rho)*S2I/(1+exp(logit_rho)),give_log); // emission density
## "
##
## # define the stepper
## sir.step<-paste0("
##   double rate[2];
##   double dN[2];
##   rate[0]=exp(log(exp(log_R0)+1) + log_mu - log(10000))*I; // Infection rate
##   rate[1]=exp(log_mu);                         // recovery rate
##   reulermultinom(1,S,&rate[0],dt,&dN[0]);      // generate the number of newly infected people
##   reulermultinom(1,I,&rate[1],dt,&dN[1]);      // generate the number of newly recovered people
##   S+=-dN[0];                                   // update the number of Susceptibles
##   I+=dN[0]-dN[1];                              // update the number of Infections
##   R+=dN[1];                                    // update the number of Recoveries
##   S2I += dN[0];
##   I2R += dN[1];
## ")
##
## # instatiate the euler stepper function
## SIR_sim <- euler.sim(step.fun = Csnippet(sir.step), delta.t = 1/7)
##
## # Define the priors
## sir.dprior <- function(params, ..., log) {
##   l <- dnorm(params["log_R0"], 0, 0.5, log = TRUE) +
##     dnorm(params["log_mu"], -0.7, 0.35, log = TRUE) +
##     dnorm(params["logit_rho"], log = T) +
##     dexp(exp(params["log_phi"]), 0.1, log = T) + params["log_phi"]
##   if(!log) l <- exp(l)
##   return(l)
## }
##
## # instatiate the pomp object
## sir_mod <- pomp(
##   data = data.frame(time = seq(1, tmax, by = 1), cases = cases),  #"cases" is the dataset, "time" is the observation time
##   times = "time",
##   t0 = 0,                      # initial time point
##   dmeasure = Csnippet(dmeas),  # evaluates the density of the measurement process
##   rmeasure = Csnippet(rmeas),  # simulates from the measurement process
##   rprocess = SIR_sim,          # simulates from the latent process
##   statenames = c("S", "I", "R", "S2I", "I2R"),  #state space variable name
##   paramnames = c("log_R0", "log_mu", "logit_rho", "log_phi"), #parameters name
##   zeronames = c("S2I", "I2R"),
##   initializer = function(params, t0, ...) {
##                         return(c(S = S0, I = I0, R = R0, S2I = 0, I2R = 0))
##                 },
##   params = c(log_R0    = log(1.5),
##              log_mu    = log(1),
##              logit_rho = logit(0.5),
##              log_phi   = log(5)),
##   dprior = sir.dprior
## )
##
## # Metropolis kernel proposal covariance matrix
## init = c(log_R0 = log(true_pars[1]), log_mu = log(true_pars[2]), logit_rho = logit(true_pars[3]), log_phi = log(true_pars[4])) + rnorm(4, 0, c(0.1, 0.1, 0.1, 0.1))
## names(init) <- c("log_R0", "log_mu", "logit_rho", "log_phi")
## cov_init <- diag(1e-3, 4); colnames(cov_init) <- rownames(cov_init) <- names(init)
##
## # adaptive phase of MCMC
## res_adapt <- pmcmc(sir_mod,
##                    Nmcmc = 1e2, # number of adaptive iterations
##                    Np = 50, # number of particles in pmmh (set to 500 in paper)
##                    start = init,
##                   proposal = mvn.rw.adaptive(rw.var = cov_init))
##
## # final sample
## sigma <- covmat(res_adapt)
## res_pomp <- pmcmc(res_adapt,
##                   Nmcmc = 2.5e2, #number of mcmc iterations
##                   Np = 50,  # number of particles in PMMH (set to 500 in paper)
##                   start = init,
##                   proposal = mvn.rw(rw.var = sigma))
##














