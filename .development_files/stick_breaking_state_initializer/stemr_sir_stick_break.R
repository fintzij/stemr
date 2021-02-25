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
    prior = c(3, 3, 0.5, 0.5),
    dist = "rsbln",
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











## ----sim_paths, echo = TRUE----------------------------------------------

# Need check storage of initial compartments
set.seed(100)
# debugonce(simulate_stem)
sim_mjp_200 <- simulate_stem(stem_object = stem_object, method = "gillespie", nsim = 20)
set.seed(100)
sim_ode_200 <- simulate_stem(stem_object = stem_object, method = "ode", nsim = 20)
set.seed(100)
sim_lna_200 <- simulate_stem(stem_object = stem_object, method = "lna", nsim = 20)


sim_mjp_200$paths
library(tidyverse)
sim_mjp_200$datasets %>% imap_dfr(~ as_tibble(.x) %>% mutate(i = .y)) %>%
  ggplot(aes(time, S2I, group = i, color = i)) +
  geom_step() +
  geom_jitter()
sim_mjp_200$failed_runs
sim_mjp_1 <- simulate_stem(stem_object = stem_object, method = "gillespie", nsim = 200)

# fais :(
sim_ode <- simulate_stem(stem_object = stem_object, method = "ode", nsim = 200)
library(tidyverse)
library(tidybayes)
sim_ode$natural_paths %>%
  imap_dfr(~as_tibble(.x) %>% mutate(sim = .y)) %>%
  pivot_longer(-c(time, sim)) %>%
  select(-sim) %>%
  group_by(time, name) %>%
  median_qi(.width = c(0.5, 0.8, 0.95)) %>%
  ggplot(aes(time, value, ymin = .lower, ymax = .upper)) +
  facet_wrap(. ~ name, scales = "free_y") +
  geom_lineribbon() +
  scale_fill_brewer()

# sim_mjp_blank <- simulate_stem(stem_object = stem_object, method = "gillespie", full_paths = T, nsim = 200)


