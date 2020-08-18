# Lump Counts -------------------------------------------------------------
## ochca_covid: sanitized ochca data from `read_santize_ochca_covid`
## time_interval_in_days: time interval to lump in days
## first day: "yyyy-mm-dd" string - filter out days before this one
## last day: "yyyy-mm-dd" string - filter out days after this one
lump_ochca_covid <- function(ochca_covid, time_interval_in_days, first_day = "0000-01-01", last_day = "9999-12-31") {
  ochca_covid %>%
    filter(posted_date >= lubridate::ymd(first_day),
           posted_date <= lubridate::ymd(last_day)) %>%
    group_by(lump = as.integer(floor((max(posted_date) - posted_date) / time_interval_in_days))) %>%
    filter(n() == time_interval_in_days) %>%
    summarize(start_date = min(posted_date),
              end_date = max(posted_date),
              new_cases = sum(new_cases),
              new_tests = sum(new_tests),
              new_deaths = sum(new_deaths)) %>%
    dplyr::select(-lump) %>%
    arrange(start_date)
}

# Create Model Objects ----------------------------------------------------
create_model_objects <- function(ochca_covid = read_rds("data/oc/ochca_covid.rds"),
                                 first_day = "0000-01-01",
                                 last_day = "9999-12-31",
                                 time_interval_in_days = 3,
                                 n_curves = 10,
                                 epi_curve_names = c("S", "E", "IE", "IP", "cumulative IE -> IP", "D", "rho incidence", "log mu death", "kappa incid", "phi death"),
                                 priors_only = F,
                                 forecast_in_days = 0,
                                 future_tests_per_day = 600,
                                 popsize = 3175692L,
                                 other_objects = list(),
                                 log_R0_normal = c(-0.2554128198465173693599, 0.6908418304107377672096),
                                 latent_dur_lognormal = c(log(7) - log(7), 0.22),
                                 early_dur_lognormal = c(log(7) - log(7), 0.22),
                                 prog_dur_lognormal = c(log(7) - log(7), 0.22),
                                 IFR_beta = c(1.5, 200),
                                 frac_carrs_beta = c(1, 1e4),
                                 frac_carrs_infec_beta = c(3, 3),
                                 frac_infec_early_beta = c(3, 3),
                                 alpha_incid_0_normal = c(4, 2),
                                 alpha_incid_1_beta = c(3, 1),
                                 sqrt_kappa_inv_incid_exponential = 1,
                                 rho_death_beta = c(8, 2),
                                 sqrt_phi_inv_death_exponential = c(1)) {

  # dat <- read_santize_ochca_covid(file_name = file_name) for old data

  lumped_ochca_covid <- lump_ochca_covid(ochca_covid,
                                         time_interval_in_days = time_interval_in_days,
                                         first_day = first_day,
                                         last_day = last_day)

  actual_first_day <- min(lumped_ochca_covid$start_date)
  actual_last_day <- max(lumped_ochca_covid$end_date)

  first_forecast_day <- min(lumped_ochca_covid$end_date)
  last_forecast_day <- actual_last_day + forecast_in_days

  n_times <- nrow(lumped_ochca_covid) + 1
  n_states <- 1
  n_obs <- n_times - 1
  num_per <- n_times

  # cases, tests, deaths, popsize, ntimes
  x_i <- matrix(c(lumped_ochca_covid$new_cases, lumped_ochca_covid$new_tests, lumped_ochca_covid$new_deaths, popsize, n_times), nrow = 1)

  dates <- seq(from = first_forecast_day, to = actual_last_day, by = time_interval_in_days)
  # dates <- lumped_ochca_covid$end_date # should match
  # 0, times since first day in weeks
  x_r <- matrix(c(0, as.integer(dates - actual_first_day + 1) / 7), nrow = n_states)


  # Posterior Predictive
  n_times_pp <- as.integer(floor((last_forecast_day - first_forecast_day) / time_interval_in_days) + 2)
  n_obs_pp <- n_times_pp - 1
  num_per_pp <- n_times_pp

  obs_inds_pp <- seq(from = n_curves, by = n_curves, length.out = n_obs_pp)
  tests <- c(lumped_ochca_covid$new_tests, rep(future_tests_per_day * time_interval_in_days, n_times_pp - n_times))

  # cases (0), tests, deaths(0), popsize, ntimes
  x_i_pp <- matrix(c(rep(0, n_obs_pp), tests, rep(0, n_obs_pp), popsize, n_times_pp), nrow = 1)
  dates_pp <- seq(from = first_forecast_day, to = last_forecast_day, by = time_interval_in_days)
  # 0, times since first day in weeks
  x_r_pp <- matrix(c(0, as.integer(dates_pp - actual_first_day + 1) / 7), nrow = n_states)


  lumped_ochca_covid_forecast <- lump_ochca_covid(ochca_covid,
                                                  time_interval_in_days = time_interval_in_days,
                                                  first_day = actual_first_day,
                                                  # last_day = dates_pp[length(dates_pp)],
                                                  last_day = ochca_covid[["posted_date"]][max(which(as.integer(dates_pp[length(dates_pp)] - ochca_covid[["posted_date"]]) %% time_interval_in_days == 0))]) %>%
    mutate(usage = ifelse(end_date <= actual_last_day, "train", "test"))

  # Build Model Objects
  model_objects <- other_objects

  model_objects$n_curves = n_curves
  model_objects$priors_only = priors_only
  model_objects$forecast_in_days = forecast_in_days
  model_objects$future_tests_per_day = future_tests_per_day
  model_objects$priors_only <- priors_only
  model_objects$time_interval_in_days <- time_interval_in_days

  model_objects$n_times <- n_times
  model_objects$n_states <- n_states
  model_objects$n_obs <- n_obs
  model_objects$num_per <- num_per
  model_objects$popsize <- popsize

  model_objects$x_i <- x_i
  model_objects$x_r <- x_r
  model_objects$dates <- dates

  model_objects$n_times_pp <- n_times_pp
  model_objects$n_obs_pp <- n_obs_pp
  model_objects$num_per_pp <- num_per_pp

  model_objects$obs_inds_pp <- obs_inds_pp
  model_objects$tests <- tests

  model_objects$x_i_pp <- x_i_pp
  model_objects$x_r_pp <- x_r_pp
  model_objects$dates_pp <- dates_pp

  model_objects$actual_first_day <- actual_first_day
  model_objects$actual_last_day <- actual_last_day

  model_objects$first_forecast_day <- first_forecast_day
  model_objects$last_forecast_day <- last_forecast_day

  model_objects$lumped_ochca_covid <- lumped_ochca_covid
  model_objects$lumped_ochca_covid_forecast <- lumped_ochca_covid_forecast

  model_objects$epi_curve_names <- epi_curve_names

  # model_objects$popsizes <- rep(popsize, n_obs_pp)

  # Priors
  model_objects$log_R0_normal <- log_R0_normal
  model_objects$latent_dur_lognormal <-   latent_dur_lognormal
  model_objects$early_dur_lognormal <-   early_dur_lognormal
  model_objects$prog_dur_lognormal <-   prog_dur_lognormal
  model_objects$IFR_beta <-   IFR_beta
  model_objects$frac_carrs_beta <-   frac_carrs_beta
  model_objects$frac_carrs_infec_beta <-   frac_carrs_infec_beta
  model_objects$frac_infec_early_beta <-   frac_infec_early_beta
  model_objects$alpha_incid_0_normal <-   alpha_incid_0_normal
  model_objects$alpha_incid_1_beta <-   alpha_incid_1_beta
  model_objects$sqrt_kappa_inv_incid_exponential <-   sqrt_kappa_inv_incid_exponential
  model_objects$rho_death_beta <-   rho_death_beta
  model_objects$sqrt_phi_inv_death_exponential <-   sqrt_phi_inv_death_exponential

  model_objects
}

#make keiths data match ours ----------------------------------------------------
ca_lump<-function(keith_data){
  keith_data<-keith_data%>%
    dplyr::select(endPt, pos, tests, death, t0)%>%
    rename(end_date=endPt, new_cases=pos, new_tests=tests, new_deaths=death)

  keith_data
}

# CA Create Model Objects (for Keith Data) ----------------------------------------------------
ca_create_model_objects <- function(lumped_data,
                                 first_day = "0000-01-01",
                                 last_day = "9999-12-31",
                                 time_interval_in_days = 4,
                                 n_curves = 8,
                                 epi_curve_names = c("S", "E", "IM", "IS", "cumulative Im -> Is", "D", "rho incidence", "log mu death"),
                                 priors_only = F,
                                 forecast_in_days = 0,
                                 future_tests_per_day = 600,
                                 popsize = 39.51e6,
                                 other_objects = list()) {

  #ochca_covid <- read_santize_ochca_covid(file_name = file_name)

  lumped_ca_covid <- lumped_data
  #popsize<-39.51e6
  #time_interval_in_days<-4
  #forecast_in_days = 0
  #n_curves = 8
  #future_tests_per_day = 600

  actual_first_day <- min(lumped_ca_covid$end_date)
  actual_last_day <- max(lumped_ca_covid$end_date)

  first_forecast_day <- min(lumped_ca_covid$end_date[lumped_ca_covid$t0!=1])
  last_forecast_day <- actual_last_day + forecast_in_days

  #now filter out start of data
  lumped_ca_covid<-lumped_ca_covid%>%
    filter(t0!=1)

  n_times <- nrow(lumped_ca_covid) + 1
  n_states <- 1
  n_obs <- n_times - 1
  num_per <- n_times

  # cases, tests, deaths, popsize, ntimes
  x_i <- matrix(c(lumped_ca_covid$new_cases, lumped_ca_covid$new_tests, lumped_ca_covid$new_deaths, popsize, n_times), nrow = 1)

  dates <- seq(from = first_forecast_day, to = actual_last_day, by = time_interval_in_days)
  # dates <- lumped_ca_covid$end_date # should match
  # 0, times since first day in weeks
  x_r <- matrix(c(0, as.integer(dates - actual_first_day + 1) / 7), nrow = n_states)


  # Posterior Predictive
  n_times_pp <- as.integer(floor((last_forecast_day - first_forecast_day) / time_interval_in_days) + time_interval_in_days-1)
  n_obs_pp <- n_times_pp - 1
  num_per_pp <- n_times_pp

  obs_inds_pp <- seq(from = n_curves, by = n_curves, length.out = n_obs_pp)
  tests <- c(lumped_ca_covid$new_tests, rep(future_tests_per_day * time_interval_in_days, n_times_pp - n_times))

  # cases (0), tests, deaths(0), popsize, ntimes
  x_i_pp <- matrix(c(rep(0, n_obs_pp), tests, rep(0, n_obs_pp), popsize, n_times_pp), nrow = 1)
  dates_pp <- seq(from = first_forecast_day, to = last_forecast_day, by = time_interval_in_days)
  # 0, times since first day in weeks
  x_r_pp <- matrix(c(0, as.integer(dates_pp - actual_first_day + 1) / 7), nrow = n_states)


  # Build Model Objects
  model_objects <- other_objects

  model_objects$n_curves = n_curves
  model_objects$priors_only = priors_only
  model_objects$forecast_in_days = forecast_in_days
  model_objects$future_tests_per_day = future_tests_per_day
  model_objects$priors_only <- priors_only
  model_objects$time_interval_in_days <- time_interval_in_days

  model_objects$n_times <- n_times
  model_objects$n_states <- n_states
  model_objects$n_obs <- n_obs
  model_objects$num_per <- num_per
  model_objects$popsize <- popsize

  model_objects$x_i <- x_i
  model_objects$x_r <- x_r
  model_objects$dates <- dates

  model_objects$n_times_pp <- n_times_pp
  model_objects$n_obs_pp <- n_obs_pp
  model_objects$num_per_pp <- num_per_pp

  model_objects$obs_inds_pp <- obs_inds_pp
  model_objects$tests <- tests

  model_objects$x_i_pp <- x_i_pp
  model_objects$x_r_pp <- x_r_pp
  model_objects$dates_pp <- dates_pp

  model_objects$actual_first_day <- actual_first_day
  model_objects$actual_last_day <- actual_last_day

  model_objects$first_forecast_day <- first_forecast_day
  model_objects$first_last_day <- last_forecast_day

  model_objects$lumped_ca_covid <- lumped_ca_covid

  model_objects$epi_curve_names <- epi_curve_names


  model_objects
}
