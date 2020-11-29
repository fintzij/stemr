library(tidyverse)
library(DescTools)
library(patchwork)
library(here)
setwd(here("/paper_V2/simulations/rev_sim_01_1wk_durs"))
source(here("/paper_V2/simulations/rev_sim_01_1wk_durs/process_fcns.R"))

true_pars =
    c(R0         = 2, # basic reproduction number
      latent_dur = 1, # latent period
      infec_dur  = 1, # infectious period duration = 1 week
      rho        = 0.5, # case detection rate
      sqrt_phi_inv = 1/sqrt(36))  # negative binomial overdispersion
 
# load results ------------------------------------------------------------

res_1 = res_4 = 
    data.frame(replication = rep(1:1e3, each = 5),
               param = c("R0", 
                         "latent_dur",
                         "infec_dur",
                         "rho",
                         "sqrt_phi_inv"),
               rel_mad = 0,
               rel_ciw = 0,
               contraction = 0,
               coverage = 0)

for(i in seq_len(1e3)[-29]) {
    
    # load fit summary
    fit_1_i = readRDS(paste0("results/seir_1_",i,".Rds"))
    fit_4_i = readRDS(paste0("results/seir_4_",i,".Rds"))
    
    # calc MADs
    res_1$rel_mad[with(res_1, which(param == "R0" & replication == i))] = 
        rel_mad("R0", fit_1_i, true_pars)
    res_1$rel_mad[with(res_1, which(param == "latent_dur" & replication == i))] = 
        rel_mad("latent_dur", fit_1_i, true_pars)
    res_1$rel_mad[with(res_1, which(param == "infec_dur" & replication == i))] = 
        rel_mad("infec_dur", fit_1_i, true_pars)
    res_1$rel_mad[with(res_1, which(param == "rho" & replication == i))] = 
        rel_mad("rho", fit_1_i, true_pars)
    res_1$rel_mad[with(res_1, which(param == "sqrt_phi_inv" & replication == i))] = 
        rel_mad("sqrt_phi_inv", fit_1_i, true_pars)
    
    # calc ciws
    res_1$rel_ciw[with(res_1, which(param == "R0" & replication == i))] = 
        rel_ciw("R0", fit_1_i, true_pars)
    res_1$rel_ciw[with(res_1, which(param == "latent_dur" & replication == i))] = 
        rel_ciw("latent_dur", fit_1_i, true_pars)
    res_1$rel_ciw[with(res_1, which(param == "infec_dur" & replication == i))] = 
        rel_ciw("infec_dur", fit_1_i, true_pars)
    res_1$rel_ciw[with(res_1, which(param == "rho" & replication == i))] = 
        rel_ciw("rho", fit_1_i, true_pars)
    res_1$rel_ciw[with(res_1, which(param == "sqrt_phi_inv" & replication == i))] = 
        rel_ciw("sqrt_phi_inv", fit_1_i, true_pars)
    
    # calc ciws
    res_1$contraction[with(res_1, which(param == "R0" & replication == i))] = 
        contraction("R0", fit_1_i)
    res_1$contraction[with(res_1, which(param == "latent_dur" & replication == i))] = 
        contraction("latent_dur", fit_1_i)
    res_1$contraction[with(res_1, which(param == "infec_dur" & replication == i))] = 
        contraction("infec_dur", fit_1_i)
    res_1$contraction[with(res_1, which(param == "rho" & replication == i))] = 
        contraction("rho", fit_1_i)
    res_1$contraction[with(res_1, which(param == "sqrt_phi_inv" & replication == i))] = 
        contraction("sqrt_phi_inv", fit_1_i)
    
    # calc coverage
    res_1$coverage[with(res_1, which(param == "R0" & replication == i))] = 
        coverage("R0", fit_1_i, true_pars)
    res_1$coverage[with(res_1, which(param == "latent_dur" & replication == i))] = 
        coverage("latent_dur", fit_1_i, true_pars)
    res_1$coverage[with(res_1, which(param == "infec_dur" & replication == i))] = 
        coverage("infec_dur", fit_1_i, true_pars)
    res_1$coverage[with(res_1, which(param == "rho" & replication == i))] = 
        coverage("rho", fit_1_i, true_pars)
    res_1$coverage[with(res_1, which(param == "sqrt_phi_inv" & replication == i))] = 
        coverage("sqrt_phi_inv", fit_1_i, true_pars)
    
    # calc MADs
    res_4$rel_mad[with(res_4, which(param == "R0" & replication == i))] = 
        rel_mad("R0", fit_4_i, true_pars)
    res_4$rel_mad[with(res_4, which(param == "latent_dur" & replication == i))] = 
        rel_mad("latent_dur", fit_4_i, true_pars)
    res_4$rel_mad[with(res_4, which(param == "infec_dur" & replication == i))] = 
        rel_mad("infec_dur", fit_4_i, true_pars)
    res_4$rel_mad[with(res_4, which(param == "rho" & replication == i))] = 
        rel_mad("rho", fit_4_i, true_pars)
    res_4$rel_mad[with(res_4, which(param == "sqrt_phi_inv" & replication == i))] = 
        rel_mad("sqrt_phi_inv", fit_4_i, true_pars)
    
    # calc ciws
    res_4$rel_ciw[with(res_4, which(param == "R0" & replication == i))] = 
        rel_ciw("R0", fit_4_i, true_pars)
    res_4$rel_ciw[with(res_4, which(param == "latent_dur" & replication == i))] = 
        rel_ciw("latent_dur", fit_4_i, true_pars)
    res_4$rel_ciw[with(res_4, which(param == "infec_dur" & replication == i))] = 
        rel_ciw("infec_dur", fit_4_i, true_pars)
    res_4$rel_ciw[with(res_4, which(param == "rho" & replication == i))] = 
        rel_ciw("rho", fit_4_i, true_pars)
    res_4$rel_ciw[with(res_4, which(param == "sqrt_phi_inv" & replication == i))] = 
        rel_ciw("sqrt_phi_inv", fit_4_i, true_pars)
    
    # calc ciws
    res_4$contraction[with(res_4, which(param == "R0" & replication == i))] = 
        contraction("R0", fit_4_i)
    res_4$contraction[with(res_4, which(param == "latent_dur" & replication == i))] = 
        contraction("latent_dur", fit_4_i)
    res_4$contraction[with(res_4, which(param == "infec_dur" & replication == i))] = 
        contraction("infec_dur", fit_4_i)
    res_4$contraction[with(res_4, which(param == "rho" & replication == i))] = 
        contraction("rho", fit_4_i)
    res_4$contraction[with(res_4, which(param == "sqrt_phi_inv" & replication == i))] = 
        contraction("sqrt_phi_inv", fit_4_i)
    
    # calc coverage
    res_4$coverage[with(res_4, which(param == "R0" & replication == i))] = 
        coverage("R0", fit_4_i, true_pars)
    res_4$coverage[with(res_4, which(param == "latent_dur" & replication == i))] = 
        coverage("latent_dur", fit_4_i, true_pars)
    res_4$coverage[with(res_4, which(param == "infec_dur" & replication == i))] = 
        coverage("infec_dur", fit_4_i, true_pars)
    res_4$coverage[with(res_4, which(param == "rho" & replication == i))] = 
        coverage("rho", fit_4_i, true_pars)
    res_4$coverage[with(res_4, which(param == "sqrt_phi_inv" & replication == i))] = 
        coverage("sqrt_phi_inv", fit_4_i, true_pars)
}

res = data.frame(obs_interval = rep(c("Weekly", "Monthly"), each = 5e3),
                 rbind(res_1, res_4))

res$obs_interval = 
    factor(res$obs_interval, levels = c("Weekly", "Monthly"))
res$param = 
    factor(res$param, 
           levels = c("R0", "latent_dur", "infec_dur", "rho", "sqrt_phi_inv"))

# generate plots ----------------------------------------------------------

mad_plot <- 
    res %>% 
    ggplot(aes(x = as.factor(obs_interval), y = rel_mad, 
               fill = as.factor(obs_interval))) + 
    geom_violin(alpha = 0.5, draw_quantiles = c(0.1, 0.5, 0.9), colour = "grey50") + 
    facet_grid(.~param, scales = "free") + 
    scale_fill_brewer(type = "qual") + 
    scale_color_brewer(type = "qual") + 
    theme_minimal() + 
    theme(legend.position = "none") + 
    labs(x = NULL, y = "Scaled MAD", title = "Scaled median absolute deviation")

ciw_plot <- 
    res %>% 
    ggplot(aes(x = as.factor(obs_interval), y = rel_ciw, 
               fill = as.factor(obs_interval))) + 
    geom_violin(alpha = 0.5, draw_quantiles = c(0.1, 0.5, 0.9), colour = "grey50") + 
    facet_grid(.~param, scales = "free") + 
    scale_fill_brewer(type = "qual") + 
    scale_color_brewer(type = "qual") + 
    theme_minimal() + 
    theme(legend.position = "none") + 
    labs(x = NULL, y = "Scaled CIW", title = "Scaled credible interval width") 

contraction_plot <- 
    res %>% 
    ggplot(aes(x = as.factor(obs_interval), y = contraction, 
               fill = as.factor(obs_interval))) + 
    geom_violin(alpha = 0.5, draw_quantiles = c(0.1, 0.5, 0.9), colour = "grey50") + 
    facet_grid(.~param, scales = "free") + 
    scale_fill_brewer(type = "qual") + 
    scale_color_brewer(type = "qual") + 
    theme_minimal() + 
    theme(legend.position = "none") + 
    labs(x = NULL, y = "Posterior contraction", title = "Posterior CI width vs. prior CI width ")

res_covg <-
    res %>% 
    group_by(param, obs_interval) %>% 
    summarise(covg = mean(coverage),
              lower = BinomCI(x = sum(coverage), 1e3, method = "wilson")[2],
              upper = BinomCI(x = sum(coverage), 1e3, method = "wilson")[3]) %>% 
    ungroup %>% as.data.frame()

covg_plot <- 
    res_covg %>% 
    ggplot(aes(x = obs_interval, y = covg, 
               ymin = lower, ymax = upper,
               colour = obs_interval)) + 
    geom_pointrange() +
    geom_hline(yintercept = 0.95, linetype = "dashed", colour = "grey50") + 
    facet_grid(.~param, scales = "free") + 
    scale_fill_brewer(type = "qual") + 
    scale_color_brewer(type = "qual") + 
    theme_minimal() + 
    theme(legend.position = "none") + 
    labs(x = NULL, y = "coverage", title = "Coverage of 95% credible intervals")

all_plots <- covg_plot / contraction_plot / mad_plot / ciw_plot

ggsave("weekly_vs_monthly_res.pdf", plot = all_plots, width = 14, height = 12)
