library(tidyverse)
library(ggthemes)
library(DescTools)
library(kableExtra)
library(knitr)
library(patchwork)
library(here)
setwd(here("paper_V2/simulations/rev_sim_02_3dis_long"))
source(here("paper_V2/simulations/rev_sim_02_3dis_long/process_fcns.R"))

true_pars_flu =
  c(R0         = 1.3, # basic reproduction number
    latent_dur = 2/7, # latent period duration = 1 week
    infec_dur  = 2.5/7, # infectious period duration = 1 week
    rho        = 0.02, # case detection rate
    sqrt_phi_inv = 1/sqrt(36))  # negative binomial overdispersion

true_pars_ebola =
  c(R0         = 1.8, # basic reproduction number
    latent_dur = 10/7, # latent period duration = 1 week
    infec_dur  = 1, # infectious period duration = 1 week
    rho        = 0.5, # case detection rate
    sqrt_phi_inv = 1/sqrt(36))  # negative binomial overdispersion

true_pars_covid =
  c(R0         = 2.5, # basic reproduction number
    latent_dur = 5/7, # latent period duration = 1 week
    infec_dur  = 10/7, # infectious period duration = 1 week
    rho        = 0.1, # case detection rate
    sqrt_phi_inv = 1/sqrt(36))  # negative binomial overdispersion
 
# load results ------------------------------------------------------------

res = data.frame(
  expand.grid(replication = 1:1e3,
              param = c("R0", 
                        "latent_dur",
                        "infec_dur",
                        "rho",
                        "sqrt_phi_inv"),
              disease = c("ebola", "flu", "covid"),
              method = c("lna", "ode", "pmmh_1hr", "pmmh_1day"),
              rel_mad = 0,
              rel_ciw = 0,
              contraction = 0,
              coverage = 0,
              ess = 0,         # total ess
              ess_per_cpu = 0, # total ess / total cpu time
              psrf = 0))

comp_res = 
  data.frame(
    expand.grid(replication = 1:1e3,
                disease = c("ebola", "flu", "covid"),
                method = c("lna", "ode", "pmmh_1hr", "pmmh_1day")),
    cpu_time = 0,
    lp_ess = 0,
    lp_ess_per_cpu = 0,
    mpsrf = 0)

combns = 
  expand.grid(replication = 1:1e3,
              disease = c("ebola", "flu", "covid"),
              method = c("lna", "ode", "pmmh_1hr", "pmmh_1day"))

pmmh_1h_files = list.files("pmmh/results/")
pmmh_1d_files = list.files("pmmh/results_1day")
lna_files  = list.files("lna/results/")
ode_files  = list.files("ode/results/")
filenames = c(lna_files, ode_files, pmmh_1h_files, pmmh_1d_files)

for(i in seq_len(nrow(combns))) {
    
  # load fit summary
  res_file = paste0("seir_", combns$disease[i], "_", combns$replication[i], ".Rds")
  file_exists = 
    if(combns$method[i] == "lna") {
      res_file %in% lna_files
    } else if(combns$method[i] == "ode") {
      res_file %in% ode_files
    } else if(combns$method[i] == "pmmh_1hr") {
      res_file %in% pmmh_1h_files
    } else if(combns$method[i] == "pmmh_1day") {
      res_file %in% pmmh_1d_files
    }
  
  if(file_exists) {
    
    # load file
    fit = 
      if(combns$method[i] == "lna") {
        readRDS(paste0("lna/results/", res_file))
      } else if(combns$method[i] == "ode") {
        readRDS(paste0("ode/results/", res_file))  
      } else if(combns$method[i] == "pmmh_1hr") {
        readRDS(paste0("pmmh/results/", res_file))
      } else if(combns$method[i] == "pmmh_1day") {
        readRDS(paste0("pmmh/results_1day/", res_file))
      }
    
    quants = fit$post_quants
    
    # indices
    ind_R0 = with(res, which(replication == combns$replication[i] & 
                               disease == combns$disease[i] & method == combns$method[i] &
                               param == "R0"))
    
    ind_lat = with(res, which(replication == combns$replication[i] & 
                               disease == combns$disease[i] & method == combns$method[i] &
                               param == "latent_dur"))
    
    ind_inf = with(res, which(replication == combns$replication[i] & 
                                disease == combns$disease[i] & method == combns$method[i] &
                                param == "infec_dur"))
    
    ind_rho = with(res, which(replication == combns$replication[i] & 
                                disease == combns$disease[i] & method == combns$method[i] &
                                param == "rho"))
    
    ind_phi = with(res, which(replication == combns$replication[i] & 
                                disease == combns$disease[i] & method == combns$method[i] &
                                param == "sqrt_phi_inv"))
    
    # grab the right true parameters
    true_pars = 
      if(combns$disease[i] == "flu") {
        true_pars_flu
      } else if(combns$disease[i] == "ebola") {
        true_pars_ebola
      } else {
        true_pars_covid
      }
    
    # calc mads
    res$rel_mad[ind_R0]  = rel_mad("R0", quants, true_pars)
    res$rel_mad[ind_lat] = rel_mad("latent_dur", quants, true_pars)
    res$rel_mad[ind_inf] = rel_mad("infec_dur", quants, true_pars)
    res$rel_mad[ind_rho] = rel_mad("rho", quants, true_pars)
    res$rel_mad[ind_phi] = rel_mad("sqrt_phi_inv", quants, true_pars)
    
    # calc ciws
    res$rel_ciw[ind_R0]  = rel_ciw("R0", quants, true_pars)
    res$rel_ciw[ind_lat] = rel_ciw("latent_dur", quants, true_pars)
    res$rel_ciw[ind_inf] = rel_ciw("infec_dur", quants, true_pars)
    res$rel_ciw[ind_rho] = rel_ciw("rho", quants, true_pars)
    res$rel_ciw[ind_phi] = rel_ciw("sqrt_phi_inv", quants, true_pars)
    
    # calc contractions
    res$contraction[ind_R0]  = contraction("R0", quants, combns$disease[i])
    res$contraction[ind_lat] = contraction("latent_dur", quants, combns$disease[i])
    res$contraction[ind_inf] = contraction("infec_dur", quants, combns$disease[i])
    res$contraction[ind_rho] = contraction("rho", quants, combns$disease[i])
    res$contraction[ind_phi] = contraction("sqrt_phi_inv", quants, combns$disease[i])
    
    # calc coverages
    res$coverage[ind_R0]  = coverage("R0", quants, true_pars)
    res$coverage[ind_lat] = coverage("latent_dur", quants, true_pars)
    res$coverage[ind_inf] = coverage("infec_dur", quants, true_pars)
    res$coverage[ind_rho] = coverage("rho", quants, true_pars)
    res$coverage[ind_phi] = coverage("sqrt_phi_inv", quants, true_pars)
    
    # computational stuff
    res$ess[ind_R0]  = sum(fit$effective_samples[,"log_R0"])
    res$ess[ind_lat] = sum(fit$effective_samples[,"log_latent_dur"])
    res$ess[ind_inf] = sum(fit$effective_samples[,"log_infec_dur"])
    res$ess[ind_rho] = sum(fit$effective_samples[,"logit_rho"])
    res$ess[ind_phi] = sum(fit$effective_samples[,"log_sqrt_phi_inv"])
    
    res$ess_per_cpu[ind_R0]  = res$ess[ind_R0] / sum(fit$times)
    res$ess_per_cpu[ind_lat] = res$ess[ind_lat] / sum(fit$times)
    res$ess_per_cpu[ind_inf] = res$ess[ind_inf] / sum(fit$times)
    res$ess_per_cpu[ind_rho] = res$ess[ind_rho] / sum(fit$times)
    res$ess_per_cpu[ind_phi] = res$ess[ind_phi] / sum(fit$times)
    
    res$psrf[ind_R0]  = fit$psrf$psrf["log_R0", 1]
    res$psrf[ind_lat] = fit$psrf$psrf["log_latent_dur", 1]
    res$psrf[ind_inf] = fit$psrf$psrf["log_infec_dur", 1]
    res$psrf[ind_rho] = fit$psrf$psrf["logit_rho", 1]
    res$psrf[ind_phi] = fit$psrf$psrf["log_sqrt_phi_inv", 1]
    
    comp_res$mpsrf[i]    = fit$psrf$mpsrf
    comp_res$cpu_time[i] = sum(fit$times)
    comp_res$lp_ess[i]   = sum(fit$effective_samples[,"logpost"])
    comp_res$lp_ess_per_cpu[i] = 
      comp_res$lp_ess[i] / comp_res$cpu_time[i]
  }
}

# relabel factors
res$param = 
      factor(res$param, 
             levels = c("R0", "latent_dur", "infec_dur", "rho", "sqrt_phi_inv"),
             labels = c("R0", "Latent period", "Infectious period", "Detection rate", "N.B. Overdisperion"))

res$disease = 
  factor(res$disease, 
         levels = c("ebola", "flu", "covid"),
         labels = c("Ebola", "Influenza", "SARS-CoV-2"))

comp_res$disease = 
  factor(comp_res$disease, 
         levels = c("ebola", "flu", "covid"),
         labels = c("Ebola", "Influenza", "SARS-CoV-2"))

res$method = 
  factor(res$method, 
         levels = c("lna", "ode", "pmmh_1hr", "pmmh_1day"),
         labels = c("LNA", "ODE", "MMTL-1 hour", "MMTL-1 day"))

comp_res$method = 
  factor(comp_res$method, 
         levels = c("lna", "ode", "pmmh_1hr", "pmmh_1day"),
         labels = c("LNA", "ODE", "MMTL-1 hour", "MMTL-1 day"))

# get rid of pmmh iterations that didn't work
res = subset(res, contraction != 0)
comp_res = subset(comp_res, cpu_time != 0)

# res <- 
#   res %>% 
#   group_by(disease, param, method) %>% 
#   mutate(out_mad = label_outs(rel_mad),
#          out_ciw = label_outs(rel_ciw),
#          out_contr = label_outs(contraction)) %>% 
#   ungroup() %>% 
#   as.data.frame()

# summaries for point-ranges
res_sum <-
  res %>% 
  group_by(disease, param, method) %>% 
  summarise(mad_X50. = quantile(rel_mad, 0.5),
            mad_X10. = quantile(rel_mad, 0.1),
            mad_X90. = quantile(rel_mad, 0.9),
            ciw_X50. = quantile(rel_ciw, 0.5),
            ciw_X10. = quantile(rel_ciw, 0.1),
            ciw_X90. = quantile(rel_ciw, 0.9),
            contr_X50. = quantile(contraction, 0.5),
            contr_X10. = quantile(contraction, 0.1),
            contr_X90. = quantile(contraction, 0.9),
            covg_est = mean(coverage),
            covg_lower = BinomCI(x = sum(coverage), 1e3, method = "wilson")[2],
            covg_upper = BinomCI(x = sum(coverage), 1e3, method = "wilson")[3]) %>% 
  ungroup() %>% 
  as.data.frame()

rep_est = function(est, ci_lower, ci_upper, digs = 2) {
  # paste0(round(est, digs), " (", round(ci_lower, digs), ", ", round(ci_upper, digs), ")")
  paste0(signif(est, digs), " (", signif(ci_lower, digs), ", ", signif(ci_upper, digs), ")")
}

comp_res_sum <-
  comp_res %>% 
  group_by(disease, method) %>% 
  summarise(cpu_time_X50. = quantile(cpu_time, 0.5),
            cpu_time_X10. = quantile(cpu_time, 0.1),
            cpu_time_X90. = quantile(cpu_time, 0.9),
            lp_ess_cpu_X50. = quantile(lp_ess_per_cpu / cpu_time, 0.5),
            lp_ess_cpu_X10. = quantile(lp_ess_per_cpu / cpu_time, 0.1),
            lp_ess_cpu_X90. = quantile(lp_ess_per_cpu / cpu_time, 0.9)) %>% 
  mutate(cpu_time = 
           rep_est(cpu_time_X50., cpu_time_X10., cpu_time_X90.),
         lp_ess_per_cpu_time = 
           rep_est(lp_ess_cpu_X50., lp_ess_cpu_X10., lp_ess_cpu_X90.)) %>% 
  select(-c(3:8)) %>% 
  ungroup() %>% 
  as.data.frame() 

comp_res_tab = 
  rbind(t(subset(comp_res_sum, disease == "Ebola"))[-1,],
        t(subset(comp_res_sum, disease == "Influenza"))[-1,],
        t(subset(comp_res_sum, disease == "SARS-CoV-2"))[-1,])

kable(comp_res_tab, 
      format = "latex",
      booktabs = TRUE) %>% 
  pack_rows("Ebola", 1, 3) %>% 
  pack_rows("Influenza", 4, 6) %>% 
  pack_rows("SARS-CoV-2", 7, 9)

# plots by disease --------------------------------------------------------
labs = c(expression(R[0]), 
         expression(1 / omega), 
         expression(1 / mu), 
         expression(rho), 
         expression(1 / sqrt(phi)))

mad_plot <-
  res_sum %>% 
  ggplot(aes(x = param, y = mad_X50., ymin = mad_X10., ymax = mad_X90.,
             colour = method, shape = method)) + 
  geom_pointrange(position = position_dodge(width = 0.8)) + 
  scale_x_discrete(labels = labs) + 
  # scale_y_continuous(limits = c(0, 2)) +
  facet_grid(.~disease, scales = "free") + 
  # scale_fill_brewer(type = "qual", palette = 7) + 
  # scale_color_brewer(type = "qual", palette = 7) + 
  scale_color_colorblind() +
  theme_minimal() + 
  theme(legend.position = "bottom",
        text = element_text(size = 20, family = "serif"),
        axis.text.x = element_text(vjust = 0.25),
        strip.text = element_text(size = 20),
        panel.spacing = unit(1.25, "lines"),
        legend.title = element_blank()) + 
  labs(x = NULL, y = expression("|" ~ hat(theta)[0.5] - theta[true] ~ "|" / theta[true]), title = "Scaled median absolute deviation")


ciw_plot <- res_sum %>% 
  ggplot(aes(x = param, y = ciw_X50., ymin = ciw_X10., ymax = ciw_X90.,
             colour = method, shape = method)) + 
  geom_pointrange(position = position_dodge(width = 0.8)) + 
  facet_grid(.~disease, scales = "free") + 
  # scale_fill_brewer(type = "qual", palette = 7) + 
  # scale_color_brewer(type = "qual", palette = 7) + 
  scale_color_colorblind() +
  scale_x_discrete(labels = labs) + 
  theme_minimal() + 
  theme(legend.position = "bottom", text = element_text(size = 20, family = "serif"),
        axis.text.x = element_text(vjust = 0.25),
        strip.text = element_text(size = 20),
        panel.spacing = unit(1.25, "lines"),
        legend.title = element_blank()) + 
  labs(x = NULL, y = expression((hat(theta)[0.975] - hat(theta)[0.025]) / theta["true"]), title = "Scaled credible interval width") 


contr_plot <- 
  res_sum %>% 
  ggplot(aes(x = param, y = contr_X50., ymin = contr_X10., ymax = contr_X90.,
             colour = method, shape = method)) + 
  geom_pointrange(position = position_dodge(width = 0.8)) + 
  facet_grid(.~disease, scales = "free") + 
  scale_x_discrete(labels = labs) + 
  # scale_fill_brewer(type = "qual", palette = 7) + 
  # scale_color_brewer(type = "qual", palette = 7) + 
  scale_color_colorblind() +
  theme_minimal() + 
  theme(legend.position = "bottom", text = element_text(size = 20, family = "serif"),
        axis.text.x = element_text(vjust = 0.25),
        strip.text = element_text(size = 20),
        panel.spacing = unit(1.25, "lines"),
        legend.title = element_blank()) + 
  labs(x = NULL, y = "Posterior contraction", title = "Posterior 80% CI width vs. prior 80% CI width ")

covg_plot <-
res_sum %>% 
  ggplot(aes(x = param, y = covg_est, ymin = covg_lower, ymax = covg_upper,
             colour = method, shape = method)) + 
  # geom_pointrange(position = position_dodge(width = 0.8)) + 
  geom_point(position = position_dodge(width = 0.8), size = 2) + 
  geom_hline(yintercept = 0.8, linetype = "dashed", colour = "grey50") + 
  # scale_y_continuous(limits = c(0.9, 1)) + 
  scale_x_discrete(labels = labs) + 
  facet_grid(.~disease, scales = "free") + 
  # scale_fill_brewer(type = "qual", palette = 7) + 
  # scale_color_brewer(type = "qual", palette = 7) +
  # scale_color_tableau(palette = "Classic Color Blind") + 
  scale_color_colorblind() +
  theme_minimal() + 
  theme(legend.position = "bottom", text = element_text(size = 20, family = "serif"),
        axis.text.x = element_text(vjust = 0.25),
        strip.text = element_text(size = 20),
        panel.spacing = unit(1.25, "lines"),
        legend.title = element_blank()) + 
  labs(x = NULL, y = "Coverage", title = "Coverage of 80% credible intervals")

all_plots <- covg_plot / contr_plot / mad_plot / ciw_plot

ggsave("threedis_res_long.pdf", plot = all_plots, width = 11, height = 13)
