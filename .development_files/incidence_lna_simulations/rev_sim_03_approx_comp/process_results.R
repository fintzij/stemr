library(tidyverse)
library(DescTools)
library(patchwork)
library(here)
setwd(here("/paper_V2/simulations/rev_sim_03_approx_comp/sir"))
source(here("/paper_V2/simulations/rev_sim_03_approx_comp/sir/process_fcns.R"))
 
# load results ------------------------------------------------------------

pmmh_files = list.files("results/pmmh/")
lna_files  = list.files("results/lna/")
ode_files  = list.files("results/ode/")
  
sim_settings = 
  expand.grid(param = c("R0", 
                        "infec_dur",
                        "rho",
                        "sqrt_phi_inv"),
              replication = 1:1e3,
              popsize = c(2e3, 1e4, 5e4))

res_ode = res_lna = res_pmmh =
    data.frame(sim_settings[,c("popsize", "param", "replication")],
               rel_mad = 0,
               rel_ciw = 0,
               contraction = 0,
               coverage = 0,
               ess = 0,         # total ess
               ess_per_cpu = 0, # total ess / total cpu time
               psrf = 0)

comp_res_ode = comp_res_lna = comp_res_pmmh = 
  data.frame(expand.grid(replication = 1:1e3,
                         popsize = c(2e3, 1e4, 5e4)),
             cpu_time = 0,
             lp_ess = 0,
             lp_ess_per_cpu = 0,
             mpsrf = 0,
             present = FALSE)

for(p in c(2e3, 1e4, 5e4)) {
  for(r in 1:1e3) {
    
    file_lna = paste0("sir_lna_popsize_", p, "_", r, ".Rds")
    file_ode = paste0("sir_ode_popsize_", p, "_", r, ".Rds")
    file_pmmh = paste0("sir_pmmh_popsize_", p, "_", r, ".Rds")
    
    if(file_lna %in% lna_files) {
      
      # load results
      fit_lna = readRDS(paste0("results/lna/",file_lna)) 
      
      # get indices
      comp_ind = with(comp_res_lna, which(popsize == p & replication == r))
      ind_R0  = with(res_lna, which(popsize == p & replication == r & param == "R0"))
      ind_inf = with(res_lna, which(popsize == p & replication == r & param == "infec_dur"))
      ind_rho = with(res_lna, which(popsize == p & replication == r & param == "rho"))
      ind_phi = with(res_lna, which(popsize == p & replication == r & param == "sqrt_phi_inv"))
      
      # record that mcmc was successful
      comp_res_lna$present[comp_ind] = TRUE
      
      ### get results
      # MADs
      res_lna$rel_mad[ind_R0]  = rel_mad("R0", fit_lna$post_quants, fit_lna$true_pars)
      res_lna$rel_mad[ind_inf] = rel_mad("infec_dur", fit_lna$post_quants, fit_lna$true_pars)
      res_lna$rel_mad[ind_rho] = rel_mad("rho", fit_lna$post_quants, fit_lna$true_pars)
      res_lna$rel_mad[ind_phi] = rel_mad("sqrt_phi_inv", fit_lna$post_quants, fit_lna$true_pars)
      
      # ciws
      res_lna$rel_ciw[ind_R0]  = rel_ciw("R0", fit_lna$post_quants, fit_lna$true_pars)
      res_lna$rel_ciw[ind_inf] = rel_ciw("infec_dur", fit_lna$post_quants, fit_lna$true_pars)
      res_lna$rel_ciw[ind_rho] = rel_ciw("rho", fit_lna$post_quants, fit_lna$true_pars)
      res_lna$rel_ciw[ind_phi] = rel_ciw("sqrt_phi_inv", fit_lna$post_quants, fit_lna$true_pars)
      
      # contractions
      res_lna$contraction[ind_R0]  = contraction("R0", fit_lna$post_quants)
      res_lna$contraction[ind_inf] = contraction("infec_dur", fit_lna$post_quants)
      res_lna$contraction[ind_rho] = contraction("rho", fit_lna$post_quants)
      res_lna$contraction[ind_phi] = contraction("sqrt_phi_inv", fit_lna$post_quants)
      
      # coverage
      res_lna$coverage[ind_R0]  = coverage("R0", fit_lna$post_quants, fit_lna$true_pars)
      res_lna$coverage[ind_inf] = coverage("infec_dur", fit_lna$post_quants, fit_lna$true_pars)
      res_lna$coverage[ind_rho] = coverage("rho", fit_lna$post_quants, fit_lna$true_pars)
      res_lna$coverage[ind_phi] = coverage("sqrt_phi_inv", fit_lna$post_quants, fit_lna$true_pars)
      
      # computational stuff
      res_lna$ess[ind_R0]  = sum(fit_lna$effective_samples[,"log_R0_m1"])
      res_lna$ess[ind_inf] = sum(fit_lna$effective_samples[,"log_infec_dur"])
      res_lna$ess[ind_rho] = sum(fit_lna$effective_samples[,"logit_rho"])
      res_lna$ess[ind_phi] = sum(fit_lna$effective_samples[,"log_sqrt_phi_inv"])
      
      res_lna$ess_per_cpu[ind_R0]  = res_lna$ess[ind_R0] / sum(fit_lna$times)
      res_lna$ess_per_cpu[ind_inf] = res_lna$ess[ind_inf] / sum(fit_lna$times)
      res_lna$ess_per_cpu[ind_rho] = res_lna$ess[ind_rho] / sum(fit_lna$times)
      res_lna$ess_per_cpu[ind_phi] = res_lna$ess[ind_phi] / sum(fit_lna$times)
      
      res_lna$psrf[ind_R0]  = fit_lna$psrf$psrf["log_R0_m1", 1]
      res_lna$psrf[ind_inf] = fit_lna$psrf$psrf["log_infec_dur", 1]
      res_lna$psrf[ind_rho] = fit_lna$psrf$psrf["logit_rho", 1]
      res_lna$psrf[ind_phi] = fit_lna$psrf$psrf["log_sqrt_phi_inv", 1]
      
      comp_res_lna$mpsrf[comp_ind]    = fit_lna$psrf$mpsrf
      comp_res_lna$cpu_time[comp_ind] = sum(fit_lna$times)
      comp_res_lna$lp_ess[comp_ind]   = sum(fit_lna$effective_samples[,"logpost"])
      comp_res_lna$lp_ess_per_cpu[comp_ind] = 
        comp_res_lna$lp_ess[comp_ind] / comp_res_lna$cpu_time[comp_ind]
    }
    
    if(file_ode %in% ode_files) {
      
      # load results
      fit_ode = readRDS(paste0("results/ode/",file_ode)) 
      
      # get indices
      comp_ind = with(comp_res_ode, which(popsize == p & replication == r))
      ind_R0  = with(res_ode, which(popsize == p & replication == r & param == "R0"))
      ind_inf = with(res_ode, which(popsize == p & replication == r & param == "infec_dur"))
      ind_rho = with(res_ode, which(popsize == p & replication == r & param == "rho"))
      ind_phi = with(res_ode, which(popsize == p & replication == r & param == "sqrt_phi_inv"))
      
      # record that mcmc was successful
      comp_res_ode$present[comp_ind] = TRUE
      
      ### get results
      # MADs
      res_ode$rel_mad[ind_R0]  = rel_mad("R0", fit_ode$post_quants, fit_ode$true_pars)
      res_ode$rel_mad[ind_inf] = rel_mad("infec_dur", fit_ode$post_quants, fit_ode$true_pars)
      res_ode$rel_mad[ind_rho] = rel_mad("rho", fit_ode$post_quants, fit_ode$true_pars)
      res_ode$rel_mad[ind_phi] = rel_mad("sqrt_phi_inv", fit_ode$post_quants, fit_ode$true_pars)
      
      # ciws
      res_ode$rel_ciw[ind_R0]  = rel_ciw("R0", fit_ode$post_quants, fit_ode$true_pars)
      res_ode$rel_ciw[ind_inf] = rel_ciw("infec_dur", fit_ode$post_quants, fit_ode$true_pars)
      res_ode$rel_ciw[ind_rho] = rel_ciw("rho", fit_ode$post_quants, fit_ode$true_pars)
      res_ode$rel_ciw[ind_phi] = rel_ciw("sqrt_phi_inv", fit_ode$post_quants, fit_ode$true_pars)
      
      # contractions
      res_ode$contraction[ind_R0]  = contraction("R0", fit_ode$post_quants)
      res_ode$contraction[ind_inf] = contraction("infec_dur", fit_ode$post_quants)
      res_ode$contraction[ind_rho] = contraction("rho", fit_ode$post_quants)
      res_ode$contraction[ind_phi] = contraction("sqrt_phi_inv", fit_ode$post_quants)
      
      # coverage
      res_ode$coverage[ind_R0]  = coverage("R0", fit_ode$post_quants, fit_ode$true_pars)
      res_ode$coverage[ind_inf] = coverage("infec_dur", fit_ode$post_quants, fit_ode$true_pars)
      res_ode$coverage[ind_rho] = coverage("rho", fit_ode$post_quants, fit_ode$true_pars)
      res_ode$coverage[ind_phi] = coverage("sqrt_phi_inv", fit_ode$post_quants, fit_ode$true_pars)
      
      # computational stuff
      res_ode$ess[ind_R0]  = sum(fit_ode$effective_samples[,"log_R0_m1"])
      res_ode$ess[ind_inf] = sum(fit_ode$effective_samples[,"log_infec_dur"])
      res_ode$ess[ind_rho] = sum(fit_ode$effective_samples[,"logit_rho"])
      res_ode$ess[ind_phi] = sum(fit_ode$effective_samples[,"log_sqrt_phi_inv"])
      
      res_ode$ess_per_cpu[ind_R0]  = res_ode$ess[ind_R0] / sum(fit_ode$times)
      res_ode$ess_per_cpu[ind_inf] = res_ode$ess[ind_inf] / sum(fit_ode$times)
      res_ode$ess_per_cpu[ind_rho] = res_ode$ess[ind_rho] / sum(fit_ode$times)
      res_ode$ess_per_cpu[ind_phi] = res_ode$ess[ind_phi] / sum(fit_ode$times)
      
      res_ode$psrf[ind_R0]  = fit_ode$psrf$psrf["log_R0_m1", 1]
      res_ode$psrf[ind_inf] = fit_ode$psrf$psrf["log_infec_dur", 1]
      res_ode$psrf[ind_rho] = fit_ode$psrf$psrf["logit_rho", 1]
      res_ode$psrf[ind_phi] = fit_ode$psrf$psrf["log_sqrt_phi_inv", 1]
      
      comp_res_ode$mpsrf[comp_ind]    = fit_ode$psrf$mpsrf
      comp_res_ode$cpu_time[comp_ind] = sum(fit_ode$times)
      comp_res_ode$lp_ess[comp_ind]   = sum(fit_ode$effective_samples[,"logpost"])
      comp_res_ode$lp_ess_per_cpu[comp_ind] = 
        comp_res_ode$lp_ess[comp_ind] / comp_res_ode$cpu_time[comp_ind]
    }
    
    if(file_pmmh %in% pmmh_files) {
      
      # load results
      fit_pmmh = readRDS(paste0("results/pmmh/",file_pmmh)) 
      
      # get indices
      comp_ind = with(comp_res_pmmh, which(popsize == p & replication == r))
      ind_R0  = with(res_pmmh, which(popsize == p & replication == r & param == "R0"))
      ind_inf = with(res_pmmh, which(popsize == p & replication == r & param == "infec_dur"))
      ind_rho = with(res_pmmh, which(popsize == p & replication == r & param == "rho"))
      ind_phi = with(res_pmmh, which(popsize == p & replication == r & param == "sqrt_phi_inv"))
      
      # record that mcmc was successful
      comp_res_pmmh$present[comp_ind] = TRUE
      
      ### get results
      # MADs
      res_pmmh$rel_mad[ind_R0]  = rel_mad("R0", fit_pmmh$post_quants, fit_pmmh$true_pars)
      res_pmmh$rel_mad[ind_inf] = rel_mad("infec_dur", fit_pmmh$post_quants, fit_pmmh$true_pars)
      res_pmmh$rel_mad[ind_rho] = rel_mad("rho", fit_pmmh$post_quants, fit_pmmh$true_pars)
      res_pmmh$rel_mad[ind_phi] = rel_mad("sqrt_phi_inv", fit_pmmh$post_quants, fit_pmmh$true_pars)
      
      # ciws
      res_pmmh$rel_ciw[ind_R0]  = rel_ciw("R0", fit_pmmh$post_quants, fit_pmmh$true_pars)
      res_pmmh$rel_ciw[ind_inf] = rel_ciw("infec_dur", fit_pmmh$post_quants, fit_pmmh$true_pars)
      res_pmmh$rel_ciw[ind_rho] = rel_ciw("rho", fit_pmmh$post_quants, fit_pmmh$true_pars)
      res_pmmh$rel_ciw[ind_phi] = rel_ciw("sqrt_phi_inv", fit_pmmh$post_quants, fit_pmmh$true_pars)
      
      # contractions
      res_pmmh$contraction[ind_R0]  = contraction("R0", fit_pmmh$post_quants)
      res_pmmh$contraction[ind_inf] = contraction("infec_dur", fit_pmmh$post_quants)
      res_pmmh$contraction[ind_rho] = contraction("rho", fit_pmmh$post_quants)
      res_pmmh$contraction[ind_phi] = contraction("sqrt_phi_inv", fit_pmmh$post_quants)
      
      # coverage
      res_pmmh$coverage[ind_R0]  = coverage("R0", fit_pmmh$post_quants, fit_pmmh$true_pars)
      res_pmmh$coverage[ind_inf] = coverage("infec_dur", fit_pmmh$post_quants, fit_pmmh$true_pars)
      res_pmmh$coverage[ind_rho] = coverage("rho", fit_pmmh$post_quants, fit_pmmh$true_pars)
      res_pmmh$coverage[ind_phi] = coverage("sqrt_phi_inv", fit_pmmh$post_quants, fit_pmmh$true_pars)
      
      # computational stuff
      res_pmmh$ess[ind_R0]  = sum(fit_pmmh$effective_samples[,"log_R0_m1"])
      res_pmmh$ess[ind_inf] = sum(fit_pmmh$effective_samples[,"log_infec_dur"])
      res_pmmh$ess[ind_rho] = sum(fit_pmmh$effective_samples[,"logit_rho"])
      res_pmmh$ess[ind_phi] = sum(fit_pmmh$effective_samples[,"log_sqrt_phi_inv"])
      
      res_pmmh$ess_per_cpu[ind_R0]  = res_pmmh$ess[ind_R0] / sum(fit_pmmh$times)
      res_pmmh$ess_per_cpu[ind_inf] = res_pmmh$ess[ind_inf] / sum(fit_pmmh$times)
      res_pmmh$ess_per_cpu[ind_rho] = res_pmmh$ess[ind_rho] / sum(fit_pmmh$times)
      res_pmmh$ess_per_cpu[ind_phi] = res_pmmh$ess[ind_phi] / sum(fit_pmmh$times)
      
      res_pmmh$psrf[ind_R0]  = fit_pmmh$psrf$psrf["log_R0_m1", 1]
      res_pmmh$psrf[ind_inf] = fit_pmmh$psrf$psrf["log_infec_dur", 1]
      res_pmmh$psrf[ind_rho] = fit_pmmh$psrf$psrf["logit_rho", 1]
      res_pmmh$psrf[ind_phi] = fit_pmmh$psrf$psrf["log_sqrt_phi_inv", 1]
      
      comp_res_pmmh$mpsrf[comp_ind]    = fit_pmmh$psrf$mpsrf
      comp_res_pmmh$cpu_time[comp_ind] = sum(fit_pmmh$times)
      comp_res_pmmh$lp_ess[comp_ind]   = sum(fit_pmmh$effective_samples[,"logpost"])
      comp_res_pmmh$lp_ess_per_cpu[comp_ind] = 
        comp_res_pmmh$lp_ess[comp_ind] / comp_res_pmmh$cpu_time[comp_ind]
    }
  }
}

res = data.frame(Method = rep(c("MMTL", "LNA", "ODE"), each = nrow(res_lna)),
                 rbind(res_pmmh, res_lna, res_ode))
comp_res = data.frame(Method = rep(c("MMTL", "LNA", "ODE"), each = nrow(comp_res_lna)),
                      rbind(comp_res_pmmh, comp_res_lna, comp_res_ode))

res$param = 
      factor(res$param, 
             levels = c("R0", "infec_dur", "rho", "sqrt_phi_inv"),
             labels = c("R0", "Infectious period", "Detection rate", "N.B. Overdispersion"))

res$popsize = 
  factor(as.character(res$popsize), 
            levels = c(2e3, 1e4, 5e4),
         labels = c("P==2000; I(t[0])==1",
                    "P==10000; I(t[0])==5",
                    "P==50000; I(t[0])==25"))

res$Method = factor(res$Method, levels = c("LNA", "ODE", "MMTL"))
comp_res$Method = factor(comp_res$Method, levels = c("LNA", "ODE", "MMTL"))

# get rid of the replications where MCMC didn't work
res = subset(res, psrf != 0)
comp_res = subset(comp_res, present == T)

# generate outlier labels
# label_outs = function(x) {
#   quants = quantile(x, c(0.1, 0.9))
#   upper_lim = quants[2] + 1.5 * abs(diff(quants))
#   lower_lim = quants[1] - 1.5 * abs(diff(quants))
#   out_lab = ifelse(x > upper_lim | x < lower_lim, 1, 0)
#   return(out_lab)
# }
# 
# res <- 
#   res %>% 
#   group_by(Method, param) %>% 
#   mutate(out_mad = label_outs(rel_mad),
#          out_ciw = label_outs(rel_ciw),
#          out_contr = label_outs(contraction)) %>% 
#   ungroup() %>% 
#   as.data.frame()

res_sum <-
  res %>% 
  group_by(popsize, param, Method) %>% 
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
  group_by(popsize, Method) %>% 
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
  
# generate plots ----------------------------------------------------------

labs = c(expression(R[0]), 
         expression(1 / mu), 
         expression(rho), 
         expression(1 / sqrt(phi)))

mad_plot <-
  res_sum %>% 
  ggplot(aes(x = param, y = mad_X50., ymin = mad_X10., ymax = mad_X90.,
             colour = Method, shape = Method)) + 
  geom_pointrange(position = position_dodge(width = 0.8)) + 
  scale_x_discrete(labels = labs) + 
  # scale_y_continuous(limits = c(0, 2)) +
  facet_grid(.~popsize, scales = "free", 
             labeller = label_parsed) + 
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
             colour = Method, shape = Method)) + 
  geom_pointrange(position = position_dodge(width = 0.8)) + 
  facet_grid(.~popsize, scales = "free",
             labeller = label_parsed) + 
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
             colour = Method, shape = Method)) + 
  geom_pointrange(position = position_dodge(width = 0.8)) + 
  facet_grid(.~popsize, scales = "free",
             labeller = label_parsed) + 
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
             colour = Method, shape = Method)) + 
  # geom_pointrange(position = position_dodge(width = 0.8)) + 
  geom_point(position = position_dodge(width = 0.8), size = 2) + 
  geom_hline(yintercept = 0.8, linetype = "dashed", colour = "grey50") + 
  # scale_y_continuous(limits = c(0.9, 1)) + 
  scale_x_discrete(labels = labs) + 
  facet_grid(.~popsize, scales = "free",
             labeller = label_parsed) + 
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

ggsave("sir_covg.pdf", plot = all_plots, width = 11, height = 13)
