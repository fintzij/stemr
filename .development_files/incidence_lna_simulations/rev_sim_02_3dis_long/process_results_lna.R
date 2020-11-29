library(tidyverse)
library(DescTools)
library(patchwork)
library(here)
setwd(here("paper_V2/simulations/rev_sim_02_3dis_diff"))
source(here("paper_V2/simulations/rev_sim_02_3dis_diff/process_fcns.R"))

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

res_flu = res_ebola = res_covid = 
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

for(i in seq_len(1e3)) {
    
  # load fit summary
  fit_flu_i   = readRDS(paste0("lna/results/seir_flu_",i,".Rds"))
  fit_covid_i = readRDS(paste0("lna/results/seir_covid_",i,".Rds"))
  fit_ebola_i = readRDS(paste0("lna/results/seir_ebola_",i,".Rds"))
  
  # get indices
  ind_R0  = with(res_flu, which(replication == i & param == "R0"))
  ind_lat = with(res_flu, which(replication == i & param == "latent_dur"))
  ind_inf = with(res_flu, which(replication == i & param == "infec_dur"))
  ind_rho = with(res_flu, which(replication == i & param == "rho"))
  ind_phi = with(res_flu, which(replication == i & param == "sqrt_phi_inv"))
  
  ### FLU results
  # calc MADs
  res_flu$rel_mad[ind_R0]  = rel_mad("R0", fit_flu_i, true_pars_flu)
  res_flu$rel_mad[ind_lat] = rel_mad("latent_dur", fit_flu_i, true_pars_flu)
  res_flu$rel_mad[ind_inf] = rel_mad("infec_dur", fit_flu_i, true_pars_flu)
  res_flu$rel_mad[ind_rho] = rel_mad("rho", fit_flu_i, true_pars_flu)
  res_flu$rel_mad[ind_phi] = rel_mad("sqrt_phi_inv", fit_flu_i, true_pars_flu)
  
  # calc ciws
  res_flu$rel_ciw[ind_R0]  = rel_ciw("R0", fit_flu_i, true_pars_flu)
  res_flu$rel_ciw[ind_lat] = rel_ciw("latent_dur", fit_flu_i, true_pars_flu)
  res_flu$rel_ciw[ind_inf] = rel_ciw("infec_dur", fit_flu_i, true_pars_flu)
  res_flu$rel_ciw[ind_rho] = rel_ciw("rho", fit_flu_i, true_pars_flu)
  res_flu$rel_ciw[ind_phi] = rel_ciw("sqrt_phi_inv", fit_flu_i, true_pars_flu)
  
  # calc contraction
  res_flu$contraction[ind_R0]  = contraction("R0", fit_flu_i, "flu")
  res_flu$contraction[ind_lat] = contraction("latent_dur", fit_flu_i, "flu")
  res_flu$contraction[ind_inf] = contraction("infec_dur", fit_flu_i, "flu")
  res_flu$contraction[ind_rho] = contraction("rho", fit_flu_i, "flu")
  res_flu$contraction[ind_phi] = contraction("sqrt_phi_inv", fit_flu_i, "flu")
  
  # calc coverage
  res_flu$coverage[ind_R0]  = coverage("R0", fit_flu_i, true_pars_flu)
  res_flu$coverage[ind_lat] = coverage("latent_dur", fit_flu_i, true_pars_flu)
  res_flu$coverage[ind_inf] = coverage("infec_dur", fit_flu_i, true_pars_flu)
  res_flu$coverage[ind_rho] = coverage("rho", fit_flu_i, true_pars_flu)
  res_flu$coverage[ind_phi] = coverage("sqrt_phi_inv", fit_flu_i, true_pars_flu)
  
  ### COVID results
  # calc MADs
  res_covid$rel_mad[ind_R0]  = rel_mad("R0", fit_covid_i, true_pars_covid)
  res_covid$rel_mad[ind_lat] = rel_mad("latent_dur", fit_covid_i, true_pars_covid)
  res_covid$rel_mad[ind_inf] = rel_mad("infec_dur", fit_covid_i, true_pars_covid)
  res_covid$rel_mad[ind_rho] = rel_mad("rho", fit_covid_i, true_pars_covid)
  res_covid$rel_mad[ind_phi] = rel_mad("sqrt_phi_inv", fit_covid_i, true_pars_covid)
  
  # calc ciws
  res_covid$rel_ciw[ind_R0]  = rel_ciw("R0", fit_covid_i, true_pars_covid)
  res_covid$rel_ciw[ind_lat] = rel_ciw("latent_dur", fit_covid_i, true_pars_covid)
  res_covid$rel_ciw[ind_inf] = rel_ciw("infec_dur", fit_covid_i, true_pars_covid)
  res_covid$rel_ciw[ind_rho] = rel_ciw("rho", fit_covid_i, true_pars_covid)
  res_covid$rel_ciw[ind_phi] = rel_ciw("sqrt_phi_inv", fit_covid_i, true_pars_covid)
  
  # calc contraction
  res_covid$contraction[ind_R0]  = contraction("R0", fit_covid_i, "covid")
  res_covid$contraction[ind_lat] = contraction("latent_dur", fit_covid_i, "covid")
  res_covid$contraction[ind_inf] = contraction("infec_dur", fit_covid_i, "covid")
  res_covid$contraction[ind_rho] = contraction("rho", fit_covid_i, "covid")
  res_covid$contraction[ind_phi] = contraction("sqrt_phi_inv", fit_covid_i, "covid")
  
  # calc coverage
  res_covid$coverage[ind_R0]  = coverage("R0", fit_covid_i, true_pars_covid)
  res_covid$coverage[ind_lat] = coverage("latent_dur", fit_covid_i, true_pars_covid)
  res_covid$coverage[ind_inf] = coverage("infec_dur", fit_covid_i, true_pars_covid)
  res_covid$coverage[ind_rho] = coverage("rho", fit_covid_i, true_pars_covid)
  res_covid$coverage[ind_phi] = coverage("sqrt_phi_inv", fit_covid_i, true_pars_covid)
  
  ### EBOLA results
  # calc MADs
  res_ebola$rel_mad[ind_R0]  = rel_mad("R0", fit_ebola_i, true_pars_ebola)
  res_ebola$rel_mad[ind_lat] = rel_mad("latent_dur", fit_ebola_i, true_pars_ebola)
  res_ebola$rel_mad[ind_inf] = rel_mad("infec_dur", fit_ebola_i, true_pars_ebola)
  res_ebola$rel_mad[ind_rho] = rel_mad("rho", fit_ebola_i, true_pars_ebola)
  res_ebola$rel_mad[ind_phi] = rel_mad("sqrt_phi_inv", fit_ebola_i, true_pars_ebola)
  
  # calc ciws
  res_ebola$rel_ciw[ind_R0]  = rel_ciw("R0", fit_ebola_i, true_pars_ebola)
  res_ebola$rel_ciw[ind_lat] = rel_ciw("latent_dur", fit_ebola_i, true_pars_ebola)
  res_ebola$rel_ciw[ind_inf] = rel_ciw("infec_dur", fit_ebola_i, true_pars_ebola)
  res_ebola$rel_ciw[ind_rho] = rel_ciw("rho", fit_ebola_i, true_pars_ebola)
  res_ebola$rel_ciw[ind_phi] = rel_ciw("sqrt_phi_inv", fit_ebola_i, true_pars_ebola)
  
  # calc contraction
  res_ebola$contraction[ind_R0]  = contraction("R0", fit_ebola_i, "ebola")
  res_ebola$contraction[ind_lat] = contraction("latent_dur", fit_ebola_i, "ebola")
  res_ebola$contraction[ind_inf] = contraction("infec_dur", fit_ebola_i, "ebola")
  res_ebola$contraction[ind_rho] = contraction("rho", fit_ebola_i, "ebola")
  res_ebola$contraction[ind_phi] = contraction("sqrt_phi_inv", fit_ebola_i, "ebola")
  
  # calc coverage
  res_ebola$coverage[ind_R0]  = coverage("R0", fit_ebola_i, true_pars_ebola)
  res_ebola$coverage[ind_lat] = coverage("latent_dur", fit_ebola_i, true_pars_ebola)
  res_ebola$coverage[ind_inf] = coverage("infec_dur", fit_ebola_i, true_pars_ebola)
  res_ebola$coverage[ind_rho] = coverage("rho", fit_ebola_i, true_pars_ebola)
  res_ebola$coverage[ind_phi] = coverage("sqrt_phi_inv", fit_ebola_i, true_pars_ebola)
}

res = data.frame(Pathogen = rep(c("Ebola", "Influenza", "SARS-CoV-2"), each = 5e3),
                 rbind(res_ebola, res_flu, res_covid))

res$param = 
      factor(res$param, 
             levels = c("R0", "latent_dur", "infec_dur", "rho", "sqrt_phi_inv"),
             labels = c("R0", "Latent period", "Infectious period", "Detection rate", "N.B. Overdisperion"))

res <- 
  res %>% 
  group_by(Pathogen, param) %>% 
  mutate(out_mad = label_outs(rel_mad),
         out_ciw = label_outs(rel_ciw),
         out_contr = label_outs(contraction)) %>% 
  ungroup() %>% 
  as.data.frame()
  
# generate plots ----------------------------------------------------------

mad_plot <-
    subset(res, out_mad == 0) %>% 
    ggplot(aes(x = as.factor(Pathogen), y = rel_mad, 
               fill = as.factor(Pathogen))) + 
    geom_violin(alpha = 0.5, draw_quantiles = c(0.1, 0.5, 0.9), colour = "grey50", trim = T, scale = "width") +
    geom_jitter(data = subset(res, out_mad == 1),
               aes(x = as.factor(Pathogen), y = rel_mad, colour = as.factor(Pathogen)),
               width = 0.1, alpha = 0.5, colour = "grey50") + 
    scale_x_discrete(labels = c("Ebola", "Influenza", "SARS-CoV-2")) + 
    scale_y_continuous(limits = c(0, 3)) +
    facet_grid(.~param, scales = "free") + 
    scale_fill_brewer(type = "qual") + 
    scale_color_brewer(type = "qual") + 
    theme_minimal() + 
    theme(legend.position = "none",
          text = element_text(size = 15, family = "serif"),
          axis.text.x = element_text(angle = 320, hjust = 0, vjust = 1),
          plot.margin = margin(10, 50, 10, 10)) + 
    labs(x = NULL, y = expression("|" ~ hat(theta)[0.5] - theta[true] ~ "|" / theta[true]), title = "Scaled median absolute deviation")

ciw_plot <-
    subset(res, out_ciw == 0) %>% 
    ggplot(aes(x = as.factor(Pathogen), y = rel_ciw, 
               fill = as.factor(Pathogen))) + 
    geom_violin(alpha = 0.5, draw_quantiles = c(0.1, 0.5, 0.9), colour = "grey50", scale = "width") +
    geom_jitter(data = subset(res, out_ciw == 1),
              aes(x = as.factor(Pathogen), y = rel_ciw, colour = as.factor(Pathogen)),
              width = 0.1, alpha = 0.5, colour = "grey50") + 
    facet_grid(.~param, scales = "free") + 
    scale_fill_brewer(type = "qual") + 
    scale_color_brewer(type = "qual") + 
    scale_x_discrete(labels = c("Ebola", "Influenza", "SARS-CoV-2")) + 
    scale_y_continuous(limits = c(0, 9)) +
    theme_minimal() + 
    theme(legend.position = "none", text = element_text(size = 15, family = "serif"),
          axis.text.x = element_text(angle = 320, hjust = 0, vjust = 1),
          plot.margin = margin(10, 50, 10, 10)) + 
    labs(x = NULL, y = expression((hat(theta)[0.975] - hat(theta)[0.025]) / theta["true"]), title = "Scaled credible interval width") 

contraction_plot <-
    subset(res, out_contr == 0) %>% 
    ggplot(aes(x = as.factor(Pathogen), y = contraction, 
               fill = as.factor(Pathogen))) + 
    geom_violin(alpha = 0.5, draw_quantiles = c(0.1, 0.5, 0.9), colour = "grey50", scale = "width") + 
    geom_jitter(data = subset(res, out_contr == 1),
              aes(x = as.factor(Pathogen), y = rel_ciw, colour = as.factor(Pathogen)),
              width = 0.1, alpha = 0.5, colour = "grey50") + 
    scale_y_continuous(limits = c(0, 1.5)) +
    scale_x_discrete(labels = c("Ebola", "Influenza", "SARS-CoV-2")) + 
    facet_grid(.~param, scales = "free") + 
    scale_fill_brewer(type = "qual") + 
    scale_color_brewer(type = "qual") + 
    theme_minimal() + 
    theme(legend.position = "none", text = element_text(size = 15, family = "serif"),
          axis.text.x = element_text(angle = 320, hjust = 0, vjust = 1),
          plot.margin = margin(10, 50, 10, 10)) + 
    labs(x = NULL, y = "Posterior contraction", title = "Posterior 95% CI width vs. prior 95% CI width ")

res_covg <-
    res %>% 
    group_by(param, Pathogen) %>% 
    summarise(covg = mean(coverage),
              lower = BinomCI(x = sum(coverage), 1e3, method = "wilson")[2],
              upper = BinomCI(x = sum(coverage), 1e3, method = "wilson")[3]) %>% 
    ungroup %>% as.data.frame()

covg_plot <-
    res_covg %>% 
    ggplot(aes(x = Pathogen, y = covg, 
               ymin = lower, ymax = upper,
               colour = Pathogen)) + 
    geom_pointrange() +
    geom_hline(yintercept = 0.95, linetype = "dashed", colour = "grey50") + 
    scale_y_continuous(limits = c(0.9, 1)) + 
    scale_x_discrete(labels = c("Ebola", "Influenza", "SARS-CoV-2")) + 
    facet_grid(.~param, scales = "free") + 
    scale_fill_brewer(type = "qual") + 
    scale_color_brewer(type = "qual") + 
    theme_minimal() + 
    theme(legend.position = "none", text = element_text(size = 15, family = "serif"),
          axis.text.x = element_text(angle = 320, hjust = 0, vjust = 1),
          plot.margin = margin(10, 50, 10, 10)) + 
    labs(x = NULL, y = "Coverage", title = "Coverage of 95% credible intervals")

all_plots <- covg_plot / contraction_plot / mad_plot / ciw_plot

ggsave("threedis_res.pdf", plot = all_plots, width = 11, height = 12)


# group by disease --------------------------------------------------------
labs = c(expression(R[0]), 
         expression(1 / omega), 
         expression(1 / mu), 
         expression(rho), 
         expression(1 / sqrt(phi)))

mad_plot <-
  subset(res, out_mad == 0) %>% 
  ggplot(aes(x = as.factor(param), y = rel_mad, 
             fill = as.factor(param))) + 
  geom_violin(alpha = 0.5, draw_quantiles = c(0.1, 0.5, 0.9), colour = "grey50", trim = T, scale = "width") +
  geom_jitter(data = subset(res, out_mad == 1),
              aes(x = as.factor(param), y = rel_mad, colour = as.factor(param)),
              width = 0.1, alpha = 0.5, colour = "grey50") + 
  # scale_x_discrete(labels = c("Ebola", "Influenza", "SARS-CoV-2")) + 
  scale_x_discrete(labels = labs) + 
  scale_y_continuous(limits = c(0, 3)) +
  facet_grid(.~Pathogen, scales = "free") + 
  scale_fill_brewer(type = "qual", palette = 7) + 
  scale_color_brewer(type = "qual", palette = 7) + 
  theme_minimal() + 
  theme(legend.position = "none",
        text = element_text(size = 15, family = "serif"),
        axis.text.x = element_text(vjust = 0.25),
        strip.text = element_text(size = 15)) + 
  labs(x = NULL, y = expression("|" ~ hat(theta)[0.5] - theta[true] ~ "|" / theta[true]), title = "Scaled median absolute deviation")

ciw_plot <-
  subset(res, out_ciw == 0) %>% 
  ggplot(aes(x = as.factor(param), y = rel_ciw, 
             fill = as.factor(param))) + 
  geom_violin(alpha = 0.5, draw_quantiles = c(0.1, 0.5, 0.9), colour = "grey50", scale = "width") +
  geom_jitter(data = subset(res, out_ciw == 1),
              aes(x = as.factor(param), y = rel_ciw, colour = as.factor(param)),
              width = 0.1, alpha = 0.5, colour = "grey50") + 
  facet_grid(.~Pathogen, scales = "free") + 
  scale_fill_brewer(type = "qual", palette = 7) + 
  scale_color_brewer(type = "qual", palette = 7) + 
  scale_x_discrete(labels = labs) + 
  scale_y_continuous(limits = c(0, 9)) +
  theme_minimal() + 
  theme(legend.position = "none", text = element_text(size = 15, family = "serif"),
        axis.text.x = element_text(vjust = 0.25),
        strip.text = element_text(size = 15)) + 
  labs(x = NULL, y = expression((hat(theta)[0.975] - hat(theta)[0.025]) / theta["true"]), title = "Scaled credible interval width") 

contraction_plot <-
  subset(res, out_contr == 0) %>% 
  ggplot(aes(x = as.factor(param), y = contraction, 
             fill = as.factor(param))) + 
  geom_violin(alpha = 0.5, draw_quantiles = c(0.1, 0.5, 0.9), colour = "grey50", scale = "width") + 
  geom_jitter(data = subset(res, out_contr == 1),
              aes(x = as.factor(param), y = rel_ciw, colour = as.factor(param)),
              width = 0.1, alpha = 0.5, colour = "grey50") + 
  scale_y_continuous(limits = c(0, 1.5)) +
  scale_x_discrete(labels = labs) + 
  facet_grid(.~Pathogen, scales = "free") + 
  scale_fill_brewer(type = "qual", palette = 7) + 
  scale_color_brewer(type = "qual", palette = 7) + 
  theme_minimal() + 
  theme(legend.position = "none", text = element_text(size = 15, family = "serif"),
        axis.text.x = element_text(vjust = 0.25),
        strip.text = element_text(size = 15)) + 
  labs(x = NULL, y = "Posterior contraction", title = "Posterior 95% CI width vs. prior 95% CI width ")

covg_plot <-
  res_covg %>% 
  ggplot(aes(x = param, y = covg, 
             ymin = lower, ymax = upper,
             colour = param)) + 
  geom_pointrange() +
  geom_hline(yintercept = 0.95, linetype = "dashed", colour = "grey50") + 
  scale_y_continuous(limits = c(0.9, 1)) + 
  scale_x_discrete(labels = labs) + 
  facet_grid(.~Pathogen, scales = "free") + 
  scale_fill_brewer(type = "qual", palette = 7) + 
  scale_color_brewer(type = "qual", palette = 7) + 
  theme_minimal() + 
  theme(legend.position = "none", text = element_text(size = 15, family = "serif"),
        axis.text.x = element_text(vjust = 0.25),
        strip.text = element_text(size = 15)) + 
  labs(x = NULL, y = "Coverage", title = "Coverage of 95% credible intervals")

all_plots <- covg_plot / contraction_plot / mad_plot / ciw_plot

ggsave("threedis_res_lna.pdf", plot = all_plots, width = 11, height = 12)
