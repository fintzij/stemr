paths_full    <- array(0.0, dim = c(53, 3, 200000))
paths_ess     <- array(0.0, dim = c(53, 3, 200000))
respaths_full <- array(0.0, dim = c(53, 3, 200000))
respaths_ess  <- array(0.0, dim = c(53, 3, 200000))

for(k in 1:200000) {
        if(k%%10000 == 0) print(k)
        path_full <- propose_lna(
                lna_times         = lna_times,
                lna_pars          = lna_parameters,
                param_update_inds = param_update_inds,
                flow_matrix       = flow_matrix,
                lna_pointer       = lna_pointer,
                set_pars_pointer  = lna_set_pars_pointer
        )

        paths_full[,,k] <- path_full$lna_path
        respaths_full[,,k] <- path_full$res_path

        path_ess <- propose_lna_ess(
                path_cur          = path_cur,
                lna_times         = lna_times,
                lna_pars          = lna_parameters,
                param_update_inds = param_update_inds,
                flow_matrix       = flow_matrix,
                lna_pointer_ess   = lna_pointer_ess,
                lna_ess_set_pars_ptr = lna_ess_set_pars_ptr
        )

        paths_ess[,,k]    <- path_ess$lna_path
        respaths_ess[,,k] <- path_ess$res_path
}

paths_comp <- data.frame(time = rep(0:52, 2), method = rep(c("ess", "full"), each = 53),
                         med = c(apply(paths_ess[,2,], 1, quantile, probs = 0.5),
                                 apply(paths_full[,2,], 1, quantile, probs = 0.5)),
                         lower25 = c(apply(paths_ess[,2,], 1, quantile, probs = 0.25),
                                      apply(paths_full[,2,], 1, quantile, probs = 0.25)),
                         upper75 = c(apply(paths_ess[,2,], 1, quantile, probs = 0.75),
                                      apply(paths_full[,2,], 1, quantile, probs = 0.75)),
                         lower125 = c(apply(paths_ess[,2,], 1, quantile, probs = 0.125),
                                   apply(paths_full[,2,], 1, quantile, probs = 0.125)),
                         upper875 = c(apply(paths_ess[,2,], 1, quantile, probs = 0.875),
                                   apply(paths_full[,2,], 1, quantile, probs = 0.875)),
                         lower025 = c(apply(paths_ess[,2,], 1, quantile, probs = 0.025),
                                      apply(paths_full[,2,], 1, quantile, probs = 0.025)),
                         upper975 = c(apply(paths_ess[,2,], 1, quantile, probs = 0.975),
                                      apply(paths_full[,2,], 1, quantile, probs = 0.975)))

res_comp <- data.frame(time = rep(0:52, 2), method = rep(c("ess", "full"), each = 53),
                       med = c(apply(respaths_ess[,2,], 1, quantile, probs = 0.5),
                               apply(respaths_full[,2,], 1, quantile, probs = 0.5)),
                       lower25 = c(apply(respaths_ess[,2,], 1, quantile, probs = 0.25),
                                   apply(respaths_full[,2,], 1, quantile, probs = 0.25)),
                       upper75 = c(apply(respaths_ess[,2,], 1, quantile, probs = 0.75),
                                   apply(respaths_full[,2,], 1, quantile, probs = 0.75)),
                       lower125 = c(apply(respaths_ess[,2,], 1, quantile, probs = 0.125),
                                    apply(respaths_full[,2,], 1, quantile, probs = 0.125)),
                       upper875 = c(apply(respaths_ess[,2,], 1, quantile, probs = 0.875),
                                    apply(respaths_full[,2,], 1, quantile, probs = 0.875)),
                       lower025 = c(apply(respaths_ess[,2,], 1, quantile, probs = 0.025),
                                    apply(respaths_full[,2,], 1, quantile, probs = 0.025)),
                       upper975 = c(apply(respaths_ess[,2,], 1, quantile, probs = 0.975),
                                    apply(respaths_full[,2,], 1, quantile, probs = 0.975)))

require(ggplot2)

pdf("LNA_full_vs_LNA_ESS.pdf")
print(ggplot(paths_comp, aes(x = time, y = med, fill = method, colour = method)) + geom_ribbon(aes(ymin = lower25, ymax = upper75), alpha = 0.1) + geom_ribbon(aes(ymin = lower125, ymax = upper875), alpha = 0.1) + geom_ribbon(aes(ymin = lower025, ymax = upper975), alpha = 0.1) + geom_line()  + labs(x = "time", y = "I", title = "Pointwise distribution of prevalence - median and 50%, 75%, and 95% bands"))

print(ggplot(res_comp, aes(x = time, y = med, fill = method, colour = method)) + geom_ribbon(aes(ymin = lower25, ymax = upper75), alpha = 0.1) + geom_ribbon(aes(ymin = lower125, ymax = upper875), alpha = 0.1) + geom_ribbon(aes(ymin = lower025, ymax = upper975), alpha = 0.1) + geom_line()  + labs(x = "time", y = "I", title = "Pointwise distribution of residual path - median and 50%, 75%, and 95% bands"))

dev.off()
