# generate outlier labels
label_outs = function(x) {
    quants = quantile(x, c(0.1, 0.9))
    upper_lim = quants[2] + 1.5 * abs(diff(quants))
    lower_lim = quants[1] - 1.5 * abs(diff(quants))
    out_lab = ifelse(x > upper_lim | x < lower_lim, 1, 0)
    return(out_lab)
}

# relative median absolute deviation
rel_mad = function(par, fit, truth) {
    
    if(par == "R0") {
        abs((exp(fit["50%", "log_R0"]) - truth["R0"]) / truth["R0"])
        
    } else if(par == "latent_dur") {
        abs((exp(fit["50%", "log_latent_dur"]) - truth["latent_dur"]) / truth["latent_dur"])
        
    } else if(par == "infec_dur") {
        abs((exp(fit["50%", "log_infecdur"]) - truth["infec_dur"]) / truth["infec_dur"])
        
    } else if(par == "rho") {
        abs((plogis(fit["50%", "logit_rho"]) - truth["rho"]) / truth["rho"])
        
    } else {
        abs((exp(fit["50%", "log_phi"]) - truth["sqrt_phi_inv"]) / truth["sqrt_phi_inv"])
    }
}

# relative credible interval width
rel_ciw = function(par, fit, truth, lower = "2.5%", upper = "97.5%") {
    
    if(par == "R0") {
        (exp(fit[upper, "log_R0"]) - exp(fit[lower, "log_R0"])) / truth["R0"]
        
    } else if(par == "latent_dur") {
        (exp(fit[upper, "log_latent_dur"]) - exp(fit[lower, "log_latent_dur"])) / truth["latent_dur"]
        
    } else if(par == "infec_dur") {
        (exp(fit[upper, "log_infecdur"]) - exp(fit[lower, "log_infecdur"])) / truth["infec_dur"]
        
    } else if(par == "rho") {
        (plogis(fit[upper, "logit_rho"]) - plogis(fit[lower, "logit_rho"])) / truth["rho"]
        
    } else {
        (exp(fit[upper, "log_phi"]) - exp(fit[lower, "log_phi"])) / truth["sqrt_phi_inv"]
    }
}

# posterior contraction
contraction = function(par, fit, lower = "2.5%", upper = "97.5%", quants = c(0.025, 0.975)) {
    
    if(par == "R0") {
        (exp(fit[upper, "log_R0"]) - exp(fit[lower, "log_R0"])) / 
            diff(exp(qnorm(quants, log(1.8), 0.5)))
        
    } else if(par == "latent_dur") {
        (exp(fit[upper, "log_latent_dur"]) - exp(fit[lower, "log_latent_dur"])) / 
            diff(exp(qnorm(quants, 0, 0.82)))
        
    } else if(par == "infec_dur") {
        (exp(fit[upper, "log_infecdur"]) - exp(fit[lower, "log_infecdur"])) /
            diff(exp(qnorm(quants, 0, 0.82)))
        
    } else if(par == "rho") {
        (plogis(fit[upper, "logit_rho"]) - plogis(fit[lower, "logit_rho"])) 
        
    } else {
        (exp(fit[upper, "log_phi"]) - exp(fit[lower, "log_phi"])) / 
            diff(qexp(quants, 2))
    }
}

# coverage
coverage = function(par, fit, truth, lower = "2.5%", upper = "97.5%") {
    
    if(par == "R0") {
        ifelse(truth["R0"] < exp(fit[upper, "log_R0"]) &  
                   truth["R0"] > exp(fit[lower, "log_R0"]), 1, 0) 
        
    } else if(par == "latent_dur") {
        ifelse(truth["latent_dur"] < exp(fit[upper, "log_latent_dur"]) &  
                   truth["latent_dur"] > exp(fit[lower, "log_latent_dur"]), 1, 0)
        
    } else if(par == "infec_dur") {
        ifelse(truth["infec_dur"] < exp(fit[upper, "log_infecdur"]) &  
                   truth["infec_dur"] > exp(fit[lower, "log_infecdur"]), 1, 0)
        
    } else if(par == "rho") {
        ifelse(truth["rho"] < plogis(fit[upper, "logit_rho"]) &  
                   truth["rho"] > plogis(fit[lower, "logit_rho"]), 1, 0)
        
    } else {
        ifelse(truth["sqrt_phi_inv"] < exp(fit[upper, "log_phi"]) &  
                   truth["sqrt_phi_inv"] > exp(fit[lower, "log_phi"]), 1, 0)
    }
}
