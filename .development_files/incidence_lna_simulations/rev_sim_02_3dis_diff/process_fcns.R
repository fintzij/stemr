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
        abs((exp(fit["50%", "log_infec_dur"]) - truth["infec_dur"]) / truth["infec_dur"])
        
    } else if(par == "rho") {
        abs((plogis(fit["50%", "logit_rho"]) - truth["rho"]) / truth["rho"])
        
    } else {
        abs((exp(fit["50%", "log_sqrt_phi_inv"]) - truth["sqrt_phi_inv"]) / truth["sqrt_phi_inv"])
    }
}

# relative credible interval width
rel_ciw = function(par, fit, truth, lower = "10%", upper = "90%") {
    
    if(par == "R0") {
        (exp(fit[upper, "log_R0"]) - exp(fit[lower, "log_R0"])) / truth["R0"]
        
    } else if(par == "latent_dur") {
        (exp(fit[upper, "log_latent_dur"]) - exp(fit[lower, "log_latent_dur"])) / truth["latent_dur"]
        
    } else if(par == "infec_dur") {
        (exp(fit[upper, "log_infec_dur"]) - exp(fit[lower, "log_infec_dur"])) / truth["infec_dur"]
        
    } else if(par == "rho") {
        (plogis(fit[upper, "logit_rho"]) - plogis(fit[lower, "logit_rho"])) / truth["rho"]
        
    } else {
        (exp(fit[upper, "log_sqrt_phi_inv"]) - exp(fit[lower, "log_sqrt_phi_inv"])) / truth["sqrt_phi_inv"]
    }
}

# posterior contraction
contraction = function(par, fit, dis, lower = "10%", upper = "90%", quants = c(0.1, 0.9)) {
    
    if(dis == "flu") {
        if(par == "R0") {
            (exp(fit[upper, "log_R0"]) - exp(fit[lower, "log_R0"])) / 
                diff(exp(qnorm(quants, log(1.5), 0.5)))
            
        } else if(par == "latent_dur") {
            (exp(fit[upper, "log_latent_dur"]) - exp(fit[lower, "log_latent_dur"])) / 
                diff(exp(qnorm(quants, log(2.25/7), 0.82)))
            
        } else if(par == "infec_dur") {
            (exp(fit[upper, "log_infec_dur"]) - exp(fit[lower, "log_infec_dur"])) /
                diff(exp(qnorm(quants, log(2.25/7), 0.82)))
            
        } else if(par == "rho") {
            (plogis(fit[upper, "logit_rho"]) - plogis(fit[lower, "logit_rho"])) 
            
        } else {
            (exp(fit[upper, "log_sqrt_phi_inv"]) - exp(fit[lower, "log_sqrt_phi_inv"])) / 
                diff(qexp(quants, 3))
        }
        
    } else if(dis == "covid") {
        if(par == "R0") {
            (exp(fit[upper, "log_R0"]) - exp(fit[lower, "log_R0"])) / 
                diff(exp(qnorm(quants, log(2.4), 0.5)))
            
        } else if(par == "latent_dur") {
            (exp(fit[upper, "log_latent_dur"]) - exp(fit[lower, "log_latent_dur"])) / 
                diff(exp(qnorm(quants, log(6/7), 0.82)))
            
        } else if(par == "infec_dur") {
            (exp(fit[upper, "log_infec_dur"]) - exp(fit[lower, "log_infec_dur"])) /
                diff(exp(qnorm(quants, log(9/7), 0.82)))
            
        } else if(par == "rho") {
            (plogis(fit[upper, "logit_rho"]) - plogis(fit[lower, "logit_rho"])) 
            
        } else {
            (exp(fit[upper, "log_sqrt_phi_inv"]) - exp(fit[lower, "log_sqrt_phi_inv"])) / 
                diff(qexp(quants, 3))
        }
        
    } else if(dis == "ebola") {
        if(par == "R0") {
            (exp(fit[upper, "log_R0"]) - exp(fit[lower, "log_R0"])) / 
                diff(exp(qnorm(quants, log(1.7), 0.5)))
            
        } else if(par == "latent_dur") {
            (exp(fit[upper, "log_latent_dur"]) - exp(fit[lower, "log_latent_dur"])) / 
                diff(exp(qnorm(quants, log(9/7), 0.82)))
            
        } else if(par == "infec_dur") {
            (exp(fit[upper, "log_infec_dur"]) - exp(fit[lower, "log_infec_dur"])) /
                diff(exp(qnorm(quants, log(9/7), 0.82)))
            
        } else if(par == "rho") {
            (plogis(fit[upper, "logit_rho"]) - plogis(fit[lower, "logit_rho"])) 
            
        } else {
            (exp(fit[upper, "log_sqrt_phi_inv"]) - exp(fit[lower, "log_sqrt_phi_inv"])) / 
                diff(qexp(quants, 3))
        }
    }
}

# coverage
coverage = function(par, fit, truth, lower = "10%", upper = "90%") {
    
    if(par == "R0") {
        ifelse(truth["R0"] < exp(fit[upper, "log_R0"]) &  
                   truth["R0"] > exp(fit[lower, "log_R0"]), 1, 0) 
        
    } else if(par == "latent_dur") {
        ifelse(truth["latent_dur"] < exp(fit[upper, "log_latent_dur"]) &  
                   truth["latent_dur"] > exp(fit[lower, "log_latent_dur"]), 1, 0)
        
    } else if(par == "infec_dur") {
        ifelse(truth["infec_dur"] < exp(fit[upper, "log_infec_dur"]) &  
                   truth["infec_dur"] > exp(fit[lower, "log_infec_dur"]), 1, 0)
        
    } else if(par == "rho") {
        ifelse(truth["rho"] < plogis(fit[upper, "logit_rho"]) &  
                   truth["rho"] > plogis(fit[lower, "logit_rho"]), 1, 0)
        
    } else {
        ifelse(truth["sqrt_phi_inv"] < exp(fit[upper, "log_sqrt_phi_inv"]) &  
                   truth["sqrt_phi_inv"] > exp(fit[lower, "log_sqrt_phi_inv"]), 1, 0)
    }
}
