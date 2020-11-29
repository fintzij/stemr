# relative median absolute deviation
rel_mad = function(par, fit, truth) {
    
    truth["R0"] = truth["R0"] - 1
    
    if(par == "R0") {
        abs(((exp(fit["50%", "log_R0_m1"])) - truth["R0"]) / truth["R0"])
        
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
    
    truth["R0"] = truth["R0"] - 1
    truth["sqrt_phi_inv"] = 1/sqrt(truth["sqrt_phi_inv"])
    
    if(par == "R0") {
        (exp(fit[upper, "log_R0_m1"]) - exp(fit[lower, "log_R0_m1"])) / truth["R0"]
        
    } else if(par == "infec_dur") {
        (exp(fit[upper, "log_infec_dur"]) - exp(fit[lower, "log_infec_dur"])) / truth["infec_dur"]
        
    } else if(par == "rho") {
        (plogis(fit[upper, "logit_rho"]) - plogis(fit[lower, "logit_rho"])) / truth["rho"]
        
    } else {
        (exp(fit[upper, "log_sqrt_phi_inv"]) - exp(fit[lower, "log_sqrt_phi_inv"])) / truth["sqrt_phi_inv"]
    }
}

# posterior contraction
contraction = function(par, fit, lower = "10%", upper = "90%", quants = c(0.1, 0.9)) {
    
    if(par == "R0") {
        (exp(fit[upper, "log_R0_m1"]) - exp(fit[lower, "log_R0_m1"])) / 
            diff(exp(qnorm(quants, 0, 0.56)))
        
    } else if(par == "infec_dur") {
        (exp(fit[upper, "log_infec_dur"]) - exp(fit[lower, "log_infec_dur"])) /
            diff(exp(qnorm(quants, 0, 0.354)))
        
    } else if(par == "rho") {
        (plogis(fit[upper, "logit_rho"]) - plogis(fit[lower, "logit_rho"])) / 
            diff(plogis(qnorm(quants, 0, 1)))
        
    } else {
        (exp(fit[upper, "log_sqrt_phi_inv"]) - exp(fit[lower, "log_sqrt_phi_inv"])) / 
            diff(qgamma(quants, 2, 4))
    }
}

# coverage
coverage = function(par, fit, truth, lower = "10%", upper = "90%") {
    
    truth["R0"] = truth["R0"] - 1
    
    if(par == "R0") {
        ifelse(truth["R0"] < exp(fit[upper, "log_R0_m1"]) &  
                   truth["R0"] > exp(fit[lower, "log_R0_m1"]), 1, 0) 
        
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
