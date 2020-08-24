bb = function(dat, pathmat, pars) {
    
    mus = pars["alpha0"] * (pathmat[-1, 2] / popsize)^pars["alpha1"] / 
          (pars["alpha0"] * (pathmat[-1, 2] / popsize)^pars["alpha1"] + 
               (1 - pathmat[-1, 2] / popsize)^pars["alpha1"])
    
    logliks = dbbinom(dat[,2], 1e3, 1 + pars["kappa"] * mus, 1 + pars["kappa"] * (1 - mus), log = T)
        
    return(logliks)
}
