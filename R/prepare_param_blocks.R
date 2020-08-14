#' Prepare the parameter blocks for MCMC
#'
#' @param param_blocks list of parameter blocks
#' @param parameters vector of parameters
#' @param param_codes vector of parameter codes
#' @param iterations number of MCMC iterations
#'
#' @return validated param_block list with mcmc bookkeeping objects
#' @export
prepare_param_blocks = function(param_blocks, parameters, param_codes, iterations) {
      
      ### first validate the parameter blocks
      # check that all parameters are in a block
      if(!identical(sort(c(sapply(param_blocks, function(x) x$pars_nat))),
                    sort(names(parameters)[!grepl("_0", names(parameters))]))) {
            stop("Not all parameters are part of a parameter block.")
      }
   
      ### get index for when MCMC adaptation is to be stopped
      adapt_iters = sapply(param_blocks, function(x) x$control$stop_adaptation)
      if(!all.equal(adapt_iters, adapt_iters[1])) {
            message("MCMC proposals for some parameter blocks are adapted longer than 
                    others. All blocks will be adapted until the maximum adaptation.")
      }
      max_adaptation = max(adapt_iters)
      
      # names of the param codes 
      code_names = names(param_codes)
         
      # functions for going to and from the estimation scale
      log_prior         = lapply(param_blocks, function(x) x$priors$logprior)
      param_inds_nat    = lapply(param_blocks, function(x) match(x$pars_nat, code_names))
      param_inds_est    = param_inds_nat
      to_est_scale      = lapply(param_blocks, function(x) x$priors$to_estimation_scale)
      from_est_scale    = lapply(param_blocks, function(x) x$priors$from_estimation_scale)
      param_initializer = lapply(param_blocks, function(x) x$initializer)
      
      # names of parameters
      param_names_nat = lapply(param_blocks, function(x) x$pars_nat)
      param_names_est = lapply(param_blocks, function(x) x$pars_est)
      
      ind_est_0 = 0
      for(s in seq_along(param_blocks)) {
            
            # get the estimation scale indices
            param_inds_est = ind_est_0 + seq_along(param_inds_nat[[s]])
            ind_est_0 = ind_est_0 + length(param_inds_nat[[s]])
            
            # initialize the parameters
            if(is.null(param_initializer[[s]])) {
                  
                  # initialize the parameters
                  param_blocks[[s]]$pars_nat = parameters[param_inds_nat[[s]]]
                  param_blocks[[s]]$pars_est = to_est_scale[[s]](param_blocks[[s]]$pars_nat)
                  
                  # check that the log-prior was not -Inf if no initializer is supplied, else initialize params
                  param_blocks[[s]]$log_pd = log_prior[[s]](param_blocks[[s]]$pars_est)
                  
                  if(is.infinite(param_blocks[[s]]$log_pd)) {
                        stop("Parameters have log prior density of negative infinity. Try another initialization.")
                  }
                  
            } else {
                  # initialize parameters
                  param_blocks[[s]]$pars_nat = param_initializer[[s]]()
                  param_blocks[[s]]$pars_est = to_est_scale[[s]](param_blocks[[s]]$pars_nat)
                  
                  # check log prior
                  param_blocks[[s]]$log_pd = log_prior[[s]](param_blocks[[s]]$pars_est)
                  
                  # keep initializing if pd is infinite
                  par_init_attempt <- 1
                  while(is.infinite(param_blocks[[s]]$log_pd) && par_init_attempt <= initialization_attempts) {
                     
                     # initialize parameters
                     param_blocks[[s]]$pars_nat = param_initializer[[s]]()
                     param_blocks[[s]]$pars_est = to_est_scale[[s]](param_blocks[[s]]$pars_nat)
                     
                     # check log prior
                     param_blocks[[s]]$log_pd = log_prior[[s]](param_blocks[[s]]$pars_est)
                        
                     par_init_attempt <- par_init_attempt + 1
                  }
                  
                  if(is.infinite(param_blocks[[s]]$log_pd)) {
                        stop("Parameters have log prior density of negative infinity. Try another initialization.")
                  }
            }
            
            # double check whether the functions for going to and from the estimation scale biject
            if(!all.equal(param_blocks[[s]]$pars_nat,
                          from_est_scale[[s]](
                             to_est_scale[[s]](
                                param_blocks[[s]]$pars_nat)))) {
               stop(paste0("Functions for going to and from the estimation scale in parameter block ", s, " do not biject."))
            }
            
            ### now set up other MCMC objects-------------
            # names and indices
            param_blocks[[s]]$param_names_nat = param_names_nat[[s]]
            param_blocks[[s]]$param_names_est = param_names_est[[s]]
            
            # parameter indices in the params_cur and params_prop matrices
            param_blocks[[s]]$param_inds_Cpp = 
               param_codes[match(param_names_nat[[s]], names(param_codes))]
            param_blocks[[s]]$param_inds_R = 
               param_blocks[[s]]$param_inds_Cpp + 1
            
            # vectors for proposals
            param_blocks[[s]]$block_size = length(param_blocks[[s]]$pars_nat)
            param_blocks[[s]]$pars_prop_nat = double(param_blocks[[s]]$block_size)
            param_blocks[[s]]$pars_prop_est = double(param_blocks[[s]]$block_size)
            copy_vec(param_blocks[[s]]$pars_prop_nat, param_blocks[[s]]$pars_nat)
            copy_vec(param_blocks[[s]]$pars_prop_est, param_blocks[[s]]$pars_est)
            param_blocks[[s]]$log_pd_prop = log_prior[[s]](param_blocks[[s]]$pars_est)
            
            # should the proposals be adapted
            param_blocks[[s]]$control$adapt = 
               param_blocks[[s]]$control$stop_adaptation != 0
            
            if(param_blocks[[s]]$control$adapt) {
               param_blocks[[s]]$control$stop_adaptation = max_adaptation
            } else {
               param_blocks[[s]]$control$stop_adaptation = 0
            }
            
            # generate sequence of gain factors
            param_blocks[[s]]$gain_factors = 
               pmin(1, param_blocks[[s]]$control$scale_constant *
                       (seq(1, iterations) * 
                           param_blocks[[s]]$control$step_size + 
                           param_blocks[[s]]$control$adaptation_offset + 1) ^ 
                       -param_blocks[[s]]$control$scale_cooling)
            
            # set nugget step size if it wasn't specified
            if(is.null(param_blocks[[s]]$control$nugget_step_size)) {
               param_blocks[[s]]$control$nugget_step_size = 
                  100 / iterations
            }
            
            # set the nugget sequence
            param_blocks[[s]]$nugget_sequence = 
               param_blocks[[s]]$control$nugget * 
               (seq(0, iterations) * 
                   param_blocks[[s]]$control$nugget_step_size + 1) ^ 
               -param_blocks[[s]]$control$nugget_cooling
            
            # zero out the adaptation and nugget sequence after adaptation ends
            param_blocks[[s]]$gain_factors[
               seq(param_blocks[[s]]$control$stop_adaptation + 1, iterations)] = 0.0
            param_blocks[[s]]$nugget_sequence[
               seq(param_blocks[[s]]$control$stop_adaptation + 1, iterations)] = 0.0
            
            # kernel_resid, _mean, and _cov for adaptation
            param_blocks[[s]]$kernel_resid = 
               double(param_blocks[[s]]$block_size)
            param_blocks[[s]]$kernel_mean  = 
               double(param_blocks[[s]]$block_size)
            param_blocks[[s]]$kernel_cov   = 
               diag(1.0, param_blocks[[s]]$block_size)
            
            # fill out the mean, covariance, and cholesky
            copy_vec(dest = param_blocks[[s]]$kernel_mean, 
                     orig = param_blocks[[s]]$pars_est)
            copy_mat(dest = param_blocks[[s]]$kernel_cov, 
                     orig = param_blocks[[s]]$sigma)
            param_blocks[[s]]$kernel_cov_chol <- 
               chol(param_blocks[[s]]$kernel_cov)
            
            # initialize the list of objects for mvnss in the block
            if(param_blocks[[s]]$alg == "mvnss") {
               
               param_blocks[[s]]$mvnss_objects = 
                  list(har_direction   = rep(0.0, param_blocks[[s]]$block_size),
                       mvn_direction   = rep(0.0, param_blocks[[s]]$block_size),
                       mvnss_propvec   = rep(0.0, param_blocks[[s]]$block_size),
                       n_expansions    = c(0.5),
                       n_contractions  = c(0.5),
                       bracket_width   = param_blocks[[s]]$control$bracket_width,
                       bracket_adaptations = 
                          sqrt(pmin(1, param_blocks[[s]]$control$scale_constant *
                                       (seq(0, iterations) * 
                                           param_blocks[[s]]$control$step_size + 
                                           param_blocks[[s]]$control$adaptation_offset + 1) ^ 
                                       -param_blocks[[s]]$control$scale_cooling)))
               
            } else if(param_blocks[[s]]$alg == "mvnmh") {
               
               param_blocks[[s]]$mvnmh_objects = 
                  list(proposal_scaling = 1.0,
                       acceptances      = 0.0)
            }
      }
      
      return(param_blocks)
}
