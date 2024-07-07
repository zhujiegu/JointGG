
# simulate data to estimate the relationship between calendar time and information level
t_info_relation <- function(N, eval_t, params, seed, up_limit, Nr.cores, recruitment_time=1, visits_age=c(seq(0,120,by=6))/12){
  
  seed=seed
  dat <- do.call(simData_seq, c(list(n = N, seed = seed, visits_age=visits_age, recruitment_time=recruitment_time, 
                                     interim_times=eval_t), params))
  
  # JM model fitting
  fit <- JM_EM(dat$interim1, init_params='two-stage', rel_tol=1e-5, steps=100, Nr.cores=Nr.cores, model_complex='normal', 
               GH_level=7, refine_GH = T)
  delta_RMST <- z_statistic(fit, up_limit = up_limit, GH_level_z=11)
  
  # using only survival data to fit a GG distribution
  var_RMST_GG <- get_var_GG(fit, up_limit)
  
  return(list(fit=fit, var_RMST=delta_RMST$var_delta_RMST, var_RMST_GG=var_RMST_GG))
}


