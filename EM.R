# dat is a list containing two sublist "longitudinal" and "survival"
# the longitudinal 
JM_EM <- function(dat, init_params='two-stage', tol=1e-3, steps=5, Nr.cores=1, model_complex, GH_level){
  if(!all.equal(unique(dat$longitudinal$ID), dat$survival$ID)) stop('please align IDs of longitudinal and survival data')
  tic <- proc.time()
  # longitudinal for two-stage initial parameter and also pseudo-adaptive GH nodes
  fit_y <- fit_longitudinal(dat)
  n_w_adj <- aGH_n_w(fit_y$reffects.individual, fit_y$var.reffects, Nr.cores=Nr.cores, GH_level=GH_level, dat, plot_nodes=F)
  
  if(init_params=='two-stage'){
    fit_t <- fit_survival(dat, model_complex, fit_y$reffects.individual)
    init_params <- list(G=fit_y$G,
                        a0=fit_y$a0,
                        a1=fit_y$a1,
                        a2=fit_y$a2,
                        sig_e2=fit_y$sig_e2,
                        beta_mu=fit_t$beta_mu,
                        beta_sigma=fit_t$beta_sigma,
                        beta_q=fit_t$beta_q)
  }
  params <- init_params
  logl <- c()
  ####################################
  # # testing
  # params_list <- list()
  ##################################
  
  IDs <- unique(dat$longitudinal$ID)
  dat_perID_t <- split(dat$survival, dat$survival$ID)
  dat_perID_t <- dat_perID_t[IDs]
  
  # for(i in 1:steps){
  #   step_outp = EM_step(dat, dat_perID_t, params, GH_level=GH_level, Nr.cores=Nr.cores, n_w_adj,  model_complex)
  #   params = step_outp$params
  #   
  #   ##################################
  #   # # testing 
  #   # params_list[[i]] <- params
  #   ##################################
  #   
  #   logl[i] = step_outp$likl_log
  #   
  #   if(i > 1){
  #     if(logl[i]-logl[i-1] < 0) print(paste('Decrease in log likelihood in step',i))
  #     if(abs(logl[i]-logl[i-1]) < tol) break
  #   }
  #   if(i %in% c(2, 1e1, 1e2, 1e3, 5e3, 1e4, 4e4)){
  #     print(data.frame(row.names = 1, steps = i, time = unname(proc.time()-tic)[3], diff = logl[i]-logl[i-1], logl = logl[i]))
  #   }
  # }
  
  # stage one, 0.9* steps, 10* tol, using pseudo-adaptive GH
  for(i in 1:floor(0.9*steps)){
    step_outp = EM_step(dat, dat_perID_t, params, GH_level=GH_level, Nr.cores=Nr.cores, n_w_adj,  model_complex)
    params = step_outp$params

    ##################################
    # # testing
    # params_list[[i]] <- params
    ##################################

    logl[i] = step_outp$likl_log

    if(i > 1){
      if(logl[i]-logl[i-1] < 0) print(paste('Decrease in log likelihood in step',i))
      if(abs(logl[i]-logl[i-1]) < 10*tol) break
    }
    if(i %in% c(2, 1e1, 1e2, 1e3, 5e3, 1e4, 4e4)){
      print(data.frame(row.names = 1, steps = i, time = unname(proc.time()-tic)[3], diff = logl[i]-logl[i-1], logl = logl[i]))
    }
  }
  print(paste('stage I of EM complete, Nr. step =', i, 'Now switch to full-adaptive GH'))
  # stage one, using adaptive GH (posterior of b conditional on full data)
  for(j in (i+1):steps){
    # update node every step
    n_w_adj2 <- aGH_n_w(step_outp$E_b_list, step_outp$E_bb_list, Nr.cores=Nr.cores, GH_level=GH_level, dat, plot_nodes=F)
    step_outp = EM_step(dat, dat_perID_t, params, GH_level=GH_level, Nr.cores=Nr.cores, n_w_adj,  model_complex)
    params = step_outp$params

    ##################################
    # # testing
    # params_list[[j]] <- params
    ##################################

    logl[j] = step_outp$likl_log


    if(logl[j]-logl[j-1] < 0) print(paste('Decrease in log likelihood in step',j))
    if(abs(logl[j]-logl[j-1]) < tol) break
    if(j %in% c(2, 1e1, 1e2, 1e3, 5e3, 1e4, 4e4)){
      print(data.frame(row.names = 1, steps = i, time = unname(proc.time()-tic)[3], diff = logl[j]-logl[j-1], logl = logl[j]))
    }
  }
  
  print('EM finished, starting computing observed information matrix')
  list_com <- aGH_common(params=params, n_w_adj=n_w_adj, dat=dat, model_complex, Nr.cores)
  # likelihood of each sample
  list_likl <- lapply(list_com, function(e){
    return(sum(e))
  })
  # Observed information matrix
  I_beta <-aGH_I_beta(Nr.cores, dat_perID_t, params, list_com, n_w_adj, list_likl,model_complex)
  I_beta <- beta_hes_transform(I_beta, model_complex, 'expand')
  
  # message("Nr steps was ", i)
  # message("Log-likelihood: ", logl[i])
  message("Nr steps was ", j)
  message("Log-likelihood: ", logl[j])
  outputt <- list(params = params, I_beta = I_beta, E_b=step_outp$E_b, logl = logl, model_complex=model_complex, 
                  GH_level=GH_level, fit_y=fit_y, fit_t=fit_t) #, debug = debug_list)
  class(outputt) <- "JMGG"
  return(outputt)
}

EM_step <- function(dat, dat_perID_t, params, GH_level, Nr.cores, n_w_adj, model_complex=c('saturated','normal')){

  # Numerical integration of the common part and store the results (all sample combined)
  list_com <- aGH_common(params=params, n_w_adj=n_w_adj, dat=dat, model_complex, Nr.cores)
  # list_nodes <- common_part$nodes
  
  # likelihood of each sample
  list_likl <- lapply(list_com, function(e){
    return(sum(e))
  })
  likl_log <- sapply(list_likl, log) %>% sum

  ##################################################################
  # Start EM for survival parameters here
  ##################################################################
  IDs <- unique(dat$longitudinal$ID)
  # Backtracking rule sequentially to find the step size for mu, sigma, q, consecutively.
  #########################################
  # mu
  #########################################
  s=0.1
  # current E(log(f(V,Delta))
  Q_beta <- aGH_Q_beta(Nr.cores=Nr.cores, update_beta=F, update_mu, update_sigma, update_q, dat_perID_t, params, 
                       list_com, n_w_adj, list_likl,model_complex)
  # Gradients
  Q_grt_beta <- aGH_Q_grt_beta(Nr.cores=Nr.cores, dat_perID_t, params, list_com, n_w_adj, list_likl,model_complex)
  Q_beta_new <- aGH_Q_beta(Nr.cores=Nr.cores, update_beta=T, update_mu=s*Q_grt_beta[1:4], update_sigma=s*Q_grt_beta[5:8], 
                           update_q=s*Q_grt_beta[9:12], dat_perID_t, params, list_com, n_w_adj, list_likl,model_complex)
  k=1 # number of search of s, to prevent infinite loop here
  while(Q_beta_new < Q_beta + 0.5*s*as.numeric(crossprod(Q_grt_beta))){
    s = 0.8*s
    # print(s)
    Q_beta_new <- aGH_Q_beta(Nr.cores=Nr.cores, update_beta=T, update_mu=s*Q_grt_beta[1:4], update_sigma=s*Q_grt_beta[5:8], 
                             update_q=s*Q_grt_beta[9:12], dat_perID_t, params, list_com, n_w_adj, list_likl,model_complex)
    k=k+1
    if(k>50){
      warning(paste0('trapped in step size search for beta'))
      s=0
      break
    }
  }
  ##################################################
  # # testing
  # del=1e-19
  # grd_num <- aGH_Q_beta(Nr.cores=1, update_beta=NULL, update_value, dat_perID_t, params, list_com, n_w_adj, list_likl,model_complex)-
  #   aGH_Q_beta(Nr.cores=1, update_beta='mu', update_value=c(del,0,0,0), dat_perID_t, params, list_com, n_w_adj, list_likl,model_complex)
  # grd_num/del
  
  # update for the next beta
  params$beta_mu <- params$beta_mu + s*Q_grt_beta[1:4]
  params$beta_sigma <- params$beta_sigma + s*Q_grt_beta[5:8]
  params$beta_q <- params$beta_q + s*Q_grt_beta[9:12]
  ##################################################################
  # Update parameters for longitudinal and random effects
  ##################################################################
  # E(b|D) and E(b^2|D)
  E_b <- aGH_b(moment='first', list_com, n_w_adj, list_likl, Nr.cores=Nr.cores)
  df_E_b <- as_tibble(do.call(rbind, E_b))
  df_E_b$ID <- names(E_b)
  colnames(df_E_b) <- c("Eb0", "Eb1", "ID")
  E_b <- lapply(E_b, function(i) data.frame(matrix(i,nrow=1)))
  
  
  E_bb <- aGH_b(moment='second', list_com, n_w_adj, list_likl, Nr.cores=Nr.cores)
  df_E_bb <- do.call(rbind, lapply(E_bb, function(e) c(Ebb0=e[1,1], Eb0b1=e[1,2], Ebb1=e[2,2]))) %>% as_tibble
  df_E_bb$ID <- names(E_b)
  # browser()
  #a0
  dat_y <- dat$longitudinal %>% left_join(df_E_b, by = "ID") %>% left_join(df_E_bb, by = "ID")
  update_y <- with(params, dat_y %>% mutate(a0_next=value-Eb0-(a1+Eb1+a2*treat)*visits_time,
                                  a1_next=value-a0-Eb0-(Eb1+a2*treat)*visits_time,
                                  a2_next=value-a0-Eb0-(a1+Eb1)*visits_time,
                                  sig_e2=(value-a0-(a1+a2*treat)*visits_time)^2 - 
                                    2*(value-a0-(a1+a2*treat)*visits_time)*(Eb0+Eb1*visits_time) + 
                                    Ebb0+2*Eb0b1*visits_time+Ebb1*visits_time^2,
                                  xt=treat*visits_time)) %>% select(-ID, -exposure) %>% colSums %>% t %>% as_tibble
  n = nrow(dat$survival)
  N = nrow(dat$longitudinal)
  params$a0 <- update_y$a0_next/N
  params$a1 <- update_y$a1_next/update_y$visits_time
  params$a2 <- update_y$a2_next/update_y$xt
  params$sig_e2 <- update_y$sig_e2/N
  params$G <- Reduce('+', E_bb)/n
  
  outp <- list(params = params, likl_log=likl_log, E_b=df_E_b, E_b_list=E_b, E_bb_list=E_bb)
  return(outp)
}


params_generate <- function(method = c("random", "specify", "method3"),params_specify=NULL){
  chosen_method <- match.arg(method)
  print(paste("Parameters generated using: ", chosen_method))
  if (chosen_method == "random"){
    params <- list(
      G=diag(runif(2,0.0001,1)),
      a0=runif(1,-10,10),
      a1=runif(1,-1,1),
      a2=runif(1,-1,1),
      sig_e2=runif(1,0.0001,1),
      beta_mu=runif(4,-0.1,0.1),
      beta_sigma=runif(4,-0.1,0.1),
      beta_q=runif(4,-0.1,0.1)
    )
  }
  if (chosen_method == "specify"){
    params = params_specify
  }
  return(params)
}


params_vec_list <- function(params_vec) {
  # Assuming the order and size of each parameter is known and fixed as per your description
  params <- list()
  
  # Extract and reshape matrix G (assuming it's 2x2 as shown)
  params$G <- matrix(params_vec[1:4], nrow = 2, byrow = TRUE)
  
  # Scalars a0, a1, a2, sig_e2
  params$a0 <- params_vec[5]
  params$a1 <- params_vec[6]
  params$a2 <- params_vec[7]
  params$sig_e2 <- params_vec[8]
  
  # Assuming beta_mu is a params_vector of length 4 (based on your unlisting)
  params$beta_mu <- params_vec[9:12]
  
  # Assuming beta_sigma is another params_vector of length 4
  params$beta_sigma <- params_vec[13:16]
  
  # Assuming beta_q is a params_vector of length 4
  params$beta_q <- params_vec[17:20]
  
  return(params)
}
