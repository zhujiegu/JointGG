# dat is a list containing two sublist "longitudinal" and "survival"
# the longitudinal 
JM_EM <- function(dat, init_params='two-stage', rel_tol=1e-7, steps=100, Nr.cores=1, model_complex='normal', 
                  GH_level=7, GH_level_I=15, refine_GH=T, full_adaptive=T, max_attp=10){
  if(!all.equal(unique(dat$longitudinal$ID), dat$survival$ID)) stop('please align IDs of longitudinal and survival data')
  # longitudinal for two-stage initial parameter and also pseudo-adaptive GH nodes
  fit_y <- fit_longitudinal(dat)
  # initial nodes based on longitudinal data alone
  n_w_adj <- aGH_n_w(fit_y$reffects.individual, fit_y$var.reffects, Nr.cores=Nr.cores, GH_level=GH_level, dat, plot_nodes=F)
  
  if(all(init_params=='two-stage')){
    fit_t <- fit_survival(dat, model_complex, fit_y$reffects.individual)
    init_params <- list(G=fit_y$G,
                        a0=fit_y$a0,
                        a1=fit_y$a1,
                        a2=fit_y$a2,
                        sig_e2=fit_y$sig_e2,
                        beta_mu=fit_t$beta_mu,
                        beta_sigma=fit_t$beta_sigma,
                        beta_q=fit_t$beta_q)
  }else{
    fit_t <- NULL
  }
  
  EM_results <- NULL
  Nr_attp=0
  while (is.null(EM_results) && Nr_attp < max_attp){
    Nr_attp <- Nr_attp + 1
    print('initial parameters:')
    print(init_params)
    EM_results <- tryCatch({
      EM_body(dat, init_params, rel_tol, steps, Nr.cores, model_complex, 
              GH_level, GH_level_I, refine_GH, full_adaptive, n_w_adj)
    }, error = function(e) {
      cat("Error occurred in attempt", Nr_attp, ":", e$message, "\n")
      dat_sub <- dat_subset(dat, proportion=0.7)
      fit_y_sub <- fit_longitudinal(dat_sub)
      fit_t_sub <- fit_survival(dat_sub, model_complex, fit_y_sub$reffects.individual)
      init_params <<- list(G=fit_y_sub$G,
                           a0=fit_y_sub$a0,
                           a1=fit_y_sub$a1,
                           a2=fit_y_sub$a2,
                           sig_e2=fit_y_sub$sig_e2,
                           beta_mu=fit_t_sub$beta_mu,
                           beta_sigma=fit_t_sub$beta_sigma,
                           beta_q=fit_t_sub$beta_q)
      
      # if(Nr_attp < max_attp/2){
      #   init_params <<- lapply(init_params, jitter)
      # }else{ # usually convergence problem is caused by beta values that are too large
      #   params_names <- names(params)
      #   init_params <<- lapply(params_names, function(name) {
      #     param <- params[[name]]
      #     if(startsWith(name, 'beta')) {
      #       param <- param / 5
      #     }else{
      #       param
      #     }
      #     return(param)
      #   })
      #   names(init_params) <<- params_names
      # }
      return(NULL)  # Return NULL to signal failure
    })
  }
  
  if(!is.null(EM_results)){
    message(paste("EM successful after", Nr_attp, "attempts"))
    message("Nr steps was ", EM_results$Nr_steps)
    message("Log-likelihood: ", EM_results$logl_final)
    EM_results$fit_y <- fit_y
    EM_results$fit_t <- fit_t
    EM_results$success <- T
    class(EM_results) <- "JMGG"
  }else{
    message(paste("EM failed after", max_attp, "attempts"))
    EM_results = list(dat=dat, params=init_params, model_complex=model_complex, 
                      GH_level=GH_level, GH_level_I=GH_level_I, n_w_adj=n_w_adj, success=F)
    class(EM_results) <- "JMGG_debug"
  }
  return(EM_results)
}

EM_body <- function(dat, init_params, rel_tol, steps, Nr.cores, model_complex, 
                    GH_level, GH_level_I, refine_GH, full_adaptive, n_w_initial){
  params <- init_params
  logl <- c()
  ####################################
  # # testing
  # params_list <- list()
  ##################################
  
  IDs <- unique(dat$longitudinal$ID)
  dat_perID_t <- split(dat$survival, dat$survival$ID)
  dat_perID_t <- dat_perID_t[IDs]
  
  tic <- proc.time()
  params_list <- list()
  step_outp_list <- list()
  params_list[[1]] <- init_params  # intial parameter [1] -> logl[1] -> step_outp_list[1]
  for(i in 1:steps){
    step_outp = EM_step(dat, dat_perID_t, params, GH_level=GH_level, Nr.cores=Nr.cores, n_w_initial,  model_complex)
    step_outp_list[[i]] <- step_outp
    params = step_outp$params
    params_list[[i+1]] <- params
    ##################################
    # # testing
    # params_list[[i]] <- params
    ##################################
    logl[i] = step_outp$likl_log
    
    if(i > 1){
      if(logl[i]-logl[i-1] < 0) print(paste('log likelihood decreased', logl[i-1]-logl[i], 'in step',i))
      if(logl[1]<0 && logl[i] < 2*logl[1]){
        error('EM algorithm may have diverged')
      }
      if(abs(logl[i]-max(logl[1:(i-1)])) < rel_tol*abs(logl[2]-logl[1])) break
      if(abs(logl[i]-logl[i-1]) < rel_tol*abs(logl[2]-logl[1])) break
    }
    if(i %in% c(2, seq(10, 1e6, by=10))){
      print(data.frame(row.names = 1, steps = i, time = unname(proc.time()-tic)[3], diff = logl[i]-logl[i-1], logl = logl[i]))
    }
  }
  # use the best parameters in the first stage
  params <- params_list[[which.max(logl)]] #note the i-th logl corresponds to the i-1 th params 
  step_outp <- step_outp_list[[which.max(logl)]]
  logl_final <- max(logl)
  
  print(paste('Stage I of EM complete, Nr. step =', i))
  
  if(refine_GH){
    print('Start second phase')
    # browser()
    if(full_adaptive){
      print(paste('Updating GH nodes in every iteration'))
    }else{
      print(paste('Updating GH nodes'))
      n_w_adj <- aGH_n_w(step_outp$E_b_list, step_outp$E_bb_var, Nr.cores=Nr.cores, GH_level=GH_level, dat, plot_nodes=F)
    }
    # stage two, using adaptive GH (posterior of b conditional on full data)
    for(j in 1:steps){
      if(full_adaptive){
        # update node every step
        n_w_adj <- aGH_n_w(step_outp$E_b_list, step_outp$E_bb_var, Nr.cores=Nr.cores, GH_level=GH_level, dat, plot_nodes=F)
      }
      step_outp = EM_step(dat, dat_perID_t, params, GH_level=GH_level, Nr.cores=Nr.cores, n_w_adj,  model_complex)
      step_outp_list[[i+j]] <- step_outp
      params = step_outp$params
      params_list[[i+j+1]] <- params
      ##################################
      # # testing
      # params_list[[j]] <- params
      ##################################
      
      logl[i+j] = step_outp$likl_log
      
      if(logl[i+j]-logl[i+j-1] < 0) print(paste('Log likelihood decreased', logl[i+j]-logl[i+j-1], 'in step',i+j))
      if(j>1){
        if(logl[1]<0 && logl[i+j] < 2*logl[1]){
          print('EM algorithm may have diverged')
          break
        }
        if(abs(logl[i+j]-max(logl[(i+1):(i+j-1)])) < rel_tol*abs(logl[2]-logl[1])) break
        if(abs(logl[i+j]-logl[i+j-1]) < rel_tol*abs(logl[2]-logl[1])) break
        
        if((i+j) %in% c(2, seq(10, 1e6, by=10))){
          print(data.frame(row.names = 1, steps = i+j, time = unname(proc.time()-tic)[3], diff = logl[i+j]-logl[i+j-1], logl = logl[i+j]))
        }
      }
    }
    # conclude with the best params in the second stage
    params <- params_list[[which.max(logl[-(1:i)]) + i]]
    step_outp <- step_outp_list[[which.max(logl[-(1:i)]) + i]]
    logl_final <- max(logl[-(1:i)])
    
  }else{
    j=0
  }
  
  print('EM finished, starting computing observed information matrix')
  I_beta <- NULL
  I_beta_attp <- 0
  while (is.null(I_beta) && I_beta_attp < 3) {
    I_beta_attp <- I_beta_attp + 1
    n_w_adj_I <- aGH_n_w(step_outp$E_b_list, step_outp$E_bb_var, Nr.cores=Nr.cores, GH_level=GH_level_I, 
                         dat, plot_nodes=F)
    I_beta <- tryCatch({
      compute_betaI(dat, dat_perID_t, params, n_w_adj_I, model_complex, Nr.cores)
    }, error = function(e) {
      cat("Error occurred in attempt", I_beta_attp, ":", e$message, "\n")
      GH_level_I <<- GH_level_I+3
      return(NULL)  # Return NULL to signal failure
    })
  }
  
  if(is.null(I_beta)) stop('computation of I_beta failed')
  # output
  EM_results <- list(params = params, I_beta = I_beta, E_b=step_outp$E_b, logl = logl,model_complex=model_complex,
                     
                     GH_level=GH_level, GH_level_I=GH_level_I, Nr_steps=i+j, logl_final=logl_final)
}

# One iteration of EM
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
  E_b <- lapply(E_b, function(i) i %>% t %>% as.data.frame)
  df_E_b <- as_tibble(do.call(rbind, E_b))
  df_E_b$ID <- names(E_b)
  colnames(df_E_b) <- c("Eb0", "Eb1", "ID")
  
  E_bb <- aGH_b(moment='second', list_com, n_w_adj, list_likl, Nr.cores=Nr.cores)
  df_E_bb <- do.call(rbind, lapply(E_bb, function(e) c(Ebb0=e[1,1], Eb0b1=e[1,2], Ebb1=e[2,2]))) %>% as_tibble
  df_E_bb$ID <- names(E_b)
  
  # variance of b for adaptive GH
  E_bb_var <- aGH_b(moment='second-center', list_com, n_w_adj, list_likl, Nr.cores=Nr.cores)
  
  # browser()
  #a0
  dat_y <- dat$longitudinal %>% left_join(df_E_b, by = "ID") %>% left_join(df_E_bb, by = "ID")
  update_y <- with(params, dat_y %>% mutate(a0_next=value-Eb0-(a1+Eb1+a2*treat)*visits_age,
                                            a1_next=value-a0-Eb0-(Eb1+a2*treat)*visits_age,
                                            a2_next=value-a0-Eb0-(a1+Eb1)*visits_age,
                                            sig_e2=(value-a0-(a1+a2*treat)*visits_age)^2 - 
                                              2*(value-a0-(a1+a2*treat)*visits_age)*(Eb0+Eb1*visits_age) + 
                                              Ebb0+2*Eb0b1*visits_age+Ebb1*visits_age^2,
                                            xt=treat*visits_age)) %>% select(-ID, -exposure) %>% colSums %>% t %>% as_tibble
  n = nrow(dat$survival)
  N = nrow(dat$longitudinal)
  params$a0 <- update_y$a0_next/N
  params$a1 <- update_y$a1_next/update_y$visits_age
  params$a2 <- update_y$a2_next/update_y$xt
  params$sig_e2 <- update_y$sig_e2/N
  params$G <- Reduce('+', E_bb)/n
  
  outp <- list(params = params, likl_log=likl_log, E_b=df_E_b, E_b_list=E_b, E_bb_list=E_bb, E_bb_var=E_bb_var)
  return(outp)
}

dat_subset <- function(dat, proportion=0.7){
  n = dat$survival %>% nrow
  n_sub <- ceiling(n*proportion)
  IDs <- dat$survival$ID
  IDs_sub <- sample(IDs, n_sub)
  dat_y <- dat$longitudinal %>% filter(ID %in% IDs_sub)
  dat_t <- dat$survival %>% filter(ID %in% IDs_sub)
  return(list(longitudinal=dat_y, survival=dat_t))
}

# compute information matrix for beta
compute_betaI <- function(dat, dat_perID_t, params, n_w_adj_I, model_complex, Nr.cores){
  list_com <- aGH_common(params=params, n_w_adj=n_w_adj_I, dat=dat, model_complex, Nr.cores)
  # likelihood of each sample
  list_likl <- lapply(list_com, function(e){
    return(sum(e))
  })
  # Observed information matrix
  I_beta <- aGH_I_beta(Nr.cores, dat_perID_t, params, list_com, n_w_adj_I, list_likl,model_complex)
  I_beta %>% ginv # test if it's invertible
  I_beta <- beta_hes_transform(I_beta, model_complex, 'expand')
  return(I_beta)
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
