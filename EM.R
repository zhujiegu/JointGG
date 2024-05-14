# dat is a list containing two sublist "longitudinal" and "survival"
# the longitudinal 
JM_EM <- function(dat, init_params, tol=1e-3, steps=5, Nr.core=1, model_complex){
  tic <- proc.time()
  params <- init_params
    params_max <- params
    logl <- c()
    ####################################
    # testing
    params_list <- list()
    ##################################
    
    
    for(i in 1:steps){
      step_outp = EM_step(dat, params, GH_level=GH_level, Nr.core=Nr.core, model_complex)
      params_next = step_outp$params
      
      ##################################
      # testing 
      params_list[[i]] <- params_next
      ##################################
      
        
      logl[i] = step_outp$likl_log

      if(i > 1){
        if(logl[i]-logl[i-1] < 0) print(paste('Decrease in log likelihood in step',i))
        if(abs(logl[i]-logl[i-1]) < tol) break
      }
      if(i %in% c(1e1, 1e2, 1e3, 5e3, 1e4, 4e4)){
        print(data.frame(row.names = 1, steps = i, time = unname(proc.time()-tic)[3], diff = logl[i]-logl[i-1], logl = logl[i]))
      }
      # solve(0)
      if(logl[i] == max(logl[1:i])) params_max <- params_next
      params = params_next
    }
  
  # latent variable values
  E_outp = E_step(dat, params, GH_level=GH_level, Nr.core=Nr.core)
  latent_var <- E_outp[names(E_outp) %in% c('b')]
  
  message("Nr steps was ", i)
  message("Log-likelihood: ", logl[i+1])
  outputt <- list(params = params, latent_var=latent_var, logl = logl[0:i+1][-1], GH_level=GH_level) #, debug = debug_list)
  class(outputt) <- "JMGG"
  return(outputt)
}

EM_step <- function(dat, params, GH_level=GH_level, Nr.core=Nr.core, model_complex=c('saturated','normal')){
  # # retrieve parameters
  # G=params$G
  # a0=params$a0
  # a1=params$a1
  # a2=params$a2
  # sig_e2=params$sig_e2
  # beta_mu=params$beta_mu
  # beta_sigma=params$beta_sigma
  # beta_q=params$beta_q
  # 
  # # define dimensions
  # n = nrow(dat_t)
  # N = nrow(dat_y)
  # b_dim = nrow(G)
  
  # Numerical integration of the common part and store the results (all sample combined)
  common_part <- GH_com(Nr.cores=Nr.cores, GH_level=GH_level, dat=dat, params=params, plot_nodes=F, model_complex)
  list_com <- common_part$list_com # values of common parts (list(length N) of list (dim^Q))
  list_nodes <- common_part$nodes
  
  # likelihood of each sample
  list_likl <- lapply(list_com, function(e){
    return(sum(e))
  })
  likl_log <- Reduce('+', list_likl) %>% log
  
  ##################################################################
  # Start ECM for survival parameters here, E and cyclic update of betas
  ##################################################################
  dat_perID_t <- split(dat$survival, dat$survival$ID)
  # Backtracking rule sequentially to find the step size for mu, sigma, q, consecutively.
  #########################################
  # mu
  #########################################
  s=1
  # current E(log(f(V,Delta))
  Q_beta <- GH_Q_beta(Nr.cores=1, update_beta=NULL, updata_value, dat_perID_t, params, common_part, list_likl)
  # Gradients
  Q_grt_beta_mu <- GH_Q_grt_beta(Nr.cores=1, beta_type='mu', dat_perID_t, params, common_part, list_likl, model_complex)
  Q_beta_new <- GH_Q_beta(Nr.cores=1, update_beta='mu', update_value=s*Q_grt_beta_mu, dat_perID_t, params, common_part, list_likl)
  k=1 # number of search of s, to prevent infinite loop here
  while(Q_beta_new < Q_beta + 0.5*s*as.numeric(crossprod(Q_grt_beta_mu))){
    s = 0.8*s
    # print(s)
    Q_beta_new <- GH_Q_beta(Nr.cores=1, update_beta='mu', update_value=s*Q_grt_beta_mu, dat_perID_t, params, common_part, list_likl)
    k=k+1
    if(k>100){
      warning('trapped in step size search for mu')
      s=0
    }
  }
  # update for the next beta
  params$beta_mu <- params$beta_mu + s*Q_grt_beta_mu
  #########################################
  # sigma
  #########################################
  s=1
  # current E(log(f(V,Delta))
  Q_beta <- Q_beta_new
  # Gradients
  Q_grt_beta_sigma <- GH_Q_grt_beta(Nr.cores=1, beta_type='sigma', dat_perID_t, params, common_part, list_likl, model_complex)
  Q_beta_new <- GH_Q_beta(Nr.cores=1, update_beta='sigma', update_value=s*Q_grt_beta_sigma, dat_perID_t, params, common_part, list_likl)
  k=1
  while(Q_beta_new < Q_beta + 0.5*s*as.numeric(crossprod(Q_grt_beta_sigma))){
    s = 0.8*s
    # print(s)
    Q_beta_new <- GH_Q_beta(Nr.cores=1, update_beta='sigma', update_value=s*Q_grt_beta_sigma, dat_perID_t, params, common_part, list_likl)
    k=k+1
    if(k>100){
      warning('trapped in step size search for sigma')
      s=0
    }
  }
  # update for the next beta
  params$beta_sigma <- params$beta_sigma + s*Q_grt_beta_sigma
  #########################################
  # q
  #########################################
  s=1
  # current E(log(f(V,Delta))
  Q_beta <- Q_beta_new
  # Gradients
  Q_grt_beta_q <- GH_Q_grt_beta(Nr.cores=1, beta_type='q', dat_perID_t, params, common_part, list_likl, model_complex)
  Q_beta_new <- GH_Q_beta(Nr.cores=1, update_beta='q', update_value=s*Q_grt_beta_q, dat_perID_t, params, common_part, list_likl)
  k=1
  while(Q_beta_new < Q_beta + 0.5*s*as.numeric(crossprod(Q_grt_beta_q))){
    s = 0.8*s
    # print(s)
    Q_beta_new <- GH_Q_beta(Nr.cores=1, update_beta='q', update_value=s*Q_grt_beta_q, dat_perID_t, params, common_part, list_likl)
    k=k+1
    if(k>100){
      warning('trapped in step size search for q')
      s=0
    }
  }
  # update for the next beta
  params$beta_q <- params$beta_q + s*Q_grt_beta_q
  ##################################################################
  # Update parameters for longitudinal and random effects
  ##################################################################
  # E(b|D) and E(b^2|D)
  E_b <- GH_b(moment='first', common_part, list_likl)
  df_E_b <- as_tibble(do.call(rbind, E_b))
  df_E_b$ID <- names(E_b)
  colnames(df_E_b) <- c("Eb0", "Eb1", "ID")
  
  E_bb <- GH_b(moment='second', common_part, list_likl)
  df_E_bb <- do.call(rbind, lapply(E_bb, function(e) c(Ebb0=e[1,1], Eb0b1=e[1,2], Ebb1=e[2,2]))) %>% as.tibble
  df_E_bb$ID <- names(E_bb)
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
  
  outp <- list(params = params, likl_log=likl_log)
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
