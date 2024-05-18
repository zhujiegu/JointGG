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
      step_outp = EM_step(dat, params, GH_level=GH_level, Nr.core=1, n_w_adj,  model_complex)
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
  E_outp = E_step(dat, params, GH_level=GH_level, Nr.cores=Nr.cores)
  latent_var <- E_outp[names(E_outp) %in% c('b')]
  
  message("Nr steps was ", i)
  message("Log-likelihood: ", logl[i+1])
  outputt <- list(params = params, latent_var=latent_var, logl = logl[0:i+1][-1], GH_level=GH_level) #, debug = debug_list)
  class(outputt) <- "JMGG"
  return(outputt)
}

EM_step <- function(dat, params, GH_level=GH_level, Nr.cores=Nr.cores, n_w_adj, model_complex=c('saturated','normal')){

  # Numerical integration of the common part and store the results (all sample combined)
  list_com <- aGH_common(params=params, n_w_adj=n_w_adj, dat=dat, model_complex, Nr.cores)
  # list_nodes <- common_part$nodes
  
  # likelihood of each sample
  list_likl <- lapply(list_com, function(e){
    return(sum(e))
  })
  likl_log <- sapply(list_likl, log) %>% sum

  ##################################################################
  # Start ECM for survival parameters here, E and cyclic update of betas
  ##################################################################
  IDs <- unique(dat$longitudinal$ID)
  dat_perID_y <- split(dat$longitudinal, dat$longitudinal$ID)
  dat_perID_y <- dat_perID_y[IDs]
  dat_perID_t <- split(dat$survival, dat$survival$ID)
  dat_perID_t <- dat_perID_t[IDs]
  # Backtracking rule sequentially to find the step size for mu, sigma, q, consecutively.
  #########################################
  # mu
  #########################################
  s=0.01
  # current E(log(f(V,Delta))
  Q_beta <- aGH_Q_beta(Nr.cores=1, update_beta=F, update_mu, update_sigma, update_q, dat_perID_t, params, 
                       list_com, n_w_adj, list_likl,model_complex)
  # Gradients
  Q_grt_beta <- aGH_Q_grt_beta(Nr.cores=1, dat_perID_t, params, list_com, n_w_adj, list_likl,model_complex)
  Q_beta_new <- aGH_Q_beta(Nr.cores=1, update_beta=T, update_mu=s*Q_grt_beta[1:4], update_sigma=s*Q_grt_beta[5:8], 
                           update_q=s*Q_grt_beta[9:12], dat_perID_t, params, list_com, n_w_adj, list_likl,model_complex)
  k=1 # number of search of s, to prevent infinite loop here
  while(Q_beta_new < Q_beta + 0.5*s*as.numeric(crossprod(Q_grt_beta))){
    s = 0.8*s
    # print(s)
    Q_beta_new <- aGH_Q_beta(Nr.cores=1, update_beta=T, update_mu=s*Q_grt_beta[1:4], update_sigma=s*Q_grt_beta[5:8], 
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
  E_b <- aGH_b(moment='first', list_com, n_w_adj, list_likl, Nr.cores=1)
  df_E_b <- as_tibble(do.call(rbind, E_b))
  df_E_b$ID <- names(E_b)
  colnames(df_E_b) <- c("Eb0", "Eb1", "ID")
  
  E_bb <- aGH_b(moment='second', list_com, n_w_adj, list_likl, Nr.cores=1)
  df_E_bb <- do.call(rbind, lapply(E_bb, function(e) c(Ebb0=e[1,1], Eb0b1=e[1,2], Ebb1=e[2,2]))) %>% as.tibble
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
