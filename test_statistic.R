# integrate out random effects
z_statistic <- function(model_fit, up_limit, GH_level_z){
  
  b_dim = ncol(model_fit$params$G)
  model_complex = model_fit$model_complex
  beta_mu <- model_fit$params$beta_mu
  beta_sigma <- model_fit$params$beta_sigma
  beta_q <- model_fit$params$beta_q
  ##############################################
  # testing
  # beta_mu <- params_true$beta_mu
  # beta_sigma <- params_true$beta_sigma
  # beta_q <- params_true$beta_q
  
  
  # standard GH rule
  GH <- gauher(GH_level_z)
  # expand grid
  nodes <- as.matrix(expand.grid(rep(list(GH$x), b_dim)))
  w <- as.matrix(expand.grid(rep(list(GH$w), b_dim))) 
  # 2^(p/2) w*exp(||nodes||^2)
  w <- 2^(b_dim/2)*apply(w, 1, prod)*exp(rowSums(nodes*nodes))
  GH_nodes = adjust_nodes(mu_b_i=as.data.frame(t(rep(0,b_dim))), var_b_i=model_fit$params$G, nodes,w)
  
  # adjust the weights by the marginal density of random effects
  f_b <- sapply(1:GH_level_z^2, function(e) Rfast::dmvnorm(as.matrix(GH_nodes$n)[e,], mu = rep(0, b_dim), 
                                                           sigma = model_fit$params$G))
  GH_nodes$w <- GH_nodes$w * f_b
  
  # treatment =1
  treat_list <- lapply(1:GH_level_z^b_dim, function(row){
    treat=1
    r=GH_nodes$n[row,] %>% as.matrix
    
    if(model_complex=='saturated'){
      z_mu =z_sigma = z_q = cbind(1,treat, r) %>% as.matrix
    }
    if(model_complex=='normal'){
      z_mu = cbind(1,treat,r) %>% as.matrix
      z_sigma = c(1,treat, 0, 0) %>% matrix(nrow=1)
      z_q = c(1, 0,0,0)%>% matrix(nrow=1)
    }
    if(model_complex=='test'){
      z_mu =z_sigma = z_q = c(1, 0,0,0)%>% matrix(nrow=1)
    }
    
    mu <- z_mu %*% beta_mu %>% as.numeric
    sigma <- exp(z_sigma %*% beta_sigma) %>% as.numeric
    q <- exp(z_q %*% beta_q) %>% as.numeric
    
    ##########################
    # mean_t
    mean_tt <- mean_t(mu, sigma, q, up_limit)
    ######################
    # gradient w.r.t beta_mu
    # mu_grt <- mean_tt
    mu_grt <- grad(function(mu) mean_t(mu, sigma, q, up_limit), mu)
    beta_mu_grt <- mu_grt*z_mu
    # gradient w.r.t beta_sigma
    # sigma_grt <- (-log(q^2)/q + q^(-1)*digamma(q^(-2)*(1+q*sigma)))*mean_tt
    sigma_grt <- grad(function(sig) mean_t(mu, sig, q, up_limit), sigma)
    beta_sigma_grt <- sigma_grt*sigma*z_sigma
    # gradient w.r.t beta_q
    q_grt <- grad(function(q) mean_t(mu, sigma, q, up_limit), q)
    beta_q_grt <- q_grt*q*z_q
    #############
    # combine
    beta_grt <- c(beta_mu_grt, beta_sigma_grt, beta_q_grt)
    beta_grt <- beta_vec_transform(beta_grt, model_complex, 'collapse')
    return(list(mean_tt = mean_tt, beta_grt = beta_grt))
  })
  
  # treatment =0
  control_list <- lapply(1:GH_level_z^b_dim, function(row){
    treat=0
    r=GH_nodes$n[row,] %>% as.matrix
    
    if(model_complex=='saturated'){
      z_mu =z_sigma = z_q = cbind(1,treat, r) %>% as.matrix
    }
    if(model_complex=='normal'){
      z_mu = cbind(1,treat,r) %>% as.matrix
      z_sigma = c(1,treat, 0, 0) %>% matrix(nrow=1)
      z_q = c(1, 0,0,0)%>% matrix(nrow=1)
    }
    if(model_complex=='test'){
      z_mu =z_sigma = z_q = c(1, 0,0,0)%>% matrix(nrow=1)
    }
    
    mu <- z_mu %*% beta_mu %>% as.numeric
    sigma <- exp(z_sigma %*% beta_sigma) %>% as.numeric
    q <- exp(z_q %*% beta_q) %>% as.numeric
    
    ##########################
    # mean_t
    mean_tt <- mean_t(mu, sigma, q, up_limit)
    ######################
    # gradient w.r.t beta_mu
    # mu_grt <- mean_tt
    mu_grt <- grad(function(mu) mean_t(mu, sigma, q, up_limit), mu)
    beta_mu_grt <- mu_grt*z_mu
    # gradient w.r.t beta_sigma
    # sigma_grt <- (-log(q^2)/q + q^(-1)*digamma(q^(-2)*(1+q*sigma)))*mean_tt
    sigma_grt <- grad(function(sig) mean_t(mu, sig, q, up_limit), sigma)
    beta_sigma_grt <- sigma_grt*sigma*z_sigma
    # gradient w.r.t beta_q
    q_grt <- grad(function(q) mean_t(mu, sigma, q, up_limit), q)
    beta_q_grt <- q_grt*q*z_q
    #############
    # combine
    beta_grt <- c(beta_mu_grt, beta_sigma_grt, beta_q_grt)
    beta_grt <- beta_vec_transform(beta_grt, model_complex, 'collapse')
    return(list(mean_tt = mean_tt, beta_grt = beta_grt))
  })
  
  mean_treat_l <- sapply(treat_list, function(e) e$mean_tt)
  grt_treat_l <- lapply(treat_list, function(e) e$beta_grt)
  mean_control_l <- sapply(control_list, function(e) e$mean_tt)
  grt_control_l <- lapply(control_list, function(e) e$beta_grt)
  
  mean_treat <- sum(mean_treat_l * GH_nodes$w)
  mean_control <- sum(mean_control_l * GH_nodes$w)
  grt_treat <- Reduce('+', Map('*', grt_treat_l, GH_nodes$w))
  grt_control <- Reduce('+', Map('*', grt_control_l, GH_nodes$w))
  
  delta_mean_t = mean_treat - mean_control
  delta_mean_grt =grt_treat-grt_control
  
  var_delta <- as.numeric(matrix(delta_mean_grt, nrow = 1) %*% 
                            solve(beta_hes_transform(model_fit$I_beta, model_complex, 'collapse')) %*%
                            matrix(delta_mean_grt, ncol = 1))
  
  return(list(delta_RMST = delta_mean_t, delta_RMST_grt =delta_mean_grt, var_delta_RMST=var_delta,
              RMST_treat= mean_treat, RMST_control=mean_control,grt_treat=grt_treat, grt_control=grt_control))  
}

# incert random effect = 0
z_statistic_naive <- function(model_fit, up_limit){
  
  b_dim = ncol(model_fit$params$G)
  model_complex = model_fit$model_complex
  beta_mu <- model_fit$params$beta_mu
  beta_sigma <- model_fit$params$beta_sigma
  beta_q <- model_fit$params$beta_q
  ##############################################
  # testing
  # beta_mu <- params_true$beta_mu
  # beta_sigma <- params_true$beta_sigma
  # beta_q <- params_true$beta_q
  
  
  # treatment =1
  treat=1
  
  if(model_complex=='saturated'){
    z_mu =z_sigma = z_q = cbind(1,treat, 0,0) %>% as.matrix
  }
  if(model_complex=='normal'){
    z_mu = cbind(1,treat,0,0) %>% as.matrix
    z_sigma = c(1,treat, 0, 0) %>% matrix(nrow=1)
    z_q = c(1, 0,0,0)%>% matrix(nrow=1)
  }
  if(model_complex=='test'){
    z_mu =z_sigma = z_q = c(1, 0,0,0)%>% matrix(nrow=1)
  }
  
  mu <- z_mu %*% beta_mu %>% as.numeric
  sigma <- exp(z_sigma %*% beta_sigma) %>% as.numeric
  q <- exp(z_q %*% beta_q) %>% as.numeric
  
  ##########################
  # mean_t
  mean_treat <- mean_t(mu, sigma, q, up_limit)
  ######################
  # gradient w.r.t beta_mu
  # mu_grt <- mean_treat
  mu_grt <- grad(function(mu) mean_t(mu, sigma, q, up_limit), mu)
  beta_mu_grt <- mu_grt*z_mu
  # gradient w.r.t beta_sigma
  #  (-log(q^2)/q + q^(-1)*digamma(q^(-2)*(1+q*sigma)))*mean_treat
  sigma_grt <-grad(function(sig) mean_t(mu, sig, q, up_limit), sigma)
  beta_sigma_grt <- sigma_grt*sigma*z_sigma
  # gradient w.r.t beta_q
  q_grt <- grad(function(q) mean_t(mu, sigma, q, up_limit), q)
  beta_q_grt <- q_grt*q*z_q
  #############
  # combine
  grt_treat <- c(beta_mu_grt, beta_sigma_grt, beta_q_grt)
  grt_treat <- beta_vec_transform(grt_treat, model_complex, 'collapse')
  
  # treatment =0
  treat=0
  if(model_complex=='saturated'){
    z_mu =z_sigma = z_q = cbind(1,treat, 0,0) %>% as.matrix
  }
  if(model_complex=='normal'){
    z_mu = cbind(1,treat,0,0) %>% as.matrix
    z_sigma = c(1,treat, 0, 0) %>% matrix(nrow=1)
    z_q = c(1, 0,0,0)%>% matrix(nrow=1)
  }
  if(model_complex=='test'){
    z_mu =z_sigma = z_q = c(1, 0,0,0)%>% matrix(nrow=1)
  }
  
  mu <- z_mu %*% beta_mu %>% as.numeric
  sigma <- exp(z_sigma %*% beta_sigma) %>% as.numeric
  q <- exp(z_q %*% beta_q) %>% as.numeric
  # browser()
  ##########################
  # mean_t
  mean_control <- mean_t(mu, sigma, q, up_limit)
  ######################
  # gradient w.r.t beta_mu
  # mu_grt <- mean_treat
  mu_grt <- grad(function(mu) mean_t(mu, sigma, q, up_limit), mu)
  beta_mu_grt <- mu_grt*z_mu
  # gradient w.r.t beta_sigma
  #  (-log(q^2)/q + q^(-1)*digamma(q^(-2)*(1+q*sigma)))*mean_treat
  sigma_grt <-grad(function(sig) mean_t(mu, sig, q, up_limit), sigma)
  beta_sigma_grt <- sigma_grt*sigma*z_sigma
  # gradient w.r.t beta_q
  q_grt <- grad(function(q) mean_t(mu, sigma, q, up_limit), q)
  beta_q_grt <- q_grt*q*z_q
  #############
  # combine
  grt_control <- c(beta_mu_grt, beta_sigma_grt, beta_q_grt)
  grt_control <- beta_vec_transform(grt_control, model_complex, 'collapse')
  
  return(list(delta_RMST = mean_treat - mean_control, delta_RMST_grt =grt_treat-grt_control,
              RMST_treat= mean_treat, RMST_control=mean_control,grt_treat=grt_treat, grt_control=grt_control))  
}

# variance of the delta RMST when fit GG distribution to survival data alone
get_var_GG <- function(model_fit, up_limit){
  beta_mu <- model_fit$fit_t$beta_noreff$beta_mu
  beta_sigma <- model_fit$fit_t$beta_noreff$beta_sigma
  beta_q <- model_fit$fit_t$beta_noreff$beta_q
  I_beta <- model_fit$fit_t$beta_noreff_hes
  
  treat=1
  if(model_fit$model_complex=='saturated'){
    z_mu =z_sigma = z_q = cbind(1,treat, 0,0) %>% as.matrix
  }
  if(model_fit$model_complex=='normal'){
    z_mu = cbind(1,treat,0,0) %>% as.matrix
    z_sigma = c(1,treat, 0, 0) %>% matrix(nrow=1)
    z_q = c(1, 0,0,0)%>% matrix(nrow=1)
  }
  if(model_fit$model_complex=='test'){
    z_mu =z_sigma = z_q = c(1, 0,0,0)%>% matrix(nrow=1)
  }
  
  mu <- z_mu %*% beta_mu %>% as.numeric
  sigma <- exp(z_sigma %*% beta_sigma) %>% as.numeric
  q <- exp(z_q %*% beta_q) %>% as.numeric
  
  # gradient w.r.t beta_mu
  mu_grt <- grad(function(mu) mean_t(mu, sigma, q, up_limit), mu)
  beta_mu_grt <- mu_grt*z_mu
  # gradient w.r.t beta_sigma
  # sigma_grt <- (-log(q^2)/q + q^(-1)*digamma(q^(-2)*(1+q*sigma)))*mean_tt
  sigma_grt <- grad(function(sig) mean_t(mu, sig, q, up_limit), sigma)
  beta_sigma_grt <- sigma_grt*sigma*z_sigma
  # gradient w.r.t beta_q
  q_grt <- grad(function(q) mean_t(mu, sigma, q, up_limit), q)
  beta_q_grt <- q_grt*q*z_q
  #############
  # combine
  grt_treat <- c(beta_mu_grt, beta_sigma_grt, beta_q_grt)
  grt_treat <- beta_vec_transform(grt_treat,model_fit$model_complex, 'collapse', rm_reff = T)
  
  #################################################
  treat=0
  if(model_fit$model_complex=='saturated'){
    z_mu =z_sigma = z_q = cbind(1,treat, 0,0) %>% as.matrix
  }
  if(model_fit$model_complex=='normal'){
    z_mu = cbind(1,treat,0,0) %>% as.matrix
    z_sigma = c(1,treat, 0, 0) %>% matrix(nrow=1)
    z_q = c(1, 0,0,0)%>% matrix(nrow=1)
  }
  if(model_fit$model_complex=='test'){
    z_mu =z_sigma = z_q = c(1, 0,0,0)%>% matrix(nrow=1)
  }
  
  mu <- z_mu %*% beta_mu %>% as.numeric
  sigma <- exp(z_sigma %*% beta_sigma) %>% as.numeric
  q <- exp(z_q %*% beta_q) %>% as.numeric
  
  # gradient w.r.t beta_mu
  mu_grt <- grad(function(mu) mean_t(mu, sigma, q, up_limit), mu)
  beta_mu_grt <- mu_grt*z_mu
  # gradient w.r.t beta_sigma
  # sigma_grt <- (-log(q^2)/q + q^(-1)*digamma(q^(-2)*(1+q*sigma)))*mean_tt
  sigma_grt <- grad(function(sig) mean_t(mu, sig, q, up_limit), sigma)
  beta_sigma_grt <- sigma_grt*sigma*z_sigma
  # gradient w.r.t beta_q
  q_grt <- grad(function(q) mean_t(mu, sigma, q, up_limit), q)
  beta_q_grt <- q_grt*q*z_q
  #############
  # combine
  grt_contl <- c(beta_mu_grt, beta_sigma_grt, beta_q_grt)
  grt_contl <- beta_vec_transform(grt_contl,model_fit$model_complex, 'collapse', rm_reff = T)
  
  
  delta_RMST_grt =grt_treat-grt_contl
  var_delta <- as.numeric(matrix(delta_RMST_grt, nrow = 1) %*% 
                            solve(I_beta) %*%
                            matrix(delta_RMST_grt, ncol = 1))
  return(var_delta)
}


mean_t <- function(mu, sigma, q, up_limit){
  error_occurred <- TRUE
  ever <- FALSE
  while (error_occurred){
    tryCatch(
      {
        t <- integrate(f = function(x) pgengamma(x,mu,sigma,q, lower.tail = F), lower = 0, upper = up_limit)
        # If no error occurs, set the flag to FALSE to exit the loop
        error_occurred <- FALSE
      },
      error = function(e) {
        # error may occur when stretch of 0 present, decrease up_limit and continue the loop
        up_limit <<- 0.8 * up_limit
        ever <<- TRUE
      }
    )
  }
  if(ever) warning('Upper limit reduced during integration. Results may not be accurate.')
  return(t$value)
}
