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
    if(model_complex=='AFT'){
      z_mu = cbind(1,treat,r) %>% as.matrix
      z_sigma = c(1,0, 0, 0) %>% matrix(nrow=1)
      z_q = c(1, 0,0,0)%>% matrix(nrow=1)
    }
    if(model_complex=='GG'){
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
    if(model_complex=='AFT'){
      z_mu = cbind(1,treat,r) %>% as.matrix
      z_sigma = c(1,0, 0, 0) %>% matrix(nrow=1)
      z_q = c(1, 0,0,0)%>% matrix(nrow=1)
    }
    if(model_complex=='GG'){
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
  
  return(list(delta_RMST = delta_mean_t, var_delta_RMST=var_delta, delta_RMST_grt =delta_mean_grt,
              RMST_treat= mean_treat, RMST_control=mean_control,grt_treat=grt_treat, grt_control=grt_control))  
}

# model_fit is the output from JM (piecewise-PH, one knot)
z_statistic_Cox <- function(model_fit, up_limit, GH_level_z, knot=2){
  
  b_dim = ncol(model_fit$coefficients$D)
  
  # get covariates and random effects of i
  treat_vec = model_fit$data.id$treat
  # get coefficients
  a0 <- model_fit$coefficients$betas['(Intercept)']
  a1 <- model_fit$coefficients$betas['visits_age']
  a2 <- model_fit$coefficients$betas['visits_age:treat']
  gam <- model_fit$coefficients$gammas['treat']
  alp <- model_fit$coefficients$alpha['alpha']
  xi1 <- model_fit$coefficients$xi['xi.1']
  xi2 <- model_fit$coefficients$xi['xi.2']
  
  # standard GH rule
  GH <- gauher(GH_level_z)
  # expand grid
  nodes <- as.matrix(expand.grid(rep(list(GH$x), b_dim)))
  w <- as.matrix(expand.grid(rep(list(GH$w), b_dim))) 
  # 2^(p/2) w*exp(||nodes||^2)
  w <- 2^(b_dim/2)*apply(w, 1, prod)*exp(rowSums(nodes*nodes))
  GH_nodes = adjust_nodes(mu_b_i=as.data.frame(t(rep(0,b_dim))), var_b_i=model_fit$coefficients$D, nodes,w)
  
  # adjust the weights by the marginal density of random effects
  f_b <- sapply(1:GH_level_z^2, function(e) Rfast::dmvnorm(as.matrix(GH_nodes$n)[e,], mu = rep(0, b_dim), 
                                                           sigma =  model_fit$coefficients$D))
  GH_nodes$w <- GH_nodes$w * f_b
  
  # treatment =1
  treat_list <- lapply(1:GH_level_z^b_dim, function(row){
    treat=1
    r=GH_nodes$n[row,] %>% as.matrix
    
    ##########################
    # mean_t
    RMST <- RMST_CoxJM_piecewise(a0,a1,a2,gam,alp,xi1,xi2, r, treat, up_limit, knot=2)
    ######################
    # gradient
    RMST_grt <- c(grad(function(x) RMST_CoxJM_piecewise(a0 = x,a1,a2,gam,alp,xi1,xi2, r, treat, up_limit, knot=2), a0),
                  grad(function(x) RMST_CoxJM_piecewise(a0,a1 = x,a2,gam,alp,xi1,xi2, r, treat, up_limit, knot=2), a1),
                  grad(function(x) RMST_CoxJM_piecewise(a0,a1,a2 = x,gam,alp,xi1,xi2, r, treat, up_limit, knot=2), a2),
                  grad(function(x) RMST_CoxJM_piecewise(a0,a1,a2,gam = x,alp,xi1,xi2, r, treat, up_limit, knot=2), gam),
                  grad(function(x) RMST_CoxJM_piecewise(a0,a1,a2,gam,alp = x,xi1,xi2, r, treat, up_limit, knot=2), alp),
                  grad(function(x) RMST_CoxJM_piecewise(a0,a1,a2,gam,alp,xi1 = x,xi2, r, treat, up_limit, knot=2), xi1),
                  grad(function(x) RMST_CoxJM_piecewise(a0,a1,a2,gam,alp,xi1,xi2 = x, r, treat, up_limit, knot=2), xi2))
    
    #############
    # combine
    return(list(RMST = RMST, RMST_grt = RMST_grt))
  })
  
  # treatment =0
  control_list <- lapply(1:GH_level_z^b_dim, function(row){
    treat=0
    r=GH_nodes$n[row,] %>% as.matrix
    
    ##########################
    # mean_t
    RMST <- RMST_CoxJM_piecewise(a0,a1,a2,gam,alp,xi1,xi2, r, treat, up_limit, knot=2)
    ######################
    # gradient
    RMST_grt <- c(grad(function(x) RMST_CoxJM_piecewise(a0 = x,a1,a2,gam,alp,xi1,xi2, r, treat, up_limit, knot=2), a0),
                  grad(function(x) RMST_CoxJM_piecewise(a0,a1 = x,a2,gam,alp,xi1,xi2, r, treat, up_limit, knot=2), a1),
                  grad(function(x) RMST_CoxJM_piecewise(a0,a1,a2 = x,gam,alp,xi1,xi2, r, treat, up_limit, knot=2), a2),
                  grad(function(x) RMST_CoxJM_piecewise(a0,a1,a2,gam = x,alp,xi1,xi2, r, treat, up_limit, knot=2), gam),
                  grad(function(x) RMST_CoxJM_piecewise(a0,a1,a2,gam,alp = x,xi1,xi2, r, treat, up_limit, knot=2), alp),
                  grad(function(x) RMST_CoxJM_piecewise(a0,a1,a2,gam,alp,xi1 = x,xi2, r, treat, up_limit, knot=2), xi1),
                  grad(function(x) RMST_CoxJM_piecewise(a0,a1,a2,gam,alp,xi1,xi2 = x, r, treat, up_limit, knot=2), xi2))
    
    #############
    # combine
    return(list(RMST = RMST, RMST_grt = RMST_grt))
  })
  
  RMST_treat_l <- sapply(treat_list, function(e) e$RMST)
  grt_treat_l <- lapply(treat_list, function(e) e$RMST_grt)
  RMST_control_l <- sapply(control_list, function(e) e$RMST)
  grt_control_l <- lapply(control_list, function(e) e$RMST_grt)
  
  RMST_treat <- sum(RMST_treat_l * GH_nodes$w)
  RMST_control <- sum(RMST_control_l * GH_nodes$w)
  grt_treat <- Reduce('+', Map('*', grt_treat_l, GH_nodes$w))
  grt_control <- Reduce('+', Map('*', grt_control_l, GH_nodes$w))
  
  delta_RMST = RMST_treat - RMST_control
  delta_RMST_grt =grt_treat-grt_control
  
  I_params <- model_fit$Hessian[c('Y.(Intercept)','Y.visits_age','Y.visits_age:treat','T.treat','T.alpha','T.xi.1','T.xi.2'),
                                c('Y.(Intercept)','Y.visits_age','Y.visits_age:treat','T.treat','T.alpha','T.xi.1','T.xi.2')]
  
  var_delta <- as.numeric(matrix(delta_RMST_grt, nrow = 1) %*% 
                            solve(I_params) %*%
                            matrix(delta_RMST_grt, ncol = 1))
  
  return(list(delta_RMST = delta_RMST, var_delta_RMST=var_delta))  
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
  if(model_complex=='AFT'){
    z_mu = cbind(1,treat,0,0) %>% as.matrix
    z_sigma = c(1,0, 0, 0) %>% matrix(nrow=1)
    z_q = c(1, 0,0,0)%>% matrix(nrow=1)
  }
  if(model_complex=='GG'){
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
  if(model_complex=='GG'){
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
get_RMST_var_GG <- function(model_fit, up_limit){
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
  if(model_fit$model_complex=='AFT'){
    z_mu = cbind(1,treat,0,0) %>% as.matrix
    z_sigma = c(1,0, 0, 0) %>% matrix(nrow=1)
    z_q = c(1, 0,0,0)%>% matrix(nrow=1)
  }
  if(model_fit$model_complex=='GG'){
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
  # RMST
  RMST_treat <- mean_t(mu, sigma, q, up_limit)
  
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
  if(model_fit$model_complex=='AFT'){
    z_mu = cbind(1,treat,0,0) %>% as.matrix
    z_sigma = c(1,0, 0, 0) %>% matrix(nrow=1)
    z_q = c(1, 0,0,0)%>% matrix(nrow=1)
  }
  if(model_fit$model_complex=='GG'){
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
  # RMST
  RMST_contl <- mean_t(mu, sigma, q, up_limit)
  
  delta_RMST = RMST_treat - RMST_contl
  delta_RMST_grt =grt_treat-grt_contl
  # variance
  var_RMST_treat <- as.numeric(matrix(grt_treat, nrow = 1) %*% 
                                 solve(I_beta) %*%
                                 matrix(grt_treat, ncol = 1))
  var_RMST_contl <- as.numeric(matrix(grt_contl, nrow = 1) %*% 
                                 solve(I_beta) %*%
                                 matrix(grt_contl, ncol = 1))
  var_delta <- as.numeric(matrix(delta_RMST_grt, nrow = 1) %*% 
                            solve(I_beta) %*%
                            matrix(delta_RMST_grt, ncol = 1))
  return(list(delta_RMST=delta_RMST, var_delta=var_delta, RMST_treat=RMST_treat,
              RMST_contl=RMST_contl, var_RMST_treat=var_RMST_treat, var_RMST_contl=var_RMST_contl))
}


RMST_CoxJM_piecewise <- function(a0,a1,a2,gam,alp,xi1,xi2, r, treat, up_limit, knot=2){
  # random effects
  r0 = r[,1]
  r1 = r[,2]
  part <- alp*(a1+r1+a2*treat)
  part_com <- exp(gam*treat+alp*(a0+r0))/part
  
  S_CoxJM_piecewise <- function(t) {
    # Define a helper function to apply the piecewise logic to a single value of t
    calculate_H <- function(ti) {
      if (ti <= knot) {
        H = xi1 * part_com * (exp(ti * part) - 1)
      } else {
        H = part_com * ((xi1 - xi2) * exp(knot * part) + xi2 * exp(ti * part) - xi1)
      }
      return(H)
    }
    H_values <- sapply(t, calculate_H)
    return(exp(-H_values))
  }
  RMST <- integrate(f = S_CoxJM_piecewise, lower = 0, upper = up_limit)
  return(RMST$value)
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
