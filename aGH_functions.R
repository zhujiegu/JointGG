# E[b|D]
aGH_b <- function(moment=c('first','second'), list_com, n_w_adj, list_likl, Nr.cores=1){
  moment = match.arg(moment)
  IDs <- names(list_com)
  k=length(n_w_adj[[1]]$w)
  if(moment=='first'){
    results <- mclapply(IDs, mc.cores =Nr.cores, function(i){
        # We multiply each element of vec with each corresponding row of nodes
        product_matrix <- sweep(n_w_adj[[i]]$n, 1, list_com[[i]], `*`)
        # Sum the resulting rows to get a single row vector (1x2)
        return(colSums(product_matrix))
      })
  }
  if(moment=='second'){
    b_crossprod <- mclapply(n_w_adj, mc.cores =Nr.cores, function(i){
      lapply(1:k, function(row) {
        i$n[row,] %>% as.matrix %>% crossprod
      })
    })
    names(b_crossprod) <- IDs
      
    results <- mclapply(IDs, mc.cores =Nr.cores, function(i) {
      # Multiply each cross_product matrix by corresponding element in vec and sum up the matrices
      summed_matrix <- Reduce(`+`, Map(function(cp, scalar) cp * scalar, b_crossprod[[i]], list_com[[i]]))
      return(summed_matrix)
    })
  }
  results <- Map('/', results, list_likl)
  names(results) <- IDs
  return(results)
}


# E[log f(V_i,Delta_i)]
aGH_Q_beta <- function(Nr.cores=1, update_beta=F, update_mu, update_sigma, update_q, dat_perID_t, params, 
                       list_com, n_w_adj, list_likl,model_complex){
  beta_mu <- params$beta_mu
  beta_sigma <- params$beta_sigma
  beta_q <- params$beta_q
  if (update_beta){
    if(length(beta_mu) != length(update_mu)) stop('please supply beta update values with the same length as beta vectors')
    beta_mu = beta_mu+update_mu
    beta_sigma = beta_sigma+update_sigma
    beta_q = beta_q+update_q
  }
  
  IDs <- names(dat_perID_t)
  k=length(n_w_adj[[1]]$w)
  list_f <- parallel::mclapply(IDs, mc.cores = Nr.cores, function(i){
    f_t <- sapply(1:k, function(row){
      r=n_w_adj[[i]]$n[row,]
      
      if(model_complex=='saturated'){
        z_mu =z_sigma = z_q = cbind(1, dat_perID_t[[i]]$treat, r) %>% as.matrix
      }
      if(model_complex=='normal'){
        z_mu = cbind(1, dat_perID_t[[i]]$treat,r) %>% as.matrix
        z_sigma = c(1, dat_perID_t[[i]]$treat, 0, 0) %>% matrix(nrow=1)
        z_q = c(1, 0,0,0)%>% matrix(nrow=1)
      }
      if(model_complex=='test'){
        z_mu =z_sigma = z_q = c(1, 0,0,0)%>% matrix(nrow=1)
      }

      mu <- z_mu %*% beta_mu %>% as.numeric
      sigma <- exp(z_sigma %*% beta_sigma) %>% as.numeric
      q <- exp(z_q %*% beta_q) %>% as.numeric
      
      ifelse(dat_perID_t[[i]]$status==1, flexsurv::dgengamma(dat_perID_t[[i]]$times, mu = mu, sigma = sigma, Q=q, log = T),
             flexsurv::pgengamma(dat_perID_t[[i]]$times, mu = mu, sigma = sigma, Q=q, lower.tail = F, log = T))
    })
  })
  names(list_f) <- IDs
  results <- lapply(IDs, function(i) sum(list_com[[i]] * list_f[[i]], na.rm = T)) # 0*Inf=NA might occur here
  # Devide by the marginal likelihood and sum over individuals to get the final expectation
  results <- Reduce('+', Map('/', results, list_likl))
  return(results)
}


# E[d log f(V_i,Delta_i)/d \beta]
aGH_Q_grt_beta <- function(Nr.cores=1, dat_perID_t, params, list_com, n_w_adj, list_likl,model_complex){
  # IDs <- names(dat_perID_t)
  beta_mu <- params$beta_mu
  beta_sigma <- params$beta_sigma
  beta_q <- params$beta_q
  k=length(n_w_adj[[1]]$w)
  IDs <- names(dat_perID_t)
  list_grt <- parallel::mclapply(IDs, mc.cores = Nr.cores, function(i){
    grt_t <- lapply(1:k, function(row){ #value of an individual on every node
      r=n_w_adj[[i]]$n[row,]
      
      if(model_complex=='saturated'){
        z_mu =z_sigma = z_q = cbind(1, dat_perID_t[[i]]$treat, r) %>% as.matrix
      }
      if(model_complex=='normal'){
        z_mu = cbind(1, dat_perID_t[[i]]$treat,r) %>% as.matrix
        z_sigma = c(1, dat_perID_t[[i]]$treat, 0, 0) %>% matrix(nrow=1)
        z_q = c(1, 0,0,0)%>% matrix(nrow=1)
      }
      if(model_complex=='test'){
        z_mu =z_sigma = z_q = c(1, 0,0,0)%>% matrix(nrow=1)
      }
      
      mu <- z_mu %*% beta_mu %>% as.numeric
      sigma <- exp(z_sigma %*% beta_sigma) %>% as.numeric
      q <- exp(z_q %*% beta_q) %>% as.numeric
      ww=(log(dat_perID_t[[i]]$times)-mu)/sigma
      
      #######################################
      # # testing
      # z <- c(1, 1, 0.1, 0.1)
      # mu <- z %*% beta_mu %>% as.numeric
      # sigma <- exp(z %*% beta_sigma) %>% as.numeric
      # q <- exp(z %*% beta_q) %>% as.numeric
      # ww=(log(3)-mu)/sigma
      #######################################
      phi_w=-q/pracma::incgam(q^-2*exp(q*ww), q^-2) *(q^-2)^(q^-2) * exp((q^-2)*(q*ww-exp(q*ww)))
      phi_q=grad(function(q) log_upper_gamma_q(q, ww=ww), q)
      #################################################
      # #tesing
      # if(is.na(phi_w))  browser()
      # stop(cat(paste0('z=',z, ', mu=',mu, ', sigma=',sigma, ', q=',q, ', w=', ww, '\n', 
      #                 'values not realistic, use realistic data or use adaptive quadrature')))
      #################################################
      if(!is.finite(phi_w)){
        warning('unrealistic gradient of Gamma(a,x) wrt w evaluated at some node, use realistic data or use adaptive quadrature')
        grd_mu = grd_sigma = z_mu*0
      }else{
        grd_mu <- if(dat_perID_t[[i]]$status==1){-z_mu/(q*sigma)*(1-exp(q*ww))} else{-z_mu*phi_w/sigma}
        grd_sigma <- if(dat_perID_t[[i]]$status==1){-z_sigma*(1+ww/q*(1-exp(q*ww)))} else{-z_sigma*ww*phi_w}
      }
      if(!is.finite(phi_q)){
        warning('unrealistic gradient of Gamma(a,x) wrt q evaluated at some node, use realistic data or use adaptive quadrature')
        grd_q <- z_q*0
      }else{
        grd_q <- if(dat_perID_t[[i]]$status==1){z_q*q^-2*(q^2+2*digamma(q^-2)+4*log(q)-2-q*ww+2*exp(q*ww)-q*ww*exp(q*ww))}else{
          z_q*q*(phi_q+2*q^-3*digamma(q^-2))
        }
      }
      return(grd=c(grd_mu,grd_sigma,grd_q))
    })
  })
  
  # multiply by the common part in GH
  # Apply this function across the corresponding elements of list_f and common_part$list_com
  results <- Map(multiply_sum_matrices, list_grt, list_com)
  # Devide by the marginal likelihood and sum over individuals to get the final expectation
  results <- Reduce('+', Map('/', results, list_likl)) %>% as.matrix %>% as.vector
  return(results)
}


# observed information matrix of beta (including mu, sigma, q)
aGH_I_beta <- function(Nr.cores=1, dat_perID_t, params, list_com, n_w_adj, list_likl,model_complex){
  # IDs <- names(dat_perID_t)
  beta_mu <- params$beta_mu
  beta_sigma <- params$beta_sigma
  beta_q <- params$beta_q
  k=length(n_w_adj[[1]]$w)
  IDs <- names(dat_perID_t)
  
  
  list_grt <- parallel::mclapply(IDs, mc.cores = Nr.cores, function(i){
    grt_t <- lapply(1:k, function(row){ #value of an individual on every node
      r=n_w_adj[[i]]$n[row,]
      
      if(model_complex=='saturated'){
        z_mu =z_sigma = z_q = cbind(1, dat_perID_t[[i]]$treat, r) %>% as.matrix
      }
      if(model_complex=='normal'){
        z_mu = cbind(1, dat_perID_t[[i]]$treat,r) %>% as.matrix
        z_sigma = c(1, dat_perID_t[[i]]$treat, 0, 0) %>% matrix(nrow=1)
        z_q = c(1, 0,0,0)%>% matrix(nrow=1)
      }
      if(model_complex=='test'){
        z_mu =z_sigma = z_q = c(1, 0,0,0)%>% matrix(nrow=1)
      }
      
      mu <- z_mu %*% beta_mu %>% as.numeric
      sigma <- exp(z_sigma %*% beta_sigma) %>% as.numeric
      q <- exp(z_q %*% beta_q) %>% as.numeric
      ww=(log(dat_perID_t[[i]]$times)-mu)/sigma
      
      phi_w=-q/pracma::incgam(q^-2*exp(q*ww), q^-2) *(q^-2)^(q^-2) * exp((q^-2)*(q*ww-exp(q*ww)))
      phi_q=grad(function(q) log_upper_gamma_q(q, ww=ww), q)

      if(!is.finite(phi_w)){
        warning('unrealistic gradient of Gamma(a,x) wrt w evaluated at some node, use realistic data or use adaptive quadrature')
        grd_mu = grd_sigma = z_mu*0
      }else{
        grd_mu <- if(dat_perID_t[[i]]$status==1){-z_mu/(q*sigma)*(1-exp(q*ww))} else{-z_mu*phi_w/sigma}
        grd_sigma <- if(dat_perID_t[[i]]$status==1){-z_sigma*(1+ww/q*(1-exp(q*ww)))} else{-z_sigma*ww*phi_w}
      }
      if(!is.finite(phi_q)){
        warning('unrealistic gradient of Gamma(a,x) wrt q evaluated at some node, use realistic data or use adaptive quadrature')
        grd_q <- z_q*0
      }else{
        grd_q <- if(dat_perID_t[[i]]$status==1){z_q*q^-2*(q^2+2*digamma(q^-2)+4*log(q)-2-q*ww+2*exp(q*ww)-q*ww*exp(q*ww))}else{
          z_q*q*(phi_q+2*q^-3*digamma(q^-2))
        }
      }
      grd <- beta_vec_transform(c(grd_mu,grd_sigma,grd_q), model_complex, 'collapse')
      return(grd)
    })
  })
  names(list_grt) <- IDs
  list_grt2 <- lapply(list_grt, function(i) lapply(i, function(e) tcrossprod(e)))
  # E[S_i(beta)]
  list_grt <- Map(multiply_sum_matrices, list_grt, list_com)
  list_grt <- Map('/', list_grt, list_likl)
  # E[S_i(beta)S_i(beta)']
  list_grt2 <- Map(multiply_sum_matrices, list_grt2, list_com)
  list_grt2 <- Map('/', list_grt2, list_likl)
  
  #############################################################
  # hessian of log f(V,Delta|b)
  vec_beta <- beta_vec_transform(c(beta_mu, beta_sigma, beta_q), model_complex, 'collapse') # reduced beta vector
  list_hessian <- parallel::mclapply(IDs, mc.cores = Nr.cores, function(i){
    hes_t <- lapply(1:k, function(row){ #value of an individual on every node
      r=n_w_adj[[i]]$n[row,]
      
      if(model_complex=='saturated'){
        z_mu =z_sigma = z_q = cbind(1, dat_perID_t[[i]]$treat, r) %>% as.matrix
      }
      if(model_complex=='normal'){
        z_mu = cbind(1, dat_perID_t[[i]]$treat,r) %>% as.matrix
        z_sigma = c(1, dat_perID_t[[i]]$treat, 0, 0) %>% matrix(nrow=1)
        z_q = c(1, 0,0,0)%>% matrix(nrow=1)
      }
      if(model_complex=='test'){
        z_mu =z_sigma = z_q = c(1, 0,0,0)%>% matrix(nrow=1)
      }
      
      hes_t <- hessian(func = function(vec){
        betas <- beta_vec_to_param(vec, model_complex)
        beta_mu <- betas$beta_mu
        beta_sigma <- betas$beta_sigma
        beta_q <- betas$beta_q
        
        mu <- z_mu %*% beta_mu %>% as.numeric
        sigma <- exp(z_sigma %*% beta_sigma) %>% as.numeric
        q <- exp(z_q %*% beta_q) %>% as.numeric
        ww=(log(dat_perID_t[[i]]$times)-mu)/sigma
        
        ifelse(dat_perID_t[[i]]$status==1, flexsurv::dgengamma(dat_perID_t[[i]]$times, mu = mu, sigma = sigma, Q=q, log = T),
               flexsurv::pgengamma(dat_perID_t[[i]]$times, mu = mu, sigma = sigma, Q=q, lower.tail = F, log = T))
      }, x=vec_beta)
      return(hes_t)
    })
  })
  names(list_hessian) <- IDs
  # # E[B_i(beta)]
  list_hessian <- Map(multiply_sum_matrices, list_hessian, list_com)
  list_hessian <- Map('/', list_hessian, list_likl)
  
  # assemble
  term1 <- - Reduce('+', list_hessian)

  term2 <- Reduce('+', list_grt2)
  
  term3 <- matrix(0, nrow = length(list_grt[[1]]), ncol = length(list_grt[[1]]))
  # Loop through each unique pair (i, j) where i != j
  for (ii in seq_along(list_grt)) {
    for (jj in seq_along(list_grt)) {
      if (ii != jj) {
        # Calculate the cross product of the two vectors
        current_crossprod <- tcrossprod(list_grt[[ii]],list_grt[[jj]])
        # Add the current cross product to the sum matrix
        term3 <- term3 + current_crossprod
      }
    }
  }
 
  results <- term1 - term2 - term3
  solve(results)
  return(results)
}


# produce adjusted notes and weights
aGH_n_w <- function(reffects.individual, var.reffects, Nr.cores=1, GH_level=9, dat, plot_nodes=F){
  # Data list per ID
  IDs <- unique(dat$longitudinal$ID)
  
  b_dim = ncol(reffects.individual[[1]])
  N = nrow(dat$longitudinal)
  n = length(IDs)
  
  # standard GH rule
  GH <- gauher(GH_level)
  # expand grid
  nodes <- as.matrix(expand.grid(rep(list(GH$x), b_dim)))
  w <- as.matrix(expand.grid(rep(list(GH$w), b_dim))) 
  # 2^(p/2) w*exp(||nodes||^2)
  w <- 2^(b_dim/2)*apply(w, 1, prod)*exp(rowSums(nodes*nodes))

  n_w_adj <- mclapply(IDs, mc.cores = Nr.cores, function(i){
    aGH_nodes = adjust_nodes(mu_b_i=reffects.individual[[i]], var_b_i=var.reffects[[i]], nodes,w)
    ####################################################################
    # # testing
    # aGH_nodes = adjust_nodes(mu_b_i=data.frame(0,0), var_b_i=params$G, nodes,w)
    
    return(aGH_nodes)
  })
  names(n_w_adj) <- IDs
  return(n_w_adj)
}

# f(Y|b)*f(V, Delta | b)*f(b|theta)* other weights for all individuals
aGH_common <- function(params, n_w_adj, dat, model_complex, Nr.cores){
  IDs <- unique(dat$longitudinal$ID)
  dat_perID_y <- split(dat$longitudinal, dat$longitudinal$ID)
  dat_perID_y <- dat_perID_y[IDs]
  dat_perID_t <- split(dat$survival, dat$survival$ID)
  dat_perID_t <- dat_perID_t[IDs]
  
  list_common <- parallel::mclapply(IDs, mc.cores = Nr.cores, function(i){
    aGH_common_i(n_w_adj[[i]], dat_y = dat_perID_y[[i]], dat_t = dat_perID_t[[i]], params = params, model_complex = model_complex)})
  names(list_common) <- IDs

  # names(list_fs) <- IDs
  return(list_common)
}

# f(Y|b)*f(V, Delta | b)*f(b|theta)* other weights for an i
aGH_common_i <- function(nw, dat_y, dat_t, params, model_complex){
  ##################################
  # testing
  # nw=n_w[[1]]
  # dat_y=dat_perID_y$ID1
  # dat_t=dat_perID_t$ID1
  ##################################
  b_dim=2
  
  fs <- sapply(1:length(nw$w), function(k){
    # assign nodes as random effects
    r <- nw$n[k, ,drop=F] %>% as.matrix() 
    b0 <- r[1,1]
    b1 <- r[1,2]
    
    ######################################
    # PRODUCT f(Y|b)
    f_y <- with(params, dat_y %>% rowwise %>% mutate(mean=a0+b0+(a1+b1+a2*treat)*visits_time,sd=sqrt(sig_e2)) %>% 
                  mutate(f=dnorm(x=value, mean, sd)) %>% ungroup() %>% summarise(fy = prod(f))) %>% as.numeric()
    
    ######################################
    # f(V, Delta | b)
    if(model_complex=='saturated'){
      z_mu =z_sigma = z_q = cbind(1, dat_t$treat,r)
    }
    if(model_complex=='normal'){
      z_mu = c(1, dat_t$treat,r)
      z_sigma = z_q = c(1, dat_t$treat, 0, 0)
    }
    if(model_complex=='test'){
      z_mu =z_sigma = z_q = c(1, 0,0,0)
    }
    mu <- z_mu %*% params$beta_mu %>% as.numeric
    sigma <- exp(z_sigma %*% params$beta_sigma) %>% as.numeric
    q <- exp(z_q %*% params$beta_q) %>% as.numeric
    f_t <- ifelse(dat_t$status==1, flexsurv::dgengamma(dat_t$times, mu = mu, sigma = sigma, Q=q, log = FALSE),
                  flexsurv::pgengamma(dat_t$times, mu = mu, sigma = sigma, Q=q, lower.tail = F, log = FALSE))
    
    ######################################
    # f(sqrt(2)*tilde(r_t)|theta)
    f_b <- Rfast::dmvnorm(r, mu = rep(0, b_dim), sigma = params$G)
    return(f_y*f_t*f_b)
  })
  return(fs*nw$w)
}

# multiplies each matrix by the corresponding scalar and sums all matrices for each ID
multiply_sum_matrices <- function(matrices, scalars) {
  Reduce(`+`, mapply(function(matrix, scalar) {
    matrix * scalar
  }, matrices, scalars, SIMPLIFY = FALSE))
}

adjust_nodes <- function(mu_b_i, var_b_i, nodes, w){
  # the lower triangle L
  L = chol(var_b_i) %>% t
  # 2^(p/2)|L|w*exp(||nodes||^2)
  w_i = det(L)*w
  # mu + \sqrt(2)L*node 
  nodes_i = slice(mu_b_i, rep(1, nrow(nodes)))  + sqrt(2)*t(L%*%t(nodes))
  ###################################################################
  # # testing
  # B_tilde_inv <- solve(chol(solve(params$G)))
  # nodes_i = sqrt(2)*t(B_tilde_inv%*%t(nodes))
  return(list(n=nodes_i, w=w_i))
}

# log upper imcomplete gamma function, for use in numerical derivative
log_upper_gamma_q <- function(q, ww){
  log(pracma::incgam(q^-2*exp(q*ww), q^-2))
}

# Standard Gauss-Hermite rule from JM
gauher <- function (n) {
  m <- trunc((n + 1)/2)
  x <- w <- rep(-1, n)
  for (i in seq_len(m)) {
    z <- if (i == 1) {
      sqrt(2*n + 1) - 1.85575 * (2*n + 1)^(-0.16667)
    } else if (i == 2) {
      z - 1.14 * n^0.426 / z
    } else if (i == 3) {
      1.86 * z - 0.86 * x[1]
    } else if (i == 4) {
      1.91 * z - 0.91 * x[2]
    } else {
      2*z - x[i - 2]
    }
    for (its in seq_len(10)) {
      p1 <- 0.751125544464943
      p2 <- 0
      for (j in seq_len(n)) {
        p3 <- p2
        p2 <- p1
        p1 <- z * sqrt(2/j) * p2 - sqrt((j - 1)/j) * p3
      }
      pp <- sqrt(2*n) * p2
      z1 <- z
      z <- z1 - p1/pp
      if (abs(z - z1) <= 3e-14) 
        break
    }
    x[i] <- z
    x[n + 1 - i] <- -z
    w[i] <- 2 / (pp * pp)
    w[n + 1 - i] <- w[i]
  }
  list(x = x, w = w)
}