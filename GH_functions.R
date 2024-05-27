# Note the functions in this file do not return the results of numerical integration, they serve as parts of GH, assembled to final results in E_step

# E[b|D]
GH_b <- function(moment=c('first','second'), common_part, list_likl){
  moment = match.arg(moment)
  if(moment=='first'){
    results <- lapply(common_part$list_com, function(vec) {
      # We multiply each element of vec with each corresponding row of nodes
      product_matrix <- sweep(common_part$nodes, 1, vec, `*`)
      # Sum the resulting rows to get a single row vector (1x2)
      return(colSums(product_matrix))
    })
  }
  if(moment=='second'){
    b_crossprod <- lapply(1:nrow(common_part$nodes), function(row) {
      common_part$nodes[row,] %>% t %>% crossprod
    })
    results <- lapply(common_part$list_com, function(vec) {
      # Multiply each cross_product matrix by corresponding element in vec and sum up the matrices
      summed_matrix <- Reduce(`+`, Map(function(cp, scalar) cp * scalar, b_crossprod, vec))
      return(summed_matrix)
    })
  }
  results <- Map('/', results, list_likl)
  return(results)
}

# E[log f(V_i,Delta_i)]
GH_Q_beta <- function(Nr.cores=1, update_beta=NULL, update_value, dat_perID_t, params, common_part, list_likl){
  beta_mu <- params$beta_mu
  beta_sigma <- params$beta_sigma
  beta_q <- params$beta_q
  if (!is.null(update_beta)){
    update_beta <- match.arg(update_beta, choices = c('mu', 'sigma', 'q'))
    if(length(beta_mu) != length(update_value)) stop('please supply an update_value with the same length as beta vectors')
    if(update_beta=='mu') beta_mu = beta_mu+update_value
    if(update_beta=='sigma') beta_sigma = beta_sigma+update_value
    if(update_beta=='q') beta_q = beta_q+update_value
  }
  
  IDs <- names(dat_perID_t)

  list_f <- parallel::mclapply(dat_perID_t, mc.cores = Nr.cores, function(i){
    f_t <- sapply(1:nrow(common_part$nodes), function(row){
      z <- c(1, i$treat, common_part$nodes[row,])
      mu <- z %*% beta_mu %>% as.numeric
      sigma <- exp(z %*% beta_sigma) %>% as.numeric
      q <- exp(z %*% beta_q) %>% as.numeric
      
      ifelse(i$status==1, flexsurv::dgengamma(i$times, mu = mu, sigma = sigma, Q=q, log = T),
             flexsurv::pgengamma(i$times, mu = mu, sigma = sigma, Q=q, lower.tail = F, log = T))
      ######################################
      # # test the pdf and cdf
      # times=10
      # mu=2.5
      # sigma=2.1
      # q=4.3
      # dgengamma(times, mu = mu, sigma = sigma, Q=q, log = T)
      # 
      # ww=(log(times)-mu)/sigma
      # log(q)-log(sigma)-log(times)-log(gamma(q^-2))-2*q^-2*log(q)+q^-2*(q*ww-exp(q*ww))
      # 
      # pgengamma(times, mu = mu, sigma = sigma, Q=q, lower.tail = T, log = F)
      # gammainc(q^-2*exp(q*ww), q^-2)
      # pgamma(q^-2*exp(q*ww), q^-2,1, lower.tail = T)
      # incgam(q^-2*exp(q*ww), q^-2)
      ######################################
    })
  })
  
  results <- lapply(IDs, function(i) sum(common_part$list_com[[i]] * list_f[[i]], na.rm = T)) # 0*Inf=NA might occur here
  # Devide by the marginal likelihood and sum over individuals to get the final expectation
  results <- Reduce('+', Map('/', results, list_likl))
  return(results)
}

# E[d log f(V_i,Delta_i)/d \beta]
GH_Q_grt_beta <- function(Nr.cores=1, beta_type=c('mu','sigma','q'), dat_perID_t, params, common_part, list_likl, model_complex){
  beta_type = match.arg(beta_type)
  # IDs <- names(dat_perID_t)
  beta_mu <- params$beta_mu
  beta_sigma <- params$beta_sigma
  beta_q <- params$beta_q
  list_f <- parallel::mclapply(dat_perID_t, mc.cores = Nr.cores, function(i){
    f_t <- lapply(1:nrow(common_part$nodes), function(row){ #value of an individual on every node
      if(model_complex=='saturated'){
        z = c(1, i$treat, common_part$nodes[row,])
      }
      if(model_complex=='normal'){
        if(beta_type=='mu') z = c(1, i$treat, common_part$nodes[row,])
        if(beta_type %in% c('sigma','q')) z = c(1, i$treat, 0, 0)
      }
      if(model_complex=='test'){
        z = c(1, 0,0,0)
      }
      mu <- z %*% beta_mu %>% as.numeric
      sigma <- exp(z %*% beta_sigma) %>% as.numeric
      q <- exp(z %*% beta_q) %>% as.numeric
      ww=(log(i$times)-mu)/sigma
      
      #######################################
      # # testing
      # z <- c(1, 1, 0.1, 0.1)
      # mu <- z %*% beta_mu %>% as.numeric
      # sigma <- exp(z %*% beta_sigma) %>% as.numeric
      # q <- exp(z %*% beta_q) %>% as.numeric
      # ww=(log(3)-mu)/sigma
      #######################################
      
      if(beta_type=='mu'){
        # use exp(log(phi)) to increase stability
        # phi_w=-exp(log(q)-log(pracma::incgam(q^-2*exp(q*ww), q^-2)) + (q^-2)*log(q^-2)+(q^-2)*(q*ww-exp(q*ww)))
        phi_w=-q/pracma::incgam(q^-2*exp(q*ww), q^-2) *(q^-2)^(q^-2) * exp((q^-2)*(q*ww-exp(q*ww)))
        #################################################
        # #tesing
        # if(is.na(phi_w))  browser()
        # stop(cat(paste0('z=',z, ', mu=',mu, ', sigma=',sigma, ', q=',q, ', w=', ww, '\n', 
        #                 'values not realistic, use realistic data or use adaptive quadrature')))
        #################################################
        if(any(is.na(phi_w),phi_w==Inf, phi_w==-Inf)){
          warning('unrealistic gradient of beta_mu evaluated at some node, use realistic data or use adaptive quadrature')
          grd <- z*0
        }else{
          grd <- if(i$status==1){-z/(q*sigma)*(1-exp(q*ww))} else{-z*phi_w/sigma} 
        }
      }
      if(beta_type=='sigma'){
        phi_w=-q/pracma::incgam(q^-2*exp(q*ww), q^-2) *(q^-2)^(q^-2) * exp((q^-2)*(q*ww-exp(q*ww)))
        if(any(is.na(phi_w),phi_w==Inf, phi_w==-Inf)){
          warning('unrealistic gradient of beta_sigma evaluated at some node, use realistic data or use adaptive quadrature')
          grd <- z*0
        }else{
          grd<- if(i$status==1){-z*(1+ww/q*(1-exp(q*ww)))} else{-z*ww*phi_w}
        }
      }
      if(beta_type=='q'){
        # phi_q=-sign(ww-2*q^-1)*exp(log(abs(ww-2*q^-1))-log(pracma::incgam(q^-2*exp(q*ww), q^-2)) + (q^-2)*log(q^-2)+(q^-2)*(q*ww-exp(q*ww)))
        phi_q=grad(function(q) log_upper_gamma_q(q, ww=ww), q)
        #################################################
        # #tesing
        # if(is.na(phi_q))  browser()
        # if(is.na(phi_q)) stop(cat(paste0('z=',z, ', mu=',mu, ', sigma=',sigma, ', q=',q, ', w=', ww, '\n', 'values not realistic, use realistic data or use adaptive quadrature')))
        #################################################
        if(any(is.na(phi_q),phi_q==Inf, phi_q==-Inf)){
          warning('unrealistic gradient of beta_q evaluated at some node, use realistic data or use adaptive quadrature')
          grd <- z*0
        }else{
          grd<- if(i$status==1){z*q^-2*(q^2+2*digamma(q^-2)+4*log(q)-2-q*ww+2*exp(q*ww)-q*ww*exp(q*ww))}else{
            z*q*(phi_q+2*q^-3*digamma(q^-2))
          }
        }
      }
      return(grd)
    })
  })
  
  # multiply by the common part in GH
  # Apply this function across the corresponding elements of list_f and common_part$list_com
  results <- Map(multiply_sum_matrices, list_f, common_part$list_com)
  # Devide by the marginal likelihood and sum over individuals to get the final expectation
  results <- Reduce('+', Map('/', results, list_likl))
  return(results)
}

# common parts in GH
GH_com <- function(Nr.cores=1, GH_level=9, dat, params, plot_nodes=F, model_complex=model_complex){
  # Data list per ID
  IDs <- unique(dat$longitudinal$ID)
  dat_perID_y <- split(dat$longitudinal, dat$longitudinal$ID)
  dat_perID_y <- dat_perID_y[IDs]
  dat_perID_t <- split(dat$survival, dat$survival$ID)
  dat_perID_t <- dat_perID_t[IDs]
  
  b_dim = nrow(params$G)
  N = nrow(dat$longitudinal)
  n = length(IDs)
  # standard GH rule
  rule <- fastGHQuad::gaussHermiteData(GH_level)
  
  # expand grid
  nodes <- as.matrix(expand.grid(lapply(apply(replicate(b_dim, rule$x), 
                                              2, list), unlist)))
  g <- as.matrix(expand.grid(lapply(apply(replicate(b_dim, rule$w), 
                                          2, list), unlist)))
  w <- apply(g, 1, prod)
  
  # adjust for mu and sigma
  mu <- rep(0,b_dim)
  sigma <- params$G
  
  nodes <- mu + t(sqrt(2)*t(chol(sigma))%*%t(nodes))
  w <- (1/sqrt(pi))^b_dim * w
  
  # visulize the nodes
  if(plot_nodes){
    plot(nodes, cex=-5/log(w), pch=19,
         xlab=expression(x[1]),
         ylab=expression(x[2]))
  }
  
  # make list
  n_w <- lapply(seq_len(nrow(nodes)), function(i) list(nodes = nodes[i, ,drop=F], w = w[i]))
  
  list_com <- parallel::mclapply(IDs, mc.cores = Nr.cores, function(i){
    sapply(n_w, fun_com, dat_y = dat_perID_y[[i]], dat_t = dat_perID_t[[i]], params = params, model_complex=model_complex)})
  
  names(list_com) <- IDs
  
  # browser()
  
  # nodes <- lapply(seq_len(nrow(nodes)), function(i) nodes[i, ,drop=F])
  return(list(list_com=list_com, nodes=nodes))
}



# common part of every integral
# PRODUCT f(Y_ij|b_i) * f(V_i,Delta_i|b_i) for each i
fun_com <- function(nw, dat_y, dat_t, params, model_complex){
  ##################################
  # testing
  # nw=n_w[[1]]
  # dat_y=dat_perID_y$ID1
  # dat_t=dat_perID_t$ID1
  ##################################
  
  # assign nodes as random effects
  b <- nw$nodes
  b0 <- b[1,1]
  b1 <- b[1,2]
  # PRODUCT f(Y|b)
  f_y <- with(params, dat_y %>% rowwise %>% mutate(mean=a0+b0+(a1+b1+a2*treat)*visits_age,sd=sqrt(sig_e2)) %>% 
                mutate(f=dnorm(x=value, mean, sd)) %>% ungroup() %>% summarise(fy = prod(f))) %>% as.numeric()
  ##################################
  # testing
  # with(params, a0+b0+(a1+b1+a2*dat_t$treat)*visits_age)
  # f_y <- with(params, dat_y %>% rowwise %>% mutate(f=dnorm(x=value, mean=a0-0.0442+(a1-1.72+a2*treat)*visits_age, sd=sqrt(sig_e2))) %>%
  #               ungroup() %>% summarise(fy = prod(f))) %>% as.numeric()
  ##################################

  
  # f(V_i,Delta_i|b_i)
  if(model_complex=='saturated'){
    z_mu =z_sigma = z_q = c(1, dat_t$treat,b)
  }
  if(model_complex=='normal'){
    z_mu = c(1, dat_t$treat,b)
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
  return(f_y*f_t*nw$w)
  # return(list(f_y=f_y, f_t=f_t, w=nw$w))
}

# multiplies each matrix by the corresponding scalar and sums all matrices for each ID
multiply_sum_matrices <- function(matrices, scalars) {
  Reduce(`+`, mapply(function(matrix, scalar) {
    matrix * scalar
  }, matrices, scalars, SIMPLIFY = FALSE))
}
