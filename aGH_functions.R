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
  # \sqrt(2)w*exp(||nodes||^2)
  w <- 2^(b_dim/2)*apply(w, 1, prod)*exp(rowSums(nodes*nodes))

  n_w_adj <- mclapply(IDs, mc.cores = Nr.cores, function(i){
    aGH_nodes = adjust_nodes(mu_b_i=reffects.individual[[i]], var_b_i=var.reffects[[i]], nodes,w)
    ####################################################################
    # # testing
    # aGH_nodes = adjust_nodes(mu_b_i=data.frame(0,0), var_b_i=params$G, nodes,w)
    
    return(aGH_nodes)
  })
  names(n_w_adj) <- IDs
  
  # # visulize the nodes
  # if(plot_nodes){
  #   plot(nodes_adj[[6]]$aGH_nodes, cex=nodes_adj[[6]]$aGH_w, pch=19,
  #        xlab=expression(x[1]),
  #        ylab=expression(x[2]))
  # }
  # 
  # browser()
  
  # nodes <- lapply(seq_len(nrow(nodes)), function(i) nodes[i, ,drop=F])
  return(n_w_adj)
}

aGH_results_ll <- function(params, n_w_adj, dat, Nr.cores, model_complex){
  IDs <- unique(dat$longitudinal$ID)
  dat_perID_y <- split(dat$longitudinal, dat$longitudinal$ID)
  dat_perID_y <- dat_perID_y[IDs]
  dat_perID_t <- split(dat$survival, dat$survival$ID)
  dat_perID_t <- dat_perID_t[IDs]
  
  list_fs <- parallel::mclapply(IDs, mc.cores = Nr.cores, function(i){
    aGH_densities_i(n_w_adj[[i]], dat_y = dat_perID_y[[i]], dat_t = dat_perID_t[[i]], params = params, model_complex = model_complex)})
  names(list_fs) <- IDs
  
  ll <- sapply(IDs, function(i){
    parts <- list_fs[[i]] %>% mutate(w=n_w_adj[[i]]$w)
    apply(parts, 1, prod) %>% sum
  })
  # names(list_fs) <- IDs
  return(sum(ll))
}

aGH_densities_i <- function(nw, dat_y, dat_t, params, model_complex){
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
    return(c(f_y=f_y, f_t=f_t, f_b=f_b))
  }) %>% t %>% as_tibble
  return(fs)
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
  # \sqrt(2)|L|w*exp(||nodes||^2)
  w_i = det(L)*w
  # mu + \sqrt(2)L*node 
  nodes_i = slice(mu_b_i, rep(1, nrow(nodes)))  + sqrt(2)*t(L%*%t(nodes))
  ###################################################################
  # # testing
  # B_tilde_inv <- solve(chol(solve(params$G)))
  # nodes_i = sqrt(2)*t(B_tilde_inv%*%t(nodes))
  return(list(n=nodes_i, w=w_i))
}

# Standard Gauss-Hermite rule from JM
gauher <-
  function (n) {
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