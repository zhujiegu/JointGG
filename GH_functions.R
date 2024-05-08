# common parts in GH
GH_com <- function(Nr.cores=1, GH_level=6, dat, params, plot_nodes=F){
  # Data list per ID
  IDs <- unique(dat$longitudinal$ID)
  dat_perID_y <- split(dat$longitudinal, dat$longitudinal$ID)
  dat_perID_y <- dat_perID_y[IDs]
  dat_perID_t <- split(dat$survival, dat$survival$ID)
  dat_perID_t <- dat_perID_t[IDs]
  
  b_dim = nrow(params$G)
  N = nrow(dat$dat)
  n = nrow(dat_t)
  N = nrow(dat_y)
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
    lapply(n_w, fun_com, dat_y=dat_perID_y$i, dat_t=dat_perID_t$i, params=params)
  })
  # browser()
  
  # lapply(n_w, fun_com, x=x[e,],y=y[e,],z=z[e], params=params)
  
  nodes <- lapply(seq_len(nrow(nodes)), function(i) nodes[i, ,drop=F])
  return(list(list_com_log=list_com_log, nodes=nodes))
}



# common part of every integral
# PRODUCT f(Y_ij|b_i) * f(V_i,Delta_i|b_i) for each i
fun_com <- function(nw, dat_y, dat_t, params){
  
  # assign nodes as random effects
  b <- nw$nodes
  b0 <- b[1,1]
  b1 <- b[1,2]
  # PRODUCT f(Y|b)
  f_y <- with(params, dat_y %>% rowwise %>% mutate(f=dnorm(x=value, mean=a0+b0+(a1+b1+a2*treat)*visits_time, sd=sqrt(sig_e2))) %>% 
         ungroup() %>% summarise(fy = prod(f))) %>% as.numeric()
  
  # f(V_i,Delta_i|b_i)
  z <- c(1,dat_t$treat,b)
  mu <- z %*% params$beta_mu %>% as.numeric
  sigma <- exp(z %*% params$beta_sigma) %>% as.numeric
  q <- exp(z %*% params$beta_q) %>% as.numeric
  
  f_t <- ifelse(dat_t$status==1, flexsurv::dgengamma(dat_t$times, mu = mu, sigma = sigma, Q=q, log = FALSE),
         flexsurv::pgengamma(dat_t$times, mu = mu, sigma = sigma, Q=q, lower.tail = F, log = FALSE))
  
  return(log(z_tu) + l_xt + l_yu + log(i$w))
}
