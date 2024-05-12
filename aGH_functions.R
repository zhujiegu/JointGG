# common parts in GH
GH_com <- function(Nr.cores=1, GH_level=9, dat, params, plot_nodes=F){
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
  GH <- gauher(GH_level)
  # expand grid
  nodes <- as.matrix(expand.grid(rep(list(GH$x), b_dim)))
  w <- as.matrix(expand.grid(rep(list(GH$w), b_dim))) 

  # adjust location for each i
  mu <- rep(0,b_dim)
  sigma <- params$G
  
  adjust_nodes <- function(mu_b_i, var_b_i, nodes){
    B_tilde_inv <- solve(chol(solve(var_b_i)))
    nodes_i = slice(mu_b_i, rep(1, nrow(nodes)))  + sqrt(2)*t(B_tilde_inv%*%t(nodes))
    return(nodes_i)
  }
  nodes_adj <- lapply(IDs, function(i) adjust_nodes(mu_b_i=reffects.individual[[i]], var_b_i=var.reffects[[i]], nodes))
  
  # 2^q/2 \pi_t exp(||b||^2)
  aGH_w <- 2^(b_dim/2) * apply(w, 1, prod) * exp(rowSums(nodes_adj * nodes_adj))
  
  # visulize the nodes
  if(plot_nodes){
    plot(nodes_adj[[6]], cex=-5/log(w), pch=19,
         xlab=expression(x[1]),
         ylab=expression(x[2]))
  }
  

  

  
  # make list
  n_w <- lapply(seq_len(nrow(nodes)), function(i) list(nodes = nodes[i, ,drop=F], w = w[i]))
  
  list_com <- parallel::mclapply(IDs, mc.cores = Nr.cores, function(i){
    sapply(n_w, fun_com, dat_y = dat_perID_y[[i]], dat_t = dat_perID_t[[i]], params = params)})
  
  names(list_com) <- IDs
  
  # browser()
  
  # nodes <- lapply(seq_len(nrow(nodes)), function(i) nodes[i, ,drop=F])
  return(list(list_com=list_com, nodes=nodes))
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