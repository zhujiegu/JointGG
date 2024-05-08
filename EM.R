# dat is a list containing two sublist "longitudinal" and "survival"
# the longitudinal 
JM_EM <- function(dat, init_params, tol, steps, Nr.core=1){
  dat_t <- dat$longtitudinal
  dat_y <- dat$survival
    params_max <- params
    for(i in 1:steps){
      E_next = E_step(dat, params, GH_level=GH_level, Nr.core=Nr.core)
      params_next = M_step(E_next, params, dat)

      if(i == 1) logl[1] = E_next$logl
      logl[i+1] = E_next$logl
      
      if(i > 1 && abs(logl[i+1]-logl[i]) < tol) break
      if(i %in% c(1, 1e1, 1e2, 1e3, 5e3, 1e4, 4e4)) {
        print(data.frame(row.names = 1, steps = i, time = unname(proc.time()-tic)[3], diff = logl[i+1]-logl[i], logl = logl[i+1]))
      }
      # solve(0)
      if(logl[i+1] > max(logl[1:i])) params_max <- params_next
      params = params_next
      # debug_list[[i]] <- list(...)
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

E_step <- function(dat, params, GH_level=GH_level, Nr.core=Nr.core){
  
  dat_t <- dat$longtitudinal
  dat_y <- dat$survival
  # retrieve parameters
  omega=params$omega
  G=params$G
  a0=params$a0
  a1=params$a1
  a2=params$a2
  sig_e2=params$sig_e2
  beta_mu=params$beta_mu
  beta_sigma=params$beta_sigma
  beta_q=params$beta_q
  
  # define dimensions
  n = nrow(dat_t)
  N = nrow(dat_y)
  b_dim = length(omega)
  
  # Numerical integration of the common part and store the results (all sample combined)
  common <- GH_com(Nr.cores=Nr.core, level=level, X,Y,Z, params, plot_nodes=F)
  list_com_log <- common$list_com_log # values of common parts (list(length N) of list (dim^Q))
  list_nodes <- common$nodes
  # divide both numerator and denominator by exp(max), to avoid "Inf"
  list_com_sub <- lapply(list_com_log, function(e){
    e <- unlist(e)
    max_e <- max(e)
    as.list(exp(e-max_e))
  }
  )
  
  # likelihood of each sample
  list_lik_log <- lapply(list_com_log, function(e){
    e <- unlist(e)
    max_e <- max(e)
    return(max_e + log(sum(exp(e-max_e))))
  })
  
  # Denominator = f(x,y,z)/exp(max)
  dnmt <- lapply(1:N, function(e) Reduce("+", list_com_sub[[e]]))
  
  # E(tu|xyz) with numerical integration
  ## different part in integrand
  int_diff <- list_nodes
  ## Numerator
  mu_TU_nu <- lapply(1:N, function(e) Reduce("+", Map("*", int_diff, list_com_sub[[e]])))
  mu_TU <- Map("/", mu_TU_nu, dnmt)
  mu_TU <- mu_TU %>% unlist %>% matrix(nrow=2) %>% t
  
  # E(Sttuu|xyz) with numerical integration
  ## different part in integrand
  int_diff <- lapply(list_nodes, crossprod)
  ## Numerator
  S_ttuu_nu <- lapply(1:N, function(e) Reduce("+", Map("*", int_diff, list_com_sub[[e]])))
  S_ttuu <- Map("/", S_ttuu_nu, dnmt)
  # 2r x 2r, Average across N samples
  # S_ttuu <- Reduce("+", S_ttuu)/N
  # Stt <- S_ttuu[1:r,1:r,drop=FALSE]
  # Suu <- S_ttuu[r+1:r,r+1:r,drop=FALSE]
  # Stu <- S_ttuu[1:r,r+1:r,drop=FALSE]
  # Chh <- Suu - t(Stu)%*%B - t(B)%*%Stu + t(B)%*%Stt%*%B
  
  # adjust (t,u) to (t,h)
  list_nodes_th <- lapply(list_nodes, function(e) {
    e[,r+1:r] <- e[,r+1:r] - e[,1:r]%*% B
    return(e)
  })
  beta <- cbind(a0,a,b)
  
  # Numerical estimation of Q_ab and gradient
  Q_ab <- GH_Q_ab(beta,list_nodes_th,Z,list_com_sub,dnmt)
  grd_ab <- GH_grd_ab(beta,list_nodes_th,Z,list_com_sub,dnmt)
  # Backtracking rule to find the step size s
  s = 1
  Q_ab_new <- GH_Q_ab(beta + s*grd_ab,list_nodes_th,Z,list_com_sub,dnmt)
  
  while(Q_ab_new < Q_ab + 0.5*s*tcrossprod(grd_ab)){
    s = 0.8*s
    # print(s)
    Q_ab_new <- GH_Q_ab(beta + s*grd_ab,list_nodes_th,Z,list_com_sub,dnmt)
  }
  
  # log likelihood
  logl <- list_lik_log %>% unlist %>% sum
  # print(list_lik %>% unlist %>% log)
  
  # output
  list(
    mu_T = mu_TU[,1:r,drop=F],
    mu_U = mu_TU[,r+1:r,drop=F],
    S_ttuu = S_ttuu,
    # Stt = Stt,
    # Suu = Suu,
    # Stu = Stu,
    # Shh = Chh,
    s = s,
    grd_ab = grd_ab,
    logl = logl,
    GH_common=common
  )
}

params_generate <- function(method = c("random", "method2", "method3")){
  chosen_method <- match.arg(method)
  print(paste("Parameters generated using: ", chosen_method))
  if (chosen_method == "random"){
    params <- list(
      G=diag(runif(2,0.0001,10)),
      a0=runif(1,-10,10),
      a1=runif(1,-10,10),
      a2=runif(1,-10,10),
      sig_e2=runif(1,0.0001,10),
      beta_mu=runif(4,-10,10),
      beta_sigma=runif(4,-1,1),
      beta_q=runif(4,-1,1)
    )
  }
  return(params)
}
