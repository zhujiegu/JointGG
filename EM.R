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
  b_dim = nrow(G)
  
  # Numerical integration of the common part and store the results (all sample combined)
  common_part <- GH_com(Nr.cores=Nr.cores, GH_level=GH_level, dat=dat, params=params, plot_nodes=F)
  list_com <- common_part$list_com # values of common parts (list(length N) of list (dim^Q))
  list_nodes <- common_part$nodes
  
  # likelihood of each sample
  list_likl <- lapply(list_com, function(e){
    return(sum(e))
  })
  
  # E(b|D) and E(b^2|D)
  E_b <- Map('/', GH_b(moment='first', common_part), list_likl)
  E_bb <- Map('/', GH_b(moment='second', common_part), list_likl)
  
  dat_perID_t <- split(dat$survival, dat$survival$ID)

  # Numerical estimation of Q_beta = E(log(f(V,Delta)) and its gradient wrt to beta
  # E(log(f(V,Delta))
  Q_beta <- Map('/', GH_Q_beta(Nr.cores=1, beta_type=c('mu','sigma','q'), dat_perID_t, params, common_part), list_likl)
  # Gradients
  Q_grt_beta_mu <- GH_Q_grt_beta(Nr.cores=1, beta_type='mu', dat_perID_t, params, common_part)
  


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

M_step <- function(E_fit, params, X, Y, Z){
  orth_x = ssq(params$Wo) > 0
  orth_y = ssq(params$Co) > 0
  orth_type = match.arg(orth_type)
  with(E_fit,{
    
    N = nrow(X)
    p = ncol(X)
    q = ncol(Y)
    r = ncol(mu_T)
    rx = ncol(params$Wo)
    ry = ncol(params$Co)
    
    # 2r x 2r, Average across N samples
    S_ttuu_avgN <- Reduce("+", S_ttuu)/N
    Stt <- S_ttuu_avgN[1:r,1:r,drop=FALSE]
    Suu <- S_ttuu_avgN[r+1:r,r+1:r,drop=FALSE]
    Stu <- S_ttuu_avgN[1:r,r+1:r,drop=FALSE]
    Shh <- Suu - t(Stu)%*%params$B - t(params$B)%*%Stu + t(params$B)%*%Stt%*%params$B
    
    # print(S_ttuu_avgN)
    # print(cbind(Stt, Suu, Shh))
    
    
    # filter out samples with negative sig2E or sig2F
    # Stt, Suu different for each sample
    tmp_sig2E = sapply(1:N, function(e){
      1/p * (sum(diag(X[e,,drop=F]%*%t(X[e,,drop=F]))) - 2*sum(diag(mu_T[e,,drop=F]%*%t(params$W)%*%t(X[e,,drop=F]))) +
               sum(diag(S_ttuu[[e]][1:r,1:r,drop=FALSE])) - sum(diag(params$SigTo)))
    })
    tmp_sig2F = sapply(1:N, function(e){
      1/q * (sum(diag(Y[e,,drop=F]%*%t(Y[e,,drop=F]))) - 2*sum(diag(mu_U[e,,drop=F]%*%t(params$C)%*%t(Y[e,,drop=F]))) +
               sum(diag(S_ttuu[[e]][r+1:r,r+1:r,drop=FALSE])) - sum(diag(params$SigUo)))
    })
    
    params$sig2E = mean(tmp_sig2E)
    params$sig2F = mean(tmp_sig2F)
    
    # if(params$sig2E < 0) params$sig2E = 1e-5
    # if(params$sig2F < 0) params$sig2F = 1e-5
    
    params$B = t(Stu) %*% MASS::ginv(Stt) * diag(1,r)
    params$SigT = Stt*diag(1,r)
    params$SigU = Suu*diag(1,r)
    params$SigH = Shh*diag(1,r)
    params$a0 = params$a0 + s*grd_ab[,1]
    params$a = params$a + s*grd_ab[,1+1:r]
    params$b = params$b + s*grd_ab[,1+r+1:r]
    
    params$W = orth(t(X/N) %*% mu_T %*% MASS::ginv(Stt),type = orth_type)
    params$C = orth(t(Y/N) %*% mu_U %*% MASS::ginv(Suu),type = orth_type)#
    
    # params$sig2E = 1/p * abs(sum(diag(X%*%t(X)))/N - 2*sum(diag(mu_T%*%t(params$W)%*%t(X)))/N +
    #                        sum(diag(params$SigT)) - sum(diag(params$SigTo)))
    # params$sig2F = 1/q * abs(sum(diag(Y%*%t(Y)))/N - 2*sum(diag(mu_U%*%t(params$C)%*%t(Y)))/N +
    #                        sum(diag(params$SigU)) - sum(diag(params$SigUo)))
    
    # SVD for high dimensional
    if(orth_x){
      dcmp_varxo <- svd_orthpart(X, params$W, params$sig2E, mu_T, Stt, rx)
      params$Wo = dcmp_varxo$V
      params$SigTo = diag(x=dcmp_varxo$D, nrow=rx)
    }
    if(orth_y){
      dcmp_varyo <- svd_orthpart(Y, params$C, params$sig2F, mu_U, Suu, ry)
      params$Co = dcmp_varyo$V
      params$SigUo = diag(x=dcmp_varyo$D, nrow=ry)
    }
    
    return(params)
  })
}


params_generate <- function(method = c("random", "method2", "method3")){
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
  return(params)
}
