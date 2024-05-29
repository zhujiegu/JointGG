

###############################################################################
# two-stage fiting
###############################################################################
# fitting longitudinal
fit_longitudinal <- function(dat){
  print('start fitting longitudinal data in two-stage')
  fit_y = suppressMessages(lme(value ~ 1 + visits_age + treat:visits_age,  
                               random=~ 1 + visits_age|ID,
                               weights=varIdent(form=~1|exposure),
                               data = dat$longitudinal,
                               control = lmeControl(maxIter=1000,msMaxIter=1000,msVerbose=TRUE,rel.tol=1e-8,
                                                    msMaxEval=1000, niterEM = 50)))
  # sometimes this does not converge, use try in pipeline
  ################################
  # extract BLUPs and variance of random effects
  ################################
  IDs <- unique(dat$longitudinal$ID)
  n <-  IDs %>% length()
  dat_perID_y <- split(dat$longitudinal, dat$longitudinal$ID)
  dat_perID_y <- dat_perID_y[IDs]
  
  reffects.individual=ranef(fit_y)
  reffects.individual <- split(reffects.individual, rownames(reffects.individual))
  reffects.individual <- reffects.individual[IDs]
  var.reffects=lapply(dat_perID_y, function(i) get_reffects_var(fit_y, i))
  return(list(
    a0=fit_y$coefficients$fixed[1],
    a1=fit_y$coefficients$fixed[2],
    a2=fit_y$coefficients$fixed[3],
    sig_e2=(fit_y$sigma)^2,
    G=Reduce('+', var.reffects)/n,
    reffects.individual=reffects.individual,
    var.reffects=var.reffects
  ))
}



################################
# extract BLUPs and variance of random effects
################################
fit_survival <- function(dat, model_complex, reffects.individual){
  print('start fitting survival data in two-stage')
  ###############################
  # define which betas need to be updated based on the complexity of the model
  if(model_complex=='normal'){
    # Initialize parameters
    # browser()
    fit_t_no_reff <- optim(rep(0.8,5),fn = function(vec) ll_t(vec, dat, model_complex, reffects.individual, rm_reff=T), 
                           gr = function(vec) ll_grt_t(vec, dat, model_complex, reffects.individual, rm_reff=T),  method = 'L-BFGS-B')
    fit_t_no_reff<- c(fit_t_no_reff$par[1:2],0,0, fit_t_no_reff$par[3:5])
    # Start optimization with error handling
    fit_t <- tryCatch({
      fit_t_tmp <- optim(fit_t_no_reff, fn = function(vec) ll_t(vec, dat, model_complex, reffects.individual), 
                         gr = function(vec) ll_grt_t(vec, dat, model_complex, reffects.individual), method = 'L-BFGS-B')
      if (!is.null(fit_t_tmp$convergence) && fit_t_tmp$convergence == 0) {
        fit_t_tmp$par
      } else {
        print("Optimization in two-stage survival failed, method: L-BFGS-B")
        fit_t_no_reff
      }
    }, error = function(e) {
      print("Optimization in two-stage survival failed, method: L-BFGS-B")
      fit_t_no_reff
    })
  }
  if(model_complex=='saturated'){
    stop('two-step to be completed')
  }
  print('finish fitting survival data in two-stage')
  betas <- beta_vec_to_param(fit_t, model_complex)
  print('initial betas_mu:') 
  print(betas$beta_mu)
  print('initial betas_sigma:') 
  print(betas$beta_sigma)
  print('initial betas_q:') 
  print(betas$beta_q)
  
  return(list(beta_mu=betas$beta_mu, beta_sigma=betas$beta_sigma, beta_q=betas$beta_q,
              beta_noreff = beta_vec_to_param(fit_t_no_reff, model_complex)
  ))
}


ll_t <- function(params_vec, dat, model_complex, reffects.individual, rm_reff=F){
  betas <- beta_vec_to_param(params_vec, model_complex, rm_reff=rm_reff)
  beta_mu <- betas$beta_mu
  beta_sigma <- betas$beta_sigma
  beta_q <- betas$beta_q
  dat_t <- dat$survival
  n=nrow(dat_t)
  f_t <- c()
  for (i in 1:n){
    if(model_complex=='saturated'){
      z_mu =z_sigma = z_q = cbind(1, dat_t[i,]$treat, reffects.individual[[i]]) %>% as.matrix
      # z_mu =z_sigma = z_q = cbind(1, dat_t[i,]$treat, dat$true_b[i,2:3]) %>% as.matrix
    }
    if(model_complex=='normal'){
      z_mu = cbind(1, dat_t[i,]$treat,reffects.individual[[i]]) %>% as.matrix
      z_sigma = cbind(1, dat_t[i,]$treat, 0, 0) %>% matrix(nrow=1)
      z_q = c(1, 0,0,0)%>% matrix(nrow=1)
    }
    if(model_complex=='test'){
      z_mu =z_sigma = z_q = c(1, 0,0,0)%>% matrix(nrow=1)
    }
    mu <- z_mu %*% beta_mu %>% as.numeric
    sigma <- exp(z_sigma %*% beta_sigma) %>% as.numeric
    q <- exp(z_q %*% beta_q) %>% as.numeric
    f_t[i] <- ifelse(dat_t$status[i]==1, flexsurv::dgengamma(dat_t$times[i], mu = mu, sigma = sigma, Q=q, log = T),
                     flexsurv::pgengamma(dat_t$times[i], mu = mu, sigma = sigma, Q=q, lower.tail = F, log = T))
  }
  return(-sum(f_t))
}

ll_grt_t <- function(params_vec, dat, model_complex, reffects.individual, rm_reff=F){
  
  betas <- beta_vec_to_param(params_vec, model_complex, rm_reff=rm_reff)
  beta_mu <- betas$beta_mu
  beta_sigma <- betas$beta_sigma
  beta_q <- betas$beta_q
  dat_t <- dat$survival
  n=nrow(dat_t)
  grd_mu <- grd_sigma <-grd_q <- matrix(NA, nrow=n, ncol=4)
  for (i in 1:n){
    if(model_complex=='saturated'){
      z_mu =z_sigma = z_q = cbind(1, dat_t[i,]$treat, reffects.individual[[i]]) %>% as.matrix
      # z_mu =z_sigma = z_q = cbind(1, dat_t[i,]$treat, dat$true_b[i,2:3]) %>% as.matrix
    }
    if(model_complex=='normal'){
      z_mu = cbind(1, dat_t[i,]$treat,reffects.individual[[i]]) %>% as.matrix
      z_sigma = cbind(1, dat_t[i,]$treat, 0, 0) %>% matrix(nrow=1)
      z_q = c(1, 0,0,0)%>% matrix(nrow=1)
    }
    if(model_complex=='test'){
      z_mu =z_sigma = z_q = c(1, 0,0,0)%>% matrix(nrow=1)
    }
    mu <- z_mu %*% beta_mu %>% as.numeric
    sigma <- exp(z_sigma %*% beta_sigma) %>% as.numeric
    q <- exp(z_q %*% beta_q) %>% as.numeric
    ww=(log(dat_t[i,]$times)-mu)/sigma
    
    phi_w=-q/pracma::incgam(q^-2*exp(q*ww), q^-2) *(q^-2)^(q^-2) * exp((q^-2)*(q*ww-exp(q*ww)))
    if(!is.finite(log_upper_gamma_q(q, ww=ww))){phi_q=NA}else{
      phi_q=grad(function(q) log_upper_gamma_q(q, ww=ww), q)
    }
    
    #################################################
    # #tesing
    # if(is.na(phi_w))  browser()
    # stop(cat(paste0('z=',z, ', mu=',mu, ', sigma=',sigma, ', q=',q, ', w=', ww, '\n', 
    #                 'values not realistic, use realistic data or use adaptive quadrature')))
    #################################################
    if(!is.finite(phi_w)){
      warning('unrealistic gradient of Gamma(a,x) wrt w evaluated at some node, use realistic data or use adaptive quadrature')
      grd_mu[i,] = grd_sigma[i,] = z_mu*0
    }else{
      grd_mu[i,] <- if(dat_t[i,]$status==1){-z_mu/(q*sigma)*(1-exp(q*ww))} else{-z_mu*phi_w/sigma} 
      grd_sigma[i,] <- if(dat_t[i,]$status==1){-z_sigma*(1+ww/q*(1-exp(q*ww)))} else{-z_sigma*ww*phi_w}
    }
    if(!is.finite(phi_q)){
      warning('unrealistic gradient of Gamma(a,x) wrt q evaluated at some node, use realistic data or use adaptive quadrature')
      grd_q[i,] <- z_q*0
    }else{
      grd_q[i,] <- if(dat_t[i,]$status==1){z_q*q^-2*(q^2+2*digamma(q^-2)+4*log(q)-2-q*ww+2*exp(q*ww)-q*ww*exp(q*ww))}else{
        z_q*q*(phi_q+2*(q^-3)*digamma(q^-2))}
    }
  }
  grd <- c(colSums(grd_mu), colSums(grd_sigma), colSums(grd_q))
  grd_reduce <- beta_vec_transform(grd, model_complex, 'collapse', rm_reff)
  return(-grd_reduce) # - correspods to minus log ll
}


get_reffects_var <- function(model_fit , dat_i){
  corr.mat=diag(1,nrow(dat_i)) 
  resid.var=(model_fit$sigma)^2
  
  x.mat=model.matrix(model_fit,data=dat_i)  #fixed effects design matrix
  coef.fixed=matrix(model_fit$coef$fixed,nrow=length(model_fit$coef$fixed),ncol=1)
  z.mat=model.matrix(formula(model_fit$modelStruct$reStr)[[1]],dat_i) #random effects design matrix
  t.z.mat=t(z.mat)
  G=unclass(getVarCov(model_fit))
  
  # # Check BLUPs
  # (G%*%t.z.mat)%*%solve((z.mat%*%G%*%t.z.mat+resid.var*corr.mat))%*%(dat_perID_y[[i]]$value-x.mat%*%coef.fixed)
  
  # variance of random effect
  var.reffects = G - (G%*%t.z.mat)%*%
    solve((z.mat%*%G%*%t.z.mat+resid.var*corr.mat))%*%(z.mat%*%t(G))
  return(var.reffects)
}

beta_vec_to_param <- function(params_vec, model_complex, rm_reff=F){
  if(model_complex=='normal'){
    if(rm_reff){
      beta_mu=c(params_vec[1], params_vec[2],0,0)
      beta_sigma=c(params_vec[3],params_vec[4],0,0)
      beta_q=c(params_vec[5],0,0,0)
    }else{
      beta_mu=c(params_vec[1], params_vec[2],params_vec[3],params_vec[4])
      beta_sigma=c(params_vec[5],params_vec[6],0,0)
      beta_q=c(params_vec[7],0,0,0) 
    }
  }
  return(list(beta_mu=beta_mu, beta_sigma=beta_sigma, beta_q=beta_q))
}

beta_vec_transform <- function(vec, model_complex, type=c('collapse', 'expand'), rm_reff=F){
  type=match.arg(type)
  if(model_complex=='normal'){
    keep <- c(1,2,3,4,5,6,9) #positions in c(beta_mu, beta_sigma, beta_q)
    if(type == 'expand'){
      result <- rep(0,12)
      result[keep] <- vec
    }else{
      result <- vec[keep]
    }
  }else{
    stop('add other model complexity')
  }
  if(rm_reff) result <- result[-(3:4)]
  return(result)
}

# expand hessian of beta to full size (incerting zeros)
beta_hes_transform <- function(hes, model_complex, type=c('collapse', 'expand')){
  type=match.arg(type)
  if(model_complex=='normal'){
    keep <- c(1,2,3,4,5,6,9)
    if(type == 'expand'){
      result <- matrix(0, nrow=12, ncol=12)
      result[keep, keep] <- hes
    }else{
      result <- hes[keep, keep]
    }
  }else{
    stop('add other model complexity')
  }
  return(result)
}
