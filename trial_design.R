library(mvtnorm)

# simulate data to estimate the relationship between calendar time and information level
t_info_relation <- function(N, eval_t, params, seed, up_limit, Nr.cores, recruitment_time=1, visits_age=c(seq(0,120,by=6))/12,
                            GH_level=7, steps=100, model_complex='normal'){
  
  seed=seed
  dat <- do.call(simData_seq, c(list(n = N, seed = seed, visits_age=visits_age, recruitment_time=recruitment_time, 
                                     interim_times=eval_t), params))
  
  # JM model fitting
  fit <- JM_EM(dat$interim1, init_params='two-stage', rel_tol=1e-5, steps=steps, Nr.cores=Nr.cores, model_complex=model_complex, 
               GH_level=GH_level, refine_GH = T)
  delta_RMST <- z_statistic(fit, up_limit = up_limit, GH_level_z=11)
  
  # using only survival data to fit a GG distribution
  # var_RMST_GG <- get_RMST_var_GG(fit, up_limit)$var_delta
  
  return(list(fit=fit, delta_RMST=delta_RMST))
}


# variance-covariance matrix of standardised test statistics given the vector of ratios I_k/I_max
generate_sig <- function(I_vec) {
  K = length(I_vec)
  sig <- matrix(0, nrow = K, ncol = K)
  for (i in 1:K) {
    for (j in 1:K) {
      sig[i, j] <- sqrt(min(I_vec[i],I_vec[j]) / max(I_vec[i],I_vec[j]))
    }
  }
  return(sig)
}

# boundaries (non-binding futility) and maximum information based on error spending functions, type I, II error, pattern of I_k
# I_vec is the vector of ratios I_k/I_max
boundaries <- function(I_vec, alpha, beta, delta, pho){
  K = length(I_vec)
  sig_z <- generate_sig(I_vec)
  a <- b <- c()
  err1 <- err2 <- c()
  ###########################
  # upper boundaries
  for(k in 1:K){
    if(k==1){
      error1 <- error_spend_f(I_vec[k], alpha, pho)
      b[k] <- -qnorm(error1)
    }else{
      error1 <- error_spend_f(I_vec[k], alpha, pho)-error_spend_f(I_vec[k-1], alpha, pho)
      b[k] <- uniroot(f = function(b_k) pmvnorm(rep(0,k), corr = sig_z[1:k,1:k], 
                                                lower=c(rep(-Inf,k-1),b_k), upper=c(b[1:k-1],Inf))-error1,
                      interval=c(0,5))$root
    }
    err1[k] <- error1
  }
  # #test
  # 1-pmvnorm(rep(0,K), corr=sig_z, lower=-Inf, upper=b, seed = 1)
  ###########################
  # lower boundaries
  # a given Imax
  a_Imax <- function(Imax){
    mean_vec <- delta*sqrt(Imax*I_vec)
    for(k in 1:K){
      if(k==1){
        error2 <- error_spend_g(I_vec[k], beta, pho)
        a[k] <- qnorm(error2, mean = mean_vec[1], sd = 1)
      }else{
        error2 <- error_spend_g(I_vec[k], beta, pho)-error_spend_g(I_vec[k-1], beta, pho)
        a[k] <- uniroot(f = function(a_k) pmvnorm(mean_vec[1:k], corr = sig_z[1:k,1:k], 
                                                  lower=c(a[1:k-1],-Inf), upper=c(b[1:k-1],a_k))-error2,
                        interval=c(-2,5))$root
      }
    }
    return(a)
  }
  
  Imax_a <- function(Imax){
    a <- a_Imax(Imax)
    return(a[K]-b[K])
  }
  I_fix <- ((-qnorm(alpha)-qnorm(beta))/delta)^2
  Imax <- uniroot(Imax_a, interval = c(I_fix, 1.5*I_fix))$root
  
  a <- a_Imax(Imax)
  for (k in 1:K) {
    err2[k] <- ifelse(k==1, error_spend_g(I_vec[k], beta, pho), 
                      error_spend_g(I_vec[k], beta, pho)-error_spend_g(I_vec[k-1], beta, pho))
  }
  # #test
  # 1-pmvnorm(delta*sqrt(Imax*I_vec), corr = sig_z, lower=a, upper=Inf)
  return(list(a=a,b=b,I_max=Imax, err1=err1,err2=err2))
}

####################################
# #test
# I_vec <- c(1/3,2/3,1)
# alpha=0.025
# beta=0.5
# delta=1
# pho=1
# bb <- boundaries(I_vec, alpha, beta=0.2, delta, pho)
# 1-pmvnorm(rep(0,K), corr=sig_z, lower=-Inf, upper=bb$b, seed = 1) # type I
# 1-pmvnorm(delta*sqrt(bb$I_max*I_vec), corr = sig_z, lower=bb$a, upper=Inf) #type II
# abi <- error_spending_func_1s(
#     K = 3,           #number of groups
#     alpha,       #required significance
#     beta=0.2,        #required power  
#     delta,       #power requirement
#     I=I_vec*bb$I_max,         #Vector of observed information levels
#     I_max=bb$I_max,       #maximum information
#     
#     #Choice of form for the error spending functions - taking input t which is the
#     #obtained fraction of the maximum information.
#     error_spend_f,
#     error_spend_g)
# 
# 1-pmvnorm(rep(0,K), corr=sig_z, lower=-Inf, upper=abi[[2]], seed = 1)
# 1-pmvnorm(delta*sqrt(bb$I_max*I_vec), corr = sig_z, lower=abi[[1]], upper=Inf)


# error spending functions
error_spend_f = function(t, alpha, pho=1){
  return(alpha * min(t^pho,1))}
error_spend_g = function(t, beta, pho=1){
  return(beta * min(t^pho,1))}
