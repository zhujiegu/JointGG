library(tibble)
library(MASS)
library(dplyr)
library(magrittr)
library(flexsurv)
library(purrr)
library(nlme)

source('dat_generate.R')
source('EM.R')
source('aGH_functions.R')
source('two_stage.R')

n=200
seed=1
visits_time=c(seq(6,60,by=6))/12
cens_time=5

Nr.cores=10
GH_level=9

# params_specify <- list(
#   G=diag(c(0.2,0.2)),
#   a0=4,
#   a1=0.5,
#   a2=-0.2,
#   sig_e2=1,
#   beta_mu=c(1,-0.5,0.02,0.01),
#   beta_sigma=c(0.5,0.4,0,0),
#   beta_q=c(0.1,-0.3,0,0)
# )

params_specify <- list(
  G=diag(c(0.5,0.1)),
  a0=4,
  a1=0.5,
  a2=-0.2,
  sig_e2=1,
  beta_mu=c(1,0.5,0.8, 0.6),
  beta_sigma=c(0.5,0.5,0,0),
  beta_q=c(0.1,0,0,0)
)

params_true <- params_generate('specify', params_specify=params_specify)


dat <- do.call(simDataJ, c(list(n = n, seed = seed, visits_time=visits_time, cens_time=cens_time), params_true))

init_params <- params_true
tol=1e-3
steps=5
Nr.core=1
model_complex='normal'

###############################################################################
# two-stage fiting
###############################################################################
fit_y = lme(value ~ 1 + visits_time + treat:visits_time,  
            random=~ 1 + visits_time|ID,
            weights=varIdent(form=~1|exposure),
            data = dat$longitudinal,
            control = lmeControl(maxIter=1000,msMaxIter=1000,msVerbose=TRUE,rel.tol=1e-8,
                                 msMaxEval=1000, niterEM = 50))
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

################################
# extract BLUPs and variance of random effects
################################
ll_param <- function(params_vec){
  beta_mu=c(params_vec[1], params_vec[2],params_vec[3],params_vec[4])
  beta_sigma=c(params_vec[5],params_vec[6],0,0)
  beta_q=c(params_vec[7],0,0,0)
  return(list(beta_mu=beta_mu, beta_sigma=beta_sigma, beta_q=beta_q))
}
param_vec_reduce <- function(full_vec){
  return(full_vec[c(1,2,3,4,5,6,9)])
}


ll_t <- function(params_vec){
  ll_param <- ll_param(params_vec)
  beta_mu <- ll_param$beta_mu
  beta_sigma <- ll_param$beta_sigma
  beta_q <- ll_param$beta_q
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

ll_grt_t <- function(params_vec){
  ll_param <- ll_param(params_vec)
  beta_mu <- ll_param$beta_mu
  beta_sigma <- ll_param$beta_sigma
  beta_q <- ll_param$beta_q
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
    phi_q=grad(function(q) log_upper_gamma_q(q, ww=ww), q)
    
    #################################################
    # #tesing
    # if(is.na(phi_w))  browser()
    # stop(cat(paste0('z=',z, ', mu=',mu, ', sigma=',sigma, ', q=',q, ', w=', ww, '\n', 
    #                 'values not realistic, use realistic data or use adaptive quadrature')))
    #################################################
    
    grd_mu[i,] <- if(dat_t[i,]$status==1){-z_mu/(q*sigma)*(1-exp(q*ww))} else{-z_mu*phi_w/sigma} 
    grd_sigma[i,] <- if(dat_t[i,]$status==1){-z_sigma*(1+ww/q*(1-exp(q*ww)))} else{-z_sigma*ww*phi_w}
    grd_q[i,] <- if(dat_t[i,]$status==1){z_q*q^-2*(q^2+2*digamma(q^-2)+4*log(q)-2-q*ww+2*exp(q*ww)-q*ww*exp(q*ww))}else{
      z_q*q*(phi_q+2*(q^-3)*digamma(q^-2))}
  }
  grd <- c(colSums(grd_mu), colSums(grd_sigma), colSums(grd_q))
  return(-param_vec_reduce(grd)) # - correspods to minus log ll
}


fit_t <- optim(rep(0.1,7),fn = ll_t, gr = ll_grt_t,  method = 'L-BFGS-B')

fit_t












reffects.individual_true=dat$true_b %>% select(-ID) %>% data.frame()
rownames(reffects.individual_true) <- IDs
colnames(reffects.individual_true) <- c('(Intercept)', 'visits_time')
reffects.individual_true <- split(reffects.individual_true, rownames(reffects.individual_true))


n_w_adj <- aGH_n_w(reffects.individual_true, var.reffects, Nr.cores=1, GH_level=GH_level, dat, plot_nodes=F)

params <- params_true



aGH_results_ll(unlist(params), n_w_adj, dat, Nr.cores, model_complex='normal')

optim_fn <- function(params_vec){
  params <- params_vec_list(params_vec)
  l <- aGH_results_ll(params =params, n_w_adj=n_w_adj, dat=dat, Nr.cores=Nr.cores, 
                 model_complex='normal')
  return(l)
}

result <- optim(par = unlist(params), fn = optim_fn, method = 'BFGS', control = list(maxit=3))
