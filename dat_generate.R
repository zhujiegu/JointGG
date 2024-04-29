library(tibble)
library(MASS)
library(dplyr)
library(magrittr)
library(flexsurv)
library(purrr)

#########################################################################
# Joint model simulation
#########################################################################
# n   : sample size
# seed : seed for simulation
# baseline: baseline survival time distribution. Generalised Gamma (GG), lambda*exp (lexp)
# Model parameters: sigmaeps2, omega, G, alpha, beta, lambda

n=100
seed=1

omega=c(0,0)
G=diag(c(0.2,0.2))
a0=0.4
a1=0.15
a2=-0.2
sig_e2=0.01
beta_mu=c(1,-0.5,2,1)
beta_sigma=c(1,-0.4,1,0.5)
beta_q=c(1,-0.8,1,1)
visits_time=c(seq(5,60,by=5))

cens_time=60

simDataJ <- function(n, seed, omega, G, a0, a1, a2, sig_e2, beta_mu, beta_sig, beta_q){
  
  #################
  # Treatment assign
  #################
  x_i = rbinom(n,1,0.5)

  #################
  # Random effects
  #################
  b <- mvrnorm(n, omega, G)
  colnames(b) <- c('b0','b1')
  #################
  # data matrix
  dat_t <- tibble(ID=paste0('ID',1:n),treat=x_i) %>%
    bind_cols(as_tibble(b))
  
  ###################
  # Survival process
  ###################
  # mu, sigma, q for GG distribution
  z <- dat_t %>% 
    mutate(intercept = 1) %>%  # Add 'intercept' column
    select(intercept, treat, b0, b1) %>% as.matrix()
  dat_t %<>% mutate(mu=c(z%*%beta_mu), sigma=c(exp(z%*%beta_sigma)), q=c(exp(z%*%beta_q)))
  
  # generate survival time
  dat_t %<>%
    mutate(surv_t = pmap(list(mu, sigma, q), function(mu, sigma, q){
    rgengamma(1, mu=mu, sigma=sigma, Q = q)
  }) %>% unlist)
  
  # event indicator
  dat_t %<>% mutate(status=ifelse(surv_t<cens_time, 1, 0), cens_t=cens_time) %>% 
    rowwise %>% mutate(times=min(surv_t,cens_t)) %>% ungroup
  
  ##############################
  # Longitudinal process  
  ##############################
  dat_y <- tibble(dat_t %>% select(ID:b1), exposure='BMI')
  dat_y %<>% cross_join(tibble(visits_time))
  dat_y %<>% mutate(value=(a0+b0+(a1+b1)*visits_time+a2*treat*visits_time +
                     rnorm(length(visits_time), 0,sqrt(sig_e2))))
  
  #---------------------------------------------------------------------
  # Creating the longitudinal and survival processes object
  #---------------------------------------------------------------------
  obj <- list(dat_y, dat_t)
  names(obj) <- c("longitudinal","survival")
  
  return(obj)
}