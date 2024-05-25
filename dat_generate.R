library(tibble)
library(MASS)
library(dplyr)
library(magrittr)
library(flexsurv)
library(purrr)
library(ggplot2)

#########################################################################
# Joint model simulation
#########################################################################
# n   : sample size
# seed : seed for simulation
# baseline: baseline survival time distribution. Generalised Gamma (GG), lambda*exp (lexp)
# Model parameters: sigmaeps2, G, alpha, beta, lambda

simData_GG <- function(n, seed, visits_time, cens_time, G, a0, a1, a2, sig_e2, beta_mu, beta_sigma, beta_q){
  
  b_dim=nrow(G)
  #################
  # Treatment assign
  #################
  x_i = rbinom(n,1,0.5)

  #################
  # Random effects
  #################
  b <- mvrnorm(n, rep(0,b_dim), G)
  colnames(b) <- c('b0','b1')
  #################
  # data matrix
  dat_t <- tibble(ID=paste0('ID',1:n),treat=x_i) %>%
    bind_cols(as_tibble(b))
  ###################
  # Survival process
  ###################
  # mu, sigma, q for GG distribution
  z <- dat_t %>% mutate(intercept = 1) %>% select(intercept, treat, b0, b1) %>% as.matrix()
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
  dat_y <- tibble(dat_t %>% select(ID:b1, times), exposure='BMI')
  dat_y %<>% cross_join(tibble(visits_time))
  dat_y %<>% filter(visits_time <= times) 
  dat_y %<>% mutate(value=(a0+b0+(a1+b1)*visits_time+a2*treat*visits_time +
                     rnorm(length(visits_time), 0,sqrt(sig_e2))))
  # plot
  dat_y$Group <- factor(dat_y$treat)
  p_y <- ggplot(dat_y, aes(x = visits_time, y = value, group = ID, color = Group)) +
    geom_line() +
    geom_point(size = 1.5) +
    labs(x = "Time (years)", y = "Systolic Blood Pressure (mmHg)") +
    theme_minimal()
  #---------------------------------------------------------------------
  # Creating the longitudinal and survival processes object
  #---------------------------------------------------------------------
  dat_b <- dat_t %>% select(ID,b0,b1)
  dat_gg <- dat_t %>% select(ID,mu,sigma,q)
  dat_y %<>% select(-b0,-b1,-Group)
  dat_t %<>% select(-b0,-b1,-mu,-sigma,-q)
  
  print(paste0('proportion of cencering: ', 1-mean(dat_t$status)))
  
  obj <- list(dat_y, dat_t, dat_b, dat_gg, p_y)
  names(obj) <- c("longitudinal","survival",'true_b','true_GG','plot_longitudinal')
  
  return(obj)
}

# flexsurvreg(formula = Surv(Time, death) ~ 1, data = prothros, dist="gengamma")
# 
# flexsurvreg(formula = Surv(Time, death) ~ 1, data = filter(prothros, treat=='prednisone'), dist="gengamma")
# flexsurvreg(formula = Surv(Time, death) ~ 1, data = filter(prothros, treat=='placebo'), dist="gengamma")

