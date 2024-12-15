library(tibble)
library(MASS)
library(dplyr)
library(magrittr)
library(flexsurv)
library(purrr)
library(ggplot2)
select <- dplyr::select

#########################################################################
# Joint model simulation
#########################################################################
# n   : sample size
# seed : seed for simulation
# baseline: baseline survival time distribution. Generalised Gamma (GG), lambda*exp (lexp)
# Model parameters: sigmaeps2, G, alpha, beta, lambda

simData_GG <- function(n, seed, visits_age, cens_time, G, a0, a1, a2, sig_e2, beta_mu, beta_sigma, beta_q){
  set.seed(seed)
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
  dat_y <- tibble(dat_t %>% select(ID:b1, times), exposure='SBP')
  dat_y %<>% cross_join(tibble(visits_age))
  dat_y %<>% filter(visits_age <= times) 
  dat_y %<>% mutate(value=(a0+b0+(a1+b1)*visits_age+a2*treat*visits_age +
                             rnorm(length(visits_age), 0,sqrt(sig_e2))))
  # plot
  dat_y$Group <- factor(dat_y$treat)
  p_y <- ggplot(dat_y, aes(x = visits_age, y = value, group = ID, color = Group)) +
    geom_line() +
    geom_point(size = 1.5) +
    scale_color_discrete(labels = c("0" = "Control", "1" = "Treatment")) +
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
  
  params_true <- list(
    G=G,
    a0=a0,
    a1=a1,
    a2=a2,
    sig_e2=sig_e2,
    beta_mu=beta_mu,
    beta_sigma=beta_sigma,
    beta_q=beta_q
  )
  
  obj <- list(dat_y, dat_t, dat_b, dat_gg, p_y, params_true)
  names(obj) <- c("longitudinal","survival",'true_b','true_GG','plot_longitudinal', 'params')
  
  return(obj)
}


simData_seq <- function(n, seed, visits_age, recruitment_time, interim_times,
                        G, a0, a1, a2, sig_e2, beta_mu, beta_sigma, beta_q, cens_interval=NULL, log_curve=F){
  set.seed(seed)
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
  
  # generate start (entering) and end time (calendar time)
  dat_t %<>% mutate(start = sort(runif(n, 0, recruitment_time))) %>% mutate(end = start+surv_t)
  
  if(!is.null(cens_interval)){
    cens_time <- runif(n, cens_interval[1], cens_interval[2])
    # event indicator
    dat_t %<>% mutate(cens_time=cens_time) %>% 
      rowwise %>% mutate(times=min(surv_t,cens_time)) %>% ungroup
  }else{
    dat_t %<>% mutate(cens_time=Inf)
  }
  
  ##############################
  # Data for interim analysis 
  ##############################
  dat <- lapply(interim_times, function(k) generate_dat_interim(dat_t, k, a0, a1, a2, sig_e2, visits_age, log_curve))
  names(dat) <- c(paste0('interim', 1:length(interim_times)))
  
  dat$true_b <- dat_t %>% select(ID,b0,b1)
  dat$true_GG <- dat_t %>% select(ID,mu,sigma,q)
  
  params_true <- list(
    G=G,
    a0=a0,
    a1=a1,
    a2=a2,
    sig_e2=sig_e2,
    beta_mu=beta_mu,
    beta_sigma=beta_sigma,
    beta_q=beta_q
  )
  dat$params <- params_true
  ##############################
  # plot
  ##############################
  # dat_y$Group <- factor(dat_y$treat)
  # p_y <- ggplot(dat_y, aes(x = visits_age, y = value, group = ID, color = Group)) +
  #   geom_line() +
  #   geom_point(size = 1.5) +
  #   labs(x = "Time (years)", y = "Systolic Blood Pressure (mmHg)") +
  #   theme_minimal()
  
  # print(paste0('proportion of cencering: ', 1-mean(dat_t$status)))
  return(dat)
}

generate_dat_interim <- function(dat_t, interim_time, a0, a1, a2, sig_e2, visits_age, log_curve){
  dat_t %<>% filter(start<=interim_time)
  dat_t %<>% rowwise %>% mutate(status=ifelse(end<min(interim_time, cens_time), 1, 0), cens_t=min(interim_time, cens_time)) %>% 
    mutate(times=min(surv_t,cens_t-start)) %>% ungroup
  # Longitudinal process
  dat_y <- tibble(dat_t %>% select(ID:b1, start, end, cens_t), exposure='SBP')
  dat_y %<>% cross_join(tibble(visits_age)) %>% mutate(visits_time = visits_age+start)
  dat_y %<>% rowwise %>% filter(visits_time <= min(end, cens_t)) %>% ungroup
  if(log_curve){
    dat_y %<>% mutate(value=b0 + b1*visits_age+a2*treat*visits_age+ 2.5*log(visits_age+1.7)+rnorm(length(visits_age), 0,sqrt(sig_e2)))
  }else{
    dat_y %<>% mutate(value=(a0+b0+(a1+b1)*visits_age+a2*treat*visits_age +
                               rnorm(length(visits_age), 0,sqrt(sig_e2))))
  }
  dat_y %<>% select(ID, treat, exposure, visits_age, value)
  dat_t %<>% select(ID, treat, status, times)
  return(list(longitudinal=dat_y, survival=dat_t))
}

# generate data from Joint model with Cox Weibull baseline
generate_dat_CoxWeibull <- function(n, seed, visits_age, recruitment_time, interim_times, G, a0, a1, a2, sig_e2, 
                                    shape, scale, association_y, association_fix, cens_interval=NULL){
  set.seed(seed)
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
  u <- runif(n)
  
  invS <- function(t,u,i,...){
    h <- function(s,...){
      XX <- cbind(s,s*x_i[i])
      ZZ <- cbind(1,s)
      fval <- as.vector(XX%*%c(a1,a2) + rowSums(ZZ*b[rep(i,nrow(ZZ)),]))
      hh <- exp(log(scale)+log(shape)+(shape-1)*log(s)+association_fix*x_i[i]+association_y*fval)
      return(hh)
    }
    return(integrate(h,lower=0,upper=t)$value + log(u))
  }
  
  trueTime <- NULL
  for(i in 1:n){
    Up <- 50
    tries <- 5
    Root <- try(uniroot(invS,interval=c(1e-05,Up),u=u[i],i=i)$root,TRUE)
    while(inherits(Root,"try-error") && tries>0){
      tries <- tries - 1 
      Up <- Up + 50
      Root <- try(uniroot(invS,interval=c(1e-05,Up),u=u[i],i=i)$root,TRUE)	
    }
    trueTime[i] <- ifelse(!inherits(Root,"try-error"), Root, NA)	
  }
  
  # True time, censoring, observed time and event indicator
  trueTime <- ifelse(is.na(trueTime),interim_times[length(interim_times)],trueTime)
  
  # mutate survival time
  dat_t %<>% mutate(surv_t = trueTime)
  
  # generate start (entering) and end time (calendar time)
  dat_t %<>% mutate(start = sort(runif(n, 0, recruitment_time))) %>% mutate(end = start+surv_t)
  
  if(!is.null(cens_interval)){
    cens_time <- runif(n, cens_interval[1], cens_interval[2])
    # event indicator
    dat_t %<>% mutate(cens_time=cens_time) %>% 
      rowwise %>% mutate(times=min(surv_t,cens_time)) %>% ungroup
  }else{
    dat_t %<>% mutate(cens_time=Inf)
  }
  
  ##############################
  # Data for interim analysis 
  ##############################
  dat <- lapply(interim_times, function(k) generate_dat_interim(dat_t, k, a0, a1, a2, sig_e2, visits_age))
  names(dat) <- c(paste0('interim', 1:length(interim_times)))
  
  dat$true_b <- dat_t %>% select(ID,b0,b1)
  
  params_true <- list(
    G=G,
    a0=a0,
    a1=a1,
    a2=a2,
    sig_e2=sig_e2,
    shape=shape, 
    scale=scale, 
    association_y=association_y, 
    association_fix=association_fix
  )
  dat$params <- params_true
  return(dat)
}

# bootstrap
generate_dat_bootstrap <- function(dat, seed, return.ID=F){
  set.seed(seed)
  # creating the bootstrap data
  IDs <- dat$survival$ID
  n=length(IDs)
  IDs_bootstr <- sample(IDs, size = n, replace = T)
  
  dat_t <- lapply(IDs_bootstr, function(i) dat$survival %>% filter(ID==i))
  dat_y <- lapply(IDs_bootstr, function(i) dat$longitudinal %>% filter(ID==i))
  for (i in 1:n){
    dat_t[[i]]$ID <- paste0('ID', i)
    dat_y[[i]]$ID <- paste0('ID', i)
  }
  dat_t <- do.call(rbind, dat_t)
  dat_y <- do.call(rbind, dat_y)
  dat = list(survival=dat_t, longitudinal=dat_y)
  if(return.ID){
    return(list(dat=dat, IDs=IDs_bootstr))
  }else{
    return(dat)
  }
}

