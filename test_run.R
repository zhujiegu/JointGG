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


fit <- JM_EM(dat, init_params='two-stage', tol=1e-3, steps=20, Nr.cores=1, model_complex='normal', GH_level=7)










reffects.individual_true=dat$true_b %>% select(-ID) %>% data.frame()
rownames(reffects.individual_true) <- IDs
colnames(reffects.individual_true) <- c('(Intercept)', 'visits_time')
reffects.individual_true <- split(reffects.individual_true, rownames(reffects.individual_true))



params <- params_true



aGH_results_ll(unlist(params), n_w_adj, dat, Nr.cores, model_complex='normal')

optim_fn <- function(params_vec){
  params <- params_vec_list(params_vec)
  l <- aGH_results_ll(params =params, n_w_adj=n_w_adj, dat=dat, Nr.cores=Nr.cores, 
                 model_complex='normal')
  return(l)
}

result <- optim(par = unlist(params), fn = optim_fn, method = 'BFGS', control = list(maxit=3))
