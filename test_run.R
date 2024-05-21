library(tibble)
library(MASS)
library(dplyr)
library(magrittr)
library(flexsurv)
library(purrr)
library(nlme)
library(parallel)
library(numDeriv)

source('dat_generate.R')
source('EM.R')
source('aGH_functions.R')
source('two_stage.R')

n=300
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


fit <- JM_EM(dat, init_params='two-stage', tol=1e-3, steps=50, Nr.cores=1, model_complex='normal', GH_level=5)

test <- z_statistic(fit, GH_level_z=9)

var_test <- matrix(test$delta_mean_grt, nrow = 1) %*% solve(beta_hes_transform(fit$I_beta, fit$model_complex, 'collapse')) %*%
  matrix(test$delta_mean_grt, ncol = 1)


#####################################################################################
# test the nodes using two-stage or JM
b_twostage <- Reduce('rbind', fit_y$reffects.individual) %>% as_tibble

test_b <- dat$true_b %>% right_join(fit$E_b, by='ID') %>% bind_cols(b_twostage)

cor(test_b$b0, test_b$Eb0)
cor(test_b$b0, test_b$`(Intercept)`)

cor(test_b$b1, test_b$Eb1)
cor(test_b$b1, test_b$visits_time)
cor(test_b$Eb1, test_b$visits_time)

test_b %>% select(-ID) %>% summarise(across(where(is.numeric), ~sd(.)))