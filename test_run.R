library(tibble)
library(MASS)
library(dplyr)
library(magrittr)
library(flexsurv)
library(purrr)
library(nlme)
library(parallel)
library(numDeriv)
select <- dplyr::select

source('dat_generate.R')
source('EM.R')
source('aGH_functions.R')
source('two_stage.R')
source('test_statistic.R')

n=100
seed=2
visits_age=c(seq(0,72,by=6))/12 # assume entering the trial at same age
cens_time=6

interim_times <- c(2,4,6)
recruitment_time = 3
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

# G=diag(c(25,1))
# a0=110
# a1=2
# a2=-4
# sig_e2=10
# beta_mu=c(1.5, 0.5, 0.01, 0.1)
# beta_sigma=c(0.4, -0.06,0,0)
# beta_q=c(-0.15,0,0,0)

params_specify <- list(
  G=diag(c(25,1)),
  a0=110,
  a1=2,
  a2=-4,
  sig_e2=10,
  beta_mu=c(1, 0.5, 0.01, 0.1),
  beta_sigma=c(0.4, -0.06,0,0),
  beta_q=c(-0.15,0,0,0)
)

params_true <- params_generate('specify', params_specify=params_specify)

dat <- do.call(simData_seq, c(list(n = n, seed = seed, visits_age=visits_age, recruitment_time=recruitment_time, 
                                   interim_times=interim_times), params_true))


fit <- JM_EM(dat, init_params='two-stage', tol=1e-3, steps=100, Nr.cores=1, model_complex='normal', 
             GH_level=5, refine_GH = T)

test <- z_statistic(fit, up_limit = 10, GH_level_z=11)
var_test <- as.numeric(matrix(test$delta_mean_grt, nrow = 1) %*% solve(beta_hes_transform(fit$I_beta, fit$model_complex, 'collapse')) %*%
  matrix(test$delta_mean_grt, ncol = 1))
z_score <- test$delta_mean_t/sqrt(var_test)
z_score
# test_n <- z_statistic_naive(fit, up_limit = 100)
# as.numeric(matrix(test_n$delta_mean_grt, nrow = 1) %*% solve(beta_hes_transform(fit$I_beta, fit$model_complex, 'collapse')) %*%
#   matrix(test_n$delta_mean_grt, ncol = 1))


#####################################################################################
# test the nodes using two-stage or JM
b_twostage <- Reduce('rbind', fit_y$reffects.individual) %>% as_tibble

test_b <- dat$true_b %>% right_join(fit$E_b, by='ID') %>% bind_cols(b_twostage)

cor(test_b$b0, test_b$Eb0)
cor(test_b$b0, test_b$`(Intercept)`)

cor(test_b$b1, test_b$Eb1)
cor(test_b$b1, test_b$visits_age)
cor(test_b$Eb1, test_b$visits_age)

test_b %>% select(-ID) %>% summarise(across(where(is.numeric), ~sd(.)))