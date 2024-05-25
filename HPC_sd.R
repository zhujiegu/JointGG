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

n=1000
seed=1
visits_time=c(seq(0,60,by=6))/12
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
  G=diag(c(25,1)),
  a0=110,
  a1=2,
  a2=-4,
  sig_e2=10,
  beta_mu=c(1.8, 0.2, 0.01, 0.1),
  beta_sigma=c(0.2, -0.03,0,0),
  beta_q=c(-0.15,0,0,0)
)

params_true <- params_generate('specify', params_specify=params_specify)


dat <- do.call(simData_GG, c(list(n = n, seed = seed, visits_time=visits_time, cens_time=cens_time), params_true))


fit <- JM_EM(dat, init_params='two-stage', tol=1e-3, steps=100, Nr.cores=10, model_complex='normal', GH_level=5)