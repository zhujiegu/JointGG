library(tibble)
library(MASS)
library(dplyr)
library(magrittr)
library(flexsurv)
library(purrr)
source('dat_generate.R')
source('EM.R')
source('GH_functions.R')

n=50
seed=1
visits_time=c(seq(6,60,by=6))/12
cens_time=5

Nr.cores=1
GH_level=13
plot_nodes=T

params_specify <- list(
  G=diag(c(0.2,0.2)),
  a0=4,
  a1=0.5,
  a2=-0.2,
  sig_e2=1,
  beta_mu=c(1,-0.5,0.02,0.01),
  beta_sigma=c(0.5,0.4,0,0),
  beta_q=c(0.1,-0.3,0,0)
)

params_specify <- list(
  G=diag(c(0.2,0.2)),
  a0=4,
  a1=0.5,
  a2=-0.2,
  sig_e2=1,
  beta_mu=c(1, 0,0,0),
  beta_sigma=c(0.5,0,0,0),
  beta_q=c(0.3,0,0,0)
)


params_true <- params_generate('specify', params_specify=params_specify)


dat <- do.call(simDataJ, c(list(n = n, seed = seed), params_true))


init_params <- params_true
tol=1e-3
steps=3
Nr.core=1
complex='test'
