library(tibble)
library(MASS)
library(dplyr)
library(magrittr)
library(flexsurv)
library(purrr)
source('dat_generate.R')
source('EM.R')
source('GH_functions.R')

n=5
seed=1
visits_time=c(seq(6,60,by=6))/12
cens_time=5

Nr.cores=1
GH_level=9
plot_nodes=T

params_true <- params_generate('random')
dat <- do.call(simDataJ, c(list(n = n, seed = seed), params_true))


init_params <- params_true
