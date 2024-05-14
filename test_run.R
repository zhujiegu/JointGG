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

n=200
seed=1
visits_time=c(seq(6,60,by=6))/12
cens_time=5

Nr.cores=1
GH_level=15
plot_nodes=T

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
model_complex='normal'

###############
#  fiting
fit_y = lme(value ~ 1 + visits_time + treat:visits_time,  
            random=~ 1 + visits_time|ID,
            weights=varIdent(form=~1|exposure),
            data = dat$longitudinal,
            control = lmeControl(maxIter=1000,msMaxIter=1000,msVerbose=TRUE,rel.tol=1e-8,
                                 msMaxEval=1000, niterEM = 50))

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
n_w_adj <- aGH_n_w(reffects.individual, var.reffects, Nr.cores=1, GH_level=9, dat, plot_nodes=F)

params <- params_true
aGH_results_ll(unlist(params), n_w_adj, dat, Nr.cores, model_complex='normal')

optim_fn <- function(params_vec){
  params <- params_vec_list(params_vec)
  l <- aGH_results_ll(params =params, n_w_adj=n_w_adj, dat=dat, Nr.cores=Nr.cores, 
                 model_complex='normal')
  return(l)
}

optim(par = unlist(params), fn = optim_fn, method = 'BFGS')
