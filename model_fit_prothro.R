library(JMbayes)
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

data(prothro)
data(prothros)

dat_t <- prothros
dat_y <- prothro

dat_y %<>% select(id, pro, time, treat) %>% mutate(exposure='pro')
colnames(dat_y) <- c('ID', 'value', 'visits_time', 'treat', 'exposure')
colnames(dat_t) <- c('ID','times','status','treat')

dat_y$ID <- paste0('ID',dat_y$ID)
dat_t$ID <- paste0('ID',dat_t$ID)

dat_t %<>% mutate(treat=ifelse(treat=='prednisone', 1,0))
dat_y %<>% mutate(treat=ifelse(treat=='prednisone', 1,0))

dat=list(longitudinal=as_tibble(dat_y), survival=as_tibble(dat_t))

all.equal(unique(dat$longitudinal$ID), dat$survival$ID)
save(dat, file='prothro.RData')



source('dat_generate.R')
source('EM.R')
source('aGH_functions.R')
source('two_stage.R')
source('test_statistic.R')

fit <- JM_EM(dat, init_params='two-stage', tol=1e-3, steps=100, Nr.cores=1, model_complex='normal', GH_level=5)

test <- z_statistic(fit, up_limit = 80, GH_level_z=11)
var_test <- as.numeric(matrix(test$delta_mean_grt, nrow = 1) %*% solve(beta_hes_transform(fit$I_beta, fit$model_complex, 'collapse')) %*%
                         matrix(test$delta_mean_grt, ncol = 1))
z_score <- test$delta_mean_t/sqrt(var_test)
z_score
