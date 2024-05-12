library(nlme)
library(fitdistrplus)

##############################################################################
# fit longitudinal alone
##############################################################################
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

get_reffects_var <- function(model_fit , dat_i){
  corr.mat=diag(1,nrow(dat_i)) 
  resid.var=(model_fit$sigma)^2
  
  x.mat=model.matrix(model_fit,data=dat_i)  #fixed effects design matrix
  coef.fixed=matrix(model_fit$coef$fixed,nrow=length(model_fit$coef$fixed),ncol=1)
  z.mat=model.matrix(formula(model_fit$modelStruct$reStr)[[1]],dat_i) #random effects design matrix
  t.z.mat=t(z.mat)
  G=as.matrix(getVarCov(model_fit))
  
  # # Check BLUPs
  # (G%*%t.z.mat)%*%solve((z.mat%*%G%*%t.z.mat+resid.var*corr.mat))%*%(dat_perID_y[[i]]$value-x.mat%*%coef.fixed)
  
  # variance of random effect
  var.reffects = G - (G%*%t.z.mat)%*%
    solve((z.mat%*%G%*%t.z.mat+resid.var*corr.mat))%*%(z.mat%*%t(G))
  return(var.reffects)
}

reffects.individual=ranef(fit_y)
reffects.individual <- split(reffects.individual, rownames(reffects.individual))
reffects.individual <- reffects.individual[IDs]
var.reffects=lapply(dat_perID_y, function(i) get_reffects_var(fit_y, i))

##############################################################################
# fit survival data
##############################################################################

#########################################
# testing
# y <- rgengamma(1000, mu=1, sigma=2, Q = 1)

dmy_gg <- function(x,mu,sig,q) dgengamma(x, mu = mu, sigma = sig, Q=q, log = F)
fit_t <- mledist(y,"my_gg",start=list(mu=1,sig=1,q=1), optim.method = 'BFGS')

