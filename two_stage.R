library(nlme)

library(fitdistrplus)

fit_rirs = lme(value ~ 1 + visits_time + treat:visits_time,  
               random=~ 1 + visits_time|ID,
               weights=varIdent(form=~1|exposure),
               data = dat$longitudinal,
               control = lmeControl(maxIter=1000,msMaxIter=1000,msVerbose=TRUE,rel.tol=1e-8,
                                    msMaxEval=1000, niterEM = 50))


y <- rgengamma(50, mu=2, sigma=1.5, Q = 1)

dmy_gg <- function(x,mu,sig,q) dgengamma(x, mu = mu, sigma = sig, Q=q, log = F)
fit_t <- mledist(y,"my_gg",start=list(mu=3,sig=1,q=2), optim.method = 'BFGS')
fit_t
fit_t$hessian %>% solve()
