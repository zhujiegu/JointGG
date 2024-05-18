library(nlme)

get_reffects_var <- function(model_fit , dat_i){
  corr.mat=diag(1,nrow(dat_i)) 
  resid.var=(model_fit$sigma)^2
  
  x.mat=model.matrix(model_fit,data=dat_i)  #fixed effects design matrix
  coef.fixed=matrix(model_fit$coef$fixed,nrow=length(model_fit$coef$fixed),ncol=1)
  z.mat=model.matrix(formula(model_fit$modelStruct$reStr)[[1]],dat_i) #random effects design matrix
  t.z.mat=t(z.mat)
  G=unclass(getVarCov(model_fit))
  
  # # Check BLUPs
  # (G%*%t.z.mat)%*%solve((z.mat%*%G%*%t.z.mat+resid.var*corr.mat))%*%(dat_perID_y[[i]]$value-x.mat%*%coef.fixed)
  
  # variance of random effect
  var.reffects = G - (G%*%t.z.mat)%*%
    solve((z.mat%*%G%*%t.z.mat+resid.var*corr.mat))%*%(z.mat%*%t(G))
  return(var.reffects)
}


##############################################################################
# fit survival data
##############################################################################

#########################################
