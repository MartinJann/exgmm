#' @importFrom stats lm pchisq pf qt sd var
#calculate externally informed GMM estimator of linear model parameters for external values
calc_beta_gmm<-function(X,y,E,mom,Omega_R_X,Omega_h,positions){
  h<-calc_h(X,y,E,mom,positions)
  solve(t(X)%*%X)%*%t(X)%*%y-Omega_R_X%*%solve(Omega_h)%*%t(h)%*%rep(1,nrow(X))
}

#calculate variance of externally informed GMM estimator of linear model parameters for external values
calc_var_gmm<-function(X,y,Omega_R_X,Omega_h,sigma){
  sigma^2*solve(t(X)%*%X)-nrow(X)*Omega_R_X%*%solve(Omega_h)%*%t(Omega_R_X)
}
