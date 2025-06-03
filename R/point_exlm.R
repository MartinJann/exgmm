calc_beta_gmm<-function(X,y,E,mom,Omega_R_X,Omega_h,positions){
  h<-calc_h(X,y,E,mom,positions)
  solve(t(X)%*%X)%*%t(X)%*%y-Omega_R_X%*%solve(Omega_h)%*%t(h)%*%rep(1,nrow(X))
}

calc_var_gmm<-function(X,y,Omega_R_X,Omega_h,sigma){
  diag(sigma^2*solve(t(X)%*%X)-nrow(X)*Omega_R_X%*%solve(Omega_h)%*%t(Omega_R_X))
}
