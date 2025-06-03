#' Confidence Union for a linear model
#'
#' @export

interval_gmm<-function(X,y,mom,I_min,I_max=I_min,grid_points=10001,intercept=TRUE){
  positions<-as.numeric(gsub("\\D","",mom))
  if(intercept==T){positions<-positions+1}
  mom[mom!="EY2"]<-gsub("[[:digit:]]+", "", mom[mom!="EY2"])
  beta_ols_h<-as.vector(lm(y~X)$coefficients)
  if(intercept==TRUE){X<-cbind(1,X)}
  int_ug<-rep(Inf,ncol(X))
  int_og<-rep(-Inf,ncol(X))
  beta_ug<-rep(Inf,ncol(X))
  beta_og<-rep(-Inf,ncol(X))
  var_ug<-rep(Inf,ncol(X))
  var_og<-rep(-Inf,ncol(X))
  g_l<-rep(list(0:(grid_points-1)),sum(I_max!=I_min))
  grid<-expand.grid(g_l)
  for (i_Ki in 1:(grid_points^sum(I_max!=I_min))) {
    E_I<-I_min
    E_I[I_max!=I_min]<-I_min[I_max!=I_min]+as.numeric(grid[i_Ki,])*(I_max-I_min)[I_max!=I_min]/(grid_points-1)
    sigma<-summary(lm(y~X[,-1]))[[6]]
    Omega_R<-calc_Omega_R_X(X,y,sigma,mom,beta_ols_h,positions)
    Omega_h<-calc_Omega_h(X,y,E_I,mom,sigma,beta_ols_h,positions)
    beta_gmm_dummy<-calc_beta_gmm(X,y,E_I,mom,Omega_R,Omega_h,positions)
    sigma<-sqrt(sum((y-X%*%beta_gmm_dummy)^2)/(nrow(X)-length(beta)))
    Omega_R<-calc_Omega_R_X(X,y,sigma,mom,beta_gmm_dummy,positions)
    Omega_h<-calc_Omega_h(X,y,E_I,mom,sigma,beta_gmm_dummy,positions)
    beta_gmm<-calc_beta_gmm(X,y,E_I,mom,Omega_R,Omega_h,positions)
    var_gmm<-calc_var_gmm(X,y,Omega_R,Omega_h,sigma)
    new_ug<-beta_gmm-qt(0.975,df=nrow(X)-ncol(X))*sqrt(var_gmm)
    new_og<-beta_gmm+qt(0.975,df=nrow(X)-ncol(X))*sqrt(var_gmm)
    if(sum(new_ug<int_ug)>0){
      int_ug[which(new_ug<int_ug)]<-new_ug[which(new_ug<int_ug)]
    }
    if(sum(new_og>int_og)>0){
      int_og[which(new_og>int_og)]<-new_og[which(new_og>int_og)]
    }
    if(sum(beta_gmm<beta_ug)>0){
      beta_ug[which(beta_gmm<beta_ug)]<-beta_gmm[which(beta_gmm<beta_ug)]
    }
    if(sum(var_gmm<var_ug)>0){
      var_ug[which(var_gmm<var_ug)]<-var_gmm[which(var_gmm<var_ug)]
    }
    if(sum(beta_gmm>beta_og)>0){
      beta_og[which(beta_gmm>beta_og)]<-beta_gmm[which(beta_gmm>beta_og)]
    }
    if(sum(var_gmm>var_og)>0){
      var_og[which(var_gmm>var_og)]<-var_gmm[which(var_gmm>var_og)]
    }
  }
  return(cbind(beta_ug,beta_og,s_ug=sqrt(var_ug),s_og=sqrt(var_og),int_ug,int_og))
}
