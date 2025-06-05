#' Test for general linear hypotheses in the externally informed linear model
#' @description Tests the general Null hypothesis: R %*% beta_0 = r for a full rank matrix R and vector r.
#'
#' @export
#' @importFrom stats lm pchisq pf qt sd var
#' @param X Numeric matrix containing independent variables in colums.
#' @param y Numeric vector containing dependent variable.
#' @param mom Character vector of moments, "EY" = expected value of y, "EXi" = expected value of ith column of X,
#'    "EXiY" = mixed moment of y and x_i, "EY2" = second moment of y, "VY" = variance of y,
#'    "CXiY" = covariance of y and x_i, "RXiY" = correlation of y and x_i, "BXiY" = simple linear regression slope
#'    of y on x_i. The index i must be replaced by a column number of X, e.g. "EX1Y" takes the first column of X.
#' @param I_min Numeric vector containing the lower bounds of external interval.
#' @param I_max Numeric vector containing the upper bounds of external interval.
#' @param var.e Numeric matrix equal to external variance.
#' @param grid_points Number of Grid points.
#' @param intercept Logical, indicates wheter intercept is used in model or not.
#' @param R Numeric matrix defining contrasts R for the nullhypothesis Rb=r.
#' @param r Numeric vector r for the nullhypothesis Rb=r.
#' @returns list containing a string describing the test result and a p-value.
#' @examples
#'
#'#overall F-test for simple linear model under external interval for EY:
#'
#'exlm_wald.test(rnorm(100),rnorm(100),c("EY"),c(-0.2),c(0.2),
#'grid_points=100,var.e=1/100,R=diag(2),r=c(0,0))
#'
#'#Simulation:
#' erg<-c()
#' for(i in 1:500){erg<-c(erg,exlm_wald.test(rnorm(100),rnorm(100),c("EY"),c(-0.2),c(0.2),
#' grid_points=100,var.e=1/100,R=diag(2),r=c(0,0))$p.value)}
#' hist(erg)

exlm_wald.test<-function(X,y,mom,I_min,I_max=I_min,grid_points=10001,intercept=TRUE,var.e,R,r){
  X<-matrix(X)
  positions<-as.numeric(gsub("\\D","",mom))
  if(intercept==T){positions<-positions+1}
  mom[mom!="EY2"]<-gsub("[[:digit:]]+", "", mom[mom!="EY2"])
  beta_ols_h<-as.vector(lm(y~X)$coefficients)
  if(intercept==TRUE){X<-cbind(1,X)}
  w_ug<-Inf
  g_l<-rep(list(0:(grid_points-1)),sum(I_max!=I_min))
  grid<-expand.grid(g_l)
  for (i_Ki in 1:(grid_points^sum(I_max!=I_min))) {
    E_I<-I_min
    E_I[I_max!=I_min]<-I_min[I_max!=I_min]+as.numeric(grid[i_Ki,])*(I_max-I_min)[I_max!=I_min]/(grid_points-1)
    sigma<-summary(lm(y~X[,-1]))[[6]]
    Omega_R<-calc_Omega_R_X(X,y,sigma,mom,beta_ols_h,positions)
    Omega_h<-calc_Omega_h(X,y,E_I,mom,sigma,beta_ols_h,positions,var.e)
    beta_gmm_dummy<-calc_beta_gmm(X,y,E_I,mom,Omega_R,Omega_h,positions)
    sigma<-sqrt(sum((y-X%*%beta_gmm_dummy)^2)/(nrow(X)-length(beta)))
    Omega_R<-calc_Omega_R_X(X,y,sigma,mom,beta_gmm_dummy,positions)
    Omega_h<-calc_Omega_h(X,y,E_I,mom,sigma,beta_gmm_dummy,positions,var.e)
    beta_gmm<-calc_beta_gmm(X,y,E_I,mom,Omega_R,Omega_h,positions)
    var_gmm<-calc_var_gmm(X,y,Omega_R,Omega_h,sigma)
    w_new<-t(R%*%beta_gmm-r)%*%MASS::ginv(R%*%var_gmm%*%t(R))%*%(R%*%beta_gmm-r)
    if(w_new<w_ug){
      w_ug<-w_new
    }
  }
  pval<-1-pchisq(as.numeric(w_ug),df=nrow(R))
  return(list(result=paste("The Gamma-maximin Wald test for linear models with",nrow(R),"degrees of freedom has p-value:",pval)
              ,p.value=pval))
}
