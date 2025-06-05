#' Externally informed Sargan-Hansen test for linear external moment functions
#'
#' @export
#' @importFrom stats lm pchisq pf qt sd var
#' @param H Numeric matrix containing external moment function values, rows are units, columns different moments.
#' @param e_min Numeric vector containing the lower bounds of external interval.
#' @param e_max Numeric vector containing the upper bounds of external interval.
#' @param var.e Numeric matrix equal to external variance.
#' @param Omega Character specifying estimator of Omega_h, "S" = sample variance matrix, "Sigma" = obvious estimator.
#' @param asymptotic Logical, if true asymptotic distribution is used, if false small sample test is used.
#' @returns list containing a string describing the test result and a p-value.
#' @examples
#' exSH(rnorm(100),c(-0.4),c(-0.3),1/30)
#'
#' exSH(cbind(rnorm(100),rnorm(100)),c(-0.2,-0.4),c(0.1,0.3),
#' Omega="S",asymptotic=FALSE,var.e=1/30*diag(2))

exSH<-function(H,e_min,e_max,var.e,Omega="Sigma",asymptotic=T){
  H<-as.matrix(H)
  if(asymptotic==T){
      pval<-switch(Omega,
             S=1-pchisq(t_I_S(H,e_min,e_max,var.e),df=ncol(H)),
             Sigma=1-pchisq(t_I_Sigma(H,e_min,e_max,var.e),df=ncol(H))
      )
  }else{
      pval<-1-pf(t_I_S(H,e_min,e_max,var.e)*(nrow(H)-ncol(H))/(ncol(H)*(nrow(H)-1)),df1=ncol(H),df2=nrow(H)-ncol(H))
  }
  if(asymptotic==T){df_rep<-ncol(H)}else{df_rep<-paste(ncol(H),"and",nrow(H)-ncol(H))}
 return(list(result=paste("The externally informed Sargan-Hansen test with",df_rep,"degrees of freedom and Omega estimator",Omega,"has p-value:",pval)
         ,p.value=pval))
}
