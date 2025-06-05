#' Externally informed Sargan-Hansen test for linear external moment functions
#'
#' @export

SH.test<-function(H,e_min,e_max,var.e,Omega="Sigma",asymptotic=T,alpha){
  if(asymptotic==T){
      switch(Omega,
             S=t_I_S(H,e_min,e_max,var.e)>qchisq(1-alpha,df=ncol(H)),
             Sigma=t_I_Sigma(H,e_min,e_max,var.e)>qchisq(1-alpha,df=ncol(H))
      )
  }else{
      t_I_S(H,e_min,e_max,var.e)*(nrow(H)-ncol(H))/(ncol(H)*(nrow(H)-1))>qf(1-alpha,df1=ncol(H),df2=nrow(H)-ncol(H))
  }
}
