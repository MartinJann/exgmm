####Teststatistics for statistics Externally informed blockinvariant GEE Estimates with external values
#' @importFrom stats lm pchisq pf qt sd var
#Wald-test
wald.test_ex<-function(data,between=NULL,within,dependent,link,Cont,r=rep(0,nrow(Cont)),eps=1e-12,ex_col=NULL,ex_val=NULL,var_ex=0){
  erg<-estimate_rmglm_ex(data,between,within,dependent,link,ex_col,ex_val,var_ex)
  var_est<-Cont%*%(erg[[2]]*erg[[3]])%*%t(Cont)
  if(rcond(var_est)>3e-16){ #check if matrix is singular
    inv_var<-solve(var_est)
  }else{
    inv_var<-MASS::ginv(var_est)
  }
  chi_stat<-erg[[3]]*t(Cont%*%erg[[1]]-r)%*%inv_var%*%(Cont%*%erg[[1]]-r)
  df<-nrow(Cont)
  p_val<-1-pchisq(chi_stat,df)
  conv=T #non-iterative => always converges
  return(data.frame(chi_stat,df,p_val,conv))
}

#Langrange-Multiplier Test

LM.test_ex<-function(data,between=NULL,within,dependent,link,Cont,r=rep(0,nrow(Cont)),eps=1e-12,ex_col=NULL,ex_val=NULL,var_ex=0){
  y<-data[,dependent]
  if(link=="ordinal"){
    y<-extend_y(y)
  }
  #Var restricted
  #get structural matrices
  erg_dum<-estimate_rmglm_ex(data,between,within,dependent,link)
  Z<-erg_dum[[4]]
  n<-erg_dum[[3]]

  ##First Step estimator
  conv=T
  if(is.null(ex_col)==F&is.null(ex_val)==F){
    ex<-rep(Inf,nrow(Z))
    ex[ex_col]<-ex_val
    beta_1<-rmglm_restricted_step(data,between,within,dependent,link,Cont,r=rep(0,nrow(Cont)),beta_start=rep(0,ncol(Z)-nrow(Cont)),eps=eps,ex_col=ex_col,ex_val=ex_val,var_ex=var_ex)
    if(beta_1[[4]]==T){conv=F} #catch first step divergence
    beta_1<-beta_1[[1]]
    Y_mu<-matrix(y-as.vector(mu(Z%*%(beta_1),link)),nrow(Z))
    H<-as.matrix(t(matrix(y,nrow(Z))-ex)[,ex!=Inf])
    cov_yh<-(Y_mu%*%H/n)%*%solve(t(H)%*%H/n+var_ex)%*%apply(t(H),mean,MARGIN=1)
  }else{
    cov_yh<-0
  }
  ##Second step estimator
  beta_est<-solve(Z)%*%mu_inv(apply(matrix(y,nrow(Z)),mean,MARGIN=1)-cov_yh,link)

  #Variance estimate
  D<-diag(as.vector(mu_deriv(Z%*%beta_est,link)))
  Y_mu<-matrix(y-as.vector(mu(Z%*%beta_est,link)),nrow(Z))
  Var_y<-Y_mu%*%t(Y_mu)/n

  if(is.null(ex_col)==F&is.null(ex_val)==F){
    var_est<-solve(D%*%Z)%*%(Var_y-(Y_mu%*%H/n)%*%solve(t(H)%*%H/n+var_ex)%*%t((Y_mu%*%H/n)))%*%solve(t(Z)%*%D)
  }else{
    var_est<-solve(D%*%Z)%*%(Var_y)%*%solve(t(Z)%*%D)
  }
  if(rcond(Cont%*%(var_est)%*%t(Cont))>3e-16){
    Invert_w<-solve(Cont%*%(var_est)%*%t(Cont))}else{
      Invert_w<-MASS::ginv(Cont%*%(var_est)%*%t(Cont))
    }
  #compute test statistic
  chi_stat<-n*t(Cont%*%beta_est-r)%*%Invert_w%*%(Cont%*%beta_est-r)
  df<-nrow(Cont)
  p_val<-1-pchisq(chi_stat,df)
  return(data.frame(chi_stat,df,p_val,conv))
}

#D_RU Test
DRU.test_ex<-function(data,between=NULL,within,dependent,link,Cont,r=rep(0,nrow(Cont)),eps=1e-12,ex_col=NULL,ex_val=NULL,var_ex=0){
  y<-data[,dependent]
  if(link=="ordinal"){
    y<-extend_y(y)
  }
  #get structural matrices
  erg_dum<-estimate_rmglm_ex(data,between,within,dependent,link)
  Z<-erg_dum[[4]]
  n<-erg_dum[[3]]
  conv=T

  ##First Step restricted estimator
  if(is.null(ex_col)==F&is.null(ex_val)==F){
    ex<-rep(Inf,nrow(Z))
    ex[ex_col]<-ex_val
    beta_1<-rmglm_restricted_step(data,between,within,dependent,link,Cont,r=r,beta_start=rep(0,ncol(Z)-nrow(Cont)),eps=eps,ex_col=ex_col,ex_val=ex_val,var_ex=var_ex)
  }else{
    beta_1<-rmglm_restricted_step(data,between,within,dependent,link,Cont,r=r,beta_start=rep(0,ncol(Z)-nrow(Cont)),eps=eps,ex_col=NULL,ex_val=NULL,var_ex=0)
  }
  if(beta_1[[4]]==T){conv=F} #catch first step divergence
  ##Second Step estimator
  if(is.null(ex_col)==F&is.null(ex_val)==F){
    beta_2<-rmglm_restricted_step(data,between,within,dependent,link,Cont,r=r,beta_start=beta_1[[2]],eps=eps,ex_col=ex_col,ex_val=ex_val,var_ex=var_ex)
    if(beta_2[[4]]==T){conv=F} #catch first second step divergence
    beta_2<-beta_2[[1]]
    Y_mu<-matrix(y-as.vector(mu(Z%*%(beta_2),link)),nrow(Z))
    H<-as.matrix(t(matrix(y,nrow(Z))-ex)[,ex!=Inf])
    cov_yh<-(Y_mu%*%H/n)%*%solve(t(H)%*%H/n+var_ex)%*%apply(t(H),mean,MARGIN=1)
    Var_y<-Y_mu%*%t(Y_mu)/n
    var_est<-Var_y-(Y_mu%*%H/n)%*%solve(t(H)%*%H/n+var_ex)%*%t((Y_mu%*%H/n))
  }else{
    beta_2<-rmglm_restricted_step(data,between,within,dependent,link,Cont,r=r,beta_start=beta_1[[2]],eps=eps,ex_col=NULL,ex_val=NULL,var_ex=0)
    if(beta_2[[4]]==T){conv=F} #catch first second step divergence
    beta_2<-beta_2[[1]]
    Y_mu<-matrix(y-as.vector(mu(Z%*%(beta_2),link)),nrow(Z))
    cov_yh<-0
    var_est<-Y_mu%*%t(Y_mu)/n
  }

  g_r<-(apply(Y_mu,mean,MARGIN=1)-cov_yh)
  #check if matrix is singular
  if(rcond(var_est)>3e-16){
    Invert_w<-solve(var_est)}else{
      Invert_w<-MASS::ginv(var_est)
    }
  #compute test stat
  if(var_ex==0 & nrow(Cont)==1 & sum(Cont!=0)==1 & sum(r^2)==0){#avoid numeric problems through ginv in case of single parameter test
    exc<-which(Cont!=0)
    chi_stat<-t(g_r)[-exc]%*%Invert_w[-exc,-exc]%*%g_r[-exc]*n}else{
      chi_stat<-t(g_r)%*%Invert_w%*%g_r*n
    }
  df<-nrow(Cont)
  p_val<-1-pchisq(chi_stat,df)
  return(data.frame(chi_stat,df,p_val,conv))

}
