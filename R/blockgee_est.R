###Estimates for statistics Externally informed blockinvariant GEE Estimates with external values
#' @importFrom stats lm pchisq pf qt sd var

#estimate and variance of (externally informed) estimate - unrestricted, single step estimator
estimate_rmglm_ex<-function(data,between=NULL,within,dependent,link,ex_col=NULL,ex_val=NULL,var_ex=0){
  y<-data[,dependent]
  Z<-create_Z(as.matrix(data[,c(between,within)]))
  n<-ncol(matrix(y,nrow(Z)))
  if(link=="ordinal"){
    k=length(unique(y))
    y<-extend_y(y)
    ##Z
    Z_old<-Z
    Z<-c()
    for(ro in 1:nrow(Z_old)){
      Z<-rbind(Z,diag(k-1)%x%t(Z_old[ro,]))
    }
  }
  ex<-rep(Inf,nrow(Z))
  ex[ex_col]<-ex_val
  Y<-matrix(y,nrow(Z))
  H<-as.matrix(t(matrix(y,nrow(Z))-ex)[,ex!=Inf])

  if(is.null(ex_col)|is.null(ex_val)){
    beta_est<-solve(Z)%*%mu_inv(apply(matrix(y,nrow(Z)),mean,MARGIN=1),link)}else{
      beta_est<-solve(Z)%*%mu_inv((apply(Y,mean,MARGIN=1)-(Y%*%H/n)%*%solve(t(H)%*%H/n+var_ex)%*%apply(t(H),mean,MARGIN=1))/as.numeric(1-apply(H,mean,MARGIN=2)%*%solve(t(H)%*%H/n+var_ex)%*%apply(t(H),mean,MARGIN=1)),link)
    }

  D<-diag(as.vector(mu_deriv(Z%*%beta_est,link)))
  Y_mu<-matrix(y-as.vector(mu(Z%*%beta_est,link)),nrow(Z))
  Var_y<-Y_mu%*%t(Y_mu)/n
  if(is.null(ex_col)|is.null(ex_val)){
    var_est<-solve(D%*%Z)%*%Var_y%*%solve(t(Z)%*%D)/n}else{
      var_est<-solve(D%*%Z)%*%(Var_y-(Y_mu%*%H/n)%*%solve(t(H)%*%H/n+var_ex)%*%t((Y_mu%*%H/n)))%*%solve(t(Z)%*%D)/n
    }
  return(list(beta_est,var_est,n,Z))
}

#restricted estimate and restricted variance estimate (under nullhypothesis)
rmglm_restricted_step<-function(data,between=NULL,within,dependent,link,Cont,r=rep(0,nrow(Cont)),beta_start,eps=1e-12,ex_col=NULL,ex_val=NULL,var_ex=0){

  A<-MASS::Null(t(Cont)) #Basis for Null-space, Matrix has to be transposed!
  sol<-t(Cont)%*%solve(Cont%*%t(Cont))%*%r #specific solution

  y<-data[,dependent]
  if(link=="ordinal"){
    y<-extend_y(y)
  }

  #get structural matrices
  erg_dum<-estimate_rmglm_ex(data,between,within,dependent,link,ex_col,ex_val,var_ex)
  Z<-erg_dum[[4]]
  n<-erg_dum[[3]]
  if(sum(beta_start!=rep(0,ncol(Z)-nrow(Cont)))==0){
    beta_start<-solve(t(A)%*%A)%*%t(A)%*%(erg_dum[[1]]-sol)
  }

  D<-diag(as.vector(mu_deriv(Z%*%(A%*%beta_start+sol),link)))
  Y_mu<-matrix(y-as.vector(mu(Z%*%(A%*%beta_start+sol),link)),nrow(Z))
  Var_y<-Y_mu%*%t(Y_mu)/n
  if(is.null(ex_col)==F&is.null(ex_val)==F){
    ex<-rep(Inf,nrow(Z))
    ex[ex_col]<-ex_val
    H<-as.matrix(t(matrix(y,nrow(Z))-ex)[,ex!=Inf])
    var_yh<-(Y_mu%*%H/n)%*%solve(t(H)%*%H/n+var_ex)%*%t((Y_mu%*%H/n))
    cov_yh<-(Y_mu%*%H/n)%*%solve(t(H)%*%H/n+var_ex)%*%apply(t(H),mean,MARGIN=1)
  }else{
    var_yh<-0
    cov_yh<-0
  }

  var_est<-solve(D%*%Z)%*%(Var_y-var_yh)%*%solve(t(Z)%*%D)

  #compute step:
  #check for singular matrices:
  if(rcond(var_est)>3e-16){
    inv_var<-solve(var_est)
  }else{
    inv_var<-MASS::ginv(var_est)
  }
  if(rcond((t(A)%*%inv_var%*%A))>3e-16){
    W<-solve((t(A)%*%inv_var%*%A))
  }else{
    W<-MASS::ginv((t(A)%*%inv_var%*%A))
  }

  #define better starting value if external information is present

  beta_erg<-beta_start
  diver=F
  if(var_ex==0 & (length(ex_col)+nrow(Cont))==nrow(Z)){ #exception: if external info and Cont specify full parameter
    beta_erg<-A%*%beta_start+sol #only this solution is possible here
  }else{
    beta_dum<-beta_start
    maxit=0
    repeat{
      beta_start<-beta_erg
      beta_erg<-beta_start+W%*%t(A)%*%inv_var%*%solve(D%*%Z)%*%(apply(matrix(y-as.vector(mu(Z%*%(A%*%beta_start+sol),link)),nrow(Z)),mean,MARGIN=1)-cov_yh)/n
      if(sum(is.nan(beta_erg))>0|maxit>10000){beta_erg<-beta_dum;diver=T;break}
      if(sum((beta_start-beta_erg)^2)<eps){break}
      maxit=maxit+1
    }
    beta_start<-beta_erg
    beta_erg<-A%*%beta_start+sol
  }
  return(list(beta_erg,beta_start,var_est,diver))
}


