t_I_S <- function(X,e_min,e_max,var.e){
  Dmat<-2*MASS:ginv(var(X)+var.e)
  Amat<-cbind(diag(-1,ncol(X)),diag(1,ncol(X)))
  bvec<-c(e_min-apply(X,FUN=mean,MARGIN=2),apply(X,FUN=mean,MARGIN=2)-e_max)
  return(quadprog::solve.QP(Dmat,rep(0,ncol(X)),Amat,bvec)$value*nrow(X))
}

t_I_Sigma <- function(X,e_min,e_max,var.e){
  Dmat<-2*var(X)*(nrow(X)-1)/nrow(X)+var.e
  Amat<-(Dmat/2)%*%cbind(diag(-1,ncol(X)),diag(1,ncol(X)))
  bvec<-c(e_min-apply(X,FUN=mean,MARGIN=2),apply(X,FUN=mean,MARGIN=2)-e_max)
  step1<-quadprog::solve.QP(Dmat,rep(0,ncol(X)),Amat,bvec)$value
  return(nrow(X)*step1/(1+step1))
}

